
#initialize and load libraries
Sys.setlocale("LC_NUMERIC", "C")
library("Seurat");
library("SCENIC");
library("loomR");
library("hdf5r");
library("SCopeLoomR");
library("AUCell");
library("RcisTarget");
library("doRNG");
library("doMC");
library("RColorBrewer");

# set variables
processori=20

# set path
setwd("/mnt/isilon/bioinformatica/pkunderfranco/data/P27_Progetto_scRNASeq_ALeo/SCENIC/Cluster_Condition_Analysis_alra/")

#  import Seurat loom file
loom <- open_loom("../Cluster_Analysis_alra/seurat_alra.loom", mode="r")
exprMat <- get_dgem(loom)
cellInfo <- get_cellAnnotation(loom)
close_loom(loom)
dim(exprMat)

# add cell info/phenodata
head(cellInfo)
cellInfo <- data.frame(cellInfo)
cellInfo<-within(cellInfo, Condizione <- paste(ExpGroup,integrated_snn_res_0_8,sep='_'))
cellTypeColumn <- "Condizione"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "Condizione"
cbind(table(cellInfo$Condizione))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables
colVars <- list(Condizione=c("T_0"="#8B0000","T_1"="#8B0000","T_2"="#8B0000","T_3"="#8B0000","T_4"="#8B0000","T_5"="#8B0000","T_6"="#8B0000","T_7"="#8B0000","T_8"="#8B0000","T_9"="#8B0000","T_10"="#8B0000",
                             "S_0"="#48D1CC", "S_1"="#48D1CC","S_2"="#48D1CC","S_3"="#48D1CC","S_4"="#48D1CC","S_5"="#48D1CC","S_6"="#48D1CC","S_7"="#48D1CC","S_8"="#48D1CC","S_9"="#48D1CC","S_10"="#48D1CC"));
colVars$Condizione <- colVars$Condizione[intersect(names(colVars$Condizione), cellInfo$Condizione)]
colVars <- list(ExpGroup=c("S"="#8B0000","T"="#48D1CC"));
colVars$ExpGroup <- colVars$ExpGroup[intersect(names(colVars$ExpGroup), cellInfo$ExpGroup)]

# Initialize SCENIC settings
org="hgnc"
dbDir="../cisTarget_databases" 
myDatasetTitle="P27"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=processori) 
scenicOptions@settings$devType="png"

# Set options
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- processori
scenicOptions@settings$seed <- 123
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Gene filter/selection
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
write.table(genesKept,file="genesKept.xls")

# run correlation
runCorrelation(exprMat_filtered, scenicOptions)

# run GENIE
runGenie3(exprMat_filtered, scenicOptions)

# Build and score the GRN
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# binarize the network activity
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
runSCENIC_4_aucell_binarize(scenicOptions)

# Regulators for clusters or known cell types
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Condizione),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, fontsize_row=6, 
                   colors = "RdYlBu", breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "Condizione", "RelativeActivity");
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),];
write.table(topRegulators, file = "regulonActivity_byCellType.xls", sep = "\t", eol = "\n", row.names = F, col.names = T, quote = F);

# Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- .5
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$Condizione), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset,fontsize_row=10, 
                   colors = "RdYlBu", breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
write.table(topRegulators, file = "regulonActivityBinary_byCellType_0.5.xls", sep = "\t", eol = "\n", row.names = F, col.names = T, quote = F);

# Export to loom/SCope
scenicOptions@fileNames$output["loomFile",] <- "output/P27_SCENIC.loom"
fileName <- getOutName(scenicOptions, "loomFile")
export2scope(scenicOptions, exprMat_filtered)

# export motifs and regulons
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
#tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Rarb"]
write.table(motifEnrichment_selfMotifs_wGenes, file = "enrichment analysis.xls", sep = "\t", eol = "\n", row.names = F, col.names = T, quote = F);

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
write.table(regulonTargetsInfo, file = "regulonTargetsInfo.xls", sep = "\t", eol = "\n", row.names = F, col.names = T, quote = F);
