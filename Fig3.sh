#!/bin/bash
#SBATCH --job-name=cellphonedb_Colangio                                 
#SBATCH --mail-type=END,FAIL                                 
#SBATCH --ntasks=12                     
#SBATCH --output=cellphonedb.log                                  
pwd; hostname; date
 
 
echo "Running cellphonedb"
module load conda/anaconda3 
source activate cellphonedb
cellphonedb method statistical_analysis ident_normal.txt matrix_normal.txt --counts-data hgnc_symbol --threshold 0.10  --debug-seed 123  --threads $SLURM_NTASKS  --iterations 500
exit
