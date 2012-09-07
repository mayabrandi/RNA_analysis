#! /bin/bash -l
#SBATCH -A a2010003
#SBATCH -p node
#SBATCH -t 01:00:00
#SBATCH -J correl
#SBATCH -e correl.err
#SBATCH -o correl.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL
cd /bubo/home/h24/mayabr/glob/RNA_analysis
R CMD BATCH '--args /bubo/home/h9/mikaelh/Rat_Baylor_2004_ensembl.bed' /correl.R
