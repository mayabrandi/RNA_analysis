#!/bin/bash
#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 1:00:00
#SBATCH -J analysis_rep
#SBATCH -e analysis_rep.err
#SBATCH -o analysis_rep.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL

project_id=$1
config_file=$2
sample_names=$3
run_name=$4

python /bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys_dev2/analysis_report.py $run_name $project_id $sample_names -c $config_file 
