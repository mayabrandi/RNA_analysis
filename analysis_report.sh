#!/bin/bash
#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 1:00:00
#SBATCH -J analysis_rep
#SBATCH -e analysis_rep.err
#SBATCH -o analysis_rep.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL


args=("$@")
fcID=${args[0]}
project_id=${args[1]}
config_file=${args[2]}
sample_names=${args[3]}

python /bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys/analysis_report.py $fcID $project_id $sample_names -c $config_file 
