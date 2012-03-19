#!/bin/bash

WP='/bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys'

config_file=/bubo/home/h24/mayabr/config/post_process.yaml
args=("$@")
fcID=${args[0]}
project_id=${args[1]}
bedfile=${args[2]}

names_counts=(`python $WP/snrc.py $fcID $project_id| tr ":" "\n"`)
names=${names_counts[0]}
name_list=(`echo ${names_counts[0]} | tr "," "\n"`)
counts=${names_counts[1]} 
path=/bubo/proj/a2010002/projects/`echo $project_id|tr '.' '_'|tr '[A-Z]' '[a-z]'`/intermediate


DEPENDENCY='afterok'

# Read Distribution
$WP/make_rd_sbatch.sh $path $names $bedfile
$WP/make_gbc_sbatch.sh $path $names $bedfile
for name in ${name_list[*]};do
	sbatch $name"_runEver_gbc.sh"
        JOBID=`sbatch $name"_runEver_rd.sh" | sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY=$DEPENDENCY:$JOBID
        rm $name"_runEver_rd.sh"
	rm $name"_runEver_gbc.sh"
done

# FPKM_PCAplot, FPKM_heatmap
JOBID=`sbatch $WP/correl.sh $path $names | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID

# Mapping Statistics
JOBID=`sbatch $WP/get_stat.sh $path $names $counts | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID

## Make report
sbatch --dependency=$DEPENDENCY $WP/analysis_report.sh $fcID $project_id $config_file $names
