#!/bin/bash

# 120125_sN188_0255_AD0H9PACXX J.Ericsson_11_01 /bubo/home/h9/mikaelh/mm9_NCBI37_Ensembl.bed
# 120127_SN1018_0062_BD0H2HACXX  M.Muurinen_11_01a /bubo/home/h9/mikaelh/hg19_GRCh37_Feb20009_Ensembl.bed

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
##K.Tammimies Prover::::::::
#names='13,14,15,16,17,18,19,20,21,22'
#counts='24640297,25693063,18565928,20280198,20172049,20782377,21250676,21958132,20883053,19533550'
#Emmas prover
#counts='21338082,15376734,18935492,75454092,15750546,20931060,20316170,18510935,17489946,18605141,23199951,17682959'
#name_list=(`echo $names | tr "," "\n"`)

path=/bubo/proj/a2010002/projects/`echo $project_id|tr '.' '_'|tr '[A-Z]' '[a-z]'`/intermediate


DEPENDENCY='afterok'

## Read Distribution
#$WP/make_rd_sbatch.sh $path $names $bedfile
#$WP/make_gbc_sbatch.sh $path $names $bedfile
#for name in ${name_list[*]};do
#	echo $name
#	sbatch $name"_runEver_gbc.sh"
#        JOBID=`sbatch $name"_runEver_rd.sh" | sed -re 's/.+\s+([0-9]+)/\1/'`
#	DEPENDENCY=$DEPENDENCY:$JOBID
#        rm $name"_runEver_rd.sh"
#	rm $name"_runEver_gbc.sh"
#done



## FPKM_PCAplot, FPKM_heatmap
JOBID=`sbatch $WP/correl.sh $path $names | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID

# Mapping Statistics
#sbatch $WP/get_stat.sh $path $names $counts
JOBID=`sbatch $WP/get_stat.sh $path $names $counts | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID

## Make report
#python $WP/analysis_report.py $fcID $project_id $names -c $config_file
sbatch --dependency=$DEPENDENCY $WP/analysis_report.sh $fcID $project_id $config_file $names
