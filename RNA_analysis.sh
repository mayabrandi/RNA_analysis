#!/bin/bash
#	This script generates some basic statistics and a report for a RNA-seq project 
#	See the README for further information

if [ $# -ne 4 ]; then
  echo "Usage:
	<flowcell id> <project id> <bed file> <gtf file>
	If this is a merging run, set <flowcell id> = 'm'
	If this is not merging run, set <gtf file> = '-'"
  exit
fi

fcID=$1
project_id=$2
bedfile=$3
gtf_file=$4

run_name=20`echo $fcID|cut -f 1 -d '_'``echo ${fcID##*_}|cut -c 1`_hiseq2000
WP=/bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys
config_file=/bubo/home/h24/mayabr/config/post_process.yaml
path=/proj/a2012043/private/nobackup/projects/`echo $project_id|tr '.' '_'|tr '[A-Z]' '[a-z]'`/intermediate
#path=/proj/a2010002/projects/`echo $project_id|tr '.' '_'|tr '[A-Z]' '[a-z]'`/intermediate

DEPENDENCY_MERGE='afterok'
DEPENDENCY_HT_cuff='afterok'
DEPENDENCY='afterok'

if [ ! -d $path ];then echo "Error: No such directory $path. 
Check that flowcell ID and project name are correctly given";exit; fi

## Fetch sample names
first_flowcell=(`ls $path|grep hiseq2000`)
name_list=`for dir in $path/$first_flowcell/tophat_out_*; do echo ${dir##*out_};done|sort -n`
names=`echo $name_list|sed -e 's/ /,/g'`

## Merge prevoius runs
if [ $fcID = "m" ];then
	## Merge prevoius runs
	JOBID=`sbatch $WP/merge.sh $path | sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY_MERGE=$DEPENDENCY_MERGE:$JOBID
	$WP/make_HTseq_cufflinks_sbatch.sh $names $path $gtf_file
	for name in ${name_list[*]};do
        	JOBID=`sbatch --dependency=$DEPENDENCY_MERGE HT_cuff_${name}.sh | sed -re 's/.+\s+([0-9]+)/\1/'`
        	DEPENDENCY_HT_cuff=$DEPENDENCY_HT_cuff:$JOBID
	done
	dep=" --dependency=$DEPENDENCY_HT_cuff"
	path=${path}/merged
else
	dep=""
	path=${path}/$run_name
fi

# Read Distribution
$WP/make_rd_sbatch.sh $path $names $bedfile
$WP/make_gbc_sbatch.sh $path $names $bedfile
for name in ${name_list[*]};do
#	sbatch $name"_runEver_gbc.sh"
        JOBID=`sbatch$dep $name"_runEver_rd.sh" | sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY=$DEPENDENCY:$JOBID
done

## FPKM_PCAplot, FPKM_heatmap
JOBID=`sbatch$dep $WP/correl.sh $path $names | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID


## Mapping Statistics
JOBID=`sbatch$dep $WP/get_stat.sh $path $names | sed -re 's/.+\s+([0-9]+)/\1/'`
DEPENDENCY=$DEPENDENCY:$JOBID

## Make report
sbatch --dependency=$DEPENDENCY $WP/analysis_report.sh $project_id $config_file $names $run_name
