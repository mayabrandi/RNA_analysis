#!/bin/bash
#	This script generates some basic statistics and a report for a RNA-seq project
#	and should be run from the intermediate directory of the project to be analysed. 
#	See the README for further information

if [ $# -ne 6 ]; then
  echo "Usage:
	stand in 'intermediate' and run

	RNA_analysis.sh <run dir> <project id> <bed file> <gtf file> <mail> <config_file>

Arguments:
        <run dir>
                - The name of the directory with the tophat_out_* -dirs.
		This is typically the same as the run name, such as
		20120323A_hiseq2000, but can be any name. The name of
		the run dir will also be the name set as the 'run name' 
		in the report. You might want to change this in the 
		rst file if your run dir name doesn't have an appropriate
		'run name'.

                - If set to 'm' the script will merge and do the analysis on 
                the merged runs

        <project id>
                - eg: M.Muurinen_11_01a

        <bed file>
                - reference gene model in bed format. Used by Ever-Seq
                to get gene body coverage and read distribution.

        <gtf fie>
                - reference annotation in gtf format, used by cufflinks and HTseq

	<mail>
		- mail adress for SLURM messages

	<config_file>
		- post_process.yaml"
  exit
fi

module unload python
module load python/2.7
export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:$PATH

run_dir=$1
project_id=$2
bedfile=$3
gtf_file=$4
mail=$5
config_file=$6
WP=/bubo/home/h24/mayabr/glob/RNA_analysis
path=`pwd`

DEPENDENCY_MERGE='afterok'
DEPENDENCY_HT='afterok'
DEPENDENCY='afterok'

if [ $run_dir = "m" ];then
	## get samplenames
	name_list=`for dir in *hiseq2000/tophat_out_*; do echo ${dir##*out_};done|sort -n|uniq`
	names=`echo $name_list|sed -e 's/ /,/g'`

	## megre old and new samples
	JOBID=`sbatch $WP/merge.sh $path | sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY_MERGE=$DEPENDENCY_MERGE:$JOBID

	## get names of samples to be merged
	rerun=`for dir in *hiseq2000/tophat_out_*; do echo ${dir##*out_};done|sort -n|uniq -d`

	## run HTseq and cufflinks on meged samples
	path=$path/merged
	for i in ${rerun[*]};do
		python $WP/make_HT_cuff.py $i $gtf_file $mail $path $config_file	
		JOBID=`sbatch --dependency=$DEPENDENCY_MERGE "HT_cuff_$i.sh"| sed -re 's/.+\s+([0-9]+)/\1/'`
		DEPENDENCY_HT=$DEPENDENCY_HT:$JOBID
		done
	if [ $DEPENDENCY_HT = 'afterok' ]; then
		dep=" --dependency=$DEPENDENCY_MERGE"
	else 
		dep=" --dependency=$DEPENDENCY_HT"
	fi
else

	## get samplenames
	name_list=`for dir in $run_dir/tophat_out_*; do echo ${dir##*out_};done|sort -n`
	names=`echo $name_list|sed -e 's/ /,/g'`
	dep=""
	path=$path/$run_dir
fi

sbatch$dep $WP/make_sbatch.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP $mail
