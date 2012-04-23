#!/bin/bash
#	This script generates some basic statistics and a report for a RNA-seq project
#	and should be run from the intermediate directory of the project to be analysed. 
#	See the README for further information

if [ $# -ne 4 ]; then
  echo "Usage:
	stand in 'intermediate' and run

	RNA_analysis.sh <run dir> <project id> <bed file> <gtf file>

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
                - reference annotation in gtf format, used by cufflinks and HTseq"
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
#run_name=20`echo $fcID|cut -f 1 -d '_'``echo ${fcID##*_}|cut -c 1`_hiseq2000
WP=/bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys
config_file=/bubo/home/h24/mayabr/config/post_process.yaml
path=`pwd`

DEPENDENCY_MERGE='afterok'
DEPENDENCY_HT='afterok'
DEPENDENCY='afterok'

if [ $run_dir = "m" ];then

	## get samplenames
	flowcells=(`ls|grep hiseq2000`)
	first_flowcell=${flowcells[0]}
	name_list=`for dir in $first_flowcell/tophat_out_*; do echo ${dir##*out_};done|sort -n`
	names=`echo $name_list|sed -e 's/ /,/g'`

	## megre old and new samples
	JOBID=`sbatch $WP/merge.sh $path | sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY_MERGE=$DEPENDENCY_MERGE:$JOBID

	## get names of samples to be merged
	rerun=""
	for dir in ${flowcells[*]};do if [ $dir != $first_flowcell ]; then rerun=$rerun" "`for i in $path/$dir/tophat_out_*; do echo ${i##*out_};done`;fi ;done
	rerun=(`echo $rerun|sed -e 's/ /\n/g'|sort -n|uniq`)

	## run HTseq and cufflinks on meged samples
	for i in ${rerun[*]};do
		
echo "#!/bin/bash -l" > HT_cuff_$i.sh
echo "#SBATCH -A a2010003" >> HT_cuff_$i.sh
echo "#SBATCH -p node" >> HT_cuff_$i.sh
echo "#SBATCH -t 20:00:00" >> HT_cuff_$i.sh
echo "#SBATCH -e HT_cuff_$i.err" >> HT_cuff_$i.sh
echo "#SBATCH -o HT_cuff_$i.out" >> HT_cuff_$i.sh
echo "#SBATCH -J HT_cuff_$i" >> HT_cuff_$i.sh
#echo "#SBATCH --qos=short" >> HT_cuff_$i.sh
echo "#SBATCH --mail-type=ALL" >> HT_cuff_$i.sh
echo "#SBATCH --mail-user=maya.brandi@scilifelab.se" >> HT_cuff_$i.sh

echo "module unload htseq
module load htseq/0.5.1
module unload cufflinks
module load cufflinks/1.2.1

samtools view $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam
python -m HTSeq.scripts.count -s no -q $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam $gtf_file > $path/merged/tophat_out_$i/$i.counts
rm $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam
samtools index $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam
cufflinks -p 8 -G $gtf_file -o $path/merged/tophat_out_$i/cufflinks_out_$i $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam" >> HT_cuff_$i.sh
	
		JOBID=`sbatch --dependency=$DEPENDENCY_MERGE "HT_cuff_$i.sh"| sed -re 's/.+\s+([0-9]+)/\1/'`
		DEPENDENCY_HT=$DEPENDENCY_HT:$JOBID
		done

	dep=" --dependency=$DEPENDENCY_HT"
	path=$path/merged
	echo $i
else

	## get samplenames
	name_list=`for dir in $run_dir/tophat_out_*; do echo ${dir##*out_};done|sort -n`
	names=`echo $name_list|sed -e 's/ /,/g'`
	dep=""
	path=$path/$run_dir
fi

sbatch$dep $WP/make_sbatch.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP

