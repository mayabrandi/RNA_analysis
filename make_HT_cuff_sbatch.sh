#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage:
        RNA_analysis.sh <sample names> <path> <gtf file>

        Arguments:
	<sample names>
		-given as a comma delimited string
	<path>
		-path to dir with tophat_out_* dirs
        <gtf fie>
                - reference annotation in gtf format, used by cufflinks and HTseq"

  exit
fi


names=(`echo $1 | tr "," "\n"`)
path=$2
gtf_file=$3


for i in ${names[*]};do

echo "#!/bin/bash -l" > HT_cuff_$i.sh
echo "#SBATCH -A a2010003" >> HT_cuff_$i.sh
echo "#SBATCH -p node" >> HT_cuff_$i.sh
echo "#SBATCH -t 5:00:00" >> HT_cuff_$i.sh
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

samtools view $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam
python -m HTSeq.scripts.count -s no -q $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam $gtf_file > $path/tophat_out_$i/$i.counts
rm $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam
samtools index $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam
cufflinks -p 8 -G $gtf_file -o $path/tophat_out_$i/cufflinks_out_$i $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam" >> HT_cuff_$i.sh

done
