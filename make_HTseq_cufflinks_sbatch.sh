#!/bin/bash -l

path=$2
samplenames=(`echo $1 | tr "," "\n"`)
gtf_path=$3
for i in ${samplenames[*]};do
        echo "#!/bin/bash -l">>"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"
        echo "#SBATCH -A a2010002">>"HT_cuff_$i.sh"
        echo "#SBATCH -p node">>"HT_cuff_$i.sh"
        echo "#SBATCH -t 15:00:00">>"HT_cuff_$i.sh"
        echo "#SBATCH -e HT_cuff_$i.err">>"HT_cuff_$i.sh"
        echo "#SBATCH -o HT_cuff_$i.out">>"HT_cuff_$i.sh"
        echo "#SBATCH -J HT_cuff_$i">>"HT_cuff_$i.sh"
        echo "#SBATCH --mail-type=ALL">>"HT_cuff_$i.sh"
        echo "#SBATCH --mail-user=maya.brandi@scilifelab.se">>"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"

        echo "module unload htseq">>"HT_cuff_$i.sh"
        echo "module load htseq/0.5.1">>"HT_cuff_$i.sh"
        echo "module unload cufflinks">>"HT_cuff_$i.sh"
        echo "module load cufflinks/1.2.1">>"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"

	echo "samtools view $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
	echo "samtools view $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
	echo "python -m HTSeq.scripts.count -s no -q accepted_hits_sorted_dupRemoved_prehtseq_$i.sam $gtf_path > $path/tophat_out_$i/$i.counts">>"HT_cuff_$i.sh"
	echo "rm accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
	echo "samtools index accepted_hits_sorted_dupRemoved_$i.bam">>"HT_cuff_$i.sh"
	echo "cufflinks -p 8 -G $gtf_path -o cufflinks_out_$i accepted_hits_sorted_dupRemoved_$i.bam">>"HT_cuff_$i.sh"

done

