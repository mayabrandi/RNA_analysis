#!/bin/bash -l
#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 15:00:00
#SBATCH -J merge
#SBATCH -e merge.err
#SBATCH -o merge.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL

path=$1
gtf_path=$2
flowcells=(`ls $path|grep hiseq2000`)
nr_fc=${#flowcells[*]}
name_list=`for dir in *hiseq2000/tophat_out_*; do echo ${dir##*out_};done|sort -n|uniq`

if [ $nr_fc  -lt 2 -o $# -lt 1 ]; then
echo "
Usage:
	merge.sh <path> [gtf file]

Arguments:
	<path> 
	- Path to the intermediate directory containing the runs to be merged. The
	script merges all directories with nameds containg 'hiseq2000'. 
       	- eg: /proj/a2010002/projects/m_muurinen_11_01a/intermediate, containing 
	the three runs: 20120127B_hiseq2000 20120228A_hiseq2000 20120328A_hiseq2000

       	[gtf fie]
	- Optional! If not given, the script will only merge, but not run cufflinks 
	and HTseq on the merged data.
        - Reference annotation in gtf format, used by cufflinks and HTseq

Output:
	- A new directory for the merged data, called merged and placed in the intermediate directory.
	"
else
	# merge first two flowcells
	mkdir $path/merged
	cd $path/${flowcells[0]}
	for samp_dir in tophat_out_*;do
		mkdir $path/merged/$samp_dir $path/merged/$samp_dir/logs
		cd $samp_dir
		bam_file_dupRem=`echo accepted_hits_sorted_dupRemoved_*.bam`
		bam_file=`echo $bam_file_dupRem|sed -e 's/sorted_dupRemoved_//g'`
		echo $samp_dir
		if [ -e $path/${flowcells[1]}/$samp_dir ]; then
			echo "merge"
			## merge read
			samtools merge $path/merged/$samp_dir/$bam_file_dupRem $path/${flowcells[0]}/$samp_dir/$bam_file_dupRem $path/${flowcells[1]}/$samp_dir/$bam_file_dupRem
			samtools merge $path/merged/$samp_dir/$bam_file $path/${flowcells[0]}/$samp_dir/$bam_file $path/${flowcells[1]}/$samp_dir/$bam_file
	
			## sum read counts
                	counts_0=(`grep 'reads have been filtered out' $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                	counts_1=(`grep 'reads have been filtered out' $path/${flowcells[1]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)                    
			counts=$((${counts_0[1]}+${counts_1[1]}))
                	sorted=$((${counts_0[0]}+${counts_1[0]}))
                	sed -e "s/${counts_0[1]}/$counts/g" $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|sed -e "s/${counts_0[0]}/$sorted/g" > $path/merged/$samp_dir/logs/prep_reads.log
		else
			echo "cp"
			## copy reads and counts
			cp -r $path/${flowcells[0]}/$samp_dir $path/merged
		fi
		cd ..
	done

	# merge all remaining flowcells with "merged"
	if [ $nr_fc -gt 2 ]; then
	echo "merge all remaining flowcells with merged"
	for ((i=2; i<$nr_fc; i++));do
		echo $i
		cd $path/merged
		pwd
		ls
 	        for samp_dir in tophat_out_*;do
	             	cd $samp_dir           
			pwd     
			bam_file_dupRem=`echo accepted_hits_sorted_dupRemoved_*.bam` 
			bam_file=`echo $bam_file_dupRem|sed -e 's/sorted_dupRemoved_//g'`
			if [ -e $path/${flowcells[$i]}/$samp_dir ]; then
				echo "merge"
				echo $samp_dir
				echo $path/${flowcells[$i]}/$samp_dir/
				ls $path/${flowcells[$i]}/$samp_dir
				echo ""
				echo $path/merged/$samp_dir
				ls $path/merged/$samp_dir/
				## merge read
				samtools merge $path/merged/$samp_dir/temp_${bam_file_dupRem} $bam_file_dupRem $path/${flowcells[$i]}/$samp_dir/$bam_file_dupRem
				samtools merge $path/merged/$samp_dir/temp_${bam_file} $bam_file $path/${flowcells[$i]}/$samp_dir/$bam_file
				mv $path/merged/$samp_dir/temp_${bam_file_dupRem} $path/merged/$samp_dir/$bam_file_dupRem
                        	mv $path/merged/$samp_dir/temp_${bam_file} $path/merged/$samp_dir/$bam_file
				## sum read counts
                        	counts_0=(`grep 'reads have been filtered out' $path/merged/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        	counts_1=(`grep 'reads have been filtered out' $path/${flowcells[$i]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        	counts=$((${counts_0[1]}+${counts_1[1]}))
                        	sorted=$((${counts_0[0]}+${counts_1[0]}))
                        	sed -e "s/${counts_0[1]}/$counts/g" $path/merged/$samp_dir/logs/prep_reads.log|sed -e "s/${counts_0[0]}/$sorted/g" > $path/merged/$samp_dir/logs/temp_prep_reads.log
                        	mv $path/merged/$samp_dir/logs/temp_prep_reads.log $path/merged/$samp_dir/logs/prep_reads.log		
                        	echo "" 
			        echo "after merge"
                                echo $samp_dir
                                echo $path/${flowcells[$i]}/$samp_dir/
                                ls $path/${flowcells[$i]}/$samp_dir
                                echo ""
                                echo $path/merged/$samp_dir
                                ls $path/merged/$samp_dir/
	
			fi
			cd ..
		done
		cd ..
        done	
	fi
	pwd
	if [ "$gtf_path" != "" ]; then
	echo "## run HTseq and cufflinks on meged samples"
	echo $gtf_path
	echo "## run HTseq and cufflinks on meged samples"
        ## get names of samples to be merged
        first_flowcell=${flowcells[0]}
        rerun=""
        for dir in ${flowcells[*]};do if [ $dir != $first_flowcell ]; then rerun=$rerun" "`for i in $path/$dir/tophat_out_*; do echo ${i##*out_};done`;fi ;done
	echo ${rerun[*]}
        rerun=(`echo $rerun|sed -e 's/ /\n/g'|sort -n|uniq`)
	echo ${rerun[*]}
	
	## run HTseq and cufflinks on meged samples
	for i in ${rerun[*]};do
        echo "#!/bin/bash -l">"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"
        echo "#SBATCH -A a2010003">>"HT_cuff_$i.sh"
        echo "#SBATCH -p core">>"HT_cuff_$i.sh"
        echo "#SBATCH -t 00:15:00">>"HT_cuff_$i.sh"
        echo "#SBATCH -e HT_cuff_$i.err">>"HT_cuff_$i.sh"
        echo "#SBATCH -o HT_cuff_$i.out">>"HT_cuff_$i.sh"
        echo "#SBATCH -J HT_cuff_$i">>"HT_cuff_$i.sh"
	echo "#SBATCH --qos=short">>"HT_cuff_$i.sh"
        echo "#SBATCH --mail-type=ALL">>"HT_cuff_$i.sh"
        echo "#SBATCH --mail-user=maya.brandi@scilifelab.se">>"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"

        echo "module unload htseq">>"HT_cuff_$i.sh"
        echo "module load htseq/0.5.1">>"HT_cuff_$i.sh"
        echo "module unload cufflinks">>"HT_cuff_$i.sh"
        echo "module load cufflinks/1.2.1">>"HT_cuff_$i.sh"
        echo "">>"HT_cuff_$i.sh"

        echo "samtools view $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
        echo "samtools view $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam |sort > $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
        echo "python -m HTSeq.scripts.count -s no -q $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam $gtf_path > $path/merged/tophat_out_$i/$i.counts">>"HT_cuff_$i.sh"
        echo "rm $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_prehtseq_$i.sam">>"HT_cuff_$i.sh"
        echo "samtools index $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam">>"HT_cuff_$i.sh"
        echo "cufflinks -p 8 -G $gtf_path -o $path/merged/tophat_out_$i/cufflinks_out_$i $path/merged/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam">>"HT_cuff_$i.sh"
	echo $i
	echo "HT_cuff_$i.sh"

        sbatch HT_cuff_$i.sh
        done
	fi
fi
