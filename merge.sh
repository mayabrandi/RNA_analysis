#! /bin/bash -l
#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 10:00:00
#SBATCH -J merge
#SBATCH -e merge.err
#SBATCH -o merge.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL

# This script takes as agrument the path to the intermediate directory containing the runs to be merged 
# Eg /proj/a2010002/projects/m_muurinen_11_01a/intermediate, containing the three runs: 
#			20120127B_hiseq2000 20120228A_hiseq2000 20120328A_hiseq2000
# It creates a new directory for the merged data, called merged 


path=$1

nr_fc=`ls $path|grep hiseq2000|wc -l`
fc=`ls $path|grep hiseq2000`
set -- $fc
flowcells=("$@")

if [ $nr_fc -lt 2 ]; then 
	echo "No files to merge"
else 
	## merge first two flowcells
	mkdir $path/merged
	for samp_dir_path in $path/${flowcells[0]}/tophat_out_*;do
		samp_dir=`echo ${samp_dir_path##*/}`
		name=`echo ${samp_dir##*out_}`
		bam_file_dupRem=accepted_hits_sorted_dupRemoved_$name.bam
		bam_file=accepted_hits_$name.bam
		mkdir $path/merged/$samp_dir $path/merged/$samp_dir/logs
		if [ -e $path/${flowcells[1]}/$samp_dir ]; then
			## merge read
			echo $name
			echo '#merge'
			samtools merge $path/merged/$samp_dir/$bam_file_dupRem $path/${flowcells[0]}/$samp_dir/$bam_file_dupRem $path/${flowcells[1]}/$samp_dir/$bam_file_dupRem
			samtools merge $path/merged/$samp_dir/$bam_file $path/${flowcells[0]}/$samp_dir/$bam_file $path/${flowcells[1]}/$samp_dir/$bam_file

                        ## sum read counts
                        counts_0=(`grep 'reads have been filtered out' $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        counts_1=(`grep 'reads have been filtered out' $path/${flowcells[1]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        counts=$((${counts_0[1]}+${counts_1[1]}))
                        sorted=$((${counts_0[0]}+${counts_1[0]}))
                        sed -e "s/${counts_0[1]}/$counts/g" $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|sed -e "s/${counts_0[0]}/$sorted/g" > $path/merged/$samp_dir/logs/prep_reads.log
		else
			## copy reads and counts
			echo $name                        
			echo '#cp'
			cp $path/${flowcells[0]}/$samp_dir/$bam_file_dupRem $path/merged/$samp_dir/$bam_file_dupRem
			cp $path/${flowcells[0]}/$samp_dir/$bam_file $path/merged/$samp_dir/$bam_file
			cp $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log $path/merged/$samp_dir/logs/
		fi
	done

	## merge all remaining flowcells with "merged"
	if [ $nr_fc -gt 2 ]; then
	for ((i=2; i<$nr_fc; i++));do
		mkdir $path/merged_temp
 	        for samp_dir_path in $path/merged/tophat_out_*;do
			samp_dir=`echo ${samp_dir_path##*/}`
                	name=`echo ${samp_dir##*out_}`
                	bam_file_dupRem=accepted_hits_sorted_dupRemoved_$name.bam
                	bam_file=accepted_hits_$name.bam
			mkdir $path/merged_temp/$samp_dir $path/merged_temp/$samp_dir/logs
			if [ -e $path/${flowcells[$i]}/$samp_dir ]; then
				## merge read
		                echo $name                        
				echo '#merge'

				samtools merge $path/merged_temp/$samp_dir/$bam_file_dupRem $samp_dir_path/$bam_file_dupRem $path/${flowcells[$i]}/$samp_dir/$bam_file_dupRem
				samtools merge $path/merged_temp/$samp_dir/$bam_file $samp_dir_path/$bam_file $path/${flowcells[$i]}/$samp_dir/$bam_file

                        	## sum read counts
                        	counts_0=(`grep 'reads have been filtered out' $path/merged/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        	counts_1=(`grep 'reads have been filtered out' $path/${flowcells[$i]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        	counts=$((${counts_0[1]}+${counts_1[1]}))
                        	sorted=$((${counts_0[0]}+${counts_1[0]}))
                        	sed -e "s/${counts_0[1]}/$counts/g" $path/merged/$samp_dir/logs/prep_reads.log|sed -e "s/${counts_0[0]}/$sorted/g" > $path/merged_temp/$samp_dir/logs/prep_reads.log
			else
				## copy reads and counts
				echo $name
				echo '#cp'
				cp $samp_dir_path/$bam_file_dupRem $path/merged_temp/$samp_dir/$bam_file_dupRem
				cp $samp_dir_path/$bam_file $path/merged_temp/$samp_dir/$bam_file
				cp $samp_dir_path/logs/prep_reads.log $path/merged_temp/$samp_dir/logs/
			fi
		done
		rm -r $path/merged
		mv $path/merged_temp/ $path/merged
        done	
	fi
fi
