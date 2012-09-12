#! /bin/bash -l

#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 50:00:00
#SBATCH -J get_stat
#SBATCH -e get_stat.err
#SBATCH -o get_stat.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL


if [ $# -ne 2 ]; then
  echo "
Usage:
        get_stat.sh <sample names> <path>

Arguments:
        <sample names>
                - given as a comma delimited string
        <path>
                - path to the directory containing the tophat_out_* directories
                from Mikaels pipeline

Output:
	stat
		- a text file containing the tot # read pairs, % mapped reads and 
		% reads left after duplicates removed for each sample. The file 
		is formated to be used by analysis_report.py
	stat.json
		- a json file with the following information about each sample:
		sample name
		tot # read pairs
		# mapped (read 1)
		# mapped (read 2)
		% mapped reads
		# mapped after duplicates removed (read 1)
		# mapped after duplicates removed (read 2)
		% of reads that mapped after duplicates removed
		# uniquely mapped after duplicates removed
		% uniquely mapped
		# spliced alignments after duplicates removed
		% spliced 
"
else




export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:$PATH



module add bioinfo-tools
module load samtools 

samplenames=(`echo $1|sed -e 's/,/ /g'`)
readpairs=();j=0; C=(); D=(); E=(); F=(); G=(); H=(); I=(); J=(); L=(); M=()
path=$2
cd $path
echo $1
echo $2
echo '{' > stat.json
for u in ${samplenames[*]}; do
	echo $u
	counts=`grep 'reads have been filtered out' tophat_out_$u/logs/prep_reads.log|cut -f 4 -d ' '`
	echo $counts
	readpairs[j]=$counts
	samplenames[j]=$u
        arg=tophat_out_${u}/accepted_hits_${u}.bam
        #CD=`samtools flagstat $arg |grep read |awk '{print $1}'`
	CD=`bam_stat.py -i $arg 2>&1|grep "Read-"|awk '{print $2}'`
	echo $CD
        set -- $CD 
        C[j]=$1
        D[j]=$2
	E[j]=`echo "100*($1+$2)/(${readpairs[j]}*2)"|bc -l|python -c "print round(float(raw_input()),2)"`
	echo ${E[j]}
        arg=tophat_out_${u}/accepted_hits_sorted_dupRemoved_${u}.bam
        #FG=`samtools flagstat $arg |grep read |awk '{print $1}'`
	FG=`bam_stat.py -i $arg 2>&1|grep "Read-"|awk '{print $2}'`
        set -- $FG
        F[j]=$1
        G[j]=$2
	H[j]=`echo "100*($1+$2)/(${C[j]}+${D[j]})"|bc -l|python -c "print round(float(raw_input()),2)"`
	I[j]=`echo "100*($1+$2)/(${readpairs[j]}*2)"|bc -l|python -c "print round(float(raw_input()),2)"`
	echo 'duprem'
	echo ${I[j]}
        J[j]=`samtools view $arg |grep -c NH:i:1$`
        K[j]=`echo "100*${J[j]}/(${F[j]}+${G[j]})"|bc -l|python -c "print round(float(raw_input()),2)"`
        L[j]=`samtools view $arg |cut -f6 |grep -c N`
        M[j]=`echo "100*${L[j]}/(${F[j]}+${G[j]})"|bc -l|python -c "print round(float(raw_input()),2)"`
	
	# Saving statistics in json format
	if [ $j -ne $(($nr_samps-1)) ]; then  
		echo "${u} :{ \"tot_#_read_pairs\" : ${readpairs[j]} , \"#_mapped_R1\" : ${C[j]} , \"#_mapped_R2\" : ${D[j]} , \"%_mapped_reads\" : ${E[j]} , \"#_mapped_(dupRem)_R1\" : ${F[j]} , \"#_mapped_(dupRem)_R2\" : ${G[j]}, \"%_of_mapped_(dupRem)\" :${I[j]} , \"#_uniquely_mapped_(dupRem)\" :${J[j]} , \"%_uniquely_mapped\" : ${K[j]}, \"#_spliced_alignments_(dupRem)\" :${L[j]} , \"%_spliced\" :${M[j]}},"   >> stat.json
	else
		echo " ${u} :{ \"tot_#_read_pairs\" : ${readpairs[j]} , \"#_mapped_R1\" : ${C[j]} , \"#_mapped_R2\" : ${D[j]} , \"%_mapped_reads\" : ${E[j]} , \"#_mapped_(dupRem)_R1\" : ${F[j]} , \"#_mapped_(dupRem)_R2\" : ${G[j]}, \"%_of_mapped_(dupRem)\" :${I[j]} , \"#_uniquely_mapped_(dupRem)\" :${J[j]} , \"%_uniquely_mapped\" : ${K[j]}, \"#_spliced_alignments_(dupRem)\" :${L[j]} , \"%_spliced\" :${M[j]}}"   >> stat.json
	fi
	j=$(($j+1))
done

echo '}' >> stat.json
     
echo 'Sample' ${samplenames[*]}>>stat
echo 'tot_#_read_pairs' ${readpairs[*]}>>stat
echo '%_mapped_reads' ${E[*]}>>stat
echo '%_reads_left_after_dup_rem' ${I[*]}>>stat

fi
