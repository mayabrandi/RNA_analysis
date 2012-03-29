#! /bin/bash -l

#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 5:00:00
#SBATCH -J get_stat
#SBATCH -e get_stat.err
#SBATCH -o get_stat.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL

module add bioinfo-tools
module load samtools

path=$1 
echo $path

samplenames=(`echo $2|sed -e 's/,/ /g'`)
echo ${samplenames[*]}
readpairs=();j=0; C=(); D=(); E=(); F=(); G=(); H=(); I=(); J=(); L=(); M=()


echo '{' > stat.json
for u in ${samplenames[*]}; do
	echo $u
	echo 'hej'
	counts=`grep 'reads have been filtered out' $path/tophat_out_$u/logs/prep_reads.log|cut -f 4 -d ' '`
	echo $path/tophat_out_$u/logs/prep_reads.log
	readpairs[j]=$counts
	echo $counts
	samplenames[j]=$u
        arg=$path/tophat_out_${u}/accepted_hits_${u}.bam
        CD=`samtools flagstat $arg |grep read |awk '{print $1}'`
        set -- $CD 
        C[j]=$1
        D[j]=$2
        E[j]=`echo "scale=1;100*($1+$2)/(${readpairs[j]}*2)"|bc`
        arg=$path/tophat_out_${u}/accepted_hits_sorted_dupRemoved_${u}.bam
        FG=`samtools flagstat $arg |grep read |awk '{print $1}'`
        set -- $FG
        F[j]=$1
        G[j]=$2
        H[j]=`echo "scale=1;100*($1+$2)/(${C[j]}+${D[j]})"|bc`
        I[j]=`echo "scale=1;100*($1+$2)/(${readpairs[j]}*2)"|bc`
        J[j]=`samtools view $arg |grep -c NH:i:1$`
        K[j]=`echo "scale=1;100*${J[j]}/(${F[j]}+${G[j]})"|bc`
        L[j]=`samtools view $arg |cut -f6 |grep -c N`
        M[j]=`echo "scale=1;100*${L[j]}/(${F[j]}+${G[j]})"|bc`

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
#echo '#_mapped_R1' ${C[*]}>>stat_more
#echo '#_mapped_R2' ${D[*]}>>stat_more
echo '%_mapped_reads' ${E[*]}>>stat
#echo '#_mapped_(dupRem)_R1' ${F[*]}>>stat_more
#echo '#_mapped_(dupRem)_R2' ${G[*]}>>stat_more
#echo '%_of_mapped_(dupRem)' ${H[*]}>>stat_more
echo '%_reads_left_after_dup_rem' ${I[*]}>>stat
#echo '#_uniquely_mapped_(dupRem)' ${J[*]}>>stat_more
#echo '%_uniquely_mapped' ${K[*]}>>stat_more
#echo '#_spliced_alignments_(dupRem)' ${L[*]}>>stat_more
#echo '%_spliced' ${M[*]}>>stat_more

