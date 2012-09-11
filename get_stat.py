import sys

name=sys.argv[1]
f=open(name+"_get_stat.sh",'w')

print >>f, """#! /bin/bash -l

#SBATCH -A a2010002
#SBATCH -p core
#SBATCH -t 2:00:00
#SBATCH -J get_stat
#SBATCH -e get_stat"""+ name + """.err
#SBATCH -o get_stat"""+ name + """.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL


export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:$PATH

module add bioinfo-tools
module load samtools

name="""+name+"""
counts=`grep 'reads have been filtered out' tophat_out_${name}/logs/prep_reads.log|cut -f 4 -d ' '`
arg=tophat_out_${name}/accepted_hits_${name}.bam
CD=`bam_stat.py -i $arg 2>&1|grep "Read-"|awk '{print $2}'`
set -- $CD
C=$1
D=$2
E=`echo "100*($1+$2)/(${counts}*2)"|bc -l|python -c "print round(float(raw_input()),2)"`
arg=tophat_out_${name}/accepted_hits_sorted_dupRemoved_${name}.bam
FG=`bam_stat.py -i $arg 2>&1|grep "Read-"|awk '{print $2}'`
set -- $FG
F=$1
G=$2
H=`echo "100*($1+$2)/(${C}+${D})"|bc -l|python -c "print round(float(raw_input()),2)"`
I=`echo "100*($1+$2)/(${counts}*2)"|bc -l|python -c "print round(float(raw_input()),2)"`
L=`samtools view $arg |cut -f6 |grep -c N`
M=`echo "100*${L}/(${F}+${G})"|bc -l|python -c "print round(float(raw_input()),2)"`

echo 'Sample' ${name}>tophat_out_${name}/stat_${name}
echo 'tot_#_read_pairs' ${counts}>>tophat_out_${name}/stat_${name}
echo '%_uniquely_mapped_reads' ${E}>>tophat_out_${name}/stat_${name}
echo '%_uniquely_mapped_reads_left_after_dup_rem' ${I}>>tophat_out_${name}/stat_${name}
echo '#_uniquely_mapped_R1' ${C}>>tophat_out_${name}/stat_${name}
echo '#_uniquely_mapped_R2' ${D}>>tophat_out_${name}/stat_${name}
echo '#_uniquely_mapped_(dupRem)_R1' ${F}>>tophat_out_${name}/stat_${name}
echo '#_uniquely_mapped_(dupRem)_R2' ${G}>>tophat_out_${name}/stat_${name}
echo '#_spliced_alignments_(dupRem)' ${L}>> tophat_out_${name}/stat_${name}
echo '%_spliced' ${M}>>tophat_out_${name}/stat_${name}

"""
