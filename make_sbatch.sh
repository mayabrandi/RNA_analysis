#!/bin/bash -l
#SBATCH -A a2010003
#SBATCH -p core
#SBATCH -t 00:15:00
#SBATCH -J make_sbatch
#SBATCH -e make_sbatch.err
#SBATCH -o make_sbatch.out
#SBATCH --qos=short
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL


names=$1
name_list=(`echo $names | tr "," "\n"`)
bedfile=$2
project_id=$3
config_file=$4
run_name=$5
path=$6
gtf_file=$7
WP=$8

cd $path

for i in ${name_list[*]};do

## read distribution
echo "#!/bin/bash -l" > ${i}_runEver_rd.sh
echo "#SBATCH -A a2010003" >> ${i}_runEver_rd.sh
echo "#SBATCH -p node" >> ${i}_runEver_rd.sh
echo "#SBATCH -t 50:00:00" >> ${i}_runEver_rd.sh
echo "#SBATCH -e Ever_rd_$i.err" >> ${i}_runEver_rd.sh
echo "#SBATCH -o Ever_rd_$i.out" >> ${i}_runEver_rd.sh
echo "#SBATCH -J Ever_rd_$i" >> ${i}_runEver_rd.sh
echo "#SBATCH --mail-type=ALL" >> ${i}_runEver_rd.sh
echo "#SBATCH --mail-user=maya.brandi@scilifelab.se" >> ${i}_runEver_rd.sh
echo "export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:"'$PYTHONPATH' >> ${i}_runEver_rd.sh
echo "export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:"'$PATH' >> ${i}_runEver_rd.sh
echo "module unload samtools" >> ${i}_runEver_rd.sh
echo "module load samtools/0.1.9" >> ${i}_runEver_rd.sh
echo "cd $path" >> ${i}_runEver_rd.sh
echo "read_distribution.py -i tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam -r $bedfile" >> ${i}_runEver_rd.sh

sbatch ${i}_runEver_rd.sh

## gene body coverage
echo "#!/bin/bash -l" > ${i}_runEver_gbc.sh
echo "#SBATCH -A a2010003" >> ${i}_runEver_gbc.sh
echo "#SBATCH -p node" >> ${i}_runEver_gbc.sh
echo "#SBATCH -t 15:00:00" >> ${i}_runEver_gbc.sh
echo "#SBATCH -e Ever_gbc_$i.err" >> ${i}_runEver_gbc.sh
echo "#SBATCH -o Ever_gbc_$i.out" >> ${i}_runEver_gbc.sh
echo "#SBATCH -J Ever_gbc_$i" >> ${i}_runEver_gbc.sh
echo "#SBATCH --mail-type=ALL" >> ${i}_runEver_gbc.sh
echo "#SBATCH --mail-user=maya.brandi@scilifelab.se" >> ${i}_runEver_gbc.sh
echo "export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:"'$PYTHONPATH' >> ${i}_runEver_gbc.sh
echo "export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:"'$PATH' >> ${i}_runEver_gbc.sh
echo "module unload samtools" >> ${i}_runEver_gbc.sh
echo "module load samtools/0.1.9" >> ${i}_runEver_gbc.sh
echo "cd $path" >> ${i}_runEver_gbc.sh
echo "samtools view /tophat_out_${i}/accepted_hits_sorted_dupRemoved_${i}.bam | geneBody_coverage.py -i - -r $bedfile -o $i" >> ${i}_runEver_gbc.sh
echo "R CMD BATCH ${i}.geneBodyCoverage_plot.r" >> ${i}_runEver_gbc.sh
echo "mv geneBody_coverage.pdf ${i}_geneBody_coverage.pdf > ${i}_runEver_gbc.sh" >> ${i}_runEver_gbc.sh

done



## quantify_rRNA
echo "#! /bin/bash -l" > quantify_rRNA.sh
echo "#SBATCH -A a2010003" >> quantify_rRNA.sh
echo "#SBATCH -p node" >> quantify_rRNA.sh
echo "#SBATCH -t 15:00:00" >> quantify_rRNA.sh
echo "#SBATCH -J quantify_rRNA" >> quantify_rRNA.sh
echo "#SBATCH -e quantify_rRNA.err" >> quantify_rRNA.sh
echo "#SBATCH -o quantify_rRNA.out" >> quantify_rRNA.sh
echo "#SBATCH --mail-user maya.brandi@scilifelab.se" >> quantify_rRNA.sh
echo "#SBATCH --mail-type=ALL" >> quantify_rRNA.sh
echo "cd $path
python /bubo/home/h7/junw/glob/quantify_rRNA.py $gtf_file" >> quantify_rRNA.sh

sbatch quantify_rRNA.sh



## correl
echo "#! /bin/bash -l" > correl.sh
echo "#SBATCH -A a2010003" >> correl.sh
echo "#SBATCH -p node" >> correl.sh
echo "#SBATCH -t 01:00:00" >> correl.sh
echo "#SBATCH -J correl" >> correl.sh
echo "#SBATCH -e correl.err" >> correl.sh
echo "#SBATCH -o correl.out" >> correl.sh
echo "#SBATCH --mail-user maya.brandi@scilifelab.se" >> correl.sh
echo "#SBATCH --mail-type=ALL" >> correl.sh
echo "cd $path
R CMD BATCH '--args ${name_list[*]}' $WP/correl.R" >> correl.sh

sbatch correl.sh


## analysis report
echo "#!/bin/bash -l" > analysis_report.sh
echo "cd $path" >> analysis_report.sh
echo "python $WP/analysis_report.py $run_name $project_id $names -c $config_file -r -s -d -f" >> analysis_report.sh


## Mapping Statistics
sbatch $WP/get_stat.sh $names $path 

chmod 777 analysis_report.sh

mkdir sbatch_scripts
mv *.sh sbatch_scripts

