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
echo names $names
bedfile=$2
echo bedfile $bedfile
project_id=$3
echo project_id $project_id
config_file=$4
echo config_file $config_file
run_name=$5
echo run_name $run_name
path=$6
echo path $path
gtf_file=$7
echo gtf_file $gtf_file
WP=$8
echo WP $WP
mail=$9
echo mail $mail
name_list=(`echo $names | tr "," "\n"`)

echo 'WP'
echo $WP
echo 'path'
echo $path
cd $path
pwd
for i in ${name_list[*]};do

## read distribution
python $WP/make_RseqQc_rd.py $i $bedfile $mail $config_file $path
sbatch RSeQC_${i}_rd.sh

## gene body coverage
python $WP/make_RseqQc_gbc.py $i $bedfile $mail $config_file $path

## statistics
python $WP/get_stat.py ${i}
sbatch ${i}_get_stat.sh

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
python $WP/quantify_rRNA.py $gtf_file" >> quantify_rRNA.sh

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

chmod 777 analysis_report.sh
pwd
mkdir sbatch_scripts
mv *.sh sbatch_scripts

