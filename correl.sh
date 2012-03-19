#! /bin/bash -l
#SBATCH -A a2010003
#SBATCH -p node
#SBATCH -t 20:00:00
#SBATCH -J correl
#SBATCH -e correl.err
#SBATCH -o correl.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL

args=("$@")
samplenames=(`echo ${args[1]} | tr "," " "`)
path=${args[0]}

echo "--args $path $samplenames"

R CMD BATCH "--args $path ${samplenames[*]}" /bubo/home/h24/mayabr/glob/useful_scripts/RNA_analys/correl.R



