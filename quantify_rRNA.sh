#! /bin/bash -l
#SBATCH -A a2010003
#SBATCH -p node
#SBATCH -t 15:00:00
#SBATCH -J quantify_rRNA
#SBATCH -e quantify_rRNA.err
#SBATCH -o quantify_rRNA.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL
cd /bubo/home/h24/mayabr/glob/RNA_analysis
python /quantify_rRNA.py maya.brandi@scilifelab.se
