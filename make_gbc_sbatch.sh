#!/bin/bash -l

path=$1
samplenames=(`echo $2 | tr "," "\n"`)
bedfile=$3

for i in ${samplenames[*]};do 
	echo "#!/bin/bash -l">>$i"_runEver_gbc.sh"
	echo "">>$i"_runEver_gbc.sh"
	echo "#SBATCH -A a2010003">>$i"_runEver_gbc.sh"
	echo "#SBATCH -p node">>$i"_runEver_gbc.sh"
	echo "#SBATCH -t 90:00:00">>$i"_runEver_gbc.sh"
	echo "#SBATCH -e Ever_gbc_$i.err">>$i"_runEver_gbc.sh" 
	echo "#SBATCH -o Ever_gbc_$i.out">>$i"_runEver_gbc.sh"
	echo "#SBATCH -J Ever_gbc_$i">>$i"_runEver_gbc.sh"
	echo "#SBATCH --mail-type=ALL">>$i"_runEver_gbc.sh"
	echo "#SBATCH --mail-user=maya.brandi@scilifelab.se">>$i"_runEver_gbc.sh"
	echo "">>$i"_runEver_gbc.sh"
	
	echo "module unload samtools">>$i"_runEver_gbc.sh"
	echo "module load samtools/0.1.9">>$i"_runEver_gbc.sh"
	echo "">>$i"_runEver_gbc.sh"

        echo "bedfile=$bedfile">>$i"_runEver_gbc.sh"
	echo 'samtools view $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam | python ~/downloads/prog/EVER-seq-1.0.5/scripts/geneBody_coverage.py -i - -r $bedfile -o '$i>>$i"_runEver_gbc.sh"
	echo "">>$i"_runEver_gbc.sh"

        echo "R CMD BATCH "$i".geneBodyCoverage_plot.r">>$i"_runEver_gbc.sh"
        echo "mv geneBody_coverage.pdf "$i"_geneBody_coverage.pdf">>$i"_runEver_gbc.sh"
done

