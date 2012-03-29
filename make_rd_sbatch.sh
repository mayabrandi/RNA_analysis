#!/bin/bash -l

path=$1
samplenames=(`echo $2 | tr "," "\n"`)
bedfile=$3

for i in ${samplenames[*]};do
        echo "#!/bin/bash -l">>$i"_runEver_rd.sh"
        echo "">>$i"_runEver_rd.sh"
        echo "#SBATCH -A a2010002">>$i"_runEver_rd.sh"
        echo "#SBATCH -p node">>$i"_runEver_rd.sh"
        echo "#SBATCH -t 15:00:00">>$i"_runEver_rd.sh"
        echo "#SBATCH -e Ever_rd_$i.err">>$i"_runEver_rd.sh"
        echo "#SBATCH -o Ever_rd_$i.out">>$i"_runEver_rd.sh"
        echo "#SBATCH -J Ever_rd_$i">>$i"_runEver_rd.sh"
        echo "#SBATCH --mail-type=ALL">>$i"_runEver_rd.sh"
        echo "#SBATCH --mail-user=maya.brandi@scilifelab.se">>$i"_runEver_rd.sh"
        echo "">>$i"_runEver_rd.sh"

        echo "module unload samtools">>$i"_runEver_rd.sh"
        echo "module load samtools/0.1.9">>$i"_runEver_rd.sh"
        echo "">>$i"_runEver_rd.sh"

        echo "bedfile=$bedfile">>$i"_runEver_rd.sh"
	echo "samtools view $path/tophat_out_$i/accepted_hits_sorted_dupRemoved_$i.bam | python ~/downloads/prog/EVER-seq-1.0.5/scripts/read_distribution.py -i - -r "'$bedfile'>>$i"_runEver_rd.sh"
        echo "">>$i"_runEver_rd.sh"
done

