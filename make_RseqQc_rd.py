import sys
from bcbio.pipeline.config_loader import load_config

"""Stand in analysisdir eg merged"""

if len(sys.argv) < 6:
        print """
Usage:

make_RseqQc_rd.py  <sample name> <bed_file> <mail> <config_file> <path>

        <sample name>           This name: /tophat_out_<sample name>
        <bed_file>      
        <mail>                  eg: maya.brandi@scilifelab.se
        <config_file>           post_process.yaml assumes that you have specified samtools 
                                version under 'custom_algorithms'/'RNA-seq analysis'
	<path>			Path to analysis dir containing the tophat_out_ directories
        """
        sys.exit()

name 		= sys.argv[1]
bed_file	= sys.argv[2]
mail            = sys.argv[3]
config_file     = sys.argv[4]
path		= sys.argv[5]
try:
        config  = load_config(config_file)
        tools   = config['custom_algorithms']['RNA-seq analysis']
        sam	= tools['sam']+'/'+tools['sam_version']
except:
        print 'ERROR: problem loading samtools version from config file'


f=open("RSeQC_"+name+"_rd.sh",'w')

print >>f, """#!/bin/bash -l
#SBATCH -A a2010003
#SBATCH -p node
#SBATCH -t 50:00:00
#SBATCH --qos=seqver
#SBATCH -e RSeQC_rd_{0}.err
#SBATCH -o RSeQC_rd_{0}.out
#SBATCH -J RSeQC_rd_{0}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={1}

export PYTHONPATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/pipelines/RSeQC-2.0.0/sw/comp/python/2.7_kalkyl/bin:$PATH

module unload samtools
module load {3}

read_distribution.py -i {4}/tophat_out_{0}/accepted_hits_sorted_dupRemoved_{0}.bam -r {2}""".format(name, mail, bed_file, sam, path)
