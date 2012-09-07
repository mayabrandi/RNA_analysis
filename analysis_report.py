#!/usr/bin/env python
"""
Usage:
     	analysis_reports.py <run name> <project id> <samle names> [Options]

Arguments:
     	<run name>
                - eg: 20120307B_hiseq2000

	<project id> 
		- eq: J.Hasmats_11_03 

	<samle names>
		- Sample names should be give as a comma delimited string.

Options:
  	-c, --config-file = <config file>   	   
		- post_process.yaml
        -r, --rRNA_table        to include an rRNA table in the report
				-requires rRNA.quantification from quantify_rRNA.sh
        -s, --Map_Stat          to include mapping statistics table in the report
				-requires stat from get_stat.sh
        -d, --Read_Dist         to include read distribution table in the report
				-requires stderr output files from read_distribution.py
        -f, --FPKM              to include fpkm-heatmap and -pca plot in the report
				-requires FPKM_PCAplot.pdf and FPKM_heatmap.pdf from correl.R	
        -c FILE, --config-file=FILE
                                FILE should be a config file. (post_process.yaml)


"""

import os
import sys
from optparse import OptionParser
from operator import itemgetter
from numpy import *
import yaml
import glob
import re
from mako.template import Template
from mako.lookup import TemplateLookup


from bcbio.log import logger, setup_logging
import bcbio.templates.mako2rst as m2r
from texttable import Texttable

from bcbio.google import bc_metrics 
from bcbio.pipeline.config_loader import load_config
from bcbio.scilifelab.google.project_metadata import ProjectMetaData

import operator


def main(run,project_id,sample_names,config_file,Map_Stat,Read_Dist,FPKM,rRNA_table):

    TEMPLATE="""\
RNA-seq analysis report for ${project_id}
=================================
   
${latex_opt}

Summary
-------------------------
**Project name:**
${project_id} (UPPMAX project ${uppnex})

**Samples:**
${samplenames}
    
**Run name:**
${runname}

**Mapping:** 
${mapping}
    
**Duplicate removal:**
${dup_rem}
    
**Read count:**
${read_count}
    
**RPKM/FPKM values:**
${quantifyer}
    
**Result directories on UPPMAX:** /proj/${uppnex}/INBOX/${project_id}/analysis/alignments (BAM files), /proj/${uppnex}/INBOX/${project_id}/analysis/quantification (FPKM files)

.. raw:: latex
       
   \clearpage    
    
Results
-------------------------"""
    
    if Map_Stat:
        TEMPLATE=TEMPLATE+"""
Mapping statistics
^^^^^^^^^^^^^^^
${Mapping_statistics}
    
Comments
~~~~~~~~
    
**tot # read pairs:** 

The total number of read pairs indicates the total number of sequenced paired-end reads. Since a paired-end read is made up of two sequenced fragments (mates), the total number of sequenced 100-bp regions is twice the number shown in this column.
    
**% mapped reads:**
    
The number of fragments that are mapped relative to the total number of sequenced fragments. 
    
**% reads left after dup rem:**
    
We remove duplicate reads (paired end reads where both mates map to the same loci as both mates in a different paired-end read because these are likely to be artifacts caused by PCR amplification or over-sequencing. Aligned files in BAM format with duplicates removed can be found in /proj/${uppnex}/INBOX/${project_id}/analysis/alignments.


.. raw:: latex
       
   \clearpage

"""
    
    TEMPLATE=TEMPLATE+"""
Expression values
^^^^^^^^^^^^^^^^^
    
The /proj/${uppnex}/INBOX/${project_id}/analysis/quantification folder contains FPKM values calculated using the Cufflinks program using ENSEMBL annotation of genes and transcripts for each sample. These files also contain the upper and lower limits of the confidence interval for the FPKM estimate. FPKM values are the paired-end equivalent of RPKM (Reads Per Kilobase per Million mapped reads; the standard measure for gene expression in RNA-seq.)
    
There is also a single fpkm_table.txt file, which contains all of the FPKM values. This can be opened in Excel or a regular text processing application.
    
For analyzing differential expression of genes or transcripts, it may be useful to have the raw read counts (the number of sequences that map to each gene/transcript) as well. These are calculated using the HTSeq software and are collected into a table called count_table.txt.


.. raw:: latex
       
   \clearpage

"""
    
    
    if FPKM:
        TEMPLATE=TEMPLATE+"""
FPKM heatmap
^^^^^^^^^^^^^^^^^
This heatmap shows the (Pearson) correlation between FPKM values of samples. 
    
${FPKM_heatmap}
    
    
.. raw:: latex
       
   \clearpage
    
FPKM PCA
^^^^^^^^^^^^^^^^^
This PCA (principal component analysis) score plot has the samples plotted according to their scores for the two principal components that explain the largest amount of variance in the FPKM data table. The number after 'expl var' in the axis labels tells you how much of the variance that is explained by each component. Similar samples should, in theory, cluster together in the PCA plot. PCA is a way to compress the information in your high-dimensional data matrix so that it can be plotted in two dimensions.
    
${FPKM_PCAplot}


.. raw:: latex
       
   \clearpage

"""
    
    if Read_Dist:
        TEMPLATE=TEMPLATE+"""
Read distribution
^^^^^^^^^^^^^^^^^
This table contain information about the extent to which sequences from each sample mapped to different structural parts of genes, like coding exons, untranslated regions, and transcription start sites. The actual number itself is less important than the relative values for the different kinds of regions. For a normal RNA-seq experiment you should have a higher value in the CDS Exon column than in the others, for example. "CDS Exon" means "coding sequence exons", "UTR" stands for "untranslated region", "TES" stands for "transcription end site", "TSS" stands for "transcription start site". "Intronic regions" should be interpreted as "intronic or intergenic regions".
Perhaps the most easily interpretable column is the final column, mRNA fraction, which gives the fraction [0-1] of sequences that mapped to ENSEMBL-annotated mRNA (including coding regions and UTRs). While this fraction is not completely accurate (because ENSEMBL doe not completely describe the transcriptome), it is a useful summary statistic which should be relatively high for an mRNA-seq experiment, typically above 0.8.

${Read_Distribution}


.. raw:: latex
       
   \clearpage

"""
    
    if rRNA_table:
        TEMPLATE=TEMPLATE+"""
Quantification of rRNA present in the samples
^^^^^^^^^^^^^^^^^^^^^
    
${rRNA_table}


.. raw:: latex
       
   \clearpage

"""

    sphinx_defs = []

    if config_file:
        config = load_config(config_file)
    else:
        config = {}

    sphinx_defs.append("('%s', '%s_analysis.tex', 'RNA-seq Analysis Report', u'SciLifeLab Stockholm', 'howto'),\n"  % (project_id, project_id))
    projectfile = "%s.mako" % (project_id) 
    fp = open(projectfile, "w")
    fp.write(TEMPLATE)
    fp.close()
    mylookup = TemplateLookup(directories=['./'])
    tmpl = Template(filename=projectfile, lookup=mylookup)

    proj_conf = {
        'id' : project_id,
	'run':run,
        'config' : config,
	'samples': sample_names.split(',')
         }

    d = generate_report(proj_conf)
    rstfile = "%s.rst" % (project_id)
    fp = open(rstfile, "w")
    fp.write(tmpl.render(**d))
    fp.close()

    sphinxconf = os.path.join(os.getcwd(), "conf.py")
    if not os.path.exists(sphinxconf):
        logger.warn("no sphinx configuration file conf.py found: you have to edit conf.py yourself!")
    else:
        fp = open(sphinxconf)
        lines = fp.readlines()
        fp.close()
        sdout = []
        modify_conf = False
        for sd in sphinx_defs:
            if not sd in lines:
                sdout.append(sd)
                modify_conf = True
        if modify_conf:
            i = lines.index("latex_documents = [\n")
            newconf = lines[:i+3] + sdout + lines[i+3:]
            fp = open("conf.py", "w")
            fp.write("".join(newconf))
            fp.close()

##-----------------------------------------------------------------------------
def generate_report(proj_conf):

    d = {
	'runname':proj_conf['run'],
	'project_id': proj_conf['id'],
        'samplenames': ' '.join(proj_conf['samples']),
        'latex_opt' : "",
        'uppnex': "",
        'mapping':"",
        'dup_rem':"",
        'read_count':"",
        'quantifyer':"",
        'gene_body_cov':"",
        'FPKM_heatmap':"",
        'FPKM_PCAplot':"",
        'Mapping_statistics': "",
        'Read_Distribution':"",
	'rRNA_table':""
        }

    ## Latex option (no of floats per page)
    floats_per_page = '.. raw:: latex\n\n   \setcounter{totalnumber}{8}'
    d['latex_opt'] = floats_per_page


    ## Metadata fetched from the 'Genomics project list' on Google Docs 
    try:
        proj_data = ProjectMetaData(proj_conf['id'], proj_conf['config'])
        uppnex_proj = proj_data.uppnex_id
    except:
        uppnex_proj = "b201YYXX"
        print "No uppnex ID fetched"
	pass
    if not uppnex_proj:
	uppnex_proj="b201YYXX"
        print "No uppnex ID fetched"
    d['uppnex'] = uppnex_proj 


    ## RNA-seq tools fetched from config file post_process.yaml
    try:
        tools      	= proj_conf['config']['custom_algorithms']['RNA-seq analysis']
        d['mapping']	= os.path.join(tools['aligner'],tools['aligner_version'])
        d['dup_rem']    = os.path.join(tools['dup_remover'],tools['dup_remover_version'])
        d['read_count'] = os.path.join(tools['counts'],tools['counts_version'])
        d['quantifyer'] = os.path.join(tools['quantifyer'],tools['quantifyer_version'])
    except:
	print "Could not fetched RNA-seq tools from config file post_process.yaml"
        d['mapping'] = "X"
        d['dup_rem'] = "X"
        d['read_count'] = "X"
        d['quantifyer'] = "X"
        pass


    ## Mapping Statistics
    tab = Texttable()
    tab.set_cols_dtype(['t','t','t','t'])
    tab.add_row(['Sample','tot_#_read_pairs','%_mapped_reads','%_reads_left_after_dup_rem'])
    try:
	for sample_name in proj_conf['samples']:
	    f=open('tophat_out_'+sample_name+'/stat_'+sample_name, 'r')
	    data = f.readlines()
	    tab.add_row([sample_name,data[1].split()[1],data[2].split()[1],data[3].split()[1]])
	    f.close()
	d['Mapping_statistics']=tab.draw()
    except:
	try:
            f=open('stat', 'r')
            data = f.readlines()
            D=dict(zip(data[0].split(),zip(data[1].split(),data[2].split(),data[3].split())))
            for sample_name in proj_conf['samples']:
	        if D.has_key(sample_name):
                    tab.add_row([sample_name,D[sample_name][0],D[sample_name][1],D[sample_name][2]])
	        else:
	            print 'kould not find '+sample_name+' in stat'
            d['Mapping_statistics']=tab.draw() 
            f.close()
        except:
	    print "Could not make Mapping Statistics table"
            pass


    ## Read Distribution 
    try:
        tab = Texttable()
        json=open('Ever_rd.json','a')
        print >> json, '{'
        Groups=["Sample:","CDS Exons:","5'UTR Exons:","3'UTR Exons:","Intronic region:","TSS up 1kb:","TES down 1kb:"]

        tab.set_cols_dtype(['t','t','t','t','t','t','t','t'])
        tab.add_row(["Sample","CDS Exon","5'UTR Exon","3'UTR Exon","Intron","TSS up 1kb","TES down 1kb","mRNA frac"])
   	
	for i in range(len(proj_conf['samples'])):
	    sample_name=proj_conf['samples'][i] 
            print >> json, sample_name+': {'
            row=[sample_name]
            Reads_counts=[]
	    try:
		f=open('RSeQC_rd_'+sample_name+'.err','r')
	    except:
            	f=open('Ever_rd_'+sample_name+'.err','r')
		pass
            for line in f:
                Group=line.split('\t')[0]
                if Group in Groups:
                    if Group=="TES down 1kb:":
                        print >> json, '"'+Group+'"'+':'+str(line.split('\t')[3].strip())
                    else:
                        print >> json, '"'+Group+'"'+':'+str(line.split('\t')[3].strip())+','
                    row.append(str(line.split('\t')[3].strip())+' ')
                    Reads_counts.append(float(line.split('\t')[2].strip()))
	    if os.path.exists('RSeQC_rd_'+sample_name+'.err'):
		t=os.popen("grep 'Total Fragments' 'RSeQC_rd_"+sample_name+".err'|sed 's/Total Fragments               //g'")
	    else:
		try:
			t=os.popen("grep 'Total Fragments' 'Ever_rd_"+sample_name+".err'|sed 's/Total Fragments               //g'")
		except:		
			pass
            tot=float(t.readline())
            frac=(Reads_counts[0]+Reads_counts[1]+Reads_counts[2])/tot
            row.append(str(round((Reads_counts[0]+Reads_counts[1]+Reads_counts[2])/tot,2)))
            tab.add_row(row)
            f.close()
            if i==(len(proj_conf['samples'])-1):
                    print >> json,'}'
            else:
                    print >> json,'},'
        print >> json, '}'
        json.close()
        d['Read_Distribution']=tab.draw()

    except:
	print "Could not make Read Distribution table"
        pass


    ## FPKM_PCAplot, FPKM_heatmap
    if os.path.exists("FPKM_PCAplot.pdf") and os.path.exists("FPKM_heatmap.pdf"):
        d['FPKM_PCAplot'] = m2r.image("FPKM_PCAplot.pdf", width="100%")
        d['FPKM_heatmap'] = m2r.image("FPKM_heatmap.pdf", width="100%")
    else:
	print "could not make FPKM PCAplot and FPKM heatmap"


    ## rRNA_table
    try:
        tab = Texttable()
        tab.set_cols_dtype(['t','t'])
        tab.add_row(["Sample","rRNA"])
	f=open('rRNA.quantification','r')
	D={}
	for line in f:
            D[str(line.split('\t')[0].strip())]=str(line.split('\t')[1].strip())
        for sample_name in proj_conf['samples']:
            if D.has_key(sample_name):
                        tab.add_row([sample_name,D[sample_name]])
        d['rRNA_table']=tab.draw()
        f.close()
    except:
	print "could not generate rRNA table"
        pass   
 
    return d


##-----------------------------------------------------------------------------
if __name__ == "__main__":
    usage = """
    analysis_reports.py <run name> <project id> <samle names> [options]
			-r, --rRNA_table	to include an rRNA table in the repport
			-s, --Map_Stat          to include mapping statistics table in the repport
  			-d, --Read_Dist         to include read distribution table in the repport
  			-f, --FPKM          	to include fpkm-heatmap and -pca plot in the repport
  			-c FILE, --config-file=FILE
						FILE should be a config file. (post_process.yaml) 

    For more extensive help type analysis_report.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",default=False)
    parser.add_option("--v1.5", dest="v1_5_fc", action="store_true", default=False)
    parser.add_option("-c", "--config-file", dest="config_file", default=None)
    parser.add_option("-r", "--rRNA_table", dest="rRNA_table",action="store_true", default=False)
    parser.add_option("-s", "--Map_Stat", dest="Map_Stat",action="store_true", default=False)
    parser.add_option("-d", "--Read_Dist", dest="Read_Dist",action="store_true", default=False)
    parser.add_option("-f", "--FPKM", dest="FPKM",action="store_true", default=False)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(
		config_file = options.config_file,
		Map_Stat = options.Map_Stat,
		Read_Dist = options.Read_Dist,
		FPKM = options.FPKM,
		rRNA_table = options.rRNA_table
		)
		
    main(*args, **kwargs)

