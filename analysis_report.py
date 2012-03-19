#!/usr/bin/env python
"""Make RNA-seq analysis for a project

Usage:
     analysis_reports.py <flowcell id> <project id> <samle names> 
                	    [--config-file=<config file>]

Sample names should be give as a comma delimited string.



The script loads the run_info.yaml and generates a report for each project.

Options:
  -c, --config-file = <config file>   	   typ post_process.yaml Vet inte vad jag ska skriva
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

from bcbio.log import create_log_handler
from bcbio.pipeline import log
import bcbio.templates.mako2rst as m2r
from texttable import Texttable

from bcbio.google import bc_metrics 
from bcbio.pipeline.config_loader import load_config
from bcbio.scilifelab.google.project_metadata import ProjectMetaData

import operator

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

**Result directories on UPPMAX:** /proj/${uppnex}/INBOX/${runname}/analysis/alignments (BAM files), /proj/${uppnex}/INBOX/${runname}/analysis/quantification (FPKM files)


Results
-------------------------

Mapping statistics
^^^^^^^^^^^^^^^
${Mapping_statistics}

Comments
~~~~~~~~~

**tot # read pairs:** 

The total number of read pairs indicates the total number of sequenced paired-end reads. Since a paired-end read is made up of two sequenced fragments (mates), the total number of sequenced 100-bp regions is twice the number shown in this column.

**% mapped reads:**

The number of fragments that are mapped relative to the total number of sequenced fragments. 

**% reads left after dup rem:**

We remove duplicate reads (paired end reads where both mates map to the same loci as both mates in a different paired-end read) because these are likely to be artifacts caused by PCR amplification or over-sequencing. Aligned files in BAM format with duplicates removed can be found in /proj/${uppnex}/INBOX/${runname}/analysis/alignments.



.. raw:: latex
   
   \clearpage

Expression values
^^^^^^^^^^^^^^^^^

The /proj/${uppnex}/INBOX/${runname}/analysis/quantification folder contains FPKM values calculated using the Cufflinks program using 
ENSEMBL annotation of genes and transcripts for each sample. These files also contain the upper and lower limits of the confidence interval for the FPKM estimate. FPKM values are the paired-end equivalent of RPKM (Reads Per Kilobase per Million mapped reads; the standard measure for gene expression in RNA-seq.)

There is also a single fpkm_table.txt file, which contains all of the FPKM values. This can be opened in Excel or a regular text processing application.

For analyzing differential expression of genes or transcripts, it may be useful to have the raw read counts (the number of sequences that map to each gene/transcript) as well. These are calculated using the HTSeq software and are collected into a table called count_table.txt.

.. raw:: latex
   
   \clearpage

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

Read distribution
^^^^^^^^^^^^^^^^^
This table contain information about the extent to which sequences from each sample mapped to different structural parts of genes, like coding exons, untranslated regions, and transcription start sites. The actual number itself is less important than the relative values for the different kinds of regions. For a normal RNA-seq experiment you should have a higher value in the CDS Exon column than in the others, for example. "CDS Exon" means "coding sequence exons", "UTR" stands for "untranslated region", "TES" stands for "transcription end site", "TSS" stands for "transcription start site". "Intronic regions" should be interpreted as "intronic or intergenic regions".


${Read_Distribution}

"""

def main(fcID,project_id,sample_names,config_file):
	#'m_ohman_11_01' ska vara k har, men stammer inte med M.Ohman_11_01 som i run_info.yaml
    run=os.listdir('/bubo/proj/a2010002/projects/'+'_'.join(project_id.lower().split('.'))+'/data')[0]
 
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
	'fcID':fcID,
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
        log.warn("no sphinx configuration file conf.py found: you have to edit conf.py yourself!")
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
        'Read_Distribution':""
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
        print "No project metadata fetched"
        pass
    d['uppnex'] = uppnex_proj 


    ## RNA-seq tools fetched from config file post_process.yaml
    try:
        tools      	= proj_conf['config']['custom_algorithms']['RNA-seq analysis']
        d['mapping']	= os.path.join(tools['aligner'],tools['aligner_version'])
        d['dup_rem']    = os.path.join(tools['dup_remover'],tools['dup_remover_version'])
        d['read_count'] = os.path.join(tools['counts'],tools['counts_version'])
        d['quantifyer'] = os.path.join(tools['quantifyer'],tools['quantifyer_version'])
    except:
        d['mapping'] = "X"
        d['dup_rem'] = "X"
        d['read_count'] = "X"
        d['quantifyer'] = "X"
        pass


    ## Mapping Statistics
    f=open('stat', 'r')
    data = f.readlines()
    max_cols=4

    for i in range((len(data[1:])-1)/max_cols+1):
       	tab = Texttable()
	tab.set_cols_dtype(['t','t','t','t'])
	data_sep = data[i*max_cols+1:max_cols*(i+1)+1]
	data_sep.insert(0,data[0])
        for j in range(len(data_sep[0].split(' '))):
	    lista=[]
     	    for col in data_sep:
                lista.append(str(col.strip().split(' ')[j].replace('_',' '))+' ')
            tab.add_row(lista)

        d['Mapping_statistics']=d['Mapping_statistics']+'\n\n'+tab.draw()

    ## Read Distribution 
    tab = Texttable()
    json=open('Ever_rd.json','a')
    print >> json, '{'
    Groups=["Samle","CDS Exons:","5'UTR Exons:","3'UTR Exons:","Intronic region:","TSS up 1kb:","TES down 1kb:"]
    tab.set_cols_dtype(['t','t','t','t','t','t','t'])
    tab.add_row(Groups)
    for i in range(len(proj_conf['samples'])):
	sample_name=proj_conf['samples'][i]
        print >> json, sample_name+': {'
	row=[sample_name]
        f=open('Ever_rd_'+sample_name+'.err','r')
        for line in f:
	    Group=line.split('\t')[0]
	    if Group in Groups:
		if Group=="TES down 1kb:":
                    print >> json, '"'+Group+'"'+':'+str(line.split('\t')[3].strip())
                else:
                    print >> json, '"'+Group+'"'+':'+str(line.split('\t')[3].strip())+','
		row.append(str(line.split('\t')[3].strip())+' ')
        tab.add_row(row)
	f.close()
	if i==(len(proj_conf['samples'])-1):
        	print >> json,'}'
	else:
		print >> json,'},'
    print >> json, '}'
    json.close()

    d['Read_Distribution']=tab.draw()

    ## FPKM_PCAplot, FPKM_heatmap
    d['FPKM_PCAplot'] = m2r.image("FPKM_PCAplot.pdf", width="100%")
    d['FPKM_heatmap'] = m2r.image("FPKM_heatmap.pdf", width="100%")

    return d


##-----------------------------------------------------------------------------
if __name__ == "__main__":
    usage = """
    analysis_reports.py <flowcell id>
                           [--archive_dir=<archive directory>]

    For more extensive help type analysis_report.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",default=False)
    parser.add_option("--v1.5", dest="v1_5_fc", action="store_true", default=False)
    parser.add_option("-c", "--config-file", dest="config_file", default=None)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(config_file = options.config_file)
    main(*args, **kwargs)
