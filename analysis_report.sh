#!/bin/bash -l
cd /bubo/home/h24/mayabr/glob/RNA_analysis
python /analysis_report.py /bubo/home/h9/mikaelh/Rattus_norvegicus.RGSC3.4.68.gtf /bubo/home/h24/mayabr/config/post_process.yaml /bubo/home/h9/mikaelh/Rat_Baylor_2004_ensembl.bed -c /proj/a2012043/private/nobackup/projects/l_olsson_12_01/intermediate/ -r -s -d -f
