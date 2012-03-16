#!/usr/bin/env python

"""
snrc.py

This script fetches samplenames and readcounts for a user specified project and flowcell run. 
For example the M.Muurinen_11_01a project run on flowcell 120127_SN1018_0062_BD0H2HACXX.

The script assumes there is a file called run_info.yaml located at (In this pearticular example)
proj/a2010002/archive/20127_SN1018_0062_BD0H2HACXX/run_info.yaml

The script prints its output to a textfile samplenames.txt in the current directory
"""

import sys
import yaml
import os
import re
from bcbio.solexa.flowcell import get_flowcell_info


if len(sys.argv) < 3:
    print "USAGE: python " + sys.argv[0] + " <flowcell ID> <project name>"
    print "Example: python " + sys.argv[0] + " 120127_SN1018_0062_BD0H2HACXX M.Muurinen_11_01a"
    sys.exit(0)


fcid=sys.argv[1]
proj=sys.argv[2]

# Fetching sample names and read counts from run_info.yaml
analysis_path = '/proj/a2010002/nobackup/illumina'
archive_path = '/proj/a2010002/archive'
fc_name, fc_date = get_flowcell_info(fcid)
fc = yaml.load(open(archive_path + '/' + fcid + '/run_info.yaml'))
i=1
dict={}
read_counts={}
for l in fc:
        if l.has_key('multiplex') & (l['description'].split(', ')[1]==proj):
		sname={}
            	samples = l['multiplex']
		for s in samples:
			sname[str(s['barcode_id'])] = re.sub('_index.+','',s['name'])
			read_counts[re.sub('_index.+','',s['name'])]=0
		dict[i]=sname
	i+=1
for l in range(1,9):
	bc_file = open(os.path.join(analysis_path,fcid,'_'.join([str(l), fc_date, fc_name, "barcode"]),'_'.join([str(l), fc_date, fc_name, "bc.metrics"])))
	for line in bc_file:
		c = line.strip().split()
		samp_name=re.sub('_index.+','',c[0])
		if (samp_name !='unmatched') and dict.has_key(l) and dict[l].has_key(samp_name):
			read_counts[dict[l][samp_name]]+=int(c[1])

# Sorting sample names if nummeric      
try:
       	vals=[]
       	nummeric_keys=[]
       	keys=[]
       	for key in read_counts.keys():
       		nummeric_keys.append(int(key))
       		nummeric_keys=sorted(nummeric_keys)
       	for key in nummeric_keys:
       		vals.append(read_counts[str(key)])
       		keys.append(str(key))
except:
       	keys=read_counts.keys()
       	vals=read_counts.values()
       	pass

f=open('samplenames.txt','a')


f.write(' '.join(keys))
f.write('\n')
for i in range(len(vals)): vals[i]=str(vals[i])
f.write(' '.join(vals))
f.close()

