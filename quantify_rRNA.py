import os
import sys
from string import *
import math
import string
import re
import commands
import operator

if len(sys.argv) < 2:
    print "USAGE: python " + sys.argv[0] + " <gff file>"
    sys.exit(0)

gffFile=sys.argv[1]

rRNAgeneList=commands.getoutput("grep 'rRNA' "+gffFile+" |awk '{print $10}'").replace('"','').replace(";","").split("\n")

DIR_tophat=commands.getoutput("ls -d tophat_out_*").split()
outList=[]
for DIR in DIR_tophat:
	try:
		countFile=commands.getoutput("ls "+DIR+"/*.counts")
		totNum=commands.getoutput("awk '{SUM+=$2} END {print SUM}' "+countFile)
		rRNAnum=0
		Lines=open(countFile).readlines()
		n=0
		for line in Lines:
			geneID=line.split()[0]
			if geneID in rRNAgeneList:
				n=n+1
				num=int(line.split()[1])
				rRNAnum=rRNAnum+num
			
		#print countFile.split('/')[-1].split('.')[0], rRNAnum, totNum
		percent=round((float(rRNAnum)/int(totNum))*100,2)
		outLine=countFile.split('/')[-1].split('.')[0]+'\t'+str(percent)+'%'+'\n'
		outList.append(outLine)
	except:
		print "could not handle " + DIR
		pass
	

outF=open("rRNA.quantification",'w')
outF.writelines(outList)
outF.close()
