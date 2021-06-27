#!/root/miniconda3/bin/python

__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser

parser = OptionParser("Usage: %prog -i merge.depth -o merge.depth.bed")
parser.add_option("-i","--input",dest = "input",
		help = "Input file; merge.depth")
parser.add_option("-o","--output",dest = "out",
		help = "Output file; bedfile.")

(options,args) = parser.parse_args()

out=open(options.out,'w')
n = 1
with open (options.input,'r') as data:
	for i in data:
		i = i.strip()
		name = i[:4]
		if i.startswith(">"):
			title = name + "_" + str(n)
			out.write(title+"\n")
			n += 1
		else:
			out.write(i + "\n")
out.close()
