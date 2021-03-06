#!/root/miniconda/bin/python
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
with open (options.input,'r') as data:
	for i in data:
		if i.startswith(">"):
			i = i.lstrip(">").split(" ")
			out.write(i[0]  + "\n")
out.close()







