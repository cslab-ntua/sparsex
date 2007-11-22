import pyspm

type_str = {
	pyspm.TYPE_DOUBLE: "double",
	pyspm.TYPE_FLOAT: "float",
}

import sys
import os
from sys import argv

if len(argv) < 4:
	print 'usage:', argv[0], 'file type index_bits (32/64)'
	sys.exit(1)

file = argv[1]
if argv[2] == 'float':
	type = pyspm.TYPE_FLOAT 
else:
	type = pyspm.TYPE_DOUBLE

if argv[3] == '32':
	crs=pyspm.CRS32
else:
	crs=pyspm.CRS64

flops = pyspm.crs_bench(file, type, crs)
print "%s pyspm_crs%s.%s %s" % (os.path.basename(file), argv[3], type_str[type], flops)
