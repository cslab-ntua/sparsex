import pyspm

test_delta1 = [357, 1, 8, 1, 1, 1, 8, 1, 1, 1, 409, 1, 3, 1, 4, 1, 1, 1]
test_delta2 = [357, 1, 357, 1, 357, 1, 357, 1, 357, 1, 357, 1, 357, 1]
test_delta3 = [1, 357, 1, 357, 1, 357, 1, 357, 1, 357, 1, 357, 1, 357]

def do_merge_on(ranges, i):
	r = ranges[i]
	ranges = ranges[:i] + ranges[i+1:]
	ranges[i-1][1] += r[1]
	ranges[i-1][2] = max(ranges[i-1][2], r[2])
	return ranges

import time
def range_merge(ci_delta):
	p = 4
	ci_size = [ None for i in xrange(len(ci_delta)) ]
	
	start = 0
	length = 1
	val = pyspm.delta_cisize(ci_delta[0])

	ranges = []
	for i in xrange(1,len(ci_delta)):
		v = pyspm.delta_cisize(ci_delta[i])
		if v != val:
			ranges.append([start, length, val])
			start += length
			length = 1
			val = v
		else:
			length += 1
	ranges.append([start, length, val])

	while True:
		#print ranges
		#time.sleep(1)
		for i in xrange(1, len(ranges)):
			#print "\t", i-1, ranges[i-1], i, ranges[i]
			# same value
			if ranges[i-1][2] == ranges[i][2]:
				ranges = do_merge_on(ranges, i)
				#print "MERGING COMMON VALS"
				break
			# small unit
			small = i-1 if ranges[i-1][2] < ranges[i][2] else i
			#print "\t small=", small
			if ranges[small][1]<p:
				ranges = do_merge_on(ranges, i)
				#print "MERGING SMALL UNITS"
				break
		else:
			break
	
	return ranges
