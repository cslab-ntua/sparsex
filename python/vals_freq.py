
from popen2 import popen2
from mmf import MMF
import string
from struct import pack, unpack

def lfreq_print(lfreq, total, desc):
	pt = 0
	print "%s (%ld/%ld => %f)" %( desc, len(lfreq), total, float(total)/float(len(lfreq)))
	for i in lfreq.__reversed__():
		p = 100*(float(i[0])/float(total))
		pt += p
		print "%-10s %10d %10.1f %10.1f"  % (hex(i[1]), i[0], p, pt)
	print

def hfreq2list(hfreq):
	rev = lambda t: (t[1],t[0])
	hlist = list([ rev(hfreq.popitem()) for i in xrange(len(hfreq)) ])
	hlist.sort()
	return hlist

def lfreq_split(lfreq, fmt_new, fmt_old):
	ret = {}
	for i in xrange(len(lfreq)):
		freq,v = lfreq.pop()
		for nv in unpack(fmt_new, pack(fmt_old, v)):
			ret[nv] = ret.get(nv, 0) + freq
	return ret
		

progress = True
if __name__ == '__main__':
	from sys import argv
	po, pi = popen2("d2ul")
	for f in argv[1:]:
		mmf = MMF(f)
		nnz = mmf.entries

		print '----- %s' %f
		hfreq = {}
		cnt = 0
		for i,j,v in mmf.iter_flstr():
			pi.write(v + "\n")
			pi.flush()
			u64 = string.atol(po.readline(), 16)
			hfreq[u64] = hfreq.get(u64, 0) + 1
			cnt += 1
			if (progress and cnt % (1024*1024) == 0):
				print cnt, nnz, float(cnt)/float(nnz)

		lfreq = hfreq2list(hfreq)
		lfreq_print(lfreq, nnz, '64-bits')

		hfreq = lfreq_split(lfreq, "II", "L")
		lfreq = hfreq2list(hfreq)
		lfreq_print(lfreq, nnz*2, '32-bits')

		hfreq = lfreq_split(lfreq, "HH", "I")
		lfreq = hfreq2list(hfreq)
		lfreq_print(lfreq, nnz*4, '16-bits')

		hfreq = lfreq_split(lfreq, "BB", "H")
		lfreq = hfreq2list(hfreq)
		lfreq_print(lfreq, nnz*8, '8-bits')
