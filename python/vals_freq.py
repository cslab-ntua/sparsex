
from popen2 import popen2
from mmf import MMF
import string
from struct import pack, unpack

def hash_freq_print(hfreq, total, desc):
	lfreq = list(hfreq.iteritems())
	lfreq.sort(cmp = lambda x,y: cmp(y[1], x[1]))
	pt = 0
	print "%s (%ld/%ld => %f)" %( desc, len(lfreq), total, float(total)/float(len(lfreq)))
	for i in lfreq:
		p = 100*(float(i[1])/float(total))
		pt += p
		print "%-10s %10d %10.1f %10.1f"  % (hex(i[0]), i[1], p, pt)
	print
	
def hash_freq_split(hfreq, fmt_new, fmt_old):
	ret = {}
	for v, freq in hfreq.iteritems():
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
		hash_freq_print(hfreq, nnz, '64-bits')

		hfreq = hash_freq_split(hfreq, "II", "L")
		hash_freq_print(hfreq, nnz*2, '32-bits')

		hfreq = hash_freq_split(hfreq, "HH", "I")
		hash_freq_print(hfreq, nnz*4, '16-bits')

		hfreq = hash_freq_split(hfreq, "BB", "H")
		hash_freq_print(hfreq, nnz*8, '8-bits')
