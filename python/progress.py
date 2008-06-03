from time import time
import sys

def cb_pstderr(cnt, itotal, prefix):
	r = (float(cnt)/float(itotal))*100
	sys.stderr.write("%s%d/%d (%5.1f %%)\n" % (prefix,cnt,itotal,r))

def progress(itotal, iprint_sh=10, cb=cb_pstderr, prefix=''):
	mask = (1<<iprint_sh) - 1
	def d(i):
		def f(*args, **kwargs):
			cnt = 0
			for item in i(*args, **kwargs):
				if cnt & mask == 0:
					cb(cnt, itotal, prefix)
				yield item
				cnt += 1
		return f
	return d

def progress_iter(iterator, itotal, iprint_sh=13, cb=cb_pstderr, prefix=''):
	cnt = 0
	mask = (1<<iprint_sh) - 1
	for item in iterator:
		if cnt & mask == 0:
			cb(cnt, itotal, prefix)
		yield item
		cnt += 1

if __name__ == '__main__':
	l = 1024*1024

	def i1():
		for i in xrange(l):
			yield i
	
	#for i in i0():
	#	pass
	
	t0 = time()
	for i in i1():
		pass
	t1 = time()
	print t1-t0

	t0 = time()
	for i in progress(l)(i1)():	
		pass
	t1 = time()
	print t1-t0

	t0 = time()
	iterator = i1()
	for i in progress_iter(iterator, l):
		pass
	t1 = time()
	print t1-t0
