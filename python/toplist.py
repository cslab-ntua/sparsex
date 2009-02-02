
# some sanity checks
assert -long(2**64) > None
assert +long(2**64) > None
assert 0 > None

class Toplist(object):
	def __init__(self, size, cmp_fn=cmp, default=None):
		self.l = [ default for i in xrange(size) ]
		self.cmp_fn = cmp_fn
	
	def add(self, element):
		self.l.append(element)
		self.l.sort(cmp=self.cmp_fn, reverse=True)
		self.l.pop()

	def __len__(self):
		return self.l.__len__() 
	
	def __iter__(self):
		return self.l.__iter__()

	def __getitem__(self, x):
		return self.l[x]
	
	def __getslice__(self, i, j):
		return self.l[i:j]
	
	def __repr__(self):
		return self.l.__repr__()

if __name__ == '__main__':
	import random
	size = 10000
	top = 10
	tl = Toplist(top)
	l = range(size)
	random.shuffle(l)
	for i in l:
		tl.add(i)
	assert tl.l == [size - i -1 for i in xrange(top)] 
