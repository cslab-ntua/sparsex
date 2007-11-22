import os

class MMF(object):
	def __init__(self, file):
		if not os.path.isfile(file):
			if os.path.isfile(file + '.mtx'):
				file = file + '.mtx'
			elif os.path.isfile(file + '.mtx.gz'):
				file = file + '.mtz.gz'
		self.file = file
		if self.file[-3:] == '.gz':
			import gzip
			self.source = gzip.open(self.file)
		else:
			self.source = open(self.file, 'r')
		
		for l in self.source:
			if not l.startswith('%'):
				break	

		l = l.split()
		assert len(l) == 3	
		self.rows,self.cols,self.entries = map(int,l)

	def iter(self):
		for l in self.source:
			if not l or l.startswith('%'):
				continue
			l = l.split()
			yield (int(l[0]) -1, int(l[1]) -1, float(l[2]))

	def iter_flstr(self):
		for l in self.source:
			if not l or l.startswith('%'):
				continue
			l = l.split()
			yield (int(l[0]) -1, int(l[1]) -1, l[2])
