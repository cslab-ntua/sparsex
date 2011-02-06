import os
import random
import math

val_max = -100
val_min = 100

def randval():
	return random.uniform(val_max, val_min)

def dumsp_rowlen(file, nz_nr, cols_nr=None, rowlen=1, rand_cols=False):
	rows_nr = nz_nr // rowlen
	if not cols_nr:
		cols_nr = rows_nr
	else:
		assert cols_nr >= row_len

	assert not os.path.lexists(file)
	f = open(file, 'w')
	s = "%d %d %d\n" % (cols_nr, rows_nr, nz_nr)

	for x in xrange(rows_nr):
		if rand_cols:
			yrange = list(random.sample(xrange(cols_nr), rowlen))
			yrange.sort()
		else:
			yrange = xrange(rowlen)
		for y in yrange:
			s += "%d %d %f\n" % ((x+1), (y+1), randval())

		if len(s) > 131072:
			f.write(s)
			s = ''

	if s:
		f.write(s)

	f.close()

def dumsp_bubble(file, cols_nr, rows_nr, bubble=4):
	assert not os.path.lexists(file)
	f = open(file, 'w')
	s = "%d %d %d\n" % (cols_nr, rows_nr, cols_nr*rows_nr/2)

	for x in xrange(0, rows_nr, bubble):
		for c in xrange(0,cols_nr):
			s += "%d %d %f\n" % ((x+1), (c+1), randval())
		f.write(s)
		s=''
	if s:
		f.write(s)
	f.close()

def dumsp_worse(file, cols_nr, rows_nr):
	assert not os.path.lexists(file)
	f = open(file, 'w')
	s = "%d %d %d\n" % (cols_nr, rows_nr, cols_nr*rows_nr/2)

	for x in xrange(rows_nr):
		for c in xrange(0, cols_nr,2):
			s += "%d %d %f\n" % ((x+1), (c+1), randval())
		f.write(s)
		s = ''

	if s:
		f.write(s)
	f.close()

def dumsp_emptyrow(file, cols_nr, rows_nr):
	assert not os.path.lexists(file)
	f = open(file, 'w')
	s = "%d %d %d\n" % (cols_nr, rows_nr, (cols_nr/2)*(rows_nr/2))

	for x in xrange(0, rows_nr, 2):
		for c in xrange(0, cols_nr,2):
			s += "%d %d %f\n" % ((x+1), (c+1), randval())
		f.write(s)
		s = ''

	if s:
		f.write(s)
	f.close()

def create1(size):
	files = []
	for rowlen in xrange(1,12):
		f = 'dumsp-%d.%d' % (size,rowlen)
		dumsp_rowlen(f, size, rowlen=rowlen)
		print f
		f += '.rand'
		dumsp_rowlen(f, size, rowlen=rowlen, rand_cols=True)
		print f

def create_cr(size):
	a = int(math.sqrt(size))
	rowlens = list(xrange(1,16))
	rowlens.extend((int(a*.05), int(a*.1), int(a*0.5), a))
	for rowlen in rowlens:
		f = 'const_row-%d.%d' % (size, rowlen)
		dumsp_rowlen(f, size, rowlen=rowlen)
		print f
		f += '.rand'
		dumsp_rowlen(f, size, rowlen=rowlen, rand_cols=True)
		print f

def create_dense(size):
	f = 'dense-%d' % size
	dumsp_rowlen(f, size*size, rowlen=size)

if __name__ == '__main__':
	print "dumsp_worse(file, cols_nr, rows_nr)"
	print "create_dense(size)"
