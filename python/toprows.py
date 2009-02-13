from mmf import MMF
from toplist import Toplist

if __name__ == '__main__':
	from sys import argv

	print 'Rows with larger number of elements'
	for f in argv[1:]:
		tl = Toplist(4, cmp_fn = lambda x,y: cmp(x[1],y[1]), default = (None, None))
		row_prev = 0
		row_elems = 0
		for row,col,val in MMF(f).iter():
			if row == row_prev:
				row_elems += 1
			else:
				tl.add((row_prev, row_elems))
				row_prev = row
				row_elems = 1
		print f, tl
