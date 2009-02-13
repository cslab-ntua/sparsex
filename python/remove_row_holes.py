from mmf import MMF

def remove_holes(file):
	mmf = MMF(file)
	print mmf.rows,  mmf.cols, mmf.entries

	pr = 0
	rholes = 0
	for r, c, v in mmf.iter():
		if r > pr + 1:
			rholes += (r - pr - 1)
		pr = r
		print r-rholes+1, c+1, v

if __name__ == '__main__':
	from sys import argv, stdout
	import os
	for f in argv[1:]:
		remove_holes(f)
		stdout.flush()

