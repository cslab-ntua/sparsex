from mmf import MMF

def detect_holes(file):
	mmf = MMF(file)
	pr = 0
	holes = []
	for r, c, v in mmf.iter():
		if r > pr + 1:
			holes.append(r)
		pr = r
	return holes

if __name__ == '__main__':
	from sys import argv, stdout
	import os
	for f in argv[1:]:
		print f, detect_holes(f)
		stdout.flush()

