from mmf import MMF

def vidx_iter(mmf_file):
	mmf = MMF(mmf_file)
	values = {}

	prev_idx = 0
	for x,y,v in mmf.iter():
		if v in values:
			yield values[v]
		else:
			idx = len(values)
			values[v] = idx
			yield prev_idx - idx
			prev_idx = idx

if __name__ == '__main__':
	from sys import argv

	for f in argv[1:]:
		for vidx in vidx_iter(f):
			print vidx

