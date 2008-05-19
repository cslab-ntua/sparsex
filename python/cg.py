from mmf import MMF

def spmv_cg(mmf_file):
	mmf = MMF(mmf_file)
	mmf_iter = mmf.iter()

	for y,x,val in mmf_iter:
		print "\ty[%s] += x[%s] * (double)%s;" % (y,x,val)

if __name__ == '__main__':
	from sys import argv
	spmv_cg(argv[1])
