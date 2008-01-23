from delta import delta_jmp_parse, delta_get_size

if __name__ == '__main__':
	from sys import argv
	import os
	seq_limit = os.getenv("SPM_DELTA_SEQ_LIMIT")
	seq_limit = int(seq_limit) if seq_limit is not None else 1000000
	for f in argv[1:]:
		d = delta_jmp_parse(f, seq_limit)
		size = delta_get_size(d)
		d.mt_check_jmp(f)
		print f, 'OK'

