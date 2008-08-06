
from vals_freq import tlist_dbfile

if __name__ == '__main__':
	from sys import argv
	for f in argv[1:]:
		tlist_dbfile(f, topsize=64)
