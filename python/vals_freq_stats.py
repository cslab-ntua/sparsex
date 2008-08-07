
from vals_freq import cache_size, hfreq_bdb_init, tlist_dbfile, entropy_dbfile

if __name__ == '__main__':
	from sys import argv
	for fname in argv[1:]:
		db = hfreq_bdb_init(fname, cache_size)
		tlist_dbfile(fname, db=db, topsize=64)
		print "Entropy: %.2f" % entropy_dbfile(fname, db)
		db.close()
