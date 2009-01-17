
import vals_freq
from vals_freq import cache_size, hfreq_bdb_init, tlist_dbfile, entropy_dbfile, hfreq_bdb_elems

vals_freq.show_progress = False

if __name__ == '__main__':
	from sys import argv
	print "%-50s %-10s %-10s %s" % ('dbfile', 'vals', 'uvals', 'entropy')
	for fname in argv[1:]:
		db = hfreq_bdb_init(fname, cache_size)
		#tlist_dbfile(fname, db=db, topsize=64)
		f = '/'.join(fname.split('/')[-2:])
		nr_elems = hfreq_bdb_elems(db)
		entropy = entropy_dbfile(fname, db, nr_elems=nr_elems)
		print "%-50s %-10d %-10d %4.2f" % (f, nr_elems, len(db), entropy)
		db.close()
