
from popen2 import popen2
from struct import pack, unpack
from itertools import imap
from math import log
import string
import sys
import os
import fcntl
import time

from mmf import MMF
from toplist import Toplist
from hfreq_bdb import hfreq_bdb_init, hfreq_bdb_add, bdb_iteritems, hfreq_bdb_get
from progress import progress_iter

def set_nonblock(fd):
	flags = fcntl.fcntl(fd, fcntl.F_GETFL)
	flags = flags | os.O_NONBLOCK
	fcntl.fcntl(fd, fcntl.F_SETFL, flags)

def lfreq_print(lfreq, total, desc):
	pt = 0
	cnt = 0
	r = float(total)/float(len(lfreq))
	print "%s (%ld/%ld => %f)" %( desc, len(lfreq), total, r)
	for i in lfreq.__reversed__():
		p = 100*(float(i[0])/float(total))
		pt += p
		print "%-10s %10d %10.1f %10.1f"  % (hex(long(i[1])), i[0], p, pt)
	print

def hfreq2list(hfreq):
	rev = lambda t: (t[1],t[0])
	hlist = [ rev(hfreq.popitem()) for i in xrange(len(hfreq)) ]
	hlist.sort()
	return hlist

def hfreq_bdb2list(hfreq_bdb):
	hlist = []
	for v, freq in bdb_iteritems(hfreq_bdb):
		v = unpack("L", v)[0]
		hlist.append((freq,v))
	hlist.sort()
	return hlist

def lfreq_split(lfreq, fmt_new, fmt_old):
	ret = {}
	for i in xrange(len(lfreq)):
		freq,v = lfreq.pop()
		for nv in unpack(fmt_new, pack(fmt_old, v)):
			ret[nv] = ret.get(nv, 0) + freq
	return ret

def hfreq_split_bdb(hfreq_bdb, fmt_new, fmt_old, new_hfreq_bdb):
	iterator = bdb_iteritems(hfreq_bdb)
	if show_progress:
		iterator = progress_iter(iterator, len(hfreq_bdb), prefix="SPLIT ")
	for v,freq in iterator:
		v = unpack("L", v)[0]
		for nv in unpack(fmt_new, pack(fmt_old, v)):
			hfreq_bdb_add(new_hfreq_bdb, pack("L", nv), freq)

def lfreq_split_bdb(lfreq, fmt_new, fmt_old, db):
	for i in xrange(len(lfreq)):
		freq,v = lfreq.pop()
		for nv in unpack(fmt_new, pack(fmt_old, v)):
			hfreq_bdb_add(db, pack("L", nv), freq)

def hfreq_mmf(mmf, po, pi):
	hfreq = {}
	cnt = 0
	for i,j,v in mmf.iter_flstr():
		pi.write(v + "\n")
		pi.flush()
		u64 = string.atol(po.readline(), 16)
		hfreq[u64] = hfreq.get(u64, 0) + 1
		cnt += 1
	return hfreq

def hfreq_mmf_xor(mmf, po, pi):
	hfreq = {}
	cnt = 0
	u64_old = long(0)
	for i,j,v in mmf.iter_flstr():
		pi.write(v + "\n")
		pi.flush()
		u64 = string.atol(po.readline(), 16)
		hfreq[u64] = hfreq.get(u64 ^ u64_old, 0) + 1
		u64_old = u64
		cnt += 1
	return hfreq

def bdb_file(mmf_file, file, new=True):
	store_dir='/home/kkourt/work/research/spmv/datafiles/stats.bdb/'
	f = store_dir + os.path.basename(mmf_file) + '/vals_freq-' + file
	if new and os.path.exists(f):
		raise RuntimeError
	return f

def hfreq_mmf_bdb(mmf, po, pi, db):
	if len(db) != 0:
		return
	iterator = mmf.iter_flstr()
	if show_progress:
		iterator = progress_iter(iterator, mmf.entries, prefix="MMF ")
	for i,j,v in iterator:
		pi.write(v + "\n")
		pi.flush()
		ul = string.atol(po.readline(), 16)
		hfreq_bdb_add(db, pack("L", ul))
	db.sync()

def hfreq_mmf_xor_bdb(mmf, po, pi, db):
	if len(db) != 0:
		return
	iterator = mmf.iter_flstr()
	if show_progress:
		iterator = progress_iter(iterator, mmf.entries, prefix="MMF ")
	ul_old = 0
	for i,j,v in iterator:
		pi.write(v + "\n")
		pi.flush()
		ul = string.atol(po.readline(), 16)
		hfreq_bdb_add(db, pack("L", ul ^ ul_old))
		ul_old = ul
	db.sync()
	
def hfreq_bdb_tlist(db, size=10):
	tl = Toplist(size)
	count = 0
	iterator = bdb_iteritems(db)
	if show_progress:
		iterator =  progress_iter(iterator, len(db), prefix="TLIST ")
	for v, freq in iterator:
		v = unpack("L", v)[0]
		tl.add((freq, v))
		count += long(freq)
	return tl, count

def tlist_repr(tl):
	return [(long(i[1]), i[0]) if i is not None else (long(0),0) for i in tl]

def tlist_dbfile(dbfile, db=None, topsize=10):
	if db is None:
		db = hfreq_bdb_init(dbfile, cache_size/2)
	tl, count = hfreq_bdb_tlist(db, size=topsize)
	count = float(count)
	print "db: %s" % dbfile
	print "%24s %7s %7s %14s " % ("value", "freq", "perc", "total_perc")
	pt = 0.0
	for v, v_freq in tlist_repr(tl):
		if v_freq == 0:
			break
		p = (v_freq/count)*100.0
		pt += p
		print "%#024x %7d    %4.1f       %4.1f" % (v, v_freq, p, pt)

def hfreq_bdb_elems(db):
	""" return the total number of elements in the set """
	iterator = bdb_iteritems(db)
	if show_progress:
		iterator = progress_iter(iterator, len(db), prefix="COUNT ")
	return reduce(lambda x,y: x+y, imap(lambda x: x[1], iterator))

def _entropy_term(freq_iter, nr_elems, base=2):
	return imap(
		lambda p : p*log(p,base),
		imap(lambda x : float(x)/float(nr_elems), freq_iter)
	)

def entropy_dbfile(dbfile, db=None, nr_elems=None):
	""" information entropy for a frequency list stored in bdb
		nr_elems: total number of elements (default: number of values in db)
	"""
	if db is None:
		db = hfreq_bdb_init(dbfile, cache_size)
	if nr_elems is None:
		nr_elems = hfreq_bdb_elems(db)
	iterator = imap(lambda x: x[1], bdb_iteritems(db))
	if show_progress:
		iterator = progress_iter(iterator, len(db), prefix="ENTROPY ")
	return -reduce(lambda x,y: x+y, _entropy_term(iterator, nr_elems))

def vals_freq_bdb(mmf_file, po, pi, mmf_bdb=hfreq_mmf_bdb, fprefix=''):
		mmf = MMF(mmf_file)
		nnz = mmf.entries
		print '%s -------' % mmf_file
		sys.stderr.write(mmf_file + "\n")
		
		f = bdb_file(mmf.file, fprefix + '64.db')
		db = hfreq_bdb_init(f, cache_size / 2)
		mmf_bdb(mmf, po, pi, db)
		#tl = hfreq_bdb_tlist(db)
		#print "%2s %s" % ("64", tlist_repr(tl))
		#sys.stdout.flush()

		args = (("32", "II", "L"), ("16", "HH", "I"), ("8",  "BB", "H"))
		for bits, new_fmt, old_fmt in  args:
			f = bdb_file(mmf.file, fprefix + bits + '.db')
			db_old = db
			db = hfreq_bdb_init(f, cache_size / 2)
			print "--SPLIT", bits
			hfreq_split_bdb(db_old, new_fmt, old_fmt, db)
			db_old.close()
			#print "--TLIST", bits
			#tl = hfreq_bdb_tlist(db)
			#print "%2s %s" % (bits, tlist_repr(tl))
			#sys.stdout.flush()
		db.close()

def vals_freq(mmf_file, po, pi):
		mmf = MMF(mmf_file)
		nnz = mmf.entries
		hfreq = hfreq_mmf(mmf, po, pi)
		lfreq = hfreq2list(hfreq)
		lfreq_print(lfreq, nnz, '64-bits')
		
		args = (("32", "II", "L"), ("16", "HH", "I"), ("8",  "BB", "H"))
		for bits, new_fmt, old_fmt in args:
			hfreq = lfreq_split(lfreq, new_fmt, old_fmt)
			lfrq = hfreq2list(hfreq)
			lfreq_print(lfreq, nnz(64/int(bits)), bits + '-bits')

show_progress = True
cache_perc = 0.85 # percentage of total memory allocated to the cache
cache_size = int((os.sysconf('SC_PHYS_PAGES')*os.sysconf('SC_PAGESIZE'))*cache_perc)
if __name__ == '__main__':
	from sys import argv
	po, pi = popen2("d2ul")
	#set_nonblock(po)
	for f in argv[1:]:
		pass
		#vals_freq_bdb(f, po, pi)
		#vals_freq_bdb(f, po, pi, hfreq_mmf_xor_bdb, fprefix="xor-")
		#vals_freq(f, po, pi)
