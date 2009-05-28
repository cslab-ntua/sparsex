from pprint import pprint

from mmf import MMF
from encoding import enc_delta, enc_rle_iter

def rle_yx_stats(xy_iter, limit=4):
	xs = list()
	rle_freqs = dict()
	rle_vals = dict()
	_e = dict()
	_e['total'] = 0
	r_prev = 1

	def update_stats():
		if len(xs) == 0:
			return
		for delta_val, freq in enc_rle_iter(enc_delta(xs)):
			_e['total'] += freq
			if freq >= limit:
				if delta_val in rle_freqs:
					rle_freqs[delta_val] += 1
					rle_vals[delta_val] += freq
				else:
					rle_freqs[delta_val] = 1
					rle_vals[delta_val] = freq

	for r,c in xy_iter:
		if r != r_prev:
			update_stats()
			r_prev = r
			del xs[:]
		xs.append(c)
	update_stats()

	total = _e['total']
	for delta in list(rle_vals.iterkeys()):
		rle_vals[delta] = 100*(float(rle_vals[delta])/float(total))
	return rle_freqs, rle_vals

def mmf_rle_stats(file):
	spm = MMF(file)
	rows = spm.rows
	# make indices 1-based for better delta encoding, so that the first
	# value is encoded as 1 and not as 0.
	data = list((1+int(y),1+int(x)) for y,x,v in spm.iter())
	stats = {}
	# horizontal
	stats['horizontal'] = rle_yx_stats(data)
	# vertical
	vdata = sorted(list((y,x) for (x,y) in data))
	stats['vertical'] = rle_yx_stats(vdata)
	del vdata
	# diagonal
	diag_f = lambda t : (t[1] - t[0] + rows, t[0] if t[0]<=t[1] else t[1])
	ddata = map(diag_f, data)
	ddata.sort()
	stats['diagonal'] = rle_yx_stats(ddata)
	del ddata

	return stats

if __name__ == '__main__':
	from sys import argv
	for file in argv[1:]:
		print file
		stats = mmf_rle_stats(file)
		for i in sorted(stats.iterkeys()):
			print "- %-15s\n" % i,
			print "\trle number      :",
			print ' '.join([ "%s:%s" % x for x in (
					sorted(stats[i][0].iteritems(), key=lambda x: x[1], reverse=True)
			) ])
			print "\trle percentage  :",
			print ' '.join([ "%s:%s" % x for x in (
					sorted(stats[i][1].iteritems(), key=lambda x: x[1], reverse=True)
			) ])
		print
