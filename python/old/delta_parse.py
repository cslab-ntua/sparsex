
from mmf import MMF
import pyspm

#from range_merge import range_merge

svn_rev = "$Rev: 371 $"

stats = {}
states = (ST_UNKNOWN, ST_MAYBE_DENSE, ST_DENSE, ST_SPARSE) = xrange(4)
events = (EV_NEW_ROW, EV_SEQ, EV_SPARSE, EV_MAXSIZE) = xrange(4)

states_str = {
	ST_UNKNOWN:     "ST_UNKNOWN",
	ST_MAYBE_DENSE: "ST_MAYBE_DENSE",
	ST_DENSE:       "ST_DENSE",
	ST_SPARSE:      "ST_SPARSE",
}

events_str = {
	EV_NEW_ROW:  "EV_NEW_ROW",
	EV_SEQ :     "EV_SEQ",
	EV_SPARSE :  "EV_SPARSE",
	EV_MAXSIZE : "EV_MAXSIZE",
}

type_str = {
	pyspm.TYPE_DOUBLE: "double",
	pyspm.TYPE_FLOAT: "float",
}

class SpmDelta_state(object):
	def __init__(self, input_iterator, delta, seq_limit=3):
		self.state = ST_UNKNOWN
		self.inp_i = input_iterator
		self.delta = delta

		inp = self.inp_i.next()
		self.x_buf = [inp[1]]
		self.x_old = 0
		self.v_buf = [inp[2]]
		self.y_cur = inp[0]
		self.y_old = 0

		self.seq_limit = seq_limit
		self.seq_cnt = 0
		self.flag_nr = False

def generate_event(inp, state):
	y,x,v = inp
	if y != state.y_cur: # XXX: This must be first
		ev = EV_NEW_ROW
	elif len(state.x_buf) >= pyspm.DELTA_CTL_SIZE_MAX:
		ev = EV_MAXSIZE
	elif x == state.x_buf[-1] + 1:
		ev = EV_SEQ
	else:
		ev = EV_SPARSE

	#print "got_event", events_str[ev], "state", states_str[state.state], 'input', inp
	return ev

def handle_new_row(event, inp, state):
	#print "handle_new_row start"
	assert len(state.x_buf) != 0
	if state.state == ST_DENSE:
		finalize_unit_de(state)
	else:
		finalize_unit_sp(state)
	state.flag_nr = True
	state.y_old = state.y_cur
	state.y_cur = inp[0]
	state.x_old = 0
	state.state = ST_UNKNOWN
	#print "handle_new_row end"

def handle_seq(event, inp, state):
	if state.state == ST_UNKNOWN:
		state.state = ST_MAYBE_DENSE
		state.seq_cnt = 1
	elif state.state == ST_MAYBE_DENSE:
		state.seq_cnt += 1
		if state.seq_cnt >= state.seq_limit:
			sp = len(state.x_buf) - state.seq_cnt
			if sp > 0:
				finalize_unit_sp(state, sp)
			state.state = ST_DENSE
	elif state.state == ST_SPARSE:
		state.seq_cnt = 1
		state.state = ST_MAYBE_DENSE
	
def handle_sparse(event, inp, state):
	if state.state != ST_SPARSE:
		if state.state == ST_DENSE:
			finalize_unit_de(state)
		state.state = ST_SPARSE

def do_finalize(state):
	if state.state == ST_DENSE:
		finalize_unit_de(state)
	else:
		finalize_unit_sp(state)
	

def handle_maxsize(event, inp, state):
	assert state.state != ST_UNKNOWN
	do_finalize(state)

event_handlers = {
	EV_NEW_ROW  : handle_new_row,
	EV_SEQ      : handle_seq,
	EV_SPARSE   : handle_sparse,
	EV_MAXSIZE  : handle_maxsize,
}

def parse(mmf_file, type=pyspm.TYPE_DOUBLE, seq_limit=3):
	mmf = MMF(mmf_file)
	delta = pyspm.delta(type=type)
	delta.prepare()
	iter = mmf.iter()
	state = SpmDelta_state(input_iterator=iter,delta=delta,seq_limit=seq_limit)

	for inp in iter:
		event = generate_event(inp, state)
		event_handlers[event](event, inp, state)
		state.x_buf.append(inp[1])
		state.v_buf.append(inp[2])
	
	do_finalize(state)
	
	delta.finalize(mmf.entries, mmf.rows, mmf.cols)
	del state
	return delta

def parse_nodelta(mmf_file, type=pyspm.TYPE_DOUBLE, seq_limit=3):
	mmf = MMF(mmf_file)
	iter = mmf.iter()
	state = SpmDelta_state(input_iterator=iter,delta=None,seq_limit=seq_limit)

	for inp in iter:
		event = generate_event(inp, state)
		event_handlers[event](event, inp, state)
		state.x_buf.append(inp[1])
		state.v_buf.append(inp[2])
	
	do_finalize(state)
	
	del state
	return None

def finalize_unit_sp_simple(state, nr=0, finlen=None):
	#print "finalize_unit_sp"
	if finlen is None:
		finlen = len(state.x_buf)
	elif finlen == 0:
		return
	
	xd_l = list()
	xd_max = 0
	xo = state.x_old
	for x in state.x_buf[:finlen]:
		xd = x - xo
		assert xd >= 0
		xd_l.append(xd)
		xd_max = max(xd_max, xd)
		xo = x
	state.x_old = xo
	del state.x_buf[:finlen]

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	v_l = state.v_buf[:finlen]
	del state.v_buf[:finlen]

	#print "\tsize=", finlen, "ci_delta_max=", xd_max, "nr_delta=", yd
	state.delta.add_sp_unit(size=finlen, ci_delta_max=xd_max, nr_delta=yd)
	#print "\tci_delta_max=", xd_max, "ci_delta_l=", xd_l
	state.delta.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
	state.delta.addvals(values_l=v_l)

##	 ranges = range_merge(xd_l)
##	 for start,length,cisize in ranges:
##			ci_delta_sl = xd_l[start:start+length]
##			ci_max = max(ci_delta_sl)
##			v_sl = v_l[start:start+length]
##			print "\tsize=", length, "ci_delta_max=", ci_max, "nr_delta=", yd
##			state.delta.add_sp_unit(size=length, ci_delta_max=ci_max, nr_delta=yd)
##			print "\tci_delta_max=", ci_max, "ci_delta_l=", ci_delta_sl
##			state.delta.add_sp_cols(ci_delta_max=ci_max, ci_delta_l=ci_delta_sl)
##			state.delta.addvals(values_l=v_sl)
##			yd = 0


def finalize_unit_sp_jmp(state, nr=0, finlen=None):
	if finlen is None:
		finlen = len(state.x_buf)
	elif finlen == 0:
		return
	
	xd_l = list()
	xd_max = 0
	xo = state.x_buf[0]
	jmp = xo - state.x_old
	for x in state.x_buf[1:finlen]:
		xd = x - xo
		assert xd >= 0
		xd_l.append(xd)
		xd_max = max(xd_max, xd)
		xo = x
	state.x_old = xo
	del state.x_buf[:finlen]

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	v_l = state.v_buf[:finlen]
	del state.v_buf[:finlen]

	#print "\t[SP-JMP] size",finlen,"jmp",jmp,"ci_delta_max",xd_max,"nr_delta",yd,"ci_delta_l",xd_l
	state.delta.add_sp_unit_jmp(size=finlen, ci_delta=jmp, ci_delta_max=xd_max, nr_delta=yd)
	state.delta.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
	state.delta.addvals(values_l=v_l)

def finalize_unit_sp_stats(state, nr=0, finlen=None):
	if finlen is None:
		finlen = len(state.x_buf)
	elif finlen == 0:
		return
	
	xd_l = list()
	xd_max = 0
	xo = state.x_old
	for x in state.x_buf[:finlen]:
		xd = x - xo
		assert xd >= 0
		xd_l.append(xd)
		xd_max = max(xd_max, xd)
		xo = x
	state.x_old = xo
	del state.x_buf[:finlen]

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	del state.v_buf[:finlen]

	if len(xd_l[1:]) == 0:
		return
	pattern = (pyspm.delta_cisize(max(xd_l[1:])), tuple(xd_l[1:]))
	if pattern in stats['sp_unit_patterns']:
		stats['sp_unit_patterns'][pattern] += 1
	else:
		stats['sp_unit_patterns'][pattern] = 1

def finalize_unit_de_simple(state, nr=0):
	assert len(state.x_buf) == len(state.v_buf)

	size = len(state.x_buf)
	xd = state.x_buf[0] - state.x_old
	state.x_old = state.x_buf[0] + size

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	#print "\t[DE] size=", size, "ci_delta=", xd, "nr_delta=", yd
	state.delta.add_de_unit(size=size, ci_delta=xd, nr_delta=yd)
	state.delta.addvals(values_l=state.v_buf)
	del state.v_buf[:]
	del state.x_buf[:]

def finalize_unit_de_stats(state, nr=0):
	assert len(state.x_buf) == len(state.v_buf)

	size = len(state.x_buf)
	xd = state.x_buf[0] - state.x_old
	state.x_old = state.x_buf[0] + size

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	del state.v_buf[:]
	del state.x_buf[:]


def check_simple(file):
	print "checking:", file
	global finalize_unit_sp
	finalize_unit_sp = finalize_unit_sp_simple
	delta = parse(file)
	delta.check(file, "spm_delta_multiply")

def check_jmp(file):
	print "checking:", file
	global finalize_unit_sp
	finalize_unit_sp = finalize_unit_sp_jmp
	delta = parse(file)
	delta.check_jmp(file, "spm_delta_multiply")

def bench_simple(file, loops=128, seq_limit=3, type=pyspm.TYPE_DOUBLE):
	global finalize_unit_sp
	finalize_unit_sp = finalize_unit_sp_simple
	delta = parse(file, seq_limit=seq_limit)
	ret = delta.bench(loops)
	del delta
	return ret

def bench_jmp(file, loops=128, seq_limit=3, type=pyspm.TYPE_DOUBLE):
	global finalize_unit_sp, finalize_unit_de
	finalize_unit_sp = finalize_unit_sp_jmp
	finalize_unit_de = finalize_unit_de_simple
	delta = parse(file, seq_limit=seq_limit)
	ret = delta.bench_jmp(loops)
	del delta
	return ret


def stat_patterns(file, seq_limit=4):
	global finalize_unit_sp, finalize_unit_de, stats
	finalize_unit_sp = finalize_unit_sp_stats
	finalize_unit_de = finalize_unit_de_stats
	stats['sp_unit_patterns'] = {}
	parse_nodelta(file, seq_limit=seq_limit)
	saved = 0
	for k, v in stats['sp_unit_patterns'].iteritems():
		init_size = (1<<k[0])*len(k[1])*v
		curr_size = (1<<k[0])*len(k[1]) + v*8
		saved += ( init_size - curr_size )
	return saved

#f = "/s/matrices/sorted/006.crystk02.mtx.sorted"
#f = "/s/matrices/sorted/006.crystk02.mtx.sorted"
if __name__ == '__main__':

	import os
	#pyspm.setaffinity(3);
	files_dir = "/s/matrices/sorted/"
	files = map( lambda x: files_dir+x, filter(
		lambda f: f.endswith(".mtx.sorted"), 
		os.listdir(files_dir)
	) )
	files.sort()

	def all(fid, type=pyspm.TYPE_DOUBLE):
		f = os.path.basename(files[fid])
		t = type_str[type]
		flops = pyspm.crs_bench(files[fid], type, pyspm.CRS64)
		print "%s pyspm_crs64.%s %s" % (f, t, flops)
		flops = pyspm.crs_bench(files[fid], type, pyspm.CRS32)
		print "%s pyspm_cr32.%s %s" % (f, t, flops)
		for sl in xrange(2,13):
			flops = bench(fid, seq_limit=sl, type=type)
			print "%s pyspm_delta.%s.%s %s" % (f, t, sl, flops)
		flops = bench(fid, seq_limit=100000, type=type)
		print "%s pyspm_delta.%s.NO_DE %s" % (f, t, flops)
	
	def do_all():
		#print "# svn_rev: %s" % svn_rev
		for type in (pyspm.TYPE_DOUBLE, pyspm.TYPE_FLOAT):
			for fid in xrange(len(files)):
				all(fid, type)
		
	def main():
		import sys
		from sys import argv
		if len(argv) < 4:
			print 'usage:', argv[0], 'file type seq_limit'
			return 1

		file = argv[1]

		if argv[2] == 'float':
			type = pyspm.TYPE_FLOAT 
		else:
			type = pyspm.TYPE_DOUBLE

		seq_l = int(argv[3])
		flops = bench_jmp(file, seq_limit=seq_l, type=type)
		print "%s pyspm_delta.%s.%s %s" % (os.path.basename(file), type_str[type], seq_l, flops)
	
	main()

	#print "check(id)"
	#print "bench(id)"
	#print "all(id)"
	#print "do_all()"
	#do_all()
