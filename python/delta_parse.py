
from mmf import MMF
import pyspm

from time import time
from sys import stdout

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
	def __init__(self, input_iter, seq_limit=3):
		self.state = ST_UNKNOWN
		self.inp_i = input_iter

		inp = self.inp_i.next()
		self.x_buf = [inp[1]]
		self.x_old = 0
		self.v_buf = [inp[2]]
		self.y_cur = inp[0]
		self.y_old = 0

		self.seq_limit = seq_limit
		self.seq_cnt = 0
		self.flag_nr = False
	
class SpmDelta_parser(object):
	def __init__(self, ll, **kwargs):
		 self.ll = ll
		 self.progress = kwargs.get('progess', False)

		 self._rows = 0
	
	def report_progerss(self):
		self._rows += 1
		if self._rows % 1024 == 0:
			current_time = time()
			ratio = self._rows/(current_time-self._t0)
			remaining =  self.mmf.rows - self._rows
			print self.state.mmf_file, '[', (current_time-self._t0)/(60*60), 'h]', 
			print 'rows remaining:', remaining,
			print 'rows/sec:', ratio,
			print 'ETA:', ((remaining/ratio))/(60*60), 'h'
			stdout.flush()
		
	def handle_new_row(self, event, inp):
		if self.progress:
			self.report_progerss()
		state = self.state
		assert len(state.x_buf) != 0
		if state.state == ST_DENSE:
			self.ll.finalize_unit_de(state)
		else:
			self.ll.finalize_unit_sp(state)
		state.flag_nr = True
		state.y_old = state.y_cur
		state.y_cur = inp[0]
		state.x_old = 0
		state.state = ST_UNKNOWN

	def handle_seq(self, event, inp):
		state = self.state
		if state.state == ST_UNKNOWN:
			state.state = ST_MAYBE_DENSE
			state.seq_cnt = 1
		elif state.state == ST_MAYBE_DENSE:
			state.seq_cnt += 1
			if state.seq_cnt >= state.seq_limit:
				sp = len(state.x_buf) - state.seq_cnt
				if sp > 0:
					self.ll.finalize_unit_sp(state, finlen=sp)
				state.state = ST_DENSE
		elif state.state == ST_SPARSE:
			state.seq_cnt = 1
			state.state = ST_MAYBE_DENSE

	def handle_sparse(self, event, inp):
		state = self.state
		if state.state != ST_SPARSE:
			if state.state == ST_DENSE:
				self.ll.finalize_unit_de(state)
			state.state = ST_SPARSE

	def do_finalize(self):
		state = self.state
		if state.state == ST_DENSE:
			self.ll.finalize_unit_de(state)
		else:
			self.ll.finalize_unit_sp(state)
	
	def handle_maxsize(self, event, inp):
		state = self.state
		assert state.state != ST_UNKNOWN
		self.do_finalize()
	
	def generate_event(self, input):
		y,x,v = input
		if y != self.state.y_cur: # XXX: This must be first
			ev = EV_NEW_ROW
		elif len(self.state.x_buf) >= pyspm.DELTA_CTL_SIZE_MAX:
			ev = EV_MAXSIZE
		elif x == self.state.x_buf[-1] + 1:
			ev = EV_SEQ
		else:
			ev = EV_SPARSE
		#print "got_event", events_str[ev], "state", states_str[self.state.state], 'input', input
		return ev

	def parse(self, mmf_file, type=pyspm.TYPE_DOUBLE, seq_limit=3,float_str=False):
		self.mmf = MMF(mmf_file)
		self._t0 = time()
		if float_str:
			mmf_iter = self.mmf.iter_flstr()
		else:
			mmf_iter = self.mmf.iter()
		self.state = SpmDelta_state(input_iter=mmf_iter,seq_limit=seq_limit)
		self.state.mmf_file = mmf_file
		self.state.nr_cols = self.mmf.cols
		self.state.nr_rows = self.mmf.rows
		self.state.nr_nz = self.mmf.entries

		for inp in mmf_iter:
			event = self.generate_event(inp)
			if event == EV_NEW_ROW:
				self.handle_new_row(event, inp)
			elif event == EV_SEQ:
				self.handle_seq(event, inp)
			elif event == EV_SPARSE:
				self.handle_sparse(event, inp)
			elif event == EV_MAXSIZE:
				self.handle_maxsize(event, inp)
			else:
				raise ValueError, 'Unknown event: %s' % ev

			self.state.x_buf.append(inp[1])
			self.state.v_buf.append(inp[2])
		self.do_finalize()

		ret = self.ll.finalize(self.state)
		
		del self.state
		return ret

class SpmDelta_ll(object):
	def __init__(self):
		pass
	
	def finalize_unit_sp(self, state, nr=0, finlen=None):
		raise NotImplementedError
	
	def finalize_unit_de(self, state, nr=0, finlen=None):
		raise NotImplementedError
	
	def finalize(self, state, **kwargs):
		raise NotImplementedError

class SpmDelta_ll_common(object):
	def __init__(self):
		self.units = []
		self.patterns = {}
		self.values = []

	def finalize_unit_de(self, state, nr=0, finlen=None):
		assert len(state.x_buf) == len(state.v_buf)
		size = len(state.x_buf)
		xd = state.x_buf[0] - state.x_old

		if state.flag_nr:
			yd = state.y_cur - state.y_old
			state.flag_nr = False
		else:
			yd = 0

		#self.units.append(('DE', size, xd))
		#self.values.extend(state.v_buf)
	
		del state.v_buf[:]
		del state.x_buf[:]

	def finalize_unit_sp(self, state, nr=0, finlen=None):
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
		ci_size = pyspm.delta_cisize(xd_max)

		#k = (pyspm.delta_cisize(max(xd_l)), tuple(xd_l))
		#if k in self.patterns:
		#	self.patterns[k] += 1
		#else:
		#	self.patterns[k] = 1
		#self.values.extend(v_l)

	def get_saved(self, state):
		saved = 0
		for k,v in self.patterns.iteritems():
			init_size = (1<<k[0])*len(k[1])*v
			curr_size = (1<<k[0])*len(k[1]) + v*4
			saved += ( init_size - curr_size )
		return saved
	
	def finalize(self, state, **kwargs):
		saved = self.get_saved(state)
		mmf = kwargs.get('mmf')
		#crs_size = float(mmf.entries*(8+4) + mmf.rows*4)
		#reduction = 100*((crs_size-saved)/crs_size)
		#print current_file, 'saved:', saved, 'crs_size', crs_size,'reduction:', reduction, '%'

import idx_list as idx_l
class SpmDelta_ll_compressed(object):
	def __init__(self):
		self.units = []
		self.sp_colinds = {
			pyspm.DELTA_CISIZE_U8  : [],
			pyspm.DELTA_CISIZE_U16 : [],
			pyspm.DELTA_CISIZE_U32 : [],
			pyspm.DELTA_CISIZE_U64 : [],
		}
		self.sp_cur_idx = {
			pyspm.DELTA_CISIZE_U8  : 0,
			pyspm.DELTA_CISIZE_U16 : 0,
			pyspm.DELTA_CISIZE_U32 : 0,
			pyspm.DELTA_CISIZE_U64 : 0,
		}
		self.sp_units = {
			pyspm.DELTA_CISIZE_U8  : 0,
			pyspm.DELTA_CISIZE_U16 : 0,
			pyspm.DELTA_CISIZE_U32 : 0,
			pyspm.DELTA_CISIZE_U64 : 0,
		}
		self.values = []
		idx_l.clear_rev_index()
	
	def _matrix_size(self):
		ws = len(self.values)*8
		ws += len(self.sp_colinds[pyspm.DELTA_CISIZE_U8])*1
		ws += len(self.sp_colinds[pyspm.DELTA_CISIZE_U16])*2
		ws += len(self.sp_colinds[pyspm.DELTA_CISIZE_U32])*4
		ws += len(self.sp_colinds[pyspm.DELTA_CISIZE_U64])*8
		ws += len(self.units)*(1+1+4) # flags, size, colid ptr
		return ws

	def _print_stats(self):
		print 'values: ', len(self.values)
		print 'sp_units: '
		print '  1 byte : ', self.sp_units[pyspm.DELTA_CISIZE_U8]
		print '  2 bytes: ', self.sp_units[pyspm.DELTA_CISIZE_U16]
		print '  4 bytes: ', self.sp_units[pyspm.DELTA_CISIZE_U32]
		print '  8 bytes: ', self.sp_units[pyspm.DELTA_CISIZE_U64]
		print 'sp_col_idxs: '
		print '  1 byte : ', len(self.sp_colinds[pyspm.DELTA_CISIZE_U8])
		print '  2 bytes: ', len(self.sp_colinds[pyspm.DELTA_CISIZE_U16])
		print '  4 bytes: ', len(self.sp_colinds[pyspm.DELTA_CISIZE_U32])
		print '  8 bytes: ', len(self.sp_colinds[pyspm.DELTA_CISIZE_U64])
		print current_file, 'MATRIX SIZE:', str(self._matrix_size())
	
	def finalize_unit_de(self, state, nr=0, finlen=None):
		assert len(state.x_buf) == len(state.v_buf)
		size = len(state.x_buf)
		xd = state.x_buf[0] - state.x_old

		if state.flag_nr:
			yd = state.y_cur - state.y_old
			state.flag_nr = False
		else:
			yd = 0

		self.units.append(('DE', size, xd))
		self.values.extend(state.v_buf)
	
		del state.v_buf[:]
		del state.x_buf[:]

	def finalize_unit_sp(self, state, nr=0, finlen=None):
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
		ci_size = pyspm.delta_cisize(xd_max)

		offset = 0
		while True:
			#limit = 8 if len(xd_l) - 8 > 8 else len(xd_l)
			if offset != len(xd_l):
				ci_idx, length = idx_l.do_add3(
					self.sp_colinds[ci_size], 
					xd_l,
					offset,
					limit = len(xd_l)
				)
			else:
				length = 0
				ci_idx = 0

			self.units.append(('SP', length, jmp, ci_size, ci_idx))
			self.sp_units[ci_size] += 1
			offset += length
			if offset == len(xd_l):
				break
			jmp = xd_l[offset]
			offset += 1
			ci_size = pyspm.delta_cisize(max(xd_l))
			
		self.values.extend(v_l)

	def finalize(self, state):
		self._print_stats()
	

current_file = ''
if __name__ == '__main__':
	try:
		import psyco
		psyco.full()
	except ImportError:
		print 'Could not load psyco'

	from sys import argv
	from os.path import basename
	for f in argv[1:]:
		current_file = basename(f)
		ll = SpmDelta_ll_common()
		parser = SpmDelta_parser(ll)
		parser.parse(f, seq_limit=8)

####################################################
#########   XXX   ##################################
####################################################

#####def finalize_unit_sp_simple(state, nr=0, finlen=None):
#####	#print "finalize_unit_sp"
#####	if finlen is None:
#####		finlen = len(state.x_buf)
#####	elif finlen == 0:
#####		return
#####	
#####	xd_l = list()
#####	xd_max = 0
#####	xo = state.x_old
#####	for x in state.x_buf[:finlen]:
#####		xd = x - xo
#####		assert xd >= 0
#####		xd_l.append(xd)
#####		xd_max = max(xd_max, xd)
#####		xo = x
#####	state.x_old = xo
#####	del state.x_buf[:finlen]
#####
#####	if state.flag_nr:
#####		yd = state.y_cur - state.y_old
#####		state.flag_nr = False
#####	else:
#####		yd = 0
#####
#####	v_l = state.v_buf[:finlen]
#####	del state.v_buf[:finlen]
#####
#####	#print "\tsize=", finlen, "ci_delta_max=", xd_max, "nr_delta=", yd
#####	state.delta.add_sp_unit(size=finlen, ci_delta_max=xd_max, nr_delta=yd)
#####	#print "\tci_delta_max=", xd_max, "ci_delta_l=", xd_l
#####	state.delta.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
#####	state.delta.addvals(values_l=v_l)
#####
#######	 ranges = range_merge(xd_l)
#######	 for start,length,cisize in ranges:
#######			ci_delta_sl = xd_l[start:start+length]
#######			ci_max = max(ci_delta_sl)
#######			v_sl = v_l[start:start+length]
#######			print "\tsize=", length, "ci_delta_max=", ci_max, "nr_delta=", yd
#######			state.delta.add_sp_unit(size=length, ci_delta_max=ci_max, nr_delta=yd)
#######			print "\tci_delta_max=", ci_max, "ci_delta_l=", ci_delta_sl
#######			state.delta.add_sp_cols(ci_delta_max=ci_max, ci_delta_l=ci_delta_sl)
#######			state.delta.addvals(values_l=v_sl)
#######			yd = 0
#####
#####
#####def finalize_unit_sp_jmp(state, nr=0, finlen=None):
#####	if finlen is None:
#####		finlen = len(state.x_buf)
#####	elif finlen == 0:
#####		return
#####	
#####	xd_l = list()
#####	xd_max = 0
#####	xo = state.x_buf[0]
#####	jmp = xo - state.x_old
#####	for x in state.x_buf[1:finlen]:
#####		xd = x - xo
#####		assert xd >= 0
#####		xd_l.append(xd)
#####		xd_max = max(xd_max, xd)
#####		xo = x
#####	state.x_old = xo
#####	del state.x_buf[:finlen]
#####
#####	if state.flag_nr:
#####		yd = state.y_cur - state.y_old
#####		state.flag_nr = False
#####	else:
#####		yd = 0
#####
#####	v_l = state.v_buf[:finlen]
#####	del state.v_buf[:finlen]
#####
#####	#print "\t[SP-JMP] size",finlen,"jmp",jmp,"ci_delta_max",xd_max,"nr_delta",yd,"ci_delta_l",xd_l
#####	state.delta.add_sp_unit_jmp(size=finlen, ci_delta=jmp, ci_delta_max=xd_max, nr_delta=yd)
#####	state.delta.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
#####	state.delta.addvals(values_l=v_l)
#####
#####def finalize_unit_sp_stats(state, nr=0, finlen=None):
#####	if finlen is None:
#####		finlen = len(state.x_buf)
#####	elif finlen == 0:
#####		return
#####	
#####	xd_l = list()
#####	xd_max = 0
#####	xo = state.x_old
#####	for x in state.x_buf[:finlen]:
#####		xd = x - xo
#####		assert xd >= 0
#####		xd_l.append(xd)
#####		xd_max = max(xd_max, xd)
#####		xo = x
#####	state.x_old = xo
#####	del state.x_buf[:finlen]
#####
#####	if state.flag_nr:
#####		yd = state.y_cur - state.y_old
#####		state.flag_nr = False
#####	else:
#####		yd = 0
#####
#####	del state.v_buf[:finlen]
#####
#####	if len(xd_l[1:]) == 0:
#####		return
#####	pattern = (pyspm.delta_cisize(max(xd_l[1:])), tuple(xd_l[1:]))
#####	if pattern in stats['sp_unit_patterns']:
#####		stats['sp_unit_patterns'][pattern] += 1
#####	else:
#####		stats['sp_unit_patterns'][pattern] = 1
#####
#####def finalize_unit_de_simple(state, nr=0):
#####	assert len(state.x_buf) == len(state.v_buf)
#####
#####	size = len(state.x_buf)
#####	xd = state.x_buf[0] - state.x_old
#####	state.x_old = state.x_buf[0] + size
#####
#####	if state.flag_nr:
#####		yd = state.y_cur - state.y_old
#####		state.flag_nr = False
#####	else:
#####		yd = 0
#####
#####	#print "\t[DE] size=", size, "ci_delta=", xd, "nr_delta=", yd
#####	state.delta.add_de_unit(size=size, ci_delta=xd, nr_delta=yd)
#####	state.delta.addvals(values_l=state.v_buf)
#####	del state.v_buf[:]
#####	del state.x_buf[:]
#####
#####def finalize_unit_de_stats(state, nr=0):
#####	assert len(state.x_buf) == len(state.v_buf)
#####
#####	size = len(state.x_buf)
#####	xd = state.x_buf[0] - state.x_old
#####	state.x_old = state.x_buf[0] + size
#####
#####	if state.flag_nr:
#####		yd = state.y_cur - state.y_old
#####		state.flag_nr = False
#####	else:
#####		yd = 0
#####
#####	del state.v_buf[:]
#####	del state.x_buf[:]
#####
#####
#####def check_simple(file):
#####	print "checking:", file
#####	global finalize_unit_sp
#####	finalize_unit_sp = finalize_unit_sp_simple
#####	delta = parse(file)
#####	delta.check(file, "spm_delta_multiply")
#####
#####def check_jmp(file):
#####	print "checking:", file
#####	global finalize_unit_sp
#####	finalize_unit_sp = finalize_unit_sp_jmp
#####	delta = parse(file)
#####	delta.check_jmp(file, "spm_delta_multiply")
#####
#####def bench_simple(file, loops=128, seq_limit=3, type=pyspm.TYPE_DOUBLE):
#####	global finalize_unit_sp
#####	finalize_unit_sp = finalize_unit_sp_simple
#####	delta = parse(file, seq_limit=seq_limit)
#####	ret = delta.bench(loops)
#####	del delta
#####	return ret
#####
#####def bench_jmp(file, loops=128, seq_limit=3, type=pyspm.TYPE_DOUBLE):
#####	global finalize_unit_sp
#####	finalize_unit_sp = finalize_unit_sp_jmp
#####	delta = parse(file, seq_limit=seq_limit)
#####	ret = delta.bench_jmp(loops)
#####	del delta
#####	return ret
#####
#####
#####def stat_patterns(file, seq_limit=4):
#####	global finalize_unit_sp, finalize_unit_de, stats
#####	finalize_unit_sp = finalize_unit_sp_stats
#####	finalize_unit_de = finalize_unit_de_stats
#####	stats['sp_unit_patterns'] = {}
#####	parse_nodelta(file, seq_limit=seq_limit)
#####	saved = 0
#####	for k, v in stats['sp_unit_patterns'].iteritems():
#####		init_size = (1<<k[0])*len(k[1])*v
#####		curr_size = (1<<k[0])*len(k[1]) + v*8
#####		saved += ( init_size - curr_size )
#####	return saved
#####
######f = "/s/matrices/sorted/006.crystk02.mtx.sorted"
######f = "/s/matrices/sorted/006.crystk02.mtx.sorted"
#####if __name__ == '__main__':
#####
#####	import os
#####	#pyspm.setaffinity(3);
#####	files_dir = "/s/matrices/sorted/"
#####	files = map( lambda x: files_dir+x, filter(
#####		lambda f: f.endswith(".mtx.sorted"), 
#####		os.listdir(files_dir)
#####	) )
#####	files.sort()
#####
#####	def all(fid, type=pyspm.TYPE_DOUBLE):
#####		f = os.path.basename(files[fid])
#####		t = type_str[type]
#####		flops = pyspm.crs_bench(files[fid], type, pyspm.CRS64)
#####		print "%s pyspm_crs64.%s %s" % (f, t, flops)
#####		flops = pyspm.crs_bench(files[fid], type, pyspm.CRS32)
#####		print "%s pyspm_cr32.%s %s" % (f, t, flops)
#####		for sl in xrange(2,13):
#####			flops = bench(fid, seq_limit=sl, type=type)
#####			print "%s pyspm_delta.%s.%s %s" % (f, t, sl, flops)
#####		flops = bench(fid, seq_limit=100000, type=type)
#####		print "%s pyspm_delta.%s.NO_DE %s" % (f, t, flops)
#####	
#####	def do_all():
#####		#print "# svn_rev: %s" % svn_rev
#####		for type in (pyspm.TYPE_DOUBLE, pyspm.TYPE_FLOAT):
#####			for fid in xrange(len(files)):
#####				all(fid, type)
#####		
#####	def main():
#####		import sys
#####		from sys import argv
#####		if len(argv) < 4:
#####			print 'usage:', argv[0], 'file type seq_limit'
#####			return 1
#####
#####		file = argv[1]
#####
#####		if argv[2] == 'float':
#####			type = pyspm.TYPE_FLOAT 
#####		else:
#####			type = pyspm.TYPE_DOUBLE
#####
#####		seq_l = int(argv[3])
#####		flops = bench_jmp(file, seq_limit=seq_l, type=type)
#####		print "%s pyspm_delta.%s.%s %s" % (os.path.basename(file), type_str[type], seq_l, flops)
#####	
#####	main()
#####
#####	#print "check(id)"
#####	#print "bench(id)"
#####	#print "all(id)"
#####	#print "do_all()"
#####	#do_all()
