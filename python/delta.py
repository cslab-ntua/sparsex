import pyspm

from delta_helpers import delta_sp, delta_sp_jump, delta_de

class Delta(object):
	def __init__(self, type=pyspm.TYPE_DOUBLE):
		self.delta = pyspm.delta(type=type)
		self.delta.prepare()

	def finalize_unit_sp(self, state, finlen=None):
		xd_l, xd_max, yd, v_l = delta_sp(state, finlen)
		#print "--->", xd_l, xd_max, yd, v_l
		self.delta.add_sp_unit(size=len(xd_l), ci_delta_max=xd_max, nr_delta=yd)
		self.delta.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
		self.delta.addvals(values_l=v_l)
	
	def finalize_unit_de(self, state, finlen=None):
		size, xd, yd, v_l = delta_de(state, finlen)
		#print "===>", size, xd, yd, v_l
		self.delta.add_de_unit(size=size, ci_delta=xd, nr_delta=yd)
		self.delta.addvals(values_l=v_l)

	def finalize(self, state):
		delta = self.delta
		del self.delta
		delta.finalize(state.nr_nz, state.nr_rows, state.nr_cols)
		return delta

class Delta_jmp(Delta):
	def finalize_unit_sp(self, state, finlen=None):
		d = self.delta
		jmp, xd_l, xd_max, yd, v_l = delta_sp_jump(state, finlen)
		d.add_sp_unit_jmp(size=len(xd_l)+1, ci_delta=jmp, ci_delta_max=xd_max, nr_delta=yd)
		d.add_sp_cols(ci_delta_max=xd_max, ci_delta_l=xd_l)
		d.addvals(values_l=v_l)

def delta_parse(f):
	from delta_parse import SpmDelta_parser
	ll = Delta()
	parser = SpmDelta_parser(ll)
	return parser.parse(f, seq_limit=8)

def delta_jmp_parse(f, seq_limit=8):
	from delta_parse import SpmDelta_parser
	ll = Delta_jmp()
	parser = SpmDelta_parser(ll)
	return parser.parse(f, seq_limit=seq_limit)

def delta_get_size(d):
	props = d.getprops()
	return props['ctl_size'] + props['nnz']*8

if __name__ == '__main__':
	from sys import argv
	import os
	seq_limit = os.getenv("SPM_DELTA_SEQ_LIMIT")
	seq_limit = int(seq_limit) if seq_limit is not None else 8
	for f in argv[1:]:
		d = delta_jmp_parse(f, seq_limit)
		size = delta_get_size(d)
		time, flops = d.bench_jmp()
		print "pyspm_delta_jmp.%s %s %d %s %s" % (seq_limit, os.path.basename(f), size, time, flops)
