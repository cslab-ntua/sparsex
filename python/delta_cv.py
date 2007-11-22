import pyspm

from delta_helpers import delta_sp_jump, delta_de

class Delta_cv(object):
	def __init__(self, type=pyspm.TYPE_DOUBLE, seq_limit=8):
		self.values = {}
		self.delta_cv = pyspm.delta_cv()
		self.delta_cv.prepare()
		self.cur_y = 0
		self.units = set()
	
	def update_vals(self, v_l):
		#print " updating vals", v_l
		vals = []
		vis = []
		idx_max = 0
		for v in v_l:
			assert isinstance(v, str)
			if v not in self.values:
				idx = self.values[v] = len(self.values)
				vals.append(v)
			else:
				idx = self.values[v]

			vis.append(idx)
			if idx > idx_max:
				idx_max=idx

		vi_s = pyspm.delta_cisize(idx_max)
		if len(vals) != 0:
			#print "   inserting new vals", list([float(v) for v in vals])
			self.delta_cv.add_vals(v_l=list([float(v) for v in vals]))
		#print "   inserting val ids", vis, "size", vi_s
		self.delta_cv.add_vi(vi_s=vi_s, vi_l=vis)

		return vi_s
		
	def finalize_unit_sp(self, state, finlen=None):
		jmp, xd_l, xd_max, yd, v_l = delta_sp_jump(state, finlen)
		self.cur_y += yd
		#print "-- finalize_unit_sp", 
		#print "jmp:", jmp, "xd_l:", xd_l, "xd_max:", xd_max, "yd:", yd, "v_l:", v_l
		vi_size = self.update_vals(v_l)
		ci_size = pyspm.delta_cisize(xd_max)

		#print " add_unit u_s:", len(xd_l), "ci_s:", ci_size, "vi_s:", vi_size, "ci_jmp:", jmp, "ri_jmp:", yd, "y:", self.cur_y
		self.delta_cv.add_unit(
			u_s=len(xd_l), 
			ci_s=ci_size, vi_s=vi_size, ci_jmp=jmp, ri_jmp=yd, 
			sp=1
		)
		self.units.add((1, ci_size, vi_size))
		#print " add_ci ci_s:", ci_size, "ci_delta_l:", xd_l
		self.delta_cv.add_ci(ci_s=ci_size, ci_delta_l=xd_l)
			
	def finalize_unit_de(self, state, finlen=None):
		#print "-- finalize_unit_de", 
		size, jmp, yd, v_l = delta_de(state, finlen)
		self.cur_y += yd
		vi_size = self.update_vals(v_l)

		#print " add_unit u_s:", size, "vi_s:", vi_size, "ci_jmp:", jmp, "ri_jmp:", yd, "y:", self.cur_y
		self.delta_cv.add_unit(
			u_s=size,
			ci_s=0, vi_s=vi_size, ci_jmp=jmp, ri_jmp=yd, 
			sp=0
		)
		self.units.add((0, vi_size))
	
	def finalize(self, state):
		delta_cv = self.delta_cv
		del self.delta_cv
		delta_cv.finalize(state.nr_nz, state.nr_rows, state.nr_cols)
		print "units:", self.units
		return delta_cv

def delta_cv_parse(f):
	from delta_parse2 import SpmDelta_parser
	ll = Delta_cv()
	parser = SpmDelta_parser(ll)
	return parser.parse(f, seq_limit=8, float_str=True)

if __name__ == '__main__':
	from sys import argv
	import os
	for f in argv[1:]:
		dcv = delta_cv_parse(f)
		flops = dcv.bench()
		print "%s pyspm_delta_cv %s" % (os.path.basename(f), flops)
