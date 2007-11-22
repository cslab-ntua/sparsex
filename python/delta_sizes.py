import pyspm
from delta_helpers import delta_sp_jump, delta_de, delta_ul_size, delta_xd_size, delta_l_size
from crs_helpers import crs_size

from os.path import basename
from sys import stdout

class Delta_size(object):
	def __init__(self, keepers):
		self.keepers = []
		for k in keepers:
			self.keepers.append(k)
		
	def finalize_unit_sp(self, state, finlen=None):
		jmp, xd_l, xd_max, yd, v_l = delta_sp_jump(state, finlen)
		for k in self.keepers:
			k.update_sp(jmp, xd_l, xd_max, yd, v_l)

	def finalize_unit_de(self, state, finlen=None):
		size, xd, yd, v_l = delta_de(state, finlen)
		for k in self.keepers:
			k.update_de(size, xd, yd, v_l)
		
	def finalize(self, state):
		v_size, i_size = crs_size(state.nr_nz, state.nr_rows)
		v_size = float(v_size)
		i_size = float(i_size)
		total_size = v_size+i_size
		f = basename(state.mmf_file)
		print f, 'crs size:', total_size, 'index percentage:', (i_size/total_size)*100
		for k in self.keepers:
			s = float(k.get_size())
			print f, k.name, 'size', s, 'reduction', ((total_size-s)/total_size)*100
		stdout.flush()
	
class Delta_jmp(object):
	name = "delta_jmp"
	
	def __init__(self, v_size=8):
		self.size = 0
		self._v_size = v_size

	def update_sp(self, jmp, xd_l, xd_max, yd, v_l):
		self.size += 1 + 1                            # flags + size
		self.size += delta_ul_size(jmp)               # jump
		self.size += delta_xd_size(xd_max)*len(xd_l)  # list of indices
		self.size += self._v_size*len(v_l)            # values
		if yd > 1:
			self.size += delta_ul_size(yd)            # yd
	
	def update_de(self, size, xd, yd, v_l):
		self.size += 1 + 1                           # flags + size
		self.size += delta_ul_size(xd)               # jump
		self.size += self._v_size*len(v_l)           # values
		if yd > 1:
			self.size += delta_ul_size(yd)           # yd
	
	def get_size(self):
		return self.size

class Delta_ci(object):
	name = "delta_jmp_ci"
	def __init__(self, v_size=8):
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
		self.sp_patterns = {
			pyspm.DELTA_CISIZE_U8  : {},
			pyspm.DELTA_CISIZE_U16 : {},
			pyspm.DELTA_CISIZE_U32 : {},
			pyspm.DELTA_CISIZE_U64 : {},
		}
		self._v_size = v_size
		self.size = 0

	def update_sp_indices(self, xd_l, xd_max):
		if len(xd_l) > 0:
			cisize = pyspm.delta_cisize(xd_max)
			sp_colinds = self.sp_colinds[cisize]
			sp_patterns = self.sp_patterns[cisize]
			pattern = (cisize, tuple(xd_l))
			if pattern in sp_patterns:
				idx = sp_patterns[pattern]
			else:
				idx = sp_patterns[pattern] = len(sp_colinds)
				sp_colinds.extend(xd_l)
				self.size += len(xd_l)*(1<<cisize)
			self.size += delta_l_size(idx - self.sp_cur_idx[cisize])
			self.sp_cur_idx[cisize] = idx

	def update_sp(self, jmp, xd_l, xd_max, yd, v_l):
		self.size += 1 + 1                            # flags + size
		self.size += delta_ul_size(jmp)               # jump

		self.size += self._v_size*len(v_l)            # values

		if yd > 1:
			self.size += delta_ul_size(yd)            # yd
		
		self.update_sp_indices(xd_l, xd_max)

	def update_de(self, size, xd, yd, v_l):
		self.size += 1 + 1                           # flags + size
		self.size += delta_ul_size(xd)               # jump
		self.size += self._v_size*len(v_l)           # values
		if yd > 1:
			self.size += delta_ul_size(yd)           # yd
	
	def get_size(self):
		return self.size

class Delta_cv(object):
	name = "delta_jmp_cv"
	def __init__(self, i_size=4, v_size=8):
		self.values = {}
		self.prev_v_idx = 0
		self.size = 0
		self._i_size = i_size
		self._v_size = v_size

	def update_vals(self, v_l):
		for v in v_l:
			assert isinstance(v, str)
			if v not in self.values:
				self.size += self._v_size
				idx = self.values[v] = len(self.values)
			else:
				idx = self.values[v]
			self.size += delta_l_size(idx - self.prev_v_idx)
			self.prev_v_idx = idx

	def update_sp(self, jmp, xd_l, xd_max, yd, v_l):
		self.size += 1 + 1                            # flags + size
		self.size += delta_ul_size(jmp)               # jump
		self.size += delta_xd_size(xd_max)*len(xd_l)  # list of indices
		if yd > 1:
			self.size += delta_ul_size(yd)            # yd
		self.update_vals(v_l)

	def update_de(self, size, xd, yd, v_l):
		self.size += 1 + 1                           # flags + size
		self.size += delta_ul_size(xd)               # jump
		if yd > 1:
			self.size += delta_ul_size(yd)           # yd
		self.update_vals(v_l)
	
	def get_size(self):
		return self.size

class Delta_cvi(Delta_ci, Delta_cv):
	name = "delta_jmp_cvi"
	def __init__(self, i_size=4, v_size=8):
		Delta_cv.__init__(self, i_size=i_size, v_size=v_size)
		Delta_ci.__init__(self, v_size=v_size)

	def update_sp(self, jmp, xd_l, xd_max, yd, v_l):
		self.size += 1 + 1
		self.size += delta_ul_size(jmp)
		if yd > 1:
			self.size += delta_ul_size(yd)

		self.update_sp_indices(xd_l, xd_max)
		self.update_vals(v_l)
	
	def update_de(self, size, xd, yd, v_l):
		self.size += 1 + 1
		self.size += delta_ul_size(xd)
		if yd > 1:
			self.size += delta_ul_size(yd)
		self.update_vals(v_l)
	
	def get_size(self):
		return self.size

class Delta_zi(object):
	name = "delta_zi"
	def __init__(self):
		pass

	def update_sp(self, jmp, xd_l, xd_max, yd, v_l):
		pass
	
	def update_de(self, size, xd, yd, v_l):
		pass
	
	def get_size(self):
		pass

from delta_parse2 import SpmDelta_parser
def do_sizes(f):
	ll = Delta_size(
		keepers = (Delta_jmp(), Delta_cv(), Delta_ci(), Delta_cvi())
	)
	parser = SpmDelta_parser(ll)
	parser.parse(f, seq_limit=8, float_str=True)

if __name__ == '__main__':
	from sys import argv
	for f in argv[1:]:
		do_sizes(f)
