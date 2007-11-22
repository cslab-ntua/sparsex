import pyspm

def delta_xd_size(xd):
	return 1<<pyspm.delta_cisize(xd)

def delta_ul_size(ul):
	i = 1
	shift = 7
	while True:
		max = (1<<shift) - 1
		if ul <= max:
			return i
		i += 1
		shift += 7

def delta_l_size(l):
	al = abs(l)
	shift = 6
	i = 1
	while True:
		max = (1<<shift) - 1
		if al <= max:
			return i
		i += 1
		shift += 7

def delta_de(state, finlen=None):
	if finlen is not None:
		raise NotImplementedError

	assert len(state.x_buf) == len(state.v_buf)
	size = len(state.x_buf)
	xd = state.x_buf[0] - state.x_old
	state.x_old = state.x_buf[0] + size

	if state.flag_nr:
		yd = state.y_cur - state.y_old
		state.flag_nr = False
	else:
		yd = 0

	#self.units.append(('DE', size, xd))
	#self.values.extend(state.v_buf)
	v_l = list(state.v_buf)

	del state.v_buf[:]
	del state.x_buf[:]

	return size, xd, yd, v_l

def delta_sp(state, finlen=None):
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

	return xd_l, xd_max, yd, v_l
	
def delta_sp_jump(state, finlen=None):
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
	
	return jmp, xd_l, xd_max, yd, v_l

