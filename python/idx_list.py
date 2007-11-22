def check_idx(l, new_l, idx):
	#print 'check_idx', l, new_l, idx
	match = 0
	if idx < 0 or idx >= len(l):
		#print idx, len(l), idx < 0, idx >= len(l)
		#print 'check_idx: (FAILED) returns', match
		return match

	for i in xrange(len(new_l)):
		if idx+i >= len(l):
			break
		if new_l[i] != l[idx+i]:
			break
		else:
			match += 1

	#print 'check_idx: returns', match
	return match
	
def check_direction(l, new_l, cur_idx, direction=1, limit=8):
	idx = cur_idx
	match = 0
	#print 'check direction', 'l', l, 'new_l', new_l, cur_idx, direction
	while True:
		if direction == 1 and idx >= len(l):
			break
		if direction == -1 and idx < 0:
			break

		match = check_idx(l, new_l, idx)
		if match >= min(limit, len(new_l)):
			#print 'FOUND MATCH: returning:', idx, match
			return idx, match
		idx += direction

	return idx, match

def do_add(l, cur_idx, new_l, limit=8):
	direction = -1 if len(l)/2 < cur_idx else 1
	#print "initial direction: ", direction

	idx, match = check_direction(l, new_l, cur_idx, direction, limit)
	if match >= min(limit, len(new_l)):
		return idx, match
	
	direction = (-1) * direction
	idx, match = check_direction(l, new_l, cur_idx, direction,limit)
	if match >= min(limit, len(new_l)):
		return idx, match

	idx = len(l)
	l.extend(new_l)
	return idx, len(new_l)

def add(l, new_l, cur_idx, limit=8):
	coords = []
	offset = 0
	while True:
		pl = new_l[offset:]
		idx,length = do_add(l, cur_idx, pl, limit)
		coords.append((idx, length))
		offset += length
		cur_idx = idx + length
		if offset == len(new_l):
			break

	return coords

def check_idx2(list, list_idx, new_list, new_list_idx):
	match = 0
	if list_idx < 0 or list_idx >= len(list):
		return match

	for i in xrange(len(new_list) - new_list_idx):
		if list_idx+i >= len(list):
			break
		if new_list[new_list_idx+i] != list[list_idx+i]:
			break
		else:
			match += 1
	return match

rev_index = {}

def clear_rev_index():
	global rev_index
	del rev_index
	rev_index = {}

def do_add2(list, new_list, new_list_idx, limit=8):
	global rev_index
	if new_list[new_list_idx] in rev_index:
		for index in rev_index[new_list[new_list_idx]]:
			match = check_idx2(list, new_list, new_list_idx, index)
			if match >= min(limit, len(new_list) - new_list_idx):
				return index, match

	idx = len(list)
	for i in xrange(new_list_idx, len(new_list)):
		elem = new_list[i]
		list.append(elem)
		if elem in rev_index:
			rev_index[elem].append(idx+i)
		else:
			rev_index[elem] = [idx+i]

	return idx, len(new_list)-new_list_idx 

def add2(list, new_l, limit=8):
	offset = 0
	coords = []
	while offset < len(new_l):
		idx,length = do_add2(list, new_l, offset, limit)
		coords.append((idx, length))
		offset += length
	return coords

rev_idx_nr = 2
def do_add3(list, new_list, new_list_idx, limit=8):
	global rev_index
	
	if len(new_list) < new_list_idx + rev_idx_nr:
		for index in xrange(len(list)):
			match = check_idx2(list, new_list, new_list_idx, index)
			if match >= min(limit, len(new_list) - new_list_idx):
				return index,match
	else:
		k = tuple( new_list[new_list_idx + i] for i in xrange(rev_idx_nr) ) 
		if k in rev_index:
			for index in rev_index[k]:
				match = check_idx2(list, new_list, new_list_idx, index)
				if match >= min(limit, len(new_list) - new_list_idx):
					return index, match

	idx = len(list)
	j = 0
	if len(list) == 0:
		for i in xrange(rev_idx_nr-1):
			list.append(new_list[new_list_idx+i])
			j += 1

	for i in xrange(new_list_idx+j, len(new_list)):
		elem = new_list[i]
		key =  tuple( 
			new_list[new_list_idx - i] 
			for i in xrange(rev_idx_nr).__reversed__()
		)
		list.append(elem)
		if key in rev_index:
			rev_index[key].append(idx+i)
		else:
			rev_index[key] = [idx+i]

	return idx, len(new_list)-new_list_idx 

def add3(list, new_l, limit=8):
	offset = 0
	coords = []
	while offset < len(new_l):
		idx,length = do_add3(list, new_l, offset, limit)
		coords.append((idx, length))
		offset += length
	return coords


def clist_gen(c_list, c_coords):
	for idx, len in c_coords:
		for i in xrange(idx, idx+len):
			yield c_list[i]

if __name__ == '__main__':
	import random
	from time import time

	print "Creating random data"
	lists = []
	orig_list = []
	total_len = 0
	for i in xrange(500):
		l = []
		for j in xrange(random.randint(5,100)):
			l.append(random.randint(0,5))
		lists.append(l)
		orig_list.extend(l)
		total_len += len(l)
	
	
	print "Executing with simple add() ...",
	c_list = []
	c_coords = []
	cur_idx = 0
	t0 = time()
	for l in lists:
		coords = add(c_list, l, cur_idx)
		c_coords.extend(coords)
		cur_idx = c_coords[-1][0] + c_coords[-1][1]
	print " (", time() - t0, " secs) ... ",

	final_list = []
	for e in clist_gen(c_list, c_coords):
		final_list.append(e)
	
	assert final_list == orig_list
	print " OK"
	print 'total_len', total_len, 'idx_len', len(c_list) + len(c_coords)

	del c_list[:]
	del c_coords[:]
	del final_list[:]

	print "Executing with simple rev_index add() ...",
	t0 = time()
	for l in lists:
		coords = add2(c_list, l)
		c_coords.extend(coords)
	print " (", time() - t0, " secs) ... ",
	
	for e in clist_gen(c_list, c_coords):
		final_list.append(e)
	
	assert final_list == orig_list
	print " OK"
	print 'total_len', total_len, 'idx_len', len(c_list) + len(c_coords)

	clear_rev_index()
	del c_list[:]
	del c_coords[:]
	del final_list[:]
	print "Executing with double rev_index add() ...",
	t0 = time()
	for l in lists:
		coords = add3(c_list, l, limit=1000000000)
		c_coords.extend(coords)
	print " (", time() - t0, " secs) ... ",
	
	for e in clist_gen(c_list, c_coords):
		final_list.append(e)
	
	assert final_list == orig_list
	print " OK"
	print 'total_len', total_len, 'idx_len', len(c_list) + len(c_coords)

	#print orig_list
	#print check_list

	#l = range(100)
	#i = len(l) - 1
	#n = range(40,101)
	#check(l, n, i)
	#check(l, n, 0)
