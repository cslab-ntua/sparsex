from os import popen2
from random import randint
from time import sleep

loops = 1000000
range = 1024

def t1():
	(chin, chout) = popen2('./phash_i')
	d = {}
	for i in xrange(loops):
		k = randint(0, range)
		v = randint(0, range)
		d[k] = v
		chin.write("I %d %d\n" % (k,v))
		chin.flush()
	tf(chin, chout, d)

def t2():
	(chin, chout) = popen2('./phash_i')
	d = {}
	for i in xrange(loops):
		k = randint(0, range)
		v = randint(0, range)
		d[k] = v + d.get(k, 0)
		chin.write("U %d %d\n" % (k,v))
		chin.flush()
	tf(chin, chout, d)


def tf(chin, chout, d):
	chin.write("S\n")
	chin.flush()
	res = chout.readline()
	if len(d) != int(res):
		print 'different size:', int(res), 'instead of', len(d)
		assert 0

	for k,v in d.iteritems():
		chin.write("G %d\n" % k)
		chin.flush()
		res = chout.readline()
		if v != int(res):
			print 'waited for', v, 'and got', int(res), '(key is', k, ')'
			assert 0

		chin.write("D %d\n" % k)
		chin.flush()

	chin.write("S\n")
	chin.flush()
	res = chout.readline()
	if int(res) != 0:
		print 'different size:', int(res), 'instead of', 0
		assert 0

	chin.close()


if __name__ == '__main__':
	t1()
	t2()
