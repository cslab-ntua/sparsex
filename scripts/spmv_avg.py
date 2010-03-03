#!/usr/bin/env python2.5

from stats import mean, deviation

if __name__ == '__main__':
	from sys import stdin, exit

	tl = []
	rl = []
	while True:
		l = stdin.readline()
		if l == '':
			break
		if not l.startswith('m:'):
			print l,
			continue

		m,f,s,t,r = l.split()
		tl.append(float(t[2:]))
		rl.append(float(r[2:]))

	if not tl:
		exit(0)

	tl.sort()
	rl.sort()
	tl = tl[1:-1]
	rl = rl[1:-1]
	print "%s %s %s t:[%.3f %.2f] r:[%.3f %.2f]" % \
	       (m, f, s, mean(tl), deviation(tl), mean(rl), deviation(rl))
