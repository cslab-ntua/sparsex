#!/usr/bin/env python2.5

def main():
	import sys
	vals = {}

	while True:
		line = sys.stdin.readline()
		#sys.stderr.write("--got: " + line)
		if not line:
			break
		key = line.rstrip("\n")
		if key in vals:
			id = vals[key]
		else:
			id = vals[key] = len(vals)
		#sys.stderr.write("--returning: " + str(id) + "\n")
		sys.stdout.write(str(id) + "\n")
		sys.stdout.flush()
		 
if __name__ == '__main__':
	main()
