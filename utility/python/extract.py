#!/bin/env python
import sys
import re

if len(sys.argv)!=7:
	print "firstrun firstburst lastburst firstevent lastevent inputfilename"
	sys.exit(0)
Frun = int(sys.argv[1])
Fburst = int(sys.argv[2])
Lburst = int(sys.argv[3])
Fevt = int(sys.argv[4])
Levt = int(sys.argv[5])
fName = sys.argv[6]

f = open(fName, 'r')
fo = open("extract.dat", 'w')

r = re.compile('[0-9]*? [0-9]*? [0-9]*?')

total = sum(1 for line in f)

f.seek(0,0)
i = 0
for line in f:
	i = i + 1
	if (i % 100)==0:
		print "Checking line %i/%i\r" % (i,total)
	if r.match(line) is not None:
		[run, burst, event] = line.split()
		if (int(run)==Frun) and (int(burst)>Fburst) and (int(burst)<Lburst):
			fo.write(line)
			
