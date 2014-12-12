#!/bin/env python
import re
import sys
import os

"""
Find events that I reject and that are selected in the passs list
"""

discarded = open(sys.argv[1], 'r')
if len(sys.argv)==3:
	outputFileName = sys.argv[2];
else:
	outputFileName = "problems.dat"

#passs = open(os.path.expanduser("~mkoval/public/data-5e6_p1.txt"), 'r')
passs = open(os.path.expanduser("~goudzovs/public/list5"), 'r')
result = open(outputFileName, 'w')

r = re.compile('[0-9]*? [0-9]*? [0-9]*? [0-9]*?')

runs = {}

# for line in passs:
# 	[run, burst, event] = [int(n) for n in line.split()]
# 	if not run in runs:
# 		runs[run] = {}
# 	if not burst in runs[run]:
# 		runs[run][burst] = set()
# 	runs[run][burst].update([event])

run = 0
for line in passs:
	spl = line.split()
 	#[run, burst, event] = [int(n) for n in spl]
 	[burst, event] = [int(n) for n in spl[:-1]]
 	#print [run,burst,event]
 	if not run in runs:
 		runs[run] = {}
 	if not burst in runs[run]:
 		runs[run][burst] = set()
 	runs[run][burst].update([event])

#total = sum(1 for line in discarded)
total = 0

passs.close()
#discarded.seek(0,0)
i=0
for line in discarded:
    i = i + 1
    if (i % 100)==0:
        sys.stdout.write("Checking line %i/%i\r" % (i,total))
        sys.stdout.flush()
    if r.match(line) is not None:
		spl = line.split()
		[run, burst, event] = [int(n) for n in spl[:-1]]
		cut = spl[-1]
		run = 0
		#print [run,burst,event]
		if run in runs:
			if burst in runs[run]:
				if event in runs[run][burst]:
					print "                                                   \r%s %s %s %s" % (run, burst, event, cut)
					result.write("%s %s %s\n" % (run, burst, event))
                
                
