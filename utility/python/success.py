#!/bin/env python
import os
import re
import sys

"""
Find events that I select and that are not selected in the ref list.
"""

select = open(sys.argv[1], 'r')
if len(sys.argv)==3:
	outputFileName = sys.argv[2];
else:
	outputFileName = "problems.dat"

#ref = open(os.path.expanduser("~mkoval/public/data-5e6_p1.txt"), 'r')
ref = open(os.path.expanduser("~goudzovs/public/list14"), 'r')
result = open(outputFileName, 'w')

r = re.compile('[0-9]*? [0-9]*? [0-9]*?')

runs = {}

#for line in ref:
#	[run, burst, event] = [int(n) for n in line.split()]
#	if not run in runs:
#		runs[run] = {}
#	if not burst in runs[run]:
#		runs[run][burst] = set()
#	runs[run][burst].update([event])

run = 0
for line in ref:
	spl = line.split()
 	#[run, burst, event] = [int(n) for n in spl]
 	[burst, event] = [int(n) for n in spl[:-1]]
 	#print [run, burst, event]
 	if not run in runs:
 		runs[run] = {}
 	if not burst in runs[run]:
 		runs[run][burst] = set()
 	runs[run][burst].update([event])
 	
#total = sum(1 for line in ref)
total = 0

ref.close()

i=0
print "Total lines " + str(total)
for line in select:
    i = i + 1
    if (i % 100)==0:
        sys.stdout.write("Checking line %i/%i\r" % (i,total))
        sys.stdout.flush()
    if r.match(line) is not None:
        [run, burst, event] = [int(n) for n in line.split()]
        good = False
        run = 0
        #print [run,burst,event]
        if run in runs:
        	if burst in runs[run]:
        		if event in runs[run][burst]:
        			good = True
        
        if good == False:
            print "                                                   \r%s %s %s" % (run, burst, event)
            result.write("%s %s %s\n" % (run, burst, event))
                
                
