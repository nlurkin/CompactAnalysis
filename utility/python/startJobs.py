#!/bin/env python
import os
import sys

script = """
export HOME=/afs/cern.ch/user/n/nlurkin
echo $HOME
cd ~/batch/output_%s
source ~/env.sh
~/Compact/compactjob -string prefix=job_%s%s:can=pi0d:nooutput=1 -i %s 

#~/Compact/outputs/fails.py job%s.txt problemsjob%s.dat
#~/Compact/outputs/success.py job%spass.txt problemsSelectjob%s.dat

#rm job%s.txt
#rm job%spass.txt
#rm job%s.root
#touch job%s.root
"""

if len(sys.argv)!=3:
	print "Missing argument: prefix listfile"
	sys.exit(0)
prefix = sys.argv[1]
listFileName = sys.argv[2]

filesList = os.path.expanduser(listFileName)
files = open(filesList, "r")
i=0
submitted=0
for line in files:
	if not os.path.exists(os.path.expanduser("~/batch/output_%s/job_%s%s.root" % (prefix,prefix,i))):
		sh = open(os.path.expanduser("~/batch/scripts_%s/job_%s%s.sh" % (prefix, prefix, i)), "w")
		sh.write(script % (prefix,prefix, i,line, i,i,i,i,i,i,i,i))
		sh.close()
		os.chmod(os.path.expanduser("~/batch/scripts_%s/job_%s%s.sh" % (prefix, prefix,i)), 0744)
		cmd = 'bsub -q 1nh -R "type=SLC5_64" ~/batch/scripts_%s/job_%s%s.sh' % (prefix,prefix,i)
		print cmd
		os.system(cmd)
		submitted+=1
	i=i+1

files.close()

print "Submitted total of " + str(submitted) + " jobs (out of " + str(i) + " files)"