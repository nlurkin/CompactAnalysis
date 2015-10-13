#!/bin/env python

import rootpy.ROOT as ROOT 
ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

from rootpy.tree import Tree
from rootpy.io import root_open
import sys
#import rootpy.compiled as C
from ROOT import ROOTBurst, ROOTRawEvent


if len(sys.argv)!=3:
    print "Not enough arguments. Expecting: input.root output.txt"
    sys.exit(0)

#C.register_file("../../userinc/exportClasses.h", ["NVtxTrack"])

f = root_open(sys.argv[1], "open")

#print list(f)

tree = f.event

with open(sys.argv[2], "w") as fd:
    for event in tree:
        fd.write("%i %i %i\n" % (event.rawBurst.nrun, event.rawBurst.time, event.rawEvent.timeStamp))