#!/bin/env python

import rootpy.ROOT as ROOT 
from signal import pause

ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

ROOT.gSystem.AddDynamicPath("/afs/cern.ch/user/n/nlurkin/Compact/ROOTComp/build/Support")
ROOT.gSystem.Load("libCuts")

from rootpy.tree import Tree
from rootpy.io import root_open
from rootpy.plotting import Hist2D, Graph2D, Canvas, Hist1D
import sys
from rootpy.interactive import wait

ROOT.gStyle.SetOptStat(0)
def printDiv(a, b, c=0, d=0, e=0, f=0):
    s = "{a:10}"
    
    dico = {"a":a, "b":"", "c":"", "d":"", "e":"", "f":""}
    if b!=0:
        dico["b"] = a*100./b
        s += "\t{b:10.3f}"
    if c!=0:
        dico["c"] = a*100./c
        s += "\t{c:10.3f}"
    if d!=0:
        dico["d"] = a*100./d
        s += "\t{d:10.3f}"
    if e!=0:
        dico["e"] = a*100./e
        s += "\t{e:10.3f}"
    if f!=0:
        dico["f"] = a*100./f
        s += "\t{f:10.3f}"
        
    return s.format(**dico);

if len(sys.argv)<2:
    print "Not enough arguments. Expecting: input.root [scanID]"
    sys.exit(0)

f = root_open(sys.argv[1], "open")

ScanIndex = 0;
if len(sys.argv)==3:
    ScanIndex = int(sys.argv[2]);

cutsSize = len(list(f.pid_pi)[0]["pid_res.good"])

class alt_pid_res():
    good=[0]*cutsSize
    bad=[0]*cutsSize

    class total:
        events = [0]*cutsSize
        class prelim:
            events = [0]*cutsSize
            class noID:
                events = [0]*cutsSize
            class manyID:
                events = [0]*cutsSize
            class associated:
                events = [0]*cutsSize
                class noID:
                    events = [0]*cutsSize
                class manyID:
                    events = [0]*cutsSize
                class ided:
                    events = [0]*cutsSize
                    class good:
                        events = [0]*cutsSize
                        class pass_t:
                            events = [0]*cutsSize
                        class noPass:
                            events = [0]*cutsSize
                    class bad:
                        events = [0]*cutsSize
                        class pass_t:
                            events = [0]*cutsSize
                        class noPass:
                            events = [0]*cutsSize
            class ided:
                events = [0]*cutsSize
                class pass_t:
                    events = [0]*cutsSize
                class noPass:
                    events = [0]*cutsSize

class alt_pid_res_mu():
    good=[0]*cutsSize
    bad=[0]*cutsSize

    class total:
        events = [0]*cutsSize
        class prelim:
            events = [0]*cutsSize
            class noID:
                events = [0]*cutsSize
            class manyID:
                events = [0]*cutsSize
            class associated:
                events = [0]*cutsSize
                class noID:
                    events = [0]*cutsSize
                class manyID:
                    events = [0]*cutsSize
                class ided:
                    events = [0]*cutsSize
                    class good:
                        events = [0]*cutsSize
                        class pass_t:
                            events = [0]*cutsSize
                        class noPass:
                            events = [0]*cutsSize
                    class bad:
                        events = [0]*cutsSize
                        class pass_t:
                            events = [0]*cutsSize
                        class noPass:
                            events = [0]*cutsSize
            class ided:
                events = [0]*cutsSize
                class pass_t:
                    events = [0]*cutsSize
                class noPass:
                    events = [0]*cutsSize
NProcessedEvents = 0
NFailedEvents = 0
NPassedEvents = 0
pid_res = alt_pid_res();
pid_res_mu = alt_pid_res_mu();

for i, h in enumerate(f.header):
#    print "Reading {0};{1};{2};{3}".format(i, h.header.NProcessedEvents, h.header.NFailedEvents, h.header.NPassedEvents)
    NProcessedEvents += h.header.NProcessedEvents
    NFailedEvents += h.header.NFailedEvents
    NPassedEvents += h.header.NPassedEvents

if hasattr(f, "pid_pi") or True:
    for i, event in enumerate(f.pid_pi):
        for index in range(cutsSize):
            pid_res.good[index] += event["pid_res.good"][index]
            pid_res.bad[index] += event["pid_res.bad"][index]
            pid_res.total.events[index] += event["pid_res.total"][index]
            pid_res.total.prelim.events[index] += event["pid_res.total.prelim"][index]
            pid_res.total.prelim.associated.events[index] += event["pid_res.total.prelim.associated"][index]
            pid_res.total.prelim.associated.noID.events[index] += event["pid_res.total.prelim.associated.noID"][index]
            pid_res.total.prelim.associated.manyID.events[index] += event["pid_res.total.prelim.associated.manyID"][index]
            pid_res.total.prelim.associated.ided.events[index] += event["pid_res.total.prelim.associated.ided"][index]
            pid_res.total.prelim.associated.ided.good.events[index] += event["pid_res.total.prelim.associated.ided.good"][index]
            pid_res.total.prelim.associated.ided.bad.events[index] += event["pid_res.total.prelim.associated.ided.bad"][index]
            pid_res.total.prelim.noID.events[index] += event["pid_res.total.prelim.noID"][index]
            pid_res.total.prelim.manyID.events[index] += event["pid_res.total.prelim.manyID"][index]
            pid_res.total.prelim.associated.ided.good.pass_t.events[index] += event["pid_res.total.prelim.associated.ided.good.pass"][index]
            pid_res.total.prelim.associated.ided.good.noPass.events[index] += event["pid_res.total.prelim.associated.ided.good.noPass"][index]
            pid_res.total.prelim.associated.ided.bad.pass_t.events[index] += event["pid_res.total.prelim.associated.ided.bad.pass"][index]
            pid_res.total.prelim.associated.ided.bad.noPass.events[index] += event["pid_res.total.prelim.associated.ided.bad.noPass"][index]
            pid_res.total.prelim.ided.events[index] += event["pid_res.total.prelim.ided"][index]
            pid_res.total.prelim.ided.pass_t.events[index] += event["pid_res.total.prelim.ided.pass"][index]
            pid_res.total.prelim.ided.noPass.events[index] += event["pid_res.total.prelim.ided.noPass"][index]
    
if hasattr(f, "pid_mu"):
    for i, event in enumerate(f.pid_mu):
        for index in range(cutsSize):
            pid_res_mu.good[index] += event["pid_res.good"][index]
            pid_res_mu.bad[index] += event["pid_res.bad"][index]
            pid_res_mu.total.events[index] += event["pid_res.total"][index]
            pid_res_mu.total.prelim.events[index] += event["pid_res.total.prelim"][index]
            pid_res_mu.total.prelim.associated.events[index] += event["pid_res.total.prelim.associated"][index]
            pid_res_mu.total.prelim.associated.noID.events[index] += event["pid_res.total.prelim.associated.noID"][index]
            pid_res_mu.total.prelim.associated.manyID.events[index] += event["pid_res.total.prelim.associated.manyID"][index]
            pid_res_mu.total.prelim.associated.ided.events[index] += event["pid_res.total.prelim.associated.ided"][index]
            pid_res_mu.total.prelim.associated.ided.good.events[index] += event["pid_res.total.prelim.associated.ided.good"][index]
            pid_res_mu.total.prelim.associated.ided.bad.events[index] += event["pid_res.total.prelim.associated.ided.bad"][index]
            pid_res_mu.total.prelim.noID.events[index] += event["pid_res.total.prelim.noID"][index]
            pid_res_mu.total.prelim.manyID.events[index] += event["pid_res.total.prelim.manyID"][index]
            pid_res_mu.total.prelim.associated.ided.good.pass_t.events[index] += event["pid_res.total.prelim.associated.ided.good.pass"][index]
            pid_res_mu.total.prelim.associated.ided.good.noPass.events[index] += event["pid_res.total.prelim.associated.ided.good.noPass"][index]
            pid_res_mu.total.prelim.associated.ided.bad.pass_t.events[index] += event["pid_res.total.prelim.associated.ided.bad.pass"][index]
            pid_res_mu.total.prelim.associated.ided.bad.noPass.events[index] += event["pid_res.total.prelim.associated.ided.bad.noPass"][index]
            pid_res_mu.total.prelim.ided.events[index] += event["pid_res.total.prelim.ided"][index]
            pid_res_mu.total.prelim.ided.pass_t.events[index] += event["pid_res.total.prelim.ided.pass"][index]
            pid_res_mu.total.prelim.ided.noPass.events[index] += event["pid_res.total.prelim.ided.noPass"][index]

print pid_res.total.prelim.ided.noPass.events[index], pid_res_mu.total.prelim.ided.noPass.events[index]
print "Processed events ->\t{0}".format(NProcessedEvents);
print "Failed events ->\t{0}".format(NFailedEvents);
print "Passed events ->\t{0}".format(NPassedEvents);

print pid_res.total.prelim.ided.pass_t.events
showIndex = ScanIndex
#for i,l in enumerate(f.cutsDefinition):
i = 0
if i==showIndex:
    #d = l.lists
    showIndex = i
#         print "\nShowing ScanIndex {0}: mpi={1}, mkmin={2}, mkmax={3}".format(i, d.maxPi0MassDiff, d.minKaonMassDiff, d.maxKaonMassDiff)

    #print "\nShowing ScanIndex {0}: cut: {1}".format(showIndex, d.unDeflectedElDist)
    print "\nMC Association (track)\n--------------";
    #print "Good = {0}\t{1}".format(pid_res.good, pid_res.good*100./(pid_res.good+pid_res.bad))
    #print "Bad = {0}\t{1}".format(pid_res.bad, pid_res.bad*100./(pid_res.good+pid_res.bad))
    print "\t\t\t{0:10}\t{1:10}\t{2:10}\t{3:10}\t{4:10}".format(0, 1, 2, 3, 4)
    print "Total: \t\t\t{0}".format(printDiv(pid_res.total.events[showIndex], pid_res.total.events[showIndex]))
    print "|-Prelim: \t\t{0}".format(printDiv(pid_res.total.prelim.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex]))
    print "  |-NoID: \t\t{0}".format(printDiv(pid_res.total.prelim.noID.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex]))
    print "  |-ManyID: \t\t{0}".format(printDiv(pid_res.total.prelim.manyID.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex]))
    print "  |-Ided: \t\t{0}".format(printDiv(pid_res.total.prelim.ided.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.ided.events[showIndex]))
    print "    |-Pass: \t\t{0}".format(printDiv(pid_res.total.prelim.ided.pass_t.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.ided.events[showIndex]))
    print "    |-NoPass: \t\t{0}".format(printDiv(pid_res.total.prelim.ided.noPass.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.ided.events[showIndex]))
    print "  |-Associated: \t{0}".format(printDiv(pid_res.total.prelim.associated.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex]))
    print "  | |-NoID: \t\t{0}".format(printDiv(pid_res.total.prelim.associated.noID.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex]))
    print "  | |-ManyID: \t\t{0}".format(printDiv(pid_res.total.prelim.associated.manyID.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex]))
    print "  | |-Ided: \t\t{0}".format(printDiv(pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex]))
    print "  |   |-Good: \t\t{0}".format(printDiv(pid_res.total.prelim.associated.ided.good.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   | |-Pass: \t{0}".format(printDiv(pid_res.total.prelim.associated.ided.good.pass_t.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   | |-NoPass: \t{0}".format(printDiv(pid_res.total.prelim.associated.ided.good.noPass.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   |-Bad: \t\t{0}".format(printDiv(pid_res.total.prelim.associated.ided.bad.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.bad.events[showIndex]))
    print "  |     |-Pass: \t{0}".format(printDiv(pid_res.total.prelim.associated.ided.bad.pass_t.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.bad.events[showIndex]))
    print "  |     |-NoPass: \t{0}".format(printDiv(pid_res.total.prelim.associated.ided.bad.noPass.events[showIndex], pid_res.total.events[showIndex], pid_res.total.prelim.events[showIndex], pid_res.total.prelim.associated.events[showIndex], pid_res.total.prelim.associated.ided.events[showIndex], pid_res.total.prelim.associated.ided.bad.events[showIndex]))
    
    print "\t\t\t{0:10}\t{1:10}\t{2:10}\t{3:10}\t{4:10}".format(0, 1, 2, 3, 4)
    print "Total: \t\t\t{0}".format(printDiv(pid_res_mu.total.events[showIndex], pid_res_mu.total.events[showIndex]))
    print "|-Prelim: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex]))
    print "  |-NoID: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.noID.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex]))
    print "  |-ManyID: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.manyID.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex]))
    print "  |-Ided: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.ided.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.ided.events[showIndex]))
    print "    |-Pass: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.ided.pass_t.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.ided.events[showIndex]))
    print "    |-NoPass: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.ided.noPass.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.ided.events[showIndex]))
    print "  |-Associated: \t{0}".format(printDiv(pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex]))
    print "  | |-NoID: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.associated.noID.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex]))
    print "  | |-ManyID: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.associated.manyID.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex]))
    print "  | |-Ided: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex]))
    print "  |   |-Good: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.good.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   | |-Pass: \t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.good.pass_t.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   | |-NoPass: \t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.good.noPass.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.good.events[showIndex]))
    print "  |   |-Bad: \t\t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.bad.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.bad.events[showIndex]))
    print "  |     |-Pass: \t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.bad.pass_t.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.bad.events[showIndex]))
    print "  |     |-NoPass: \t{0}".format(printDiv(pid_res_mu.total.prelim.associated.ided.bad.noPass.events[showIndex], pid_res_mu.total.events[showIndex], pid_res_mu.total.prelim.events[showIndex], pid_res_mu.total.prelim.associated.events[showIndex], pid_res_mu.total.prelim.associated.ided.events[showIndex], pid_res_mu.total.prelim.associated.ided.bad.events[showIndex]))

if not hasattr(f, "cutsDefinition"):
    sys.exit(0)
    
#cut = [0]
#for val in f.cutsDefinition:
#    cut.append(val.lists.unDeflectedElDist)
    
#passed = Hist1D(320, 10, 70)
#print len(pid_res.total.prelim.ided.pass_t.events)
#for i, val  in enumerate(f.cutsDefinition):
#    if i==len(pid_res.total.prelim.ided.pass_t.events):
#        continue
#    print i, val.lists.unDeflectedElDist
    #passed.Fill(val.lists.unDeflectedElDist, pid_res.total.prelim.ided.pass_t.events[i])

#c1 = Canvas()
#passed.Draw()
#wait()

#sys.exit(0)
cuts1 = [0]
cuts2 = [0]
for i, val in enumerate(f.cutsDefinition):
    if i>= cutsSize-1:
        break

    cuts1.append((val.lists.k2pi.maxPi0MassDiff-val.lists.k2pi.minPi0MassDiff)/2.)
    cuts2.append((val.lists.k2pi.maxKaonMassDiff-val.lists.k2pi.minKaonMassDiff)/2.)
    
cuts1 = list(set(cuts1))
cuts1.sort()

cuts2 = list(set(cuts2))
cuts2.sort()

lastCuts1 = cuts1[-1]
lastCuts2 = cuts2[-1]

print cuts1
print cuts2

bins1 = []
bins2 = []
for i,el in enumerate(cuts1[0:-1]):
    bins1.append(el + (cuts1[i+1]-el)/2.)
    #bins1 = [f+(s-f)/2. for f,s in zip(cuts1[0:-1], cuts1[1:])]
for i,el in enumerate(cuts2[0:-1]):
    bins2.append(el + (cuts2[i+1]-el)/2.)
#bins2 = [f+(s-f)/2. for f,s in zip(cuts2[0:-1], cuts2[1:])]
bins1.append(lastCuts1 + (lastCuts1-bins1[-1]))
bins2.append(lastCuts2 + (lastCuts2-bins2[-1]))

#[0, 0.005, 0.008, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3]
#[0.4, 0.43, 0.46, 0.47, 0.475, 0.485]

#bins1 = [0.00255, 0.0066, 0.0141, 0.0301, 0.0501, 0.0701, 0.0901, 0.1251, 0.175, 0.2251, 0.2751, 1]
#bins1 = [0, 1]
#test = Hist1D(bins1)
print bins1
print bins2
total = Hist2D(bins1, bins2)
total.SetTitle("Total")
#passed = Graph2D(ScanIndex)
passed = Hist2D(bins1, bins2)
#passed = Hist1D(bins1)
passed.SetTitle("Passed")
#ided = Graph2D(ScanIndex)
ided = Hist2D(bins1, bins2)
ided.SetTitle("ided")
#noid = Graph2D(ScanIndex)
noid = Hist2D(bins1, bins2)
noid.SetTitle("noid")
#manyid = Graph2D(ScanIndex)
manyid = Hist2D(bins1, bins2)
manyid.SetTitle("manyid")
badid = Hist2D(bins1, bins2)
badid.SetTitle("bad")
for i, val  in enumerate(f.cutsDefinition):
    if i>= cutsSize-1:
        break
    if i==0:
        continue
    pi0diff = (val.lists.k2pi.maxPi0MassDiff-val.lists.k2pi.minPi0MassDiff)/2.
    kdiff = (val.lists.k2pi.maxKaonMassDiff-val.lists.k2pi.minKaonMassDiff)/2.
    
    total_events = pid_res.total.prelim.associated.ided.events[i] + pid_res.total.prelim.associated.manyID.events[i] + pid_res.total.prelim.associated.noID.events[i] 
    total.Fill(pi0diff, kdiff, total_events)
    
    passed.Fill(pi0diff, kdiff, pid_res.total.prelim.associated.ided.events[i])
    
    manyid.Fill(pi0diff, kdiff, pid_res.total.prelim.associated.manyID.events[i])
    print pid_res.total.prelim.associated.manyID.events[i]
    noid.Fill(pi0diff, kdiff, pid_res.total.prelim.associated.noID.events[i])
    badid.Fill(pi0diff, kdiff, pid_res.total.prelim.associated.ided.bad.events[i])
    
    
#    if i>= cutsSize:
#        break
#    if i==0:
#        continue
    #if val.lists.k2pi.minKaonMassDiff==0.43 and val.lists.k2pi.maxKaonMassDiff==0.51:
    #    continue
    #print val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, val.lists.maxKaonMassDiff
    #passed.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    #passed.Fill(val.lists.k2pi.maxPi0MassDiff, val.lists.k2pi.maxKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    #passed.Fill(val.lists.maxPi0MassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    
    #if(val.lists.maxKaonMassDiff>0.05 and val.lists.maxKaonMassDiff<0.12):
    #    print val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i]
    #ided.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.ided.events[i])
    #ided.Fill(val.lists.k2pi.maxPi0MassDiff, val.lists.k2pi.maxKaonMassDiff, pid_res.total.prelim.ided.events[i])
    #noid.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.noID.events[i])
    #noid.Fill(val.lists.k2pi.maxPi0MassDiff, val.lists.k2pi.maxKaonMassDiff, pid_res.total.prelim.noID.events[i])
    #manyid.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.manyID.events[i])
    #manyid.Fill(val.lists.k2pi.maxPi0MassDiff, val.lists.k2pi.maxKaonMassDiff, pid_res.total.prelim.manyID.events[i])
#c0 = Canvas()
#total.Draw("colTEXT")
#c1 = Canvas()
#passed.Draw("colTEXT")
#c2 = Canvas()
#ided.Draw("colTEXT")
#c3 = Canvas()
#noid.Draw("colTEXT")
c3 = Canvas()
manyid.Draw("colTEXT")

badid.Divide(passed)
passed.Divide(total)
manyid.Divide(total)
noid.Divide(total)

c4 = Canvas()
manyid.Draw("colTEXT")
wait()
def minmax(hist2d):
    vals = []
    for x in range(0, hist2d.GetNbinsX()):
        for y in range(0, hist2d.GetNbinsY()):
            vals.append(hist2d.GetBinContent(hist2d.GetBin(x,y)))
    vals = [v for v in vals if v!=0]
    if len(vals)>0:
        return (min(vals), max(vals))
    else:
        return (0,0)

print "Passed MinMax=", minmax(passed)
print "ManyID MinMax=", minmax(manyid)
print "NoID MinMax=", minmax(noid)
print "BadID MinMax=", minmax(badid)

wait()