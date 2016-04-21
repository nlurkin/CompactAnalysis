#!/bin/env python

import rootpy.ROOT as ROOT 
from signal import pause

ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

ROOT.gSystem.AddDynamicPath("/afs/cern.ch/user/n/nlurkin/Compact/ROOTComp/build/Support")
ROOT.gSystem.Load("libCuts")

from rootpy.tree import Tree
from rootpy.io import root_open
from rootpy.plotting import Hist2D, Graph2D, Canvas, Hist1D, Graph1D
import sys
from rootpy.interactive import wait

ROOT.gStyle.SetOptStat(0)


if len(sys.argv)<1:
    print "Not enough arguments. Expecting: input.root [scanID]"
    sys.exit(0)

f = root_open(sys.argv[1], "open")

runmap = {}
i = 0
for evt in f.event:
    if not runmap.has_key(evt.rawBurst.nrun):
        runmap[evt.rawBurst.nrun] = [[evt.rawBurst.abcog_params.alpha], [evt.rawBurst.abcog_params.beta]]
    else:
        if runmap[evt.rawBurst.nrun][0][-1] != evt.rawBurst.abcog_params.alpha:
            runmap[evt.rawBurst.nrun][0].append(evt.rawBurst.abcog_params.alpha)
        if runmap[evt.rawBurst.nrun][1][-1] != evt.rawBurst.abcog_params.alpha:
            runmap[evt.rawBurst.nrun][1].append(evt.rawBurst.abcog_params.alpha)
    
#     i = i +1
#     if i>1e4:
#         break
    
alphaPlot = Graph1D()
alphaPlot.SetName("Alpha")
alphaPlot.SetTitle("Alpha")
betaPlot = Graph1D()
betaPlot.SetName("Beta")
betaPlot.SetTitle("Beta")

i_alpha = 0
i_beta = 0
for run in sorted(runmap.keys()):
    npoints = len(runmap[run][0])
    dist = 1./npoints
    for x in range(0, npoints):
        alphaPlot.SetPoint(i_alpha, run+(x*dist), runmap[run][0][x])
        i_alpha = i_alpha+1
    npoints = len(runmap[run][1])
    dist = 1./npoints
    for x in range(0, npoints):
        betaPlot.SetPoint(i_beta, run+(x*dist), runmap[run][1][x])
        i_beta = i_beta+1

c1 = Canvas()
alphaPlot.Draw("")
c2 = Canvas()
betaPlot.Draw("")

wait()






if not hasattr(f, "cutsDefinition"):
    sys.exit(0)
    
cut = [0]
for val in f.cutsDefinition:
    cut.append(val.lists.unDeflectedElDist)
    
passed = Hist1D(320, 10, 70)
for i, val  in enumerate(f.cutsDefinition):
    passed.Fill(val.lists.unDeflectedElDist, pid_res.total.prelim.ided.pass_t.events[i])

c1 = Canvas()
passed.Draw()
wait()

sys.exit(0)
cuts1 = [0]
cuts2 = [0]
cuts3 = [0]
for val in f.cutsDefinition:
    cuts1.append(val.lists.maxPi0MassDiff)
    cuts2.append(val.lists.maxKaonMassDiff)
    cuts3.append(val.lists.minKaonMassDiff)
    
cuts1 = list(set(cuts1))
cuts1.sort()

cuts2 = list(set(cuts2))
cuts2.sort()

cuts3 = list(set(cuts3))
cuts3.sort()

lastCuts1 = cuts1[-1]
lastCuts2 = cuts2[-1]

print cuts1
print cuts2
print cuts3

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
for i, val  in enumerate(f.cutsDefinition):
    if i>= cutsSize:
        break
#    if i==0:
#        continue
    if val.lists.minKaonMassDiff==0.43 and val.lists.maxKaonMassDiff==0.51:
        continue
    #print val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, val.lists.maxKaonMassDiff
    #passed.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    passed.Fill(val.lists.maxPi0MassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    #passed.Fill(val.lists.maxPi0MassDiff, pid_res.total.prelim.ided.pass_t.events[i])
    
    #if(val.lists.maxKaonMassDiff>0.05 and val.lists.maxKaonMassDiff<0.12):
    #    print val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.ided.pass_t.events[i]
    #ided.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.ided.events[i])
    ided.Fill(val.lists.maxPi0MassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.ided.events[i])
    #noid.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.noID.events[i])
    noid.Fill(val.lists.maxPi0MassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.noID.events[i])
    #manyid.SetPoint(i, val.lists.maxPi0MassDiff, val.lists.minKaonMassDiff, pid_res.total.prelim.manyID.events[i])
    manyid.Fill(val.lists.maxPi0MassDiff, val.lists.maxKaonMassDiff, pid_res.total.prelim.manyID.events[i])
c1 = Canvas()
passed.Draw("colTEXT")
c2 = Canvas()
ided.Draw("colTEXT")
c3 = Canvas()
noid.Draw("colTEXT")
c4 = Canvas()
#manyid.Draw("surf1")
manyid.Draw("colTEXT")
wait()