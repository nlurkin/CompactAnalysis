#!/bin/env python

import rootpy.ROOT as ROOT 
from signal import pause

ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

ROOT.gSystem.AddDynamicPath("/afs/cern.ch/user/n/nlurkin/Compact/ROOTComp/build/Support")
ROOT.gSystem.Load("libCuts")

def propagateAfter(zplane, pt, rawEvent):
    t = rawEvent.track[pt.trackID];
    return t.aDetPos + t.aMomentum*((zplane-t.aDetPos.Z())/t.aMomentum.Z());

def pi0d_L3Trigger(plane, t, rawEvent, pbWall):
    propPos = propagateAfter(plane, t, rawEvent)
    lkrAcceptance = t.lkr_acc
    goodAcceptance = False
    if lkrAcceptance==0 and t.p>=5.5 and t.E/t.p>0.8 and rawEvent.track[t.trackID].dDeadCell>2:
        goodAcceptance=True

    goodPBWall = True
    if pbWall:
        if propPos.Y()>-33.575 and propPos.Y() < -11.850:
            goodPBWall = False

    if goodAcceptance and goodPBWall:
        return True
    return False

from rootpy.tree import Tree
from rootpy.io import root_open
from rootpy.plotting import Hist2D, Graph2D, Canvas, Hist1D, Graph1D
import sys
from rootpy.interactive import wait

ROOT.gStyle.SetOptStat(0)

def min(a,b):
    if a<=b:
        return a
    else:
        return b
    
def max(a,b):
    if a<=b:
        return b
    else:
        return a
    
if len(sys.argv)<3:
    print "Not enough arguments. Expecting: input.root output.dat"
    sys.exit(0)

f = root_open(sys.argv[1], "open")

for h in f.header:
    lkrz = h.geom.Lkr.z
    break

listSort = {}
total = len(f.event)
xDistrib_tot = Hist1D(100, 0, 1)
for i, evt in enumerate(f.event):
    if (i % 1000) ==0:
        print "{0}/{1}\r".format(i, total),
        sys.stdout.flush()
    ELKr_ep = 0
    ELKr_em = 0
    t_ep = evt.corrEvent.pTrack[evt.pi0dEvent.ep.parentTrack];
    t_em = evt.corrEvent.pTrack[evt.pi0dEvent.em.parentTrack];
    
    propPos = propagateAfter(lkrz, t_ep, evt.rawEvent)
    
    goodPBWall = True
    if evt.rawBurst.pbWall and (propPos.Y()>-33.575 and propPos.Y() < -11.850):
        goodPBWall = False
    if t_ep.lkr_acc==0 and goodPBWall and evt.rawEvent.track[t_ep.trackID].dDeadCell>2.:
        ELKr_ep = t_ep.p;

    propPos = propagateAfter(lkrz, t_em, evt.rawEvent)
    goodPBWall = True
    if evt.rawBurst.pbWall and (propPos.Y()>-33.575 and propPos.Y() < -11.850):
        goodPBWall = False
    if t_em.lkr_acc==0 and goodPBWall and evt.rawEvent.track[t_em.trackID].dDeadCell>2.:
        ELKr_em = t_em.p
 
    E_lkr = ELKr_ep + ELKr_em + evt.pi0dEvent.gamma.P.E()
    
    #E_lkr = min(t_ep.p, t_em.p)
    #E_lkr = max(t_ep.E/t_ep.p, t_em.E/t_ep.p)
    
    if not E_lkr in listSort:
        listSort[E_lkr] = []
    listSort[E_lkr].append({"run":evt.rawBurst.nrun, "time":evt.rawBurst.time, "timestamp":evt.rawEvent.timeStamp, "x":evt.pi0dEvent.x})
    xDistrib_tot.Fill(evt.pi0dEvent.x)
    

n_rejected = 0
frac_to_reject = 0.03
nreject = frac_to_reject*total

xDistrib = Hist1D(100, 0, 1)
with open(sys.argv[2], "write") as fd:
    for E in sorted(listSort):
        for evt in listSort[E]:
            fd.write("{run} {time} {timestamp}\n".format(**evt))
            n_rejected += 1
            xDistrib.Fill(evt["x"])
    
        if n_rejected>=nreject:
            break

xDistrib_tot.Scale((1/xDistrib_tot.integral()) * xDistrib.integral())
xDistrib_tot.SetMarkerColor("red")
xDistrib_tot.SetMarkerSize(0.5)
xDistrib.SetMarkerSize(0.5)
xDistrib.Draw()
xDistrib_tot.Draw("SAME")

wait()
