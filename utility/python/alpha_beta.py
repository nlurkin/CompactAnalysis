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
for i, evt in enumerate(f.event):
    #if i>20000:
    #    break
    if not runmap.has_key(evt.rawBurst.nrun):
        runmap[evt.rawBurst.nrun] = [[evt.rawBurst.abcog_params.alpha], [evt.rawBurst.abcog_params.beta], [evt.rawBurst.alpha]]
    else:
        if runmap[evt.rawBurst.nrun][0][-1] != evt.rawBurst.abcog_params.alpha:
            runmap[evt.rawBurst.nrun][0].append(evt.rawBurst.abcog_params.alpha)
        if runmap[evt.rawBurst.nrun][1][-1] != evt.rawBurst.abcog_params.alpha:
            runmap[evt.rawBurst.nrun][1].append(evt.rawBurst.abcog_params.alpha)
        if runmap[evt.rawBurst.nrun][1][-1] != evt.rawBurst.alpha:
            runmap[evt.rawBurst.nrun][1].append(evt.rawBurst.alpha)
    
#     i = i +1
#     if i>1e4:
#         break
    
alphaPlot = Graph1D()
alphaPlot.SetName("Alpha")
alphaPlot.SetTitle("")
alphaPlot.markercolor = ROOT.gStyle.GetColorPalette(180)
alphaPlot.linecolor = 0
alphaPlot.markerstyle = 7
betaPlot = Graph1D()
betaPlot.SetName("Beta")
betaPlot.SetTitle("")
betaPlot.markercolor = ROOT.gStyle.GetColorPalette(180)
betaPlot.linecolor = 0
betaPlot.markerstyle = 7
lambdaPlot = Graph1D()
lambdaPlot.SetName("Lambda")
lambdaPlot.SetTitle("")
lambdaPlot.markercolor = ROOT.gStyle.GetColorPalette(180)
lambdaPlot.linecolor = 0
lambdaPlot.markerstyle = 7

i_alpha = 0
i_beta = 0
i_lambda = 0
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
    npoints = len(runmap[run][2])
    dist = 1./npoints
    for x in range(0, npoints):
        lambdaPlot.SetPoint(i_lambda, run+(x*dist), runmap[run][2][x])
        i_lambda = i_lambda+1

c1 = Canvas()
alphaPlot.Draw("")
c2 = Canvas()
betaPlot.Draw("")
c3 = Canvas()
lambdaPlot.Draw("")

alphaPlot.GetXaxis().SetTitle("run #")
alphaPlot.GetXaxis().SetTitleColor(1)
alphaPlot.GetYaxis().SetTitle("#alpha factor")
alphaPlot.GetYaxis().SetTitleColor(1)
betaPlot.GetXaxis().SetTitle("run #")
betaPlot.GetXaxis().SetTitleColor(1)
betaPlot.GetYaxis().SetTitle("#beta factor")
betaPlot.GetYaxis().SetTitleColor(1)
lambdaPlot.GetXaxis().SetTitle("run #")
lambdaPlot.GetXaxis().SetTitleColor(1)
lambdaPlot.GetYaxis().SetTitle("#lambda_{r}")
lambdaPlot.GetYaxis().SetTitleColor(1)
wait()
