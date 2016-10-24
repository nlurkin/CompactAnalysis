#!/bin/env python

import rootpy.ROOT as ROOT 
from signal import pause

ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

ROOT.gSystem.AddDynamicPath("/afs/cern.ch/user/n/nlurkin/Compact/ROOTComp/build/Support")
ROOT.gSystem.Load("libCuts")

from rootpy.tree import Tree
from rootpy.io import root_open
from rootpy.plotting import Hist1D, Efficiency, Canvas
import sys
from rootpy.interactive import wait
from root_numpy import tree2array
import os

fRatioMap = {}
fAverageRatio = None

def loadWeights(fileName):
    global fRatioMap
    global fAverageRatio
    
    print "Load Weights"

    if not os.path.exists(fileName):
        print "[ERROR] Weight file {0} does not exist".format(fileName)
        return false

    with open(fileName, "r") as fd:
        burst =-1
        sumRatio = 0
        number = 0

        for line in fd:
            fields = line.split()

            sumRatio += float(fields[1])
            number += 1
            fRatioMap[int(fields[0])] = float(fields[1])

    fAverageRatio = sumRatio / float(number)
    print "Average {0} {1}".format(sumRatio, fAverageRatio)
    return True

def applyWeights(run, wpk, usewpk):
    if not run in fRatioMap:
        print "[ERROR] No weight found for run {0}".format(run)
        return -1

    ratio = fRatioMap[run]
    weight = fAverageRatio / ratio
    if not usewpk:
        wpk = 1
    return weight*wpk


ROOT.gStyle.SetOptStat(0)

if len(sys.argv)<3:
    print "Not enough arguments. Expecting: input.root [scanID]"
    sys.exit(0)

loadWeights("pi0dalitz_weights.dat")

ftmp = root_open(sys.argv[2], "recreate")
fullHist = Hist1D(1000, 0, 1, name="FullHisto")
finalHist = Hist1D(1000, 0, 1, name="FinalHisto")

if len(sys.argv)==3:
    with open(sys.argv[1], "r") as fd:
        for line in fd:
            print "Reading {0}".format(line.strip())
            f = root_open("xroot://eosna62.cern.ch//{0}".format(line.strip()), "open")
            
            rec = tree2array(f.event, branches=["mc.xTrue", "rawBurst.nrun", "corrEvent.weight"])
            total = len(rec)
            for i, evt in enumerate(rec):
                print "\r{0}/{1} {2}".format(i, total, i/float(total)),
                fullHist.Fill(evt[0], applyWeights(evt[1], evt[2], False))
            print 
    
    f.close()
    
    ftmp.cd()
    f = root_open("/afs/cern.ch/user/n/nlurkin/work/MinCompactExtRange/Dflt/pi.root", "open")
    rec = tree2array(f.event, branches=["mc.xTrue", "rawBurst.nrun", "corrEvent.weight"])
    total = len(rec)
    for i, evt in enumerate(rec):
        print "\r{0}/{1} {2}".format(i, total, i/float(total)),
        finalHist.Fill(evt[0], applyWeights(evt[1], evt[2], False))
    
    f.close()
    ftmp.cd()
    fullHist.write()
    finalHist.write()
else:
    with open(sys.argv[1], "r") as fd:
        for line in fd:
            print "Reading {0}".format(line.strip())
            f = root_open(line.strip(), "open")
            
            fullHist.add(f.FullHisto, 1.)
    
    finalHist.add(f.FinalHisto, 1.)
    
    f.close()
    ftmp.cd()
    c1 = Canvas()
    c1.cd()
    fullHist.Draw()
    c = Canvas()
    c.cd()
    finalHist.draw()
    c = Canvas()
    accHist = Hist1D(1000, 0, 1, name="accHisto")
    for i in range(0,1001):
        #print finalHist.GetBinCenter(i), finalHist.GetBinContent(i), fullHist.GetBinContent(i)
        if fullHist.GetBinContent(i)>0:
            accHist.Fill(finalHist.GetBinCenter(i), (finalHist.GetBinContent(i)/fullHist.GetBinContent(i))*100.)
    c.cd()
    #teff = Efficiency(finalHist, fullHist)
    #teff.Draw()
    accHist.draw("HIST")
    wait()
