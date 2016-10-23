#!/bin/env python

from signal import pause
import sys
import argparse

import rootpy.ROOT as ROOT 

from rootpy.interactive import wait
from rootpy.io import root_open, Directory
from rootpy.plotting import Hist2D, Graph2D, Canvas, Hist1D, HistStack, Legend
from rootpy.tree import Tree, TreeModel, IntCol
from rootpy.plotting.style import set_style
from rootpy.plotting.utils import draw
from rootpy.io.pickler import dump,load

from root_numpy import tree2array
import os
import yaml

tfd = None
maxfetch = 1000000
#maxfetch = 100
limit = None

class Description(object):
    def __init__(self):
        self.path = ""
        self.shortName = ""
        self.nbins = -1
        self.xmin = -1
        self.xmax = -1
        self.ytitle = ""
        self.xtitle = ""
        self.arrowLoc = None
        self.force = False
        self.samples = ["data", "pi", "mu", "e"]
        
        self.data = None
        self.pi = None
        self.mu = None
        self.e = None
        
        self.hdata = None
        self.hpi = None
        self.hmu = None
        self.he = None
    
    def testSample(self, sample):
        return sample in self.samples
    
    def createHisto(self):
        if self.testSample("data"):
            self.data = Hist1D(self.nbins, self.xmin, self.xmax, name="{0}_{1}".format(self.shortName, "data"))
        if self.testSample("pi"):
            self.pi = Hist1D(self.nbins, self.xmin, self.xmax, name="{0}_{1}".format(self.shortName, "pi"))
        if self.testSample("mu"):
            self.mu = Hist1D(self.nbins, self.xmin, self.xmax, name="{0}_{1}".format(self.shortName, "mu"))
        if self.testSample("e"):
            self.e  = Hist1D(self.nbins, self.xmin, self.xmax, name="{0}_{1}".format(self.shortName, "e"))

    def setHeader(self, header, sample):
        setattr(self, "h{0}".format(sample), header)

    def setPlot(self, plot, sample):
        getattr(self, sample).Add(plot, 1.0)
        
    def write(self):
        global tfd
        if self.testSample("data"):
            self.data.write()
        if self.testSample("pi"):
            self.pi.write()
        if self.testSample("mu"):
            self.mu.write()
        if self.testSample("e"):
            self.e.write()
        dump([self.hdata, self.hpi, self.hmu, self.he], tfd, key="head_{0}".format(self.shortName))
    
    def read(self):
        global tdf
        if self.testSample("data"):
            self.data = tfd.Get("{0}/{0}_data".format(self.shortName))
        if self.testSample("pi"):
            self.pi = tfd.Get("{0}/{0}_pi".format(self.shortName))
        if self.testSample("mu"):
            self.mu = tfd.Get("{0}/{0}_mu".format(self.shortName))
        if self.testSample("e"):
            self.e = tfd.Get("{0}/{0}_e".format(self.shortName))
        [self.hdata, self.hpi, self.hmu, self.he] = load(tfd, key="head_{0}".format(self.shortName))
        
    def setStyle(self):
        if self.testSample("data"):
            self.data.linecolor='red'
            self.data.markercolor='red'
            self.data.title = "Data"
            self.data.markersize = 0.5
            self.data.SetStats(False)

        if self.testSample("pi"):
            self.pi.title = "Pi"
            self.pi.fillcolor = ROOT.TColor.GetColor(137,116,232)#ROOT.gStyle.GetColorPalette(180)
            self.pi.linecolor = ROOT.TColor.GetColor(137,116,232)#ROOT.gStyle.GetColorPalette(180)
            self.pi.fillstyle = "solid"
            self.pi.SetStats(False)
    
        if self.testSample("mu"):
            self.mu.title = "Mu"
            self.mu.fillcolor = ROOT.TColor.GetColor(116,232,136)#ROOT.gStyle.GetColorPalette(480)
            self.mu.linecolor = ROOT.TColor.GetColor(116,232,136)#ROOT.gStyle.GetColorPalette(480)
            self.mu.fillstyle = "solid"
            self.mu.SetStats(False)
        
        if self.testSample("e"):
            self.e.title = "e"
            self.e.fillcolor = ROOT.TColor.GetColor(232,126,116)#ROOT.gStyle.GetColorPalette(400)
            self.e.linecolor = ROOT.TColor.GetColor(232,126,116)#ROOT.gStyle.GetColorPalette(400)
            self.e.fillstyle = "solid"
            self.e.SetStats(False)

    def scaleSample(self):
        br_e = 5.07E-2
        br_mu = 3.353E-2
        br_pi = 2.066E-1
        if self.pi:
            scaleSample(self.hpi, br_pi, self.pi)
        if self.mu:
            scaleSample(self.hmu, br_mu, self.mu)
        if self.e:
            scaleSample(self.he, br_e, self.e)
    
    def scaleToData(self):
        totalMC = 0
        if self.pi:
            totalMC += self.pi.Integral()
        if self.mu:
            totalMC += self.mu.Integral()
        if self.e:
            totalMC += self.e.Integral()
            
        if self.pi:
            scaleHisto(totalMC, self.data.Integral(), self.pi)
        if self.mu:
            scaleHisto(totalMC, self.data.Integral(), self.mu)
        if self.e:
            scaleHisto(totalMC, self.data.Integral(), self.e)
    
    def makeLegend(self):
        nentry = 0
        if self.data:
            nentry += 1
        if self.pi:
            nentry += 1
        if self.mu:
            nentry += 1
        if self.e:
            nentry += 1
        legend = Legend(nentry, entryheight=0.04, entrysep=0.01, textsize=20, textfont=43, leftmargin=0.6, rightmargin=-0.07)
        if self.data:
            legend.AddEntry(self.data, "Data", style="lp")
        if self.pi:
            legend.AddEntry(self.pi, "K^{#pm}#rightarrow#pi^{#pm}#pi^{0}_{D}", style="lf")
        if self.mu:
            legend.AddEntry(self.mu, "K^{#pm}#rightarrow#pi^{0}_{D}#mu^{#pm}#nu", style="lf")
        if self.e:
            legend.AddEntry(self.e, "K^{#pm}#rightarrow#pi^{0}_{D}e^{#pm}#nu", style="lf")
        return legend

    def drawArrow(self):
        print self.arrowLoc
        if not self.arrowLoc is None:
            for arr in self.arrowLoc:
                binID = self.data.FindBin(arr)
                val1 = self.data.GetMaximum()*0.90
                factor1 = 1.2
                factor2 = 0.05
                if self.log:
                    factor1 = 1.02
                    factor2 = 0.005
                val2 = max(self.data.GetBinContent(binID)*factor1, self.data.GetMaximum()*factor2)
                a = ROOT.TArrow(arr, max(val1, val2), arr, min(val1, val2), 0.01, "|>")
                a.Draw()

class LineDescription(Description):
    def __init__(self):
        super(LineDescription, self).__init__()
        self.extract = False
        self.varName = ""
        self.cut = None
        self.usewpk = True
        self.rmin = 0.89
        self.rmax = 1.11
        self.doRatio = True
        self.log = False
        
        self.branchIndex = -1
    
    def fromLine(self, line):
        self.__dict__.update(line)
    
    def fillTestCut(self, event, sample):
        if self.cut is None or (hasattr(event[1], '__iter__') and event[1][self.cut]):
            if sample=="data":
                getattr(self, sample).Fill(event[self.branchIndex])
            else:
                getattr(self, sample).Fill(event[self.branchIndex], applyWeights(event[0], event[2], self.usewpk, event[3], event[4]))

class ScanDescription(Description):
    def __init__(self):
        super(ScanDescription, self).__init__()
        
class Header(object):
    def __init__(self, cheader=None):
        if cheader is None:
            self.NProcessedEvents = 0
            self.NFailedEvents = 0
            self.NPassedEvents = 0
        else:
            self.NProcessedEvents = cheader.NProcessedEvents
            self.NFailedEvents = cheader.NFailedEvents
            self.NPassedEvents = cheader.NPassedEvents
        
    def __add__(self, other):
        self.NProcessedEvents += other.NProcessedEvents
        self.NFailedEvents += other.NFailedEvents
        self.NPassedEvents += other.NPassedEvents
        return self
    
    def __repr__(self):
        return "{0} {1} {2}".format(self.NProcessedEvents, self.NFailedEvents, self.NPassedEvents)

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

def secondOrder_Pk(nrun, pkgen, mcPartType):
  # 4 Feb 2010 "crude" quadratic weights
  slope = 0.0
  if (nrun>=20154 and nrun<=20256):
    slope = -0.05
  elif (nrun>=20268 and nrun<=20303):
    slope = -0.08
  elif (nrun>=20304 and nrun<=20324): 
    slope = -0.02
  elif (nrun>=20332 and nrun<=20371):
    slope = +0.08
  elif (nrun>=20387 and nrun<=20486): 
    slope = +0.07
  elif (nrun>=20487 and nrun<=20521): 
    slope = +0.06
  elif (nrun>=20522 and nrun<=20612): 
    slope = +0.01
  elif (nrun>=20613 and nrun<=20695): 
    slope = +0.10

  wt_pk = 1 + slope*pow(pkgen-74.0, 2)

  # Dec 2010 fine corrections to the weights (K+)
  if (mcPartType==1):

    # period 1
    if (nrun>=20114 and nrun<=20154):
      wt_pk *= (1.0 - 0.005*(pkgen-74.0))
      if (pkgen>77.5): 
 		wt_pk *= (1.0+4.0*(pkgen-77.5))
    if (nrun>=20155 and nrun<=20174):
      wt_pk *= (1.0 + 0.065*(pkgen-74.0))
      if (pkgen>77.5): 
 		wt_pk *= (1.0+1.0*(pkgen-77.5))
    if (nrun>=20175 and nrun<=20203):
      wt_pk *= (1.0 + 0.10*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+0.5*(pkgen-77.0))
      if (pkgen>70.0 and pkgen<71.5): 
 		wt_pk *= (1.0-0.5*(71.5-pkgen))

    # periods 2+3
    if (nrun>=20209 and nrun<=20226):
      wt_pk *= (1.0 + 0.05*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+0.5*(pkgen-77.0))
      if (pkgen>70.0 and pkgen<71.5): 
 		wt_pk *= (1.0-0.5*(71.5-pkgen))
    if (nrun>=20228 and nrun<=20256):
      wt_pk *= (1.0 + 0.035*(pkgen-74.0))
      if (pkgen>76.8): 
 		wt_pk *= (1.0+0.7*(pkgen-76.8))
      if (pkgen>70.0 and pkgen<71.6): 
 		wt_pk *= (1.0-0.5*(71.6-pkgen))
    if (nrun>=20268 and nrun<=20291):
      wt_pk *= (1.0 + 0.04*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+6.0*(pkgen-77.0))
    if (nrun>=20296 and nrun<=20303):
      wt_pk *= (1.0 + 0.015*(pkgen-74.0))
      if (pkgen>76.9): 
 		wt_pk *= (1.0+5.0*(pkgen-76.9))
    if (nrun>=20304 and nrun<=20324):
      wt_pk *= (1.0 + 0.02*(pkgen-74.0))
      if (pkgen>77.2): 
 		wt_pk *= (1.0+1.3*(pkgen-77.2))
      if (pkgen>70.0 and pkgen<71.5): 
 		wt_pk *= (1.0-0.25*(71.5-pkgen))

    # period 4
    if (nrun>=20332 and nrun<=20351):
      if (pkgen>77.3): 
 		wt_pk *= (1.0+4.0*(pkgen-77.3))
      if (pkgen>70.0 and pkgen<71.0): 
 		wt_pk *= (1.0+1.0*(71.0-pkgen))
    if (nrun>=20352 and nrun<=20371):
      wt_pk *= (1.0 - 0.01*(pkgen-74.0))
      if (pkgen>77.3): 
 		wt_pk *= (1.0+6.0*(pkgen-77.3))
      if (pkgen>70.0 and pkgen<71.0): 
 		wt_pk *= (1.0+2.0*(71.0-pkgen))
    if (nrun>=20387 and nrun<=20404):
      wt_pk *= (1.0 - 0.015*(pkgen-74.0))
      if (pkgen>77.3): 
 		wt_pk *= (1.0+4.5*(pkgen-77.3))
      if (pkgen>70.0 and pkgen<71.0): 
 		wt_pk *= (1.0+4.0*(71.0-pkgen))

    # period 5
    if (nrun>=20410 and nrun<=20424):
      wt_pk *= (1.0 - 0.01*(pkgen-74.0))
      if (pkgen>77.5): 
 		wt_pk *= (1.0+12.0*(pkgen-77.5))
      if (pkgen<71.2): 
 		wt_pk *= (1.0+ 2.5*(71.2-pkgen))
    if (nrun>=20438 and nrun<=20453):
      wt_pk *= (1.0 - 0.015*(pkgen-74.0))
      if (pkgen>77.5): 
 		wt_pk *= (1.0+16.0*(pkgen-77.5))
      if (pkgen<71.1): 
 		wt_pk *= (1.0+ 5.0*(71.1-pkgen))
    if (nrun>=20459 and nrun<=20478):
      wt_pk *= (1.0 - 0.01*(pkgen-74.0))
      if (pkgen>77.4): 
 		wt_pk *= (1.0+10.0*(pkgen-77.4))
      if (pkgen<71.0): 
 		wt_pk *= (1.0+ 3.0*(71.0-pkgen))
    if (nrun>=20482 and nrun<=20485):
      wt_pk *= (1.0 + 0.00*(pkgen-74.0))
      if (pkgen>77.3): 
 		wt_pk *= (1.0+6.0*(pkgen-77.3))
      if (pkgen<71.0): 
 		wt_pk *= (1.0+5.0*(71.0-pkgen))

  # Jan-Feb 2011 fine corrections to the weights (K-) 
  else:

    # period 1
    if (nrun>=20114 and nrun<=20154):
      wt_pk *= (1.0 - 0.005*(pkgen-74.0))
      if (pkgen>77.5): 
 		wt_pk *= (1.0+4.0*(pkgen-77.5))
    if (nrun>=20155 and nrun<=20174):
      if (pkgen>77.5): 
 		wt_pk *= (1.0+1.0*(pkgen-77.5))
    if (nrun>=20175 and nrun<=20203):
      wt_pk *= (1.0 + 0.05*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+0.5*(pkgen-77.0))
      if (pkgen>70.5 and pkgen<71.5): 
 		wt_pk *= (1.0-0.5*(71.5-pkgen))

    # periods 2+3
    if (nrun>=20209 and nrun<=20226):
      wt_pk *= (1.0 + 0.03*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+0.5*(pkgen-77.0))
      if (pkgen>70.0 and pkgen<71.5): 
 		wt_pk *= (1.0-0.5*(71.5-pkgen))
    if (nrun>=20228 and nrun<=20256):
      wt_pk *= (1.0 + 0.035*(pkgen-74.0))
      if (pkgen>76.8): 
 		wt_pk *= (1.0+0.7*(pkgen-76.8))
      if (pkgen>70.0 and pkgen<71.6): 
 		wt_pk *= (1.0-0.5*(71.6-pkgen))
    if (nrun>=20268 and nrun<=20291):
      wt_pk *= (1.0 + 0.04*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+6.0*(pkgen-77.0))
    if (nrun>=20296 and nrun<=20303):
      wt_pk *= (1.0 + 0.05*(pkgen-74.0))
      if (pkgen>76.9): 
 		wt_pk *= (1.0+5.0*(pkgen-76.9))

    # period 6
    if (nrun>=20487 and nrun<=20531):
      wt_pk *= (1.0 + 0.02*(pkgen-74.0))
      if (pkgen>77.0): 
 		wt_pk *= (1.0+4.0*(pkgen-77.0))
      if (pkgen>70.0 and pkgen<71.5): 
 		wt_pk *= (1.0+0.4*(71.5-pkgen))

  if (wt_pk<0): 
 		wt_pk = 0.01

  return wt_pk


def applyWeights(run, wpk, usewpk, pkgen, mcPartType):
    if not run in fRatioMap:
        print "[ERROR] No weight found for run {0}".format(run)
        return -1

    ratio = fRatioMap[run]
    weight = fAverageRatio / ratio
    if not usewpk:
        wpk = 1
    elif usewpk==2:
        wpk = secondOrder_Pk(run, pkgen, mcPartType)
    
    return weight*wpk

ROOTWAY = False
def readData(path, fList):
    print "Reading : {0}".format(path)
    f = root_open(path, "open")

    header = Header()
    for i, h in enumerate(f.header):
        header += Header(h.header)
        
    if ROOTWAY:
        for el in fList:
            if el.cut is None:
                cutCondition = ""
            else:
                cutCondition = "cutsResult[{0}]==1".format(el.cut)

            tfd.cd(el.shortName)
            f.event.Draw("{varName}>>{0}_data({1},{2},{3})".format(el.shortName, el.nbins, el.xmin, el.xmax, varName=el.varName), cutCondition)
            el.data = ROOT.gROOT.FindObject("{0}_data".format(el.shortName))
            el.hdata = header 
    else:
        total = len(f.event)
    
        listBranch = ["rawBurst.nrun"]
        if "cutsResult" in [x.GetName() for x in f.event.iterbranches()]:
            listBranch.append("cutsResult")
        else:
            listBranch.append("rawBurst.time")
        
        atLeastOne = False
        for el in fList:
            tfd.cd(el.shortName)
            el.setHeader(header, "data")
            if not el.testSample("data"):
                continue
            atLeastOne = True
            if not el.varName in listBranch:
                el.branchIndex = len(listBranch)
                listBranch.append(el.varName)
            else:
                el.branchIndex = listBranch.index(el.varName)
        tfd.purge()        
        tfd.cd()
        if not atLeastOne:
            f.close()
            return

        i=0
        factor = 0
        while i<total:
            rec = tree2array(f.event, branches=listBranch, start=factor*maxfetch, stop=(factor+1)*maxfetch)
            factor += 1
            for event in rec:
                print "\r{0}/{1} {2}".format(i, total, i/float(total)),
                for el in fList:
                    el.fillTestCut(event, "data")
                i += 1
            if not limit is None and i>limit:
                break
    
    print
    f.close()

def printKey(name):
    global tfd
    if hasattr(tfd, name):
        print getattr(tfd, name).keys(False)

def extractPlot(path, component, fList):
    global tfd
    print "Reading: {0}".format(path)
    f = root_open(path, "open")

    for el in fList:
        if el.cut is None:
            cutCondition = 0
        else:
            cutCondition = el.cut
        header = Header()
        fs = f.get("{0}/fitStruct".format(cutCondition))
        fs.Draw("totEvents>>hTotEvents")
        temp = ROOT.gROOT.FindObject("hTotEvents")
        header.NProcessedEvents = temp.GetBinLowEdge(temp.GetMaximumBin())
        fs.Draw("selEvents>>hSelEvents")
        temp = ROOT.gROOT.FindObject("hSelEvents")
        header.NPassedEvents = temp.GetBinLowEdge(temp.GetMaximumBin())
        header.NFailedEvents = 0

        tfd.cd(el.shortName)
        plot = f.get("{0}/{1}".format(cutCondition, el.varName))
        el.setHeader(header, component)
        el.setPlot(plot, component)
    
    print
    f.close()
    
def readMC(path, fList, sample):
    atLeastOne = False
    
    for el in fList:
        if not el.testSample(sample):
            continue
        atLeastOne = True
    
    if not atLeastOne:
        return
    
    print "Reading: {0}".format(path)
    f = root_open(path, "open")
    tfd.cd()
    header = Header()
    for i, h in enumerate(f.header):
        header += Header(h.header)

    total = len(f.event)
    
    listBranch = ["rawBurst.nrun"]
    if "cutsResult" in [x.GetName() for x in f.event.iterbranches()]:
        listBranch.append("cutsResult")
    else:
        listBranch.append("rawBurst.time")
    listBranch.append("corrEvent.weight")
    listBranch.append("mc.k.P.Vect().Mag()")
    listBranch.append("mc.k.pdgID")
    
    for el in fList:
        tfd.cd(el.shortName)
        el.setHeader(header, sample)
        if not el.varName in listBranch:
            el.branchIndex = len(listBranch)
            listBranch.append(el.varName)
        else:
            el.branchIndex = listBranch.index(el.varName)
    tfd.purge()
    tfd.cd()
    i=0
    factor = 0
    while i<total:
        rec = tree2array(f.event, branches=listBranch, start=factor*maxfetch, stop=(factor+1)*maxfetch)
        factor += 1
        for event in rec:
            print "\r{0}/{1} {2}".format(i, total, i/float(total)),
            for el in fList:
                if el.testSample(sample):
                    el.fillTestCut(event, sample)
            i += 1
        if not limit is None and i>limit:
            break
    print 
    f.close()

def scaleHisto(totalMCNorm, ndata, hist):
    factor = (float(ndata)) / totalMCNorm
    hist.Scale(factor)

def scaleSample(header, BR, hist, scaleFactor=1.):
    hist.Scale(BR / (header.NProcessedEvents * scaleFactor))

def buildRatio(mc, data):
    mcSum = Hist1D(data.nbins(), data.lowerbound(), data.upperbound())
    for c in mc:
        mcSum.Add(c, 1)
    
    r = Hist1D(data.nbins(), data.lowerbound(), data.upperbound())
    mcSum.Sumw2()
    r.Sumw2()
    return r.divide(data, mcSum, 1, 1, "B")

def plotRatio(pad, stack, data, ratio, legend, plotProperties):
    pad.Divide(1, 2)
    pad.cd(1)
    pad1 = pad.GetPad(1)
    pad1.SetPad(0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.SetRightMargin(0.05)
    #pad1.SetGrid()
    pad1.cd()
    #pad1->SetLogy()
    #fStack->SetMinimum(5)
    stack.Draw("HIST")
    #fStack->GetXaxis()->SetRangeUser(0.46, 0.52)
    data.Draw("SAME E P")
    #fStack2->GetXaxis()->SetRangeUser(0.46, 0.52)
    plotProperties.drawArrow()
    legend.Draw()
    if plotProperties.log:
        pad.SetLogy()

    #Y axis mc plot settings
    stack.GetYaxis().SetTitle(plotProperties.ytitle)
    stack.GetYaxis().SetTitleSize(25)
    stack.GetYaxis().SetTitleFont(43)
    stack.GetYaxis().SetTitleOffset(1.70)
    stack.GetXaxis().SetLabelOffset(999)

    pad2 = pad.GetPad(2)
    pad2.SetPad(0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetRightMargin(0.05)
    pad2.SetGrid()
    pad2.cd()
    ratio.SetStats(0)
    ratio.Draw("ep")
    ratio.markercolor = "red"
    ratio.markersize = 0.5
    ratio.GetYaxis().SetRangeUser(plotProperties.rmin, plotProperties.rmax)
    #fSecondary->GetXaxis()->SetRangeUser(0.46, 0.52)
    
    
    # Ratio plot (ratio) settings
    ratio.SetTitle("")
    
    # Y axis ratio plot settings
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(25)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.70)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(15)
    
    # X axis ratio plot settings
    ratio.GetXaxis().SetTitle(plotProperties.xtitle)
    ratio.GetXaxis().SetTitleSize(25)
    ratio.GetXaxis().SetTitleColor(1)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(4.)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(15)

def plotSimple(pad, stack, data, legend, plotProperties):
    pad.cd()
    pad.SetRightMargin(0.05)
    if plotProperties.log:
        pad.SetLogy()
        stack.SetMinimum(2)
    stack.Draw("HIST")
    data.Draw("SAME E P")
    plotProperties.drawArrow()
    legend.Draw()

    #Y axis mc plot settings
    stack.GetYaxis().SetTitle(plotProperties.ytitle)
    stack.GetYaxis().SetTitleSize(25)
    stack.GetYaxis().SetTitleFont(43)
    stack.GetYaxis().SetTitleOffset(1.50)

    # X axis ratio plot settings
    stack.GetXaxis().SetTitle(plotProperties.xtitle)
    stack.GetXaxis().SetTitleSize(25)
    stack.GetXaxis().SetTitleColor(1)
    stack.GetXaxis().SetTitleFont(43)
    stack.GetXaxis().SetLabelFont(43)
    stack.GetXaxis().SetLabelSize(15)
 
def cleanall_el(name):
    global tfd
    print "Passing by cleanall"
    if hasattr(tfd, name):
        tdir = getattr(tfd, name) 
        for keys in tdir.keys():
            if not str(keys) == "None":
                print "Deleting {0};*".format(keys.name)
                tdir.delete("{0};*".format(keys.name))
    if name in [x.name[5:] for x in tfd.keys() if "head" in x.name]:
        print "Deleting head_{0};*".format(name)
        tfd.delete("head_{0};*".format(name))

def makeScans(path, fLists):
    global tfd
    print "Generating plots: {0}".format(",".join([x.shortName for x in fLists]))
    for el in fLists:
        if el.force:
            cleanall_el(el.shortName)
        if not hasattr(tfd, el.shortName):
            tfd.mkdir(el.shortName)
        tfd.cd(el.shortName)
        el.createHisto()
        tfd.cd()

    extractPlot("{0}/data.root".format(path), "data", fLists)
    extractPlot("{0}/pi.root".format(path), "pi", fLists)
    extractPlot("{0}/mu.root".format(path), "mu", fLists)
    extractPlot("{0}/e.root".format(path), "e", fLists)

def makePlots(path, fLists):
    global tfd
    print "Generating plots: {0}".format(",".join([x.shortName for x in fLists]))
    for el in fLists:
        if el.force:
            cleanall_el(el.shortName)
        if not hasattr(tfd, el.shortName):
            tfd.mkdir(el.shortName)
        tfd.cd(el.shortName)
        el.createHisto()
        tfd.cd()
    
    if fLists[0].extract:
        extractPlot("{0}/data.root".format(path), "data", fLists)
        extractPlot("{0}/pi.root".format(path), "pi", fLists)
        extractPlot("{0}/mu.root".format(path), "mu", fLists)
        extractPlot("{0}/e.root".format(path), "e", fLists)
    else:
        readData("{0}/data.root".format(path), fLists)
        readMC("{0}/pi.root".format(path), fLists, "pi")
        readMC("{0}/mu.root".format(path), fLists, "mu")
        readMC("{0}/e.root".format(path),  fLists, "e")

    for el in fLists:
        tfd.cd(el.shortName)
        el.write()
        tfd.cd()
    tfd.cd()
    tfd.purge()
    
def drawPlots(plotProperties):
    global tfd
    print plotProperties.__dict__
    ROOT.gStyle.SetOptStat(0)
    
    shortName = plotProperties.shortName
    
    plotProperties.read()
    plotProperties.setStyle()
    
    plotProperties.scaleSample()
    plotProperties.scaleToData()
    

    if plotProperties.doRatio:
        c = Canvas(600,700)
    else:
        c = Canvas(600,600)
    legend = plotProperties.makeLegend()
    stack = HistStack()
    
    if plotProperties.e:
        stack.Add(plotProperties.e)
    if plotProperties.mu:
        stack.Add(plotProperties.mu)
    if plotProperties.pi:
        stack.Add(plotProperties.pi)
    
    if plotProperties.doRatio:
        r = buildRatio([plotProperties.pi, plotProperties.mu, plotProperties.e], plotProperties.data)
        plotRatio(c, stack, plotProperties.data, r, legend, plotProperties)
    else:
        plotSimple(c, stack, plotProperties.data, legend, plotProperties)
        
    c.SaveAs("{0}.png".format(shortName))
    c.SaveAs("{0}.pdf".format(shortName))

def processFile(fileName, doBuild, rootFile):
    global tfd
    if rootFile is None:
        rootFile = "tmp.root"
        
    if doBuild:
        if os.path.exists(rootFile):
            tfd = root_open(rootFile, "UPDATE")
        else:
            tfd = root_open(rootFile, "CREATE")
    else:
        tfd = root_open(rootFile, "open")
    
    loadWeights("pi0dalitz_weights.dat")
    plotsInFile = [x.name[5:] for x in tfd.keys() if "head" in x.name]
    fLists = {}
    with open(fileName, "r") as fd:
        elements = yaml.load(fd)
        for line in elements:
            ld = LineDescription()
            ld.fromLine(line)
            if doBuild and (ld.shortName in plotsInFile) and not ld.force:
                print "{0} already in the file ... skipping".format(ld.shortName)
                continue
            elif not doBuild and (not ld.shortName in plotsInFile):
                print "Error: missing requested plot in file: {0}".format(ld.shortName)
            if ld.path in fLists:
                fLists[ld.path].append(ld)
            else:
                fLists[ld.path] = [ld]
    
    print "Going to process... "
    for path in fLists:
        print path + ": " + ",".join([x.shortName for x in fLists[path]])
        

    for path in fLists:
        if doBuild:
            makePlots(path, fLists[path])
        else:
            for el in fLists[path]:
                drawPlots(el)
                wait()
            
    tfd.purge()
    tfd.close()

def processScan(fileName, doBuild, rootFile):
    global tfd
    if rootFile is None:
        rootFile = "tmp_scan.root"
        
    if doBuild:
        if os.path.exists(rootFile):
            tfd = root_open(rootFile, "UPDATE")
        else:
            tfd = root_open(rootFile, "CREATE")
    else:
        tfd = root_open(rootFile, "open")
    
    loadWeights("pi0dalitz_weights.dat")
    plotsInFile = [x.name[5:] for x in tfd.keys() if "head" in x.name]
    fLists = {}
    with open(fileName, "r") as fd:
        elements = yaml.load(fd)
        for line in elements:
            ld = ScanDescription()
            ld.fromLine(line)
            if doBuild and (ld.shortName in plotsInFile) and not ld.force:
                print "{0} already in the file ... skipping".format(ld.shortName)
                continue
            elif not doBuild and (not ld.shortName in plotsInFile):
                print "Error: missing requested plot in file: {0}".format(ld.shortName)
            if ld.path in fLists:
                fLists[ld.path].append(ld)
            else:
                fLists[ld.path] = [ld]
                
    for path in fLists:
        if doBuild:
            makeScans(path, fLists[path])
            pass
        else:
            for el in fLists[path]:
                #drawPlots(el)
                pass
            
    tfd.purge()
    tfd.close()
if __name__=="__main__":
    set_style("MyStyle")
    
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    
    ratio = subparser.add_parser("ratio")
    ratio.set_defaults(func=processFile)
    ratio.add_argument("-c", "--config", required=True)
    ratio.add_argument("-b", "--build", action="store_true", default=False)
    ratio.add_argument("-r", "--root")
    scan = subparser.add_parser("scan")
    scan.set_defaults(func=processScan)
    scan.add_argument("-c", "--config", required=True)
    scan.add_argument("-b", "--build", action="store_true", default=False)
    scan.add_argument("-r", "--root")
    
    args = parser.parse_args()
    print args

    ROOT.gROOT.SetBatch(True)
    if args.build:
        ROOT.gSystem.AddDynamicPath("../../obj")
        ROOT.gSystem.Load("libexportClasses")

    args.func(args.config, args.build, args.root)
