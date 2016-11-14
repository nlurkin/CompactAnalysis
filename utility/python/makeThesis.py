#!/bin/env python

from plotLib import loadWeights, applyWeights, plotRatio, plotSimple
from scanLib import extractScan
from pyconfig import LineDescription, ScanDescription, Header
import pyconfig

#from signal import pause
#import sys
import argparse

import rootpy.ROOT as ROOT 

from rootpy.interactive import wait
from rootpy.io import root_open, Directory
#from rootpy.plotting import Hist2D, Graph2D, Canvas, Hist1D, HistStack, Legend
from rootpy.plotting import Canvas, HistStack, Hist1D
#from rootpy.tree import Tree, TreeModel, IntCol
from rootpy.plotting.style import set_style
#from rootpy.plotting.utils import draw

#from root_numpy import tree2array
import os
import yaml

def buildRatio(mc, data):
    mcSum = Hist1D(data.nbins(), data.lowerbound(), data.upperbound())
    for c in mc:
        mcSum.Add(c, 1)
    
    r = Hist1D(data.nbins(), data.lowerbound(), data.upperbound())
    mcSum.Sumw2()
    r.Sumw2()
    return r.divide(data, mcSum, 1, 1, "B")

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
    print "Generating plots: {0}".format(",".join([x.shortName for x in fLists]))
    for el in fLists:
        if el.force:
            cleanall_el(el.shortName)
        if not hasattr(pyconfig.tfd, el.shortName):
            pyconfig.tfd.mkdir(el.shortName)
        pyconfig.tfd.cd(el.shortName)
        pyconfig.tfd.cd()

    extractScan("{0}/data.root".format(path), "data", fLists)
    extractScan("{0}/pi.root".format(path), "pi", fLists)
    extractScan("{0}/mu.root".format(path), "mu", fLists)
    extractScan("{0}/e.root".format(path), "e", fLists)
    
    for el in fLists:
        pyconfig.tfd.cd(el.shortName)
        el.write()
        pyconfig.tfd.cd()
    
    pyconfig.tfd.cd()
    pyconfig.tfd.purge()
    

def makePlots(path, fLists):
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
    ROOT.gStyle.SetOptStat(0)
    
    shortName = plotProperties.shortName
    plotProperties.read(tfd)
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

def drawScans(plotProperties):
    ROOT.gStyle.SetOptStat(0)
    
    shortName = plotProperties.shortName
    
    plotProperties.read()
    plotProperties.setStyle()
    
    plotProperties.scaleSample()
    plotProperties.scaleToData()
    

    c = Canvas(600,400)
    
    
    
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
    if rootFile is None:
        rootFile = "tmp_scan.root"
        
    if doBuild:
        if os.path.exists(rootFile):
            pyconfig.tfd = root_open(rootFile, "UPDATE")
        else:
            pyconfig.tfd = root_open(rootFile, "CREATE")
    else:
        pyconfig.tfd = root_open(rootFile, "open")
    
    loadWeights("pi0dalitz_weights.dat")
    plotsInFile = [x.name[5:] for x in pyconfig.tfd.keys() if "head" in x.name]
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
                drawScans(el)
                pass
            
    pyconfig.tfd.purge()
    pyconfig.tfd.close()
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
