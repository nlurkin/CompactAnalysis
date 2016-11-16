import os

maxfetch = 1000000
#maxfetch = 100
limit = None

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

def plotRatio(pad, stack, data, ratio, legend, plotProperties):
    pad.Divide(1, 2)
    pad.cd(1)
    pad1 = pad.GetPad(1)
    pad1.SetPad(0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.SetRightMargin(0.05)
    pad1.SetLeftMargin(0.16)
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
    stack.GetYaxis().SetTitleOffset(2.2)
    stack.GetYaxis().SetLabelFont(43)
    stack.GetYaxis().SetLabelSize(20)
    stack.GetXaxis().SetLabelOffset(999)
    stack.GetXaxis().SetLabelSize(0)

    pad2 = pad.GetPad(2)
    pad2.SetPad(0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetRightMargin(0.05)
    pad2.SetLeftMargin(0.16)
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
    ratio.GetYaxis().SetTitleOffset(2.2)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(20)
    
    # X axis ratio plot settings
    ratio.GetXaxis().SetTitle(plotProperties.xtitle)
    ratio.GetXaxis().SetTitleSize(25)
    ratio.GetXaxis().SetTitleColor(1)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(4.)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(20)

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
    stack.GetYaxis().SetLabelFont(43)
    stack.GetYaxis().SetLabelSize(20)


    # X axis ratio plot settings
    stack.GetXaxis().SetTitle(plotProperties.xtitle)
    stack.GetXaxis().SetTitleSize(25)
    stack.GetXaxis().SetTitleColor(1)
    stack.GetXaxis().SetTitleFont(43)
    stack.GetXaxis().SetLabelFont(43)
    stack.GetXaxis().SetLabelSize(20)
