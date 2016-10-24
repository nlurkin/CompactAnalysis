from rootpy.io import root_open
from rootpy.io import Directory
from rootpy.plotting import Graph
import pyconfig

import rootpy.ROOT as ROOT 


def extractScan(path, component, fList):
    print "Reading: {0}".format(path)
    f = root_open(path, "open")
    
    varName = "sig" if component=="data" else "dNew"

    for el in fList:
        header = pyconfig.Header()
        fs = f.get("0/fitStruct")
        fs.Draw("totEvents>>hTotEvents")
        temp = ROOT.gROOT.FindObject("hTotEvents")
        header.NProcessedEvents = temp.GetBinLowEdge(temp.GetMaximumBin())

        pyconfig.tfd.cd(el.shortName)
        for dir in [x.name for x in f.objects(Directory)]:
            plot = f.get("{0}/{1}".format(dir, varName))
            el.setHeader(header, component)
            el.setPlot(plot, component)
    
    print
    f.close()

def generateResult(fList):
    #Creation
    scanWErr = Graph()
    scanwUncorr = Graph()
    scanDefault = Graph()
    scanWErr.SetName("ScanWErr")
    scanwUncorr.SetName("scanwUncorr")
    scanDefault.SetName("scanDefault")

    #Points
    for i in fList.genCuts():
        scanWErr.SetPoint(i, fScanValues[i], fResultValues[i])
        scanwUncorr.SetPoint(i, fScanValues[i], fResultValues[i])
        scanWErr.SetPointError(i, 0, fResultErrors[i])
        scanwUncorr.SetPointError(i, 0, fUncorrErrors[i])

    if (*min_element(fScanValues.begin(), fScanValues.end()) == 0) {
        int size = fScanValues.size()
        scanWErr.SetPoint(size, -0.5, 0)
    }

    scanDefault.SetPoint(fDefaultCutValue, fScanValues[fDefaultCutValue],
            0.04)
    scanDefault.SetPointError(fDefaultCutValue, 0, 10)

    //Style
    scanWErr.SetMarkerStyle(23)
    scanWErr.SetMarkerColor(4)
    scanWErr.SetLineColor(4)
    scanWErr.SetFillStyle(0)

    scanwUncorr.SetMarkerStyle(0)
    scanwUncorr.SetMarkerColor(4)
    scanwUncorr.SetLineColor(6)
    scanwUncorr.SetFillStyle(0)

    scanDefault.SetMarkerStyle(0)
    scanDefault.SetLineStyle(7)
    scanDefault.SetLineWidth(2)
    scanDefault.SetLineColor(8)
    scanDefault.SetFillStyle(0)

    THStack* mcStack = getStackFromFile("mcStack")
    THStack* dataStack = getStackFromFile("sigStack")
    TLegend* legend = getLegendFromFile()

    //Plotting
    if (mcStack || dataStack) {
        pad.Divide(1, 2)
        pad.cd(1)
        TPad* pad1 = static_cast<TPad*>(pad.GetPad(1))
        pad1.SetPad(0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(3)
        pad1.SetGrid()
        pad1.cd()               // pad1 becomes the current pad

        if (mcStack)
            mcStack.Draw("HIST")
        if (dataStack)
            dataStack.Draw("SAME E P")
        if (legend)
            legend.Draw()
        if (mcStack)
            mcStack.GetXaxis().SetRangeUser(scanWErr.GetXaxis().GetXmin(),
                    scanWErr.GetXaxis().GetXmax())

        TPad* pad2 = static_cast<TPad*>(pad.GetPad(2))
        pad2.SetPad(0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.2)
        pad2.SetBottomMargin(0.2)
        pad2.SetGrid() // vertical grid
        pad2.cd()
    } else {
        pad.SetTopMargin(0.2)
        pad.SetBottomMargin(0.2)
        pad.SetGrid()
        pad.cd()
    }

    scanWErr.SetTitle("FF Slope fit result")
    scanWErr.GetXaxis().SetTitle("Cut value")
    scanWErr.GetYaxis().SetTitle("FF Slope a")
    scanWErr.GetYaxis().SetTitleOffset(2.2)

    scanWErr.SetTitle("")
    scanWErr.SetLineWidth(2)
    scanWErr.GetYaxis().SetNdivisions(505)
    scanWErr.GetYaxis().SetTitleSize(15)
    scanWErr.GetYaxis().SetTitleFont(43)
    scanWErr.GetYaxis().SetTitleOffset(1.4)
    scanWErr.GetYaxis().SetLabelFont(43) // Absolute font size in pixel (precision 3)
    scanWErr.GetYaxis().SetLabelSize(15)
    scanWErr.GetXaxis().SetLabelFont(43)
    scanWErr.GetXaxis().SetLabelSize(15)
    scanWErr.GetYaxis().SetRangeUser(0.0249, 0.051)
    scanWErr.Draw("AP")

    scanwUncorr.SetTitle("")
    scanwUncorr.SetLineWidth(3)
    scanwUncorr.GetYaxis().SetRangeUser(0.0249, 0.051)
    scanwUncorr.Draw("PSAME")
    scanDefault.Draw("PSAME")