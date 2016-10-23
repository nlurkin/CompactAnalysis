#!/bin/env python

import rootpy.ROOT as ROOT 
from rootpy.interactive import wait
from rootpy.io import root_open, Directory
from rootpy.plotting import Hist1D, Hist2D, Graph2D, Canvas, Hist1D, HistStack, Legend
from rootpy.tree import Tree, TreeModel, IntCol
from rootpy.plotting.style import set_style
from rootpy.plotting.utils import draw
from rootpy.io.pickler import dump,load
from libPyROOT import gROOT

fdtmp = None
hdflt = None
hdflt_pi = None

def prepare_dflt():
    global fdtmp, hdflt, hdflt_pi
    fdflt = root_open("/afs/cern.ch/user/n/nlurkin/work/MinCompactExtRange/Dflt/data.root", "open")
    fdtmp.cd()
    fdflt.event.Draw("pi0dEvent.x>>hdflt(100, 0, 1)")
    fdflt.Close()
    hdflt = Hist1D(ROOT.gROOT.FindObject("hdflt"))
    fdflt = root_open("/afs/cern.ch/user/n/nlurkin/work/MinCompactExtRange/Dflt/pi.root", "open")
    fdtmp.cd()
    fdflt.event.Draw("pi0dEvent.x>>hdflt_pi(100, 0, 1)")
    fdflt.Close()
    hdflt_pi = Hist1D(ROOT.gROOT.FindObject("hdflt_pi"))

def do_it(file, output, alsoPi):
    global fdtmp, hdflt, hdflt_pi
    fextra = root_open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "open")
    fdtmp.cd()
    fextra.event.Draw("pi0dEvent.x>>hextra(100, 0, 1)")
    fextra.Close()

    hextra =  ROOT.gROOT.FindObject("hextra")
    n = Hist1D(100, 0, 1, name="diff")
    n.Add(hextra, 1)
    n.Add(hdflt, -1)

    n_pi = None
    if alsoPi:
        fextra = root_open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/pi.root", "open")
        textra = fextra.Get("event")
        fdtmp.cd()
        textra.Draw("pi0dEvent.x>>hextra_pi(100, 0, 1)")
        fextra.Close()

        hextra_pi = ROOT.gROOT.FindObject("hextra_pi")
        n_pi = Hist1D(100, 0, 1, name="diff_pi")
        n_pi.Add(hextra_pi, 1)
        n_pi.Add(hdflt_pi, -1)

    c = Canvas(600, 600)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.15)
    n.markercolor = 'red'
    n.linecolor = 'red'
    n.markerstyle = 8
    n.markersize = 0.8
    n.fillcolor = 'red'
    n.fillstyle = 'solid'
    n.SetStats(False)
    if alsoPi:
        n_pi.Scale(n.Integral()/n_pi.Integral())
        n_pi.fillcolor = ROOT.TColor.GetColor(137,116,232)
        n_pi.linecolor = ROOT.TColor.GetColor(137,116,232)
        n_pi.SetStats(False)

        l = Legend(2, entryheight=0.04, entrysep=0.01, textsize=20, textfont=43, leftmargin=0.6, rightmargin=0.00)
        l.AddEntry(n, "Data", style="lp")
        l.AddEntry(n_pi, "K^{#pm}#rightarrow#pi^{#pm}#pi^{0}_{D}", style="lf")
        
        c.SetRightMargin(0.05)
        n.Draw("HIST P")
        n_pi.Draw("HIST SAME")
        l.Draw()
        f = n
    else:
        n.Draw("HIST")
        f = n

    f.GetYaxis().SetTitle("Events/(0.01)")
    f.GetYaxis().SetTitleSize(25)
    f.GetYaxis().SetTitleFont(43)
    f.GetYaxis().SetTitleOffset(1.70)

    # X axis ratio plot settings
    f.GetXaxis().SetTitle("x")
    f.GetXaxis().SetTitleSize(25)
    f.GetXaxis().SetTitleColor(1)
    f.GetXaxis().SetTitleFont(43)
    f.GetXaxis().SetLabelFont(43)
    f.GetXaxis().SetLabelSize(15)

    c.SaveAs(output + ".png")
    c.SaveAs(output + ".pdf")

def do_time_track(file, name):
    fdtmp.cd()
    htime = Hist1D(160, -80, 80, name="htime")

    ftime = root_open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "open")
    ttime = ftime.Get("event")
    fdtmp.cd()
    ttime.Draw("rawEvent.track.time-rawBurst.tOffst.Dch>>htimetmp(160, -80, 80)")
    ftime.Close()

    htimetmp = ROOT.gROOT.FindObject("htimetmp")
    htime.Add(htimetmp, 1)

    c = Canvas(600, 600)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.15)
    c.SetLogy(True)
    htime.SetTitle("")
    htime.markercolor ='red'
    htime.linecolor = 'red'
    htime.markerstyle = 8
    htime.markersize = 0.8
    htime.fillcolor = 'red'
    htime.fillstyle = 'solid'
    htime.SetStats(False)
    htime.Draw("HIST")

    htime.GetYaxis().SetTitle("Events/(0.01)")
    htime.GetYaxis().SetTitleSize(25)
    htime.GetYaxis().SetTitleFont(43)
    htime.GetYaxis().SetTitleOffset(1.70)
    htime.GetYaxis().SetRangeUser(10, 5e5)

    # X axis ratio plot settings
    htime.GetXaxis().SetTitle("t_{track}-TOffset_{DCH}")
    htime.GetXaxis().SetTitleSize(25)
    htime.GetXaxis().SetTitleColor(1)
    htime.GetXaxis().SetTitleFont(43)
    htime.GetXaxis().SetLabelFont(43)
    htime.GetXaxis().SetLabelSize(15)

    c.SaveAs(name + ".png")
    c.SaveAs(name + ".pdf")

def do_time_vtx(file, name):
    ftime = root_open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "open")
    fdtmp.cd()
    ftime.event.Draw("rawEvent.vtx.time-rawBurst.tOffst.Dch>>htime(160, -80, 80)")
    ftime.Close()

    htime = ROOT.gROOT.FindObject("htime")

    c = Canvas(600, 600)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.15)
    c.SetLogy(True)
    htime.SetTitle("")
    htime.markercolor = 'red'
    htime.linecolor = 'red'
    htime.markerstyle = 8
    htime.markersize = 0.8
    htime.fillcolor = 'red'
    htime.fillstyle = 'solid'
    htime.SetStats(False)
    htime.Draw("HIST")
    a1 = ROOT.TArrow(25, 15000, 25, 1000, 0.01, "|>")
    a1.Draw()
    a2 = ROOT.TArrow(-25, 15000, -25, 1000, 0.01, "|>")
    a2.Draw()

    htime.GetYaxis().SetTitle("Events/(0.01)")
    htime.GetYaxis().SetTitleSize(25)
    htime.GetYaxis().SetTitleFont(43)
    htime.GetYaxis().SetTitleOffset(1.70)
    htime.GetYaxis().SetRangeUser(0.5, 1e5)

    # X axis ratio plot settings
    htime.GetXaxis().SetTitle("t_{vtx} - TOffset_{DCH}")
    htime.GetXaxis().SetTitleSize(25)
    htime.GetXaxis().SetTitleColor(1)
    htime.GetXaxis().SetTitleFont(43)
    htime.GetXaxis().SetLabelFont(43)
    htime.GetXaxis().SetLabelSize(15)

    c.SaveAs(name + ".png")
    c.SaveAs(name + ".pdf")

if __name__=="__main__":
    set_style("MyStyle")
    
    ROOT.gROOT.SetBatch(True)
    fdtmp = root_open("tmp_diff.root", "RECREATE")
    prepare_dflt()

    #do_it("Stage2_MinCompactExtRange_correct/extratrack", "x_extra_tracks", False)
    #do_it("Stage2_MinCompactExtRange_correct/vtx_charge", "x_vertex_qvtx", False)
    #do_it("Stage2_MinCompactExtRange_correct/track_time", "x_out_of_time", False)
    #do_it("Stage2_MinCompactExtRange/rndpid", "two_hypo", True)
    #do_time_track("Stage2_MinCompactExtRange_correct/track_time", "track_time")
    do_time_vtx("Stage2_MinCompactExtRange_correct/track_time", "vtx_time")

