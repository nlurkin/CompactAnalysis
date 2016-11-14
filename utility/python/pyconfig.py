from rootpy.io.pickler import dump,load
from rootpy.plotting import Legend
import rootpy.ROOT as ROOT

def scaleHisto(totalMCNorm, ndata, hist):
    factor = (float(ndata)) / totalMCNorm
    hist.Scale(factor)

def scaleSample(header, BR, hist, scaleFactor=1.):
    hist.Scale(BR / (header.NProcessedEvents * scaleFactor))

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
        setattr(self, sample, plot.clone("{0}_{1}".format(plot.name, sample)))
        #getattr(self, sample).Add(plot, 1.0)
        
    def write(self, tfd):
        if self.testSample("data"):
            self.data.write()
        if self.testSample("pi"):
            self.pi.write()
        if self.testSample("mu"):
            self.mu.write()
        if self.testSample("e"):
            self.e.write()
        dump([self.hdata, self.hpi, self.hmu, self.he], tfd, key="head_{0}".format(self.shortName))
    
    def read(self, tfd):
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
        
    def fromLine(self, line):
        self.__dict__.update(line)
        
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

