#!/bin/env python

import rootpy.ROOT as ROOT 
from signal import pause
import sys

ROOT.gSystem.AddDynamicPath("../../obj")
ROOT.gSystem.Load("libexportClasses")

ROOT.gSystem.AddDynamicPath("/afs/cern.ch/user/n/nlurkin/Compact/ROOTComp/build/Support")
ROOT.gSystem.Load("libCuts")

from rootpy.tree import Tree
from rootpy.io import root_open

def checkColor(val, listChanged):
    if val in listChanged:
        return "\033[31m"
    else:
        return ""

def printCuts(left,right, index, changed):
    print "----------------------  Cut index {0} ------------------------------".format(index)
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("varName", "Def. Value", "Cut Value", color=checkColor("varName", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("triggerMask", left.triggerMask, right.triggerMask, color=checkColor("triggerMask", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("numVertex3", left.numVertex3, right.numVertex3, color=checkColor("numVertex3", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("minZVertex", left.minZVertex, right.minZVertex, color=checkColor("minZVertex", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("maxZVertex", left.maxZVertex, right.maxZVertex, color=checkColor("maxZVertex", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("maxChi2Vertex", left.maxChi2Vertex, right.maxChi2Vertex, color=checkColor("maxChi2Vertex", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("maxExtraTracks", left.maxExtraTracks, right.maxExtraTracks, color=checkColor("maxExtraTracks", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("maxTrackTime", left.maxTrackTime, right.maxTrackTime, color=checkColor("maxTrackTime", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("boolBadTrack", left.boolBadTrack, right.boolBadTrack, color=checkColor("boolBadTrack", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("numBadTrackCombi", left.numBadTrackCombi, right.numBadTrackCombi, color=checkColor("numBadTrackCombi", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("numXCandidates", left.numXCandidates, right.numXCandidates, color=checkColor("numXCandidates", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("boolBadECandidates", left.boolBadECandidates, right.boolBadECandidates, color=checkColor("boolBadECandidates", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("minTrackMomentum", left.minTrackMomentum, right.minTrackMomentum, color=checkColor("minTrackMomentum", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("maxTrackMomentum", left.maxTrackMomentum, right.maxTrackMomentum, color=checkColor("maxTrackMomentum", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("numAddGoodCluster", left.numAddGoodCluster, right.numAddGoodCluster, color=checkColor("numAddGoodCluster", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("lkrAcceptance", left.lkrAcceptance, right.lkrAcceptance, color=checkColor("lkrAcceptance", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("minGammaEnergy", left.minGammaEnergy, right.minGammaEnergy, color=checkColor("minGammaEnergy", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("minDeadCellDist", left.minDeadCellDist, right.minDeadCellDist, color=checkColor("minDeadCellDist", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("minGammaDCHRadius", left.minGammaDCHRadius, right.minGammaDCHRadius, color=checkColor("minGammaDCHRadius", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("deflectedElDist", left.deflectedElDist, right.deflectedElDist, color=checkColor("deflectedElDist", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("unDeflectedElDist", left.unDeflectedElDist, right.unDeflectedElDist, color=checkColor("unDeflectedElDist", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.minTotalMomentum", left.k2pi.minTotalMomentum, right.k2pi.minTotalMomentum, color=checkColor("k2pi.minTotalMomentum", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.maxTotalMomentum", left.k2pi.maxTotalMomentum, right.k2pi.maxTotalMomentum, color=checkColor("k2pi.maxTotalMomentum", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.maxPt", left.k2pi.maxPt, right.k2pi.maxPt, color=checkColor("k2pi.maxPt", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.minPi0MassDiff", left.k2pi.minPi0MassDiff, right.k2pi.minPi0MassDiff, color=checkColor("k2pi.minPi0MassDiff", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.maxPi0MassDiff", left.k2pi.maxPi0MassDiff, right.k2pi.maxPi0MassDiff, color=checkColor("k2pi.maxPi0MassDiff", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.minKaonMassDiff", left.k2pi.minKaonMassDiff, right.k2pi.minKaonMassDiff, color=checkColor("k2pi.minKaonMassDiff", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.maxKaonMassDiff", left.k2pi.maxKaonMassDiff, right.k2pi.maxKaonMassDiff, color=checkColor("k2pi.maxKaonMassDiff", changed))
    print "{color}{0:25} {1:>10} {2:>10}\033[0m".format("k2pi.pi0Mass2DiffCoarse", left.k2pi.pi0Mass2DiffCoarse, right.k2pi.pi0Mass2DiffCoarse, color=checkColor("k2pi.pi0Mass2DiffCoarse", changed))
    
def compareCuts(left, right):
    listDiff = []
    if left.triggerMask != right.triggerMask:
        listDiff.append("triggerMask")
    if left.numVertex3              != right.numVertex3:
        listDiff.append("numVertex3")
    if left.minZVertex              != right.minZVertex:
        listDiff.append("minZVertex")
    if left.maxZVertex              != right.maxZVertex:
        listDiff.append("maxZVertex")
    if left.maxChi2Vertex           != right.maxChi2Vertex:
        listDiff.append("maxChi2Vertex")
    if left.maxExtraTracks          != right.maxExtraTracks:
        listDiff.append("maxExtraTracks")
    if left.maxTrackTime            != right.maxTrackTime:
        listDiff.append("maxTrackTime")
    if left.numBadTrackCombi        != right.numBadTrackCombi:
        listDiff.append("numBadTrackCombi")
    if left.numXCandidates          != right.numXCandidates:
        listDiff.append("numXCandidates")
    if left.minTrackMomentum        != right.minTrackMomentum:
        listDiff.append("minTrackMomentum")
    if left.maxTrackMomentum        != right.maxTrackMomentum:
        listDiff.append("maxTrackMomentum")
    if left.numAddGoodCluster       != right.numAddGoodCluster:
        listDiff.append("numAddGoodCluster")
    if left.lkrAcceptance           != right.lkrAcceptance:
        listDiff.append("lkrAcceptance")
    if left.minGammaEnergy          != right.minGammaEnergy:
        listDiff.append("minGammaEnergy")
    if left.minDeadCellDist         != right.minDeadCellDist:
        listDiff.append("minDeadCellDist")
    if left.minGammaDCHRadius       != right.minGammaDCHRadius:
        listDiff.append("minGammaDCHRadius")
    if left.deflectedElDist       != right.deflectedElDist:
        listDiff.append("deflectedElDist")
    if left.unDeflectedElDist       != right.unDeflectedElDist:
        listDiff.append("unDeflectedElDist")
    if left.k2pi.minTotalMomentum   != right.k2pi.minTotalMomentum:
        listDiff.append("k2pi.minTotalMomentum")
    if left.k2pi.maxTotalMomentum   != right.k2pi.maxTotalMomentum:
        listDiff.append("k2pi.maxTotalMomentum")
    if left.k2pi.maxPt              != right.k2pi.maxPt:
        listDiff.append("k2pi.maxPt")
    if left.k2pi.minPi0MassDiff     != right.k2pi.minPi0MassDiff:
        listDiff.append("k2pi.minPi0MassDiff")
    if left.k2pi.maxPi0MassDiff     != right.k2pi.maxPi0MassDiff:
        listDiff.append("k2pi.maxPi0MassDiff")
    if left.k2pi.minKaonMassDiff    != right.k2pi.minKaonMassDiff:
        listDiff.append("k2pi.minKaonMassDiff")
    if left.k2pi.maxKaonMassDiff    != right.k2pi.maxKaonMassDiff:
        listDiff.append("k2pi.maxKaonMassDiff")
    if left.k2pi.pi0Mass2DiffCoarse != right.k2pi.pi0Mass2DiffCoarse:
        listDiff.append("k2pi.pi0Mass2DiffCoarse")
    if left.kmu3.maxTotalMomentum   != right.kmu3.maxTotalMomentum:
        listDiff.append("kmu3.maxTotalMomentum")
    if left.kmu3.minPt              != right.kmu3.minPt:
        listDiff.append("kmu3.minPt")
    if left.kmu3.maxPt              != right.kmu3.maxPt:
        listDiff.append("kmu3.maxPt")
    if left.kmu3.maxPi0MassDiff     != right.kmu3.maxPi0MassDiff:
        listDiff.append("kmu3.maxPi0MassDiff")
    if left.kmu3.maxMissMassSq      != right.kmu3.maxMissMassSq:
        listDiff.append("kmu3.maxMissMassSq")
        
    return listDiff


if len(sys.argv)<2:
    print "Not enough arguments. Expecting: input.root [scanID]"
    sys.exit(0)

f = root_open(sys.argv[1], "open")

cutsSize = len(list(f.pid_pi)[0]["pid_res.good"])

default_cut = None
cutsLists = []
for i,cut in enumerate(f.cutsDefinition):
    default = cut.lists.getDefaultIndex()
    if i==default:
        default_cut = ROOT.ScanCuts(cut.lists)
    cutsLists.append(ROOT.ScanCuts(cut.lists))
    if i>=cutsSize:
        break

def filterVars(val):
    if "_" in val:
        return False
    if val[0]=="k":
        return False
    if val[0].isupper():
        return False
    if "operator" in val:
        return False
    
    return True

listVars = dir(ROOT.ScanCuts())
listVars = filter(filterVars, listVars)

for i, cuts in enumerate(cutsLists):
    if i==default:
        print "--------------- Default cut -------------------"
    printCuts(default_cut, cuts, i, compareCuts(default_cut, cuts))
    
    
