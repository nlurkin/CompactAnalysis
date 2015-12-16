#!/bin/env python
import subprocess
import sys
import threading
import time

import fcntl
import os
import shutil
import re


filePiTemplate = """mcfiles= /afs/cern.ch/user/n/nlurkin/work/L2_L3/{sampleName}/listpi{fileID}.lst 
mcIndex=0
mcout=pi{fileID}_{scanV}.root
mccolors=630 740 850
brs=2.066E-1
mclegends=K^{{+}}->#pi^{{+}}#pi^{{0}}_{{d}}

binsfile=/afs/cern.ch/user/n/nlurkin/work/bins/bins_x0.1.dat
equalbin=true

{comment}scanid={scanV}"""

fileGlobalTemplate = """mcfiles= /afs/cern.ch/user/n/nlurkin/work/L2_L3/{sampleName}/listmu.lst 
mcIndex=0
mcout=mu_{scanV}.root
mccolors=330 440 550
brs=3.353E-2
mclegends=K^{{+}}->#pi^{{0}}_{{d}}#mu^{{+}}#nu

datafiles= /afs/cern.ch/user/n/nlurkin/work/L2_L3/{sampleName}/listdata.lst
dataout=data_{scanV}.root
datacolors=360
datalegends=Data

binsfile=/afs/cern.ch/user/n/nlurkin/work/bins/bins_x0.1.dat
equalbin=true

{comment}scanid={scanV}"""

fileAllTemplate = """
mcIndex=0 1
mcout=mu_{scanV}.root pi_{scanV}.root
mccolors=330 440 550 630 740 850
brs=3.353E-2 2.066E-1
mclegends=K^{{+}}->#pi^{{0}}_{{d}}#mu^{{+}}#nu K^{{+}}->#pi^{{+}}#pi^{{0}}_{{d}}

dataout=data_{scanV}.root
datacolors=360
datalegends=Data
datafactor=1

binsfile=/afs/cern.ch/user/n/nlurkin/work/bins/bins_x0.1.dat
equalbin=true"""

class myThread (threading.Thread):
    def __init__(self, tID, scanID, listFile, final):
        threading.Thread.__init__(self)
        self.listFile = listFile
        self.final = final
        self.scanID = scanID
        self.tID = tID
    
    def run(self):
        printTo(self.tID+10, "\rRunning " + self.listFile)
        if self.final:
            resultLine = runSSH("lxplus", "source env32.sh; cd {1}; ~/Compact/utility/build/fit {0} -g;".format(self.listFile, os.path.abspath(".")), self.tID, True)
            with open("result.dat", "a") as fd:
                fd.write("{0}={1}\n".format(self.scanID, resultLine))
        else:
            runSSH("lxplus", "source env32.sh; cd {1}; ~/Compact/utility/build/fit {0};".format(self.listFile, os.path.abspath(".")), self.tID)
        
        printTo(self.tID+10, "\rThread {0} finished for scan {1}".format(self.listFile, self.scanID))
        

def generateFiles(sample, scanV, noScan):
    comment = ""
    if noScan:
        comment = "#"
    for fileID in range(1,7):
        with open("listFile{0}.cfg".format(fileID), "w") as fd:
            fd.write(filePiTemplate.format(sampleName=sample, scanV=scanV, fileID=fileID, comment=comment))

    with open("listFile_global.cfg".format(fileID), "w") as fd:
        fd.write(fileGlobalTemplate.format(sampleName=sample, scanV=scanV, fileID=fileID, comment=comment))

    with open("listFile.cfg".format(fileID), "w") as fd:
        fd.write(fileAllTemplate.format(sampleName=sample, scanV=scanV, fileID=fileID, comment=comment))

def setNonBlocking(fd):
    """
    Set the file description of the given file descriptor to non-blocking.
    """
    flags = fcntl.fcntl(fd, fcntl.F_GETFL)
    flags = flags | os.O_NONBLOCK
    fcntl.fcntl(fd, fcntl.F_SETFL, flags)

def printTo(line, message):
    print "\033[{0};0H{1}                                                  ".format(5+line, message)
    
def clearScr():
    #print "\033[2J"
    print "\033[5:0H\033[J"

def runSSH(host, command, tID, ignErr=False):
    cmd = ["ssh", host, command]
    #print " ".join(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    #stdout = p.communicate()[0]
    #setNonBlocking(p.stdout)
    
    line = ""
    resultLine = ""
    while True:
        while True:
            try:
                line = p.stdout.readline()
            except IOError:
                continue
            else:
                break
        line = line.rstrip('\n')
        if line!="":
            pass
            #print line
        if "--Done--" in line or "brute" in line:
            printTo(tID, "\r[{0}] Complete ...".format(tID))
            time.sleep(1);
            break
        if "[ERROR]" in line and not ignErr:
            printTo(tID, "\r[{0}]File already exists... fail".format(tID))
            time.sleep(1);
            break;
        if "*** Break ***" in line:
            printTo(tID, "\r[{0}] stderr fired ... abort".format(tID))
            break
        if "RESULTLINE:" in line:
            resultLine = line[11:]
        if "Processing file" in line:
            m = re.findall("([0-9]+)/([0-9]+)", line)
            if m:
                printTo(tID, "\r[{0}] File {1}/{2}".format(tID, m[0][0], m[0][1]))
        
    p.terminate()
    printTo(tID+10, "\r[{0}] runSSH done".format(tID))
    return resultLine
    
    
def runParts(scanID):
    threadsL = []
    threadsL.append(myThread(0, scanID, "listFile1.cfg", False))
    threadsL.append(myThread(1, scanID, "listFile2.cfg", False))
    threadsL.append(myThread(2, scanID, "listFile3.cfg", False))
    threadsL.append(myThread(3, scanID, "listFile4.cfg", False))
    threadsL.append(myThread(4, scanID, "listFile5.cfg", False))
    threadsL.append(myThread(5, scanID, "listFile6.cfg", False))
    threadsL.append(myThread(6, scanID, "listFile_global.cfg", False))
    
    clearScr()
    
    for t in threadsL:
        t.daemon = True
        t.start()
        time.sleep(1)
            
    while threading.active_count()>1:
        time.sleep(0.1)
    #for t in threadsL:
    #    t.join()

def runFit(scanID):
    os.system("hadd pi_{0}.root pi?_{0}.root".format(scanID))
    t = myThread(0, scanID, "listFile.cfg", True)
    t.daemon = 1
    t.start()
    while threading.active_count()>1:
        time.sleep(0.1)

def testSample(sample):
    sampleFolder = "/afs/cern.ch/user/n/nlurkin/work/L2_L3/{0}".format(sample)
    if not os.path.exists(sampleFolder):
        print "Sample folder not found: {0}".format(sampleFolder)
        return False
    if not os.path.exists("{0}/listdata.lst".format(sampleFolder)):
        print "List files don't seem to exist for sample"
        return False
    return True

def mvFiles(scanID):
    myPath = "{0}".format(scanID) 
    if not os.path.exists(myPath):
        os.mkdir(myPath)
    shutil.move("pi1_{0}.root".format(myPath), "{0}/pi1.root".format(myPath))
    shutil.move("pi2_{0}.root".format(myPath), "{0}/pi2.root".format(myPath))
    shutil.move("pi3_{0}.root".format(myPath), "{0}/pi3.root".format(myPath))
    shutil.move("pi4_{0}.root".format(myPath), "{0}/pi4.root".format(myPath))
    shutil.move("pi5_{0}.root".format(myPath), "{0}/pi5.root".format(myPath))
    shutil.move("pi6_{0}.root".format(myPath), "{0}/pi6.root".format(myPath))
    shutil.move("pi_{0}.root".format(myPath), "{0}/pi.root".format(myPath))
    shutil.move("mu_{0}.root".format(myPath), "{0}/mu.root".format(myPath))
    shutil.move("data_{0}.root".format(myPath), "{0}/data.root".format(myPath))

if __name__=="__main__":
    if len(sys.argv)<3:
        print """Usage
        doScan.py sampleName nScan"""
        sys.exit(0)

    sample = sys.argv[1]
    nScan = int(sys.argv[2])
    startScan = 0
    if len(sys.argv)==4:
        startScan = int(sys.argv[3])
        
    if not testSample(sample):
        sys.exit(0)
    
    if nScan==0:
        generateFiles(sample, 0, True)
        runParts(0)
        runFit(0)
    for scanV in range(startScan, nScan):
        generateFiles(sample, scanV, False)
        runParts(scanV)
        runFit(scanV)
        mvFiles(scanV)
        printTo(-4, "Scan for {0} finished".format(scanV))
