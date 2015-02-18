/*
 * ROOTOutput.h
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

#ifndef ROOTOUTPUT_H_
#define ROOTOUTPUT_H_

// Compact includes
#include "exportClasses.h"

// Std includes
#include <fstream>

// Local includes
#include "CompactImport.h"

class TTree;
class TFile;
class ScanCuts;

class ROOTOutput {
public:
	ROOTOutput();
	virtual ~ROOTOutput();

	bool openOutput(bool doMC, bool doOutput, bool doScan, ROOTBurst &rootBurst,
			ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent,
			NGeom &rootGeom, ROOTMCEvent &rootMC, ROOTPhysicsEvent &rootPhysics,
			ROOTFileHeader &outputFileHeader, ScanCuts &cutsDefinition);
	void fillEvent();
	void fillCuts();
	void close();

	std::ofstream& f1() {
		return fprt;
	};
	std::ofstream& f2() {
		return fprt;
	};

	void newResult(bool pass){ scanPass.push_back(pass);};
	void resetResult(){ scanPass.clear();};
private:
	std::string generateROOTName();
	std::string generatePassName();
	std::string generateFailName();
	std::string expandPrefix();

	std::string prefix;
	std::ofstream fprt, fprt2;
	TTree *outTree, *outHeaderTree, *outCuts;
	TFile *outFile;
	std::vector<bool> scanPass;
};

#endif /* ROOTOUTPUT_H_ */
