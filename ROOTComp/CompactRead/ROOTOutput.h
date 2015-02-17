/*
 * ROOTOutput.h
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

#ifndef ROOTOUTPUT_H_
#define ROOTOUTPUT_H_

#include <fstream>
#include "exportClasses.h"
#include "CompactImport.h"

class TTree;
class TFile;

class ROOTOutput {
public:
	ROOTOutput();
	virtual ~ROOTOutput();

	bool openOutput(bool doMC, bool doOutput, ROOTBurst &rootBurst,
			ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent,
			NGeom &rootGeom, ROOTMCEvent &rootMC, ROOTPhysicsEvent &rootPhysics,
			ROOTFileHeader &outputFileHeader);
	void fill();
	void close();

	std::ofstream& f1() {
		return fprt;
	}
	;
	std::ofstream& f2() {
		return fprt;
	}
	;
private:
	std::string prefix;
	std::ofstream fprt, fprt2;
	TTree *outTree, *outHeaderTree;
	TFile *outFile;
};

#endif /* ROOTOUTPUT_H_ */
