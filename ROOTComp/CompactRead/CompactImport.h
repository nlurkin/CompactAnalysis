/*
 * CompactImport.h
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

#ifndef COMPACTIMPORT_H_
#define COMPACTIMPORT_H_

#include <TString.h>

class TChain;
class ROOTRawEvent;
class ROOTCorrectedEvent;
class ROOTBurst;
class ROOTFileHeader;
class NGeom;
class ROOTMCEvent;
class ROOTOutput;

class CompactImport {
public:
	CompactImport();
	virtual ~CompactImport();

	bool readInput(TString fName, bool isList);

	int getNEvents();
	int getNHeaders();

	int nextEvent(ROOTFileHeader &outputHeader);
	int firstEvent(ROOTFileHeader &outputHeader, int first=0);
	bool eof();

	void associateTrees(ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent,
			ROOTBurst &rootBurst, ROOTFileHeader &rootFileHeader,
			NGeom &rootGeom, ROOTMCEvent &rootMC);

	bool getHasMC() {
		return hasMC;
	}
	;
private:
	int currentEvent;
	int currentFile;
	int nEvents;
	int nFiles;
	bool hasMC;
	TString currentFileName;
	TChain *inTree;
	TChain *headerTree;

	//### Ptr for TTree
	ROOTRawEvent *rawEvent_ptr;
	ROOTCorrectedEvent *corrEvent_ptr;
	ROOTBurst *rootBurst_ptr;
	ROOTFileHeader *rootFileHeader_ptr;
	NGeom *Geom;
	ROOTMCEvent *rootMC_ptr;
};

#endif /* COMPACTIMPORT_H_ */
