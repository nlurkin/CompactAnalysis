/*
 * ROOTOutput.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

// ROOT includes
#include <TFile.h>
#include <TTree.h>

// Local includes
#include "ROOTOutput.h"
#include "ScanCuts.h"

#include <iostream>


#ifdef __MAKECINT__
#pragma link C++ class vector<bool>+;
#endif

ROOTOutput::ROOTOutput() :
		prefix("output"), outTree(nullptr), outHeaderTree(nullptr), outCuts(nullptr), outFile(nullptr) {

}

ROOTOutput::~ROOTOutput() {
	close();
}

bool ROOTOutput::openOutput(bool doMC, bool doOutput, bool doScan, ROOTBurst &rootBurst,
		ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent, NGeom &rootGeom,
		ROOTMCEvent &rootMC, ROOTPhysicsEvent &rootPhysics,
		ROOTFileHeader &outputFileHeader, ScanCuts &cutsDefinition) {
	outFile = TFile::Open(generateROOTName().c_str(), "RECREATE");
	outTree = new TTree("event", "Event");
	outHeaderTree = new TTree("header", "Header");

	outTree->Branch("rawBurst", "ROOTBurst", &rootBurst);
	outTree->Branch("rawEvent", "ROOTRawEvent", &rawEvent);
	outTree->Branch("corrEvent", "ROOTCorrectedEvent", &corrEvent);
	if (doMC)
		outTree->Branch("mc", "ROOTMCEvent", &rootMC);

	outTree->Branch("pi0dEvent", "ROOTPhysicsEvent", &rootPhysics);
	outHeaderTree->Branch("header", "ROOTFileHeader", &outputFileHeader);
	outHeaderTree->Branch("geom", "NGeom", &rootGeom);

	if (doOutput) {
		fprt.open(generateFailName().c_str(), ofstream::out);
		fprt2.open(generatePassName().c_str(), ofstream::out);
	}
	if(doScan) {
		outCuts = new TTree("cutsDefinition", "CutsDefinition");
		outCuts->Branch("lists", "ScanCuts", &cutsDefinition);
		outTree->Branch("cutsResult", &scanPass);
	}
	return true;
}

void ROOTOutput::fillEvent() {
	outTree->Fill();
}

void ROOTOutput::fillCuts() {
	if(!outCuts) return;
	outCuts->Fill();
}

void ROOTOutput::close() {
	if (outFile) {
		if (outFile->IsOpen()) {
			outTree->Write();
			outHeaderTree->Fill();
			outHeaderTree->Write();
			if(outCuts) outCuts->Write();

			outFile->Close();
		}
	}
}

std::string ROOTOutput::generateROOTName() {
	return expandPrefix() + ".root";
}

std::string ROOTOutput::generatePassName() {
	return expandPrefix() + "pass.txt";
}

std::string ROOTOutput::generateFailName() {
	return expandPrefix() + ".txt";
}

std::string ROOTOutput::expandPrefix() {
	if (prefix.find('~') != std::string::npos)
		prefix = prefix.replace(prefix.find('~'), 1,
				std::string("/afs/cern.ch/user/n/nlurkin"));
	return prefix;
}
