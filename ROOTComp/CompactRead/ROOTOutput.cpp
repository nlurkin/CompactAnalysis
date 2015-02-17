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


ROOTOutput::ROOTOutput() :
		prefix("output"), outTree(nullptr), outHeaderTree(nullptr), outFile(nullptr) {

}

ROOTOutput::~ROOTOutput() {
	close();
}

bool ROOTOutput::openOutput(bool doMC, bool doOutput, ROOTBurst &rootBurst,
		ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent, NGeom &rootGeom,
		ROOTMCEvent &rootMC, ROOTPhysicsEvent &rootPhysics,
		ROOTFileHeader &outputFileHeader) {
	outFile = TFile::Open(generateROOTName().c_str(), "RECREATE");
	outTree = new TTree("event", "Event");
	outHeaderTree = new TTree("header", "Header");

	outTree->Branch("rawBurst", "ROOTBurst", &rootBurst);
	outTree->Branch("rawEvent", "ROOTRawEvent", &rawEvent);
	outTree->Branch("corrEvent", "ROOTCorrectedEvent", &corrEvent);
	outTree->Branch("geom", "NGeom", &rootGeom);
	if (doMC)
		outTree->Branch("mc", "ROOTMCEvent", &rootMC);

	outTree->Branch("pi0dEvent", "ROOTPhysicsEvent", &rootPhysics);
	outHeaderTree->Branch("header", "ROOTFileHeader", &outputFileHeader);

	if (doOutput) {
		fprt.open(generateFailName().c_str(), ofstream::out);
		fprt2.open(generatePassName().c_str(), ofstream::out);
	}
	return true;
}

void ROOTOutput::fill() {
	outTree->Fill();
}

void ROOTOutput::close() {
	if (outFile) {
		if (outFile->IsOpen()) {
			outTree->Write();
			outHeaderTree->Fill();
			outHeaderTree->Write();

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
