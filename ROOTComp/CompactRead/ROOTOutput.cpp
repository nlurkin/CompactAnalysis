/*
 * ROOTOutput.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

#include "ROOTOutput.h"
#include <TFile.h>
#include <TTree.h>
#include "exportClasses.h"
#include <iostream>

ROOTOutput::ROOTOutput() {
	// TODO Auto-generated constructor stub

}

ROOTOutput::~ROOTOutput() {
	// TODO Auto-generated destructor stub
	close();
}

bool ROOTOutput::openOutput(bool doMC, bool doOutput, ROOTBurst &rootBurst,
		ROOTRawEvent &rawEvent, ROOTCorrectedEvent &corrEvent, NGeom &rootGeom,
		ROOTMCEvent &rootMC, ROOTPhysicsEvent &rootPhysics,
		ROOTFileHeader &outputFileHeader) {
	outFile = gFile;
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
		fprt.open((prefix + ".txt").c_str(), ofstream::out);
		fprt2.open((prefix + "pass.txt").c_str(), ofstream::out);
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
