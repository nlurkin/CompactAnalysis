/*
 * CompactImport.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: ncl
 */

#include "CompactImport.h"
#include <iostream>
#include <TChain.h>
#include <fstream>
#include "exportClasses.h"
#include <TFile.h>
#include "ROOTOutput.h"

CompactImport::CompactImport() :
		currentEvent(-1), currentFile(-1), nEvents(0), nFiles(0), hasMC(false), currentFileName(), inTree(
				new TChain("event")), headerTree(new TChain("header")), rawEvent_ptr(
				nullptr), corrEvent_ptr(nullptr), rootBurst_ptr(nullptr), rootFileHeader_ptr(
				nullptr), Geom(nullptr), rootMC_ptr(nullptr) {
}

CompactImport::~CompactImport() {
	// TODO Auto-generated destructor stub
}

bool CompactImport::readInput(TString fName, bool isList) {

	if (isList) {
		if (fName.Contains(".root")) {
			std::cerr << "ROOT file cannot be a list" << std::endl;
			return false;
		}
		TString inputFileName;
		std::ifstream inputList(fName.Data());
		while (inputFileName.ReadLine(inputList)) {
			if (inputFileName.Contains("/castor/")
					&& !inputFileName.Contains(
							"root://castorpublic.cern.ch//")) {
				TString svcClass = getenv("STAGE_SVCCLASS");
				if (svcClass == "")
					svcClass = "na62";
				inputFileName = "root://castorpublic.cern.ch//" + inputFileName
						+ "?svcClass=" + svcClass;
			}
			if (inputFileName.Contains("/eos/")
					&& !inputFileName.Contains("root://eosna62.cern.ch//")) {
				inputFileName = "root://eosna62.cern.ch//" + inputFileName;
			}
			inTree->AddFile(inputFileName);
			headerTree->AddFile(inputFileName);
		}
	} else {
		inTree->AddFile(fName);
		headerTree->AddFile(fName);
	}

	inTree->SetBranchAddress("rawBurst", &rootBurst_ptr);
	inTree->SetBranchAddress("rawEvent", &rawEvent_ptr);
	inTree->SetBranchAddress("corrEvent", &corrEvent_ptr);
	headerTree->SetBranchAddress("header", &rootFileHeader_ptr);
	headerTree->SetBranchAddress("geom", &Geom);
	if (inTree->GetListOfBranches()->Contains("mc")) {
		inTree->SetBranchAddress("mc", &rootMC_ptr);
		hasMC = true;
	}

	nEvents = inTree->GetEntries();
	nFiles = headerTree->GetEntries();

	return true;
}

int CompactImport::getNEvents() {
	return nEvents;
}

int CompactImport::getNHeaders() {
	return nFiles;
}

int CompactImport::firstEvent(ROOTFileHeader &outputHeader) {
	currentEvent = -1;
	currentFile = -1;
	currentFileName = "";
	return nextEvent(outputHeader);
}
int CompactImport::nextEvent(ROOTFileHeader &outputHeader) {
	currentEvent++;
	inTree->GetEntry(currentEvent);
	if (currentFileName != inTree->GetCurrentFile()->GetName()) {
		currentFileName = inTree->GetCurrentFile()->GetName();
		++currentFile;
		headerTree->GetEntry(currentFile);
		outputHeader.NProcessedEvents += rootFileHeader_ptr->NProcessedEvents;
		outputHeader.NFailedEvents += rootFileHeader_ptr->NFailedEvents;
	}

	return currentEvent;
}

bool CompactImport::eof() {
	if (currentEvent == nEvents)
		return true;
	else
		return false;
}

void CompactImport::associateTrees(ROOTRawEvent& rawEvent,
		ROOTCorrectedEvent& corrEvent, ROOTBurst& rootBurst,
		ROOTFileHeader& rootFileHeader, NGeom& rootGeom, ROOTMCEvent& rootMC) {

	rawEvent_ptr = &rawEvent;
	corrEvent_ptr = &corrEvent;
	rootBurst_ptr = &rootBurst;
	rootFileHeader_ptr = &rootFileHeader;
	Geom = &rootGeom;
	rootMC_ptr = &rootMC;
}
