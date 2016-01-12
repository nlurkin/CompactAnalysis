/*
 * Sample.cpp
 *
 *  Created on: Dec 21, 2015
 *      Author: nlurkin
 */

#include "Sample.h"
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "ScanCuts.h"

using namespace std;

Sample::Sample() :
		fIndex(-1), fBr(1), fOutputFD(nullptr), fCfg(nullptr), fWeights(
				nullptr), fMainSubSample(0) {
}

Sample::Sample(int index, ConfigFile *cfg) :
		fIndex(index), fBr(1), fOutputFD(nullptr), fCfg(cfg), fWeights(nullptr), fMainSubSample(
				0) {
}

Sample::~Sample() {
	// TODO Auto-generated destructor stub
}

bool Sample::addFile(string fileName) {
	if (!TString(fileName).Contains(".root")) {
		cout << "List file detected..." << endl;
		ifstream listFile(fileName);
		string buffer;
		while (getline(listFile, buffer)) {
			fListFiles.push_back(buffer);
		}
	} else {
		fListFiles.push_back(fileName);
	}

	return true;
}

void Sample::fill(TFile* tempFD, int nbins, double* bins) {
	initOutput();
	for (auto ss : fSubSamples)
		ss->initHisto(nbins, bins, fCfg);
	int iFile = 0;
	double totalFiles = fListFiles.size();
	TFile *inFileFD;
	cout << "Sample " << fIndex << " with BR " << fBr << endl;
	for (auto files : fListFiles) {
		cout << "Processing file " << files << " " << setprecision(2)
				<< std::fixed << iFile * 100. / totalFiles << "% " << iFile
				<< "/" << totalFiles << endl;
		inFileFD = TFile::Open(files.c_str());

		//Request the TTree reading function
		doFill(inFileFD, tempFD);

		//Close the input file
		inFileFD->Close();
		iFile++;
	}
	closeOutput(tempFD);
}

void Sample::get(TFile* tempFD) {
	TFile *inFileFD;
	cout << "Sample " << fIndex << " with BR " << fBr << endl;
	cout << "Opening file " << fOutputFile << endl;
	inFileFD = TFile::Open(fOutputFile.c_str());

	//Request the Histo reading function
	doGet(inFileFD, tempFD);

	//Close the input file
	inFileFD->Close();
}

void Sample::initOutput() {
	fOutputFD = TFile::Open(fOutputFile.c_str(), "RECREATE");
	int index = 0;
	for (auto ss : fSubSamples) {
		fOutputFD->mkdir(Form("%i", index));
		fOutputFD->cd(Form("%i", index));
		ss->initOutput();
	}
}

void Sample::closeOutput(TFile* tempFD) {
	fOutputFD->cd();

	int index = 0;
	for (auto ss : fSubSamples) {
		fOutputFD->cd(Form("%i", index));
		ss->writeTree();
		ss->doWrite();
		fOutputFD->cd();
		++index;
	}

	fOutputFD->Close();

	tempFD->cd();
	for (auto ss : fSubSamples)
		ss->doSetName();
}

Sample* Sample::Add(const Sample* other) {
	for (unsigned int i=0; i < fSubSamples.size(); i++)
		fSubSamples[i]->Add(other->fSubSamples[i]);
	return this;
}

void Sample::doFill(TFile* inputFD, TFile* tempFD) {
	//Get the TTree
	//Input
	ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
	ROOTBurst *burstBrch = new ROOTBurst();
	ROOTRawEvent *rawBrch = new ROOTRawEvent();
	ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
	ROOTFileHeader *headerBrch = new ROOTFileHeader();
	ROOTMCEvent *mcEvent = 0;
	NGeom *geomBrch = new NGeom();
	vector<bool> *cutsPass = 0;

	TTree *t = (TTree*) inputFD->Get("event");
	TTree *th = (TTree*) inputFD->Get("header");
	if (t->GetListOfBranches()->Contains("mc"))
		mcEvent = new ROOTMCEvent();
	if (t->GetListOfBranches()->Contains("cutsResult")) {
		cutsPass = new vector<bool>;
	}

	t->SetBranchAddress("pi0dEvent", &eventBrch);
	t->SetBranchAddress("rawBurst", &burstBrch);
	t->SetBranchAddress("rawEvent", &rawBrch);
	t->SetBranchAddress("corrEvent", &corrBrch);
	th->SetBranchAddress("header", &headerBrch);
	th->SetBranchAddress("geom", &geomBrch);
	if (mcEvent)
		t->SetBranchAddress("mc", &mcEvent);
	if (cutsPass) {
		t->SetBranchAddress("cutsResult", &cutsPass);
		if (fCfg->getScanId() == -1) {
			TTree *tc = (TTree*) inputFD->Get("cutsDefinition");
			ScanCuts *cutsLists = new ScanCuts();
			tc->SetBranchAddress("lists", &cutsLists);
			tc->GetEntry(0);
			fMainSubSample = cutsLists->getDefaultIndex();
		} else
			fMainSubSample = fCfg->getScanId();
	}

	tempFD->cd();
	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = 0;
	for (int i = 0; i < th->GetEntries(); i++) {
		th->GetEntry(i);
		totalChanEvents += headerBrch->NProcessedEvents;
	}

	//Read events and fill histo
	int i = 0;

	for (auto ss : fSubSamples)
		ss->initNewFile(totalChanEvents, nevt);

	cout << "Filling " << nevt << endl;
	for (; i < nevt; ++i) {
		if (i % 10000 == 0)
			cout << setprecision(2) << i * 100. / (double) nevt << "% " << i
					<< "/" << nevt << "\r";
		cout.flush();
		t->GetEntry(i);
		for (auto ss : fSubSamples)
			ss->processEvent(eventBrch, burstBrch, rawBrch, corrBrch,
					headerBrch, mcEvent, geomBrch, cutsPass, fCfg, fWeights);
	}
}

void Sample::doGet(TFile* inputFD, TFile* tempFD) {
	int index = 0;
	for (auto ss : fSubSamples) {
		inputFD->cd(Form("%i", index));
		ss->doGet(gDirectory, tempFD);
	}
}
