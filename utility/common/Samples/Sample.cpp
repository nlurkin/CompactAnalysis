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

using namespace std;

void initFitStruct(fitStruct &s){
	s.n1 = 0;
	s.nx = 0;
	s.nxx = 0;
	s.selEvents = 0;
	s.totEvents = 0;
}

void sumTreeFitStruct(fitStruct &in, TTree *t, fitStruct &out, double factor){
	for(int i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		out.n1 			+= factor*in.n1;
		out.nx 			+= factor*in.nx;
		out.nxx 		+= factor*in.nxx;
		out.selEvents 	+= factor*in.selEvents;
		out.totEvents 	+= factor*in.totEvents;
	}

	std::cout << "Total events: \t" << out.totEvents << std::endl;
	std::cout << "Sel.  events: \t" << out.selEvents << std::endl;
	std::cout << "n1    events: \t" << out.n1 << std::endl;
	std::cout << "nx    events: \t" << out.nx << std::endl;
	std::cout << "nxx   events: \t" << out.nxx << std::endl;
}


Sample::Sample() :
		fIndex(-1), fBr(1), fFitTree(nullptr),
		fOutputFD(nullptr), fCfg(nullptr), fWeights(nullptr)
{
	initFitStruct(fFitBrch);
}

Sample::Sample(int index, ConfigFile *cfg) :
		fIndex(index), fBr(1), fFitTree(nullptr),
		fOutputFD(nullptr), fCfg(cfg), fWeights(nullptr)
{
	initFitStruct(fFitBrch);
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
	initHisto(nbins, bins);
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
		//Input::getInputMCFill(ffd, fdo, cfg.getBrs()[prevIndex], prevIndex);

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
	fFitTree = new TTree("fitStruct", "fitStruct tree");
	fFitTree->Branch("fitStruct", &fFitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	fFitBrch.n1 = 0;
	fFitBrch.nx = 0;
	fFitBrch.nxx = 0;
	fFitBrch.selEvents = 0;
	fFitBrch.totEvents = 0;
}

void Sample::closeOutput(TFile* tempFD) {
	fOutputFD->cd();
	fFitTree->Fill();
	fFitTree->Write();

	doWrite();

	fOutputFD->Close();

	tempFD->cd();
	doSetName();
}

void Sample::scale(TH1 *histo, double scaleFactor){
	cout << "Scaling " << histo->GetName() << " " << histo->GetEntries() << " " << histo->Integral() << " " << fBr << " " << fFitBrch.totEvents << " " << scaleFactor << endl;
	histo->Scale(fBr/(fFitBrch.totEvents*scaleFactor));
}

Sample& operator +=(Sample &first, const Sample* other) {
	first.fFitBrch += other->fFitBrch;
	return first;
}

