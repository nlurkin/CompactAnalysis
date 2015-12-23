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

Sample::Sample(int index, ConfigFile &cfg) :
		fIndex(index),
		fCfg(cfg)
{
	// TODO Auto-generated constructor stub
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
