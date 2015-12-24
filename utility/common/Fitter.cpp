/*
 * Fitter.cpp
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#include "Fitter.h"
#include <TString.h>
#include <TFile.h>
#include <iostream>

#include "FitMCSample.h"

using namespace std;

Fitter::Fitter() :
		fDataSample(nullptr),
		fBinning(nullptr),
		fRunWeights(nullptr),
		fNBins(0)
{
	struct stat buffer;
	srand(time(0));

	TString tempFileName = ".tempComb";
	tempFileName += rand() % 99999;
	tempFileName += ".root";
	while (stat(tempFileName.Data(), &buffer) == 0) {
		tempFileName = ".tempComb";
		tempFileName += rand() % 99999;
		tempFileName += ".root";
	}
	fTempFile = TFile::Open(tempFileName, "RECREATE");
}

Fitter::~Fitter() {
	TString tempFileName = fTempFile->GetName();
	fTempFile->Close();
	remove(tempFileName.Data());
}

void Fitter::prepareSamples(ConfigFile &cfg) {
	int prevIndex = -1;
	int newIndex;

	FitMCSample *tempSample;
	//Getting MC
	for (unsigned int i = 0; i < cfg.getMcFileNames().size(); ++i) {
		//TODO do this at ConfigFile level (align all vectors)
		//If not enough index in index vector -> assume same as previous
		if (i >= cfg.getMcIndexes().size())
			newIndex = prevIndex;
		else
			newIndex = cfg.getMcIndexes()[i];
		//If index is different from previous -> new sample
		if (prevIndex != newIndex) {
			//New sample
			tempSample = new FitMCSample(newIndex, cfg);
			tempSample->setBr(cfg.getBrs()[newIndex]);
			tempSample->setWeights(fRunWeights);
			tempSample->setOutputFile(cfg.getMcOutputFiles()[newIndex]);
			prevIndex = newIndex;
		}

		//Open new input file
		cout << cfg.getMcFileNames()[i] << endl;
		tempSample->addFile(cfg.getMcFileNames()[i]);

		fMCSamples.push_back(tempSample);
	}

	if (cfg.getDataFileNames().size() > 0){
		fDataSample = new FitDataSample(0, cfg);
		fDataSample->setBr(1);
		fDataSample->setTestA(cfg.getTestA());
		fDataSample->setWeights(fRunWeights);
		fDataSample->setOutputFile(cfg.getDataOutputFiles()[0]);
	}

	for (unsigned int i = 0; i < cfg.getDataFileNames().size(); ++i) {
		fDataSample->addFile(cfg.getDataFileNames()[i]);
	}
}

void Fitter::fillSamples() {
	for (auto mcSample : fMCSamples) {
		mcSample->fill(fTempFile, fNBins, fBinning);
	}

	if(fDataSample) fDataSample->fill(fTempFile, fNBins, fBinning);
}

void Fitter::getSamples(){
	for (auto mcSample : fMCSamples) {
		mcSample->get(fTempFile);
	}

	if(fDataSample) fDataSample->get(fTempFile);
}
