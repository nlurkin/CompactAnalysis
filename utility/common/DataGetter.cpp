/*
 * Combiner.cpp
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#include "DataGetter.h"

#include <sys/stat.h>
#include <TFile.h>
#include <TString.h>

#include "ConfigFile.h"

using namespace std;

DataGetter::DataGetter() :
		fBinning(nullptr), fRunWeights(nullptr), fNBins(0) {
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

DataGetter::~DataGetter() {
	TString tempFileName = fTempFile->GetName();
	fTempFile->Close();
	remove(tempFileName.Data());
}

void DataGetter::fillSamples() {
	for (auto mcSample : fMCSamples)
		mcSample->fill(fTempFile, fNBins, fBinning);

	for (auto dataSample : fDataSamples)
		dataSample->fill(fTempFile, fNBins, fBinning);
}

void DataGetter::getSamples() {
	for (auto mcSample : fMCSamples)
		mcSample->get(fTempFile);

	for (auto dataSample : fDataSamples)
		dataSample->get(fTempFile);
}
