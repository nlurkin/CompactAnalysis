/*
 * Combiner.h
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DATAGETTER_H_
#define COMMON_DATAGETTER_H_

#include "Interface/MCSample.h"
#include "Interface/DataSample.h"
#include "Interface/Sample.h"
#include "ConfigFile.h"
#include <iostream>

class RunWeights;
class TFile;

class DataGetter {
public:
	DataGetter();
	virtual ~DataGetter();

	template <class TMCSample, class TDataSample>
	void prepareSamples(ConfigFile &cfg);
	void fillSamples();
	void getSamples();
	template <class TMCSample, class TDataSample>
	void mergeSamples();

	void setBinning(int nbins, double* binning) {
		fNBins = nbins;
		fBinning = binning;
	}

	void setRunWeights(const RunWeights* runWeights) {
		fRunWeights = runWeights;
	}

protected:
	TFile *fTempFile;
	std::vector<Sample*> fMCSamples;
	std::vector<Sample*> fDataSamples;
	double *fBinning;
	const RunWeights *fRunWeights;
	int fNBins;

	Sample *fFinalMCSample;
	Sample *fFinalDataSample;
};

template <class TMCSample, class TDataSample>
void DataGetter::prepareSamples(ConfigFile& cfg) {
	int prevIndex = -1;
	int newIndex;

	fFinalMCSample = new Sample(998, &cfg);
	fFinalDataSample = new Sample(999, &cfg);

	fFinalMCSample->setCfg(&cfg);
	fFinalDataSample->setCfg(&cfg);

	Sample *tempSample;
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
			tempSample = new Sample(newIndex, &cfg);
			tempSample->setBr(cfg.getBrs()[newIndex]);
			tempSample->setWeights(fRunWeights);
			tempSample->setOutputFile(cfg.getMcOutputFiles()[newIndex]);
			tempSample->setLegend(cfg.getMcLegendTitle()[newIndex]);
			prevIndex = newIndex;
		}

		//Open new input file
		tempSample->addFile(cfg.getMcFileNames()[i]);

		fMCSamples.push_back(tempSample);
	}

	Sample *dataSample;
	prevIndex = -1;
	for (unsigned int i = 0; i < cfg.getDataFileNames().size(); ++i) {
		//TODO do this at ConfigFile level (align all vectors)
		//If not enough index in index vector -> assume same as previous
		if (i >= cfg.getDataIndexes().size())
			newIndex = prevIndex;
		else
			newIndex = cfg.getDataIndexes()[i];
		//If index is different from previous -> new sample
		if (prevIndex != newIndex) {
			dataSample = new Sample(newIndex, &cfg);
			dataSample->setBr(1);
			dataSample->setTestA(cfg.getTestA());
			dataSample->setFactor(cfg.getDataFactor()[newIndex]);
			dataSample->setWeights(fRunWeights);
			dataSample->setOutputFile(cfg.getDataOutputFiles()[newIndex]);
			dataSample->setLegend(cfg.getDataLegendTitle()[newIndex]);
			prevIndex = newIndex;
		}

		dataSample->addFile(cfg.getDataFileNames()[i]);
		fDataSamples.push_back(dataSample);
	}

	fFinalDataSample->prepareNSubSamples<TDataSample>(cfg.getNScan());
	fFinalDataSample->prepareNSubSamples<TMCSample>(cfg.getNScan());

	for(auto sample : fDataSamples) sample->prepareNSubSamples<TDataSample>(cfg.getNScan());
	for(auto sample : fMCSamples) sample->prepareNSubSamples<TMCSample>(cfg.getNScan());
}

template <class TMCSample, class TDataSample>
void DataGetter::mergeSamples() {
	fFinalDataSample->initHisto(fNBins, fBinning);
	fFinalDataSample->renameHisto();
	fFinalMCSample->initHisto(fNBins, fBinning);
	fFinalMCSample->renameHisto();
	for (auto sample : fDataSamples)
		fFinalDataSample->Add(sample);

	std::vector<typename TMCSample::bContent> b;
	for (auto sample : fMCSamples) {
		std::vector<typename TMCSample::bContent> nVec = sample->getIntegrals<TMCSample>();
		for(unsigned int i=0; i<nVec.size(); i++) b[i] += nVec[i];
		fFinalMCSample->Add(sample);
	}

	for (auto sample : fMCSamples) {
		sample->scaleToData<TMCSample>(b, fFinalDataSample->getTotalSize());
	}

	fFinalMCSample->scaleToData<TMCSample>(b, fFinalDataSample->getTotalSize());
}


#endif /* COMMON_DATAGETTER_H_ */
