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

	fFinalMCSample = new TMCSample();
	fFinalDataSample = new TDataSample();

	static_cast<TMCSample*>(fFinalMCSample)->setCfg(&cfg);
	static_cast<TDataSample*>(fFinalDataSample)->setCfg(&cfg);

	TMCSample *tempSample;
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
			tempSample = new TMCSample(newIndex, &cfg);
			tempSample->setBr(cfg.getBrs()[newIndex]);
			tempSample->setWeights(fRunWeights);
			tempSample->setOutputFile(cfg.getMcOutputFiles()[newIndex]);
			prevIndex = newIndex;
		}

		//Open new input file
		tempSample->addFile(cfg.getMcFileNames()[i]);

		fMCSamples.push_back(tempSample);
	}

	TDataSample *dataSample;
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
			dataSample = new TDataSample(newIndex, &cfg);
			dataSample->setBr(1);
			dataSample->setTestA(cfg.getTestA());
			dataSample->setFactor(cfg.getDataFactor()[newIndex]);
			dataSample->setWeights(fRunWeights);
			dataSample->setOutputFile(cfg.getDataOutputFiles()[newIndex]);
			prevIndex = newIndex;
		}

		dataSample->addFile(cfg.getDataFileNames()[i]);
		fDataSamples.push_back(dataSample);
	}
}

template <class TMCSample, class TDataSample>
void DataGetter::mergeSamples() {
	static_cast<TDataSample*>(fFinalDataSample)->initHisto(fNBins, fBinning);
	static_cast<TMCSample*>(fFinalMCSample)->initHisto(fNBins, fBinning);
	for (auto sample : fDataSamples)
		static_cast<TDataSample*>(fFinalDataSample)->Add(static_cast<TDataSample*>(sample));

	typename TMCSample::bContent b;
	for (auto sample : fMCSamples) {
		b += static_cast<TMCSample*>(sample)->getIntegrals();
		static_cast<TMCSample*>(fFinalMCSample)->Add(static_cast<TMCSample*>(sample));
	}

	for (auto sample : fMCSamples) {
		static_cast<TMCSample*>(sample)->scaleToData(b, static_cast<TDataSample*>(fFinalDataSample)->getTotalSize());
	}

	static_cast<TMCSample*>(fFinalMCSample)->scaleToData(b, static_cast<TDataSample*>(fFinalDataSample)->getTotalSize());
}


#endif /* COMMON_DATAGETTER_H_ */
