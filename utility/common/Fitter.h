/*
 * Fitter.h
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITTER_H_
#define COMMON_FITTER_H_

#include <sys/stat.h>
#include "FitMCSample.h"
#include "FitDataSample.h"

class RunWeights;
class TFile;

class Fitter {
public:
	Fitter();
	virtual ~Fitter();

	void prepareSamples(ConfigFile &cfg);
	void fillSamples();
	void getSamples();
	void mergeSamples();

	void setBinning(int nbins, double* binning) {
		fNBins = nbins;
		fBinning = binning;
	}

	void setRunWeights(const RunWeights* runWeights) {
		fRunWeights = runWeights;
	}

private:
	TFile *fTempFile;
	std::vector<FitMCSample*> fMCSamples;
	std::vector<FitDataSample*> fDataSamples;
	double *fBinning;
	const RunWeights *fRunWeights;
	int fNBins;

	FitMCSample   fFinalMCSample;
	FitDataSample fFinalDataSample;
};

#endif /* COMMON_FITTER_H_ */
