/*
 * Fitter.h
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITTER_H_
#define COMMON_FITTER_H_

#include <sys/stat.h>
#include "../Samples/FitMCSample.h"
#include "../Samples/FitDataSample.h"

class RunWeights;
class TFile;

class Fitter{
public:
	Fitter();
	virtual ~Fitter();

	void prepareSamples(ConfigFile &cfg);
	void fillSamples();
	void getSamples();
	void mergeSamples();

	void PrepareHistos(std::vector<int> allColors, std::vector<int> dataColors);

	void fit(bool useNew, bool userROOT);

	void setBinning(int nbins, double* binning) {
		fNBins = nbins;
		fBinning = binning;
	}

	void setRunWeights(const RunWeights* runWeights) {
		fRunWeights = runWeights;
	}

	static void minFunctionStatic(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,	Int_t flag);

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
