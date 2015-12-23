/*
 * FitSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITMCSAMPLE_H_
#define COMMON_FITMCSAMPLE_H_

#include "Sample.h"

class FitMCSample: public Sample {
public:
	FitMCSample(int index, ConfigFile &cfg);
	virtual ~FitMCSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);

private:
	TH1D *d1, *d2, *d3;
	TH1D *dNew;
	TH1D *dAlpha, *dBeta, *dGamma;
};

#endif /* COMMON_FITMCSAMPLE_H_ */
