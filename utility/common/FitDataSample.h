/*
 * FitDataSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITDATASAMPLE_H_
#define COMMON_FITDATASAMPLE_H_

#include "Sample.h"

class FitDataSample: public Sample {
public:
	FitDataSample(int index, ConfigFile &cfg);
	virtual ~FitDataSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);

	double getTestA() const {
		return fTestA;
	}

	void setTestA(double testA) {
		fTestA = testA;
	}

private:
	TH1D *dSig;
	double fTestA;
};

#endif /* COMMON_FITDATASAMPLE_H_ */
