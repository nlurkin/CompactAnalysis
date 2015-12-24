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
	FitDataSample();
	FitDataSample(int index, ConfigFile &cfg);
	virtual ~FitDataSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);
	virtual void scaleToData(double nData) {};

	void scale();

	double getTestA() const {
		return fTestA;
	}

	void setTestA(double testA) {
		fTestA = testA;
	}

	double getFactor() const {
		return fFactor;
	}

	void setFactor(double factor) {
		fFactor = factor;
	}

	friend FitDataSample& operator+=(FitDataSample &first, const FitDataSample* other);
private:
	TH1D *dSig;
	double fTestA;
	double fFactor;
};

#endif /* COMMON_FITDATASAMPLE_H_ */
