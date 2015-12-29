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
	typedef struct bContent_t{
		double dSig;
	} bContent;

	FitDataSample();
	FitDataSample(int index, ConfigFile *cfg);
	virtual ~FitDataSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);
	virtual void scaleToData(double) {};
	virtual void setPlotStyle(std::vector<int> color);
	virtual void populateStack(InputFitDrawer &drawer);
	virtual void populateFit(FitResultDrawer &drawer, double norm, double a);

	virtual TH1D* getMainHisto() { return dSig; }

	void scale();

	bContent getBinContent(int bin);

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
