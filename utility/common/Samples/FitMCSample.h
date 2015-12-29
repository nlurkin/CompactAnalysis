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
	typedef struct bContent_t{
		double d1,d2,d3,dNew,dAlpha,dBeta,dGamma;
	} bContent;

	FitMCSample();
	FitMCSample(int index, ConfigFile *cfg);
	virtual ~FitMCSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);
	virtual void scaleToData(double nData);
	virtual void setPlotStyle(std::vector<int> color);
	virtual void populateStack(InputFitDrawer &drawer);
	virtual void populateFit(FitResultDrawer &drawer, double norm, double a);
	virtual TH1D* getMainHisto() { return dAlpha; }

	void scale();

	bContent getBinContent(int bin);
	friend FitMCSample& operator+=(FitMCSample &first, const FitMCSample* other);
private:
	TH1D *d1, *d2, *d3;
	TH1D *dNew;
	TH1D *dAlpha, *dBeta, *dGamma;
};

#endif /* COMMON_FITMCSAMPLE_H_ */
