/*
 * FitDataSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITDATASAMPLE_H_
#define COMMON_FITDATASAMPLE_H_

#include "../Interface/DataSample.h"
#include "../Interface/Sample.h"

class FitDataSample: public Sample, public DataSample{
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
	virtual void populateStack(HistoDrawer *drawer);
	virtual void populateFit(HistoDrawer *drawer, double norm, double a);
	virtual double getFFIntegral(double a);
	virtual void renameHisto() {};

	virtual TH1D* getMainHisto() { return dSig; }

	void scale();

	bContent getBinContent(int bin);

	FitDataSample* Add(const FitDataSample* other);
private:
	TH1D *dSig;
};

#endif /* COMMON_FITDATASAMPLE_H_ */
