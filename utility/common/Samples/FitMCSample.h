/*
 * FitSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITMCSAMPLE_H_
#define COMMON_FITMCSAMPLE_H_

#include "../Interface/Sample.h"
#include "../Interface/MCSample.h"

class FitMCSample: public Sample, public MCSample {
public:
	typedef struct bContent_t{
		double d1,d2,d3,dNew,dAlpha,dBeta,dGamma;

		struct bContent_t &operator +=(struct bContent_t other){
			d1 += other.d1;
			d2 += other.d2;
			d3 += other.d3;
			dNew += other.dNew;
			dAlpha += other.dAlpha;
			dBeta += other.dBeta;
			dGamma += other.dGamma;

			return *this;
		}
	} bContent;

	FitMCSample();
	FitMCSample(int index, ConfigFile *cfg);
	virtual ~FitMCSample();

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins);
	virtual void scaleToData(bContent totalMC, double nData);
	virtual void setPlotStyle(std::vector<int> color);
	virtual void populateStack(HistoDrawer *drawer);
	virtual void populateFit(HistoDrawer *drawer, double norm, double a);
	virtual TH1D* getMainHisto() { return dAlpha; }
	virtual double getFFIntegral(double a);

	bContent getIntegrals();

	void scale();

	bContent getBinContent(int bin);
	FitMCSample* Add(const FitMCSample* other);
private:
	TH1D *d1, *d2, *d3;
	TH1D *dNew;
	TH1D *dAlpha, *dBeta, *dGamma;
};

#endif /* COMMON_FITMCSAMPLE_H_ */
