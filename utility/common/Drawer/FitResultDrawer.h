/*
 * FitResultDrawer.h
 *
 *  Created on: Dec 29, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_FITRESULTDRAWER_H_
#define COMMON_DRAWER_FITRESULTDRAWER_H_

#include "../Interface/HistoDrawer.h"

class FitResultDrawer: public HistoDrawer {
public:
	FitResultDrawer();
	virtual ~FitResultDrawer();

	virtual void draw();

	void setMc(TH1D* mc, TH1D* data) {
		fMC = mc;
		fSig = data;
	}

	void setBinning(int nbins, double* binning) {
		fNBins = nbins;
		fBinning = binning;
	}

	void setResult(double norm, double a, double aerr) {
		fNorm= norm;
		fA = a;
		fAErr = aerr;
	}
private:
	TH1D *fMC, *fSig;
	THStack *fFit;
	THStack *fFitSig;
	TLegend *fLegFit;
	TLegend *fLegFitSig;

	double *fBinning;
	int fNBins;
	double fNorm;
	double fA, fAErr;
};

#endif /* COMMON_DRAWER_FITRESULTDRAWER_H_ */
