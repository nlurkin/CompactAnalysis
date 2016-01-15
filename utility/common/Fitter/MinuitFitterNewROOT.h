/*
 * MinuitFitter.h
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_MINUITFITTERNEWROOT_H_
#define COMMON_MINUITFITTERNEWROOT_H_

#include <TH1D.h>
#include "../Fitter/MinuitFitter.h"

class MinuitFitterNewROOT: public MinuitFitter {
public:
	MinuitFitterNewROOT(int bins);
	virtual ~MinuitFitterNewROOT();

	void init(double*binning);

	void fit();
	double minFunction(double n, double a);

	double function(double G, double a, double a_i, double b_i, double g_i);

private:
	int fNBins;

	TH1D *fComp, *fSigComp;
};

#endif /* COMMON_MINUITFITTERNEWROOT_H_ */
