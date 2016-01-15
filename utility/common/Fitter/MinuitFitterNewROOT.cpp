/*
 * MinuitFitter.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#include "../Fitter/MinuitFitterNewROOT.h"

#include <iostream>

using namespace std;

MinuitFitterNewROOT::MinuitFitterNewROOT(int bins) :
		fNBins(bins), fComp(nullptr), fSigComp(nullptr) {
}

void MinuitFitterNewROOT::init(double* binning) {
	fComp = new TH1D("comp", "comp", fNBins-1, binning);
	fSigComp = new TH1D("sigcomp", "sigcomp", fNBins-1, binning);
}

MinuitFitterNewROOT::~MinuitFitterNewROOT() {
	// TODO Auto-generated destructor stub
}

void MinuitFitterNewROOT::fit() {
	mnfree(0);
	FixParameter(0);

	MinuitFitter::fit();
}

double MinuitFitterNewROOT::minFunction(double, double a) {
	double chi2 = 0.;
	//double M_i, D_i, m_i;
	//double a_i, b_i, g_i;
	//double err;

	fComp->Add(fMCSamples->getAlpha(), 1);
	fComp->Add(fMCSamples->getBeta(), 2*a);
	fComp->Add(fMCSamples->getGamma(), a*a);
	fSigComp->Add(fDataSamples->getSig(), 1);
	/*for (int i = 0; i <= fNBins; ++i) {
		M_i = 0;
		a_i = 0;
		b_i = 0;
		g_i = 0;
		FitMCSample::bContent bins = fMCSamples->getBinContent(i);
		M_i += bins.dNew;
		a_i += bins.dAlpha;
		b_i += bins.dBeta;
		g_i += bins.dGamma;
		FitMCSample::bContent binsErr = fMCSamples->getBinError(i);
		err = 0;
		//err += binsErr.dNew;
		err += binsErr.dAlpha;
		err += binsErr.dBeta;
		err += binsErr.dGamma;


		D_i = 0;
		D_i += fDataSamples->getBinContent(i).dSig;

		m_i = function(1, a, a_i, b_i, g_i);

		fComp->Fill(fComp->GetBinCenter(i), m_i);
		fComp->SetBinError(i, err);
		fSigComp->SetBinContent(i, D_i);
	}*/

	chi2 = fSigComp->Chi2Test(fComp, "UW CHI2 P");
	fComp->Reset("M");
	fSigComp->Reset("M");
	return chi2;
}

double MinuitFitterNewROOT::function(double G, double a, double a_i, double b_i,
		double g_i) {
	return G * (a_i + 2. * a * b_i + a * a * g_i);
}
