/*
 * MinuitFitter.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#include "../Fitter/MinuitFitter.h"

#include <iostream>
#include <TMath.h>
#include <TGraph.h>

using namespace std;

MinuitFitter *gTemp;

MinuitFitter::MinuitFitter() :
		fNorm(0), fNormErr(0), fFormFactor(0), fFormFactorErr(0), fMinimum(0) {
	int flag;
	double args[1];

	args[0] = 0;
	mnexcm("SET PRINTOUT", args, 1, flag);
	args[0] = 1;
	mnexcm("SET ERROR", args, 1, flag);
	args[0] = 1;
	mnexcm("SET STRATEGY", args, 1, flag);
	mnparm(0, "G", 1., 100, 0, 0, flag);
	mnparm(1, "a", 0.05, 0.001, 0, 0, flag);

	// Chi2: 1.
	// -logl: 0.5
	SetErrorDef(0.5);
	SetFCN(minFunctionStatic);
}

MinuitFitter::~MinuitFitter() {
}

void MinuitFitter::fit() {
	int flag;
	double args[1];

	//Run MINUIT
	gTemp = this;
	args[0] = 100000;
	mnexcm("MIGRAD", args, 1, flag);

	//Get MINUIT results
	double edm, errdef;
	int nvpar, nparx, icstat;

	GetParameter(0, fNorm, fNormErr);
	GetParameter(1, fFormFactor, fFormFactorErr);
	mnstat(fMinimum, edm, errdef, nvpar, nparx, icstat);

	cout << fMinimum << " " << edm << " " << errdef << " " << nvpar << " "
			<< nparx << " " << icstat << endl;
}

void MinuitFitter::minFunctionStatic(Int_t&, Double_t*, Double_t& f,
		Double_t* par, Int_t) {
	f = gTemp->minFunction(par[0], par[1]);
}

void MinuitFitter::printResult() {
	double chi2Prob = TMath::Prob(fMinimum, 50 - 2);
	cout << "######## Procedure result #########" << endl
			<< "-------------------------------------" << endl;
	cout << "Global normalization : " << fNorm << "+-" << fNormErr << endl;
	cout << "Slope a : " << fFormFactor << "+-" << fFormFactorErr << endl;
	cout << "Chi2 : " << fMinimum << " prob : " << chi2Prob << " p-value : "
			<< endl;
}

void MinuitFitter::drawResult(vector<Sample*> mcSamples, int nbins, double *binning) {
	FitResultDrawer drawer;
	drawer.setTitle(fName);

	for (auto sample : mcSamples) {
		static_cast<FitMCSample*>(sample)->populateFit(&drawer, fNorm, fFormFactor);
	}

	fDataSamples->populateFit(&drawer, fNorm, fFormFactor);

	drawer.setBinning(nbins, binning);
	drawer.setMc(fMCSamples->getMainHisto(), fDataSamples->getMainHisto());
	drawer.setResult(fNorm, fFormFactor, fFormFactorErr);
	drawer.draw();
}

double MinuitFitter::chi2pValue(double chi2, int ndof) {
	double step = 0.0001;
	double currChi2 = 0;
	double integral = 0;
	//cout << currChi2 << " " << chi2 << endl;
	while (currChi2 < chi2) {
		integral += (1 - TMath::Prob(currChi2, ndof)) * step;
		currChi2 += step;
	}

	return 1 - integral;
}


void MinuitFitter::chi2Profile(TString) {
//	double chi2Min, chi2Max;
//	double step;
//	double chi2 = fFormFactor;
//	chi2Min = chi2 * 0.8;
//	chi2Max = chi2 * 1.2;
//	step = (chi2Max - chi2Min) / 100.;
//
//	TGraph* profile = new TGraph();
//	profile->SetTitle("Chi2 Profile");
//	profile->SetName("chi2Profile" + name);
//
//	int npar = 3;
//	double params[3] = { fNorm, 0, 1 };
//	double chi2Value;
//	for (int i = 0; i < 100; i++) {
//		params[1] = chi2Min + i * step;
//		//MinuitFitter::minFct(npar, NULL, chi2Value, params, 0);
//		profile->SetPoint(i, chi2Min + i * step, chi2Value);
//	}
//
////	new TCanvas("chi2Profile" + name, "chi2Profile");
//	profile->Draw("A*");
}
