/*
 * MinuitFitterNewChi2.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#include "../Fitter/MinuitFitterNewChi2.h"
#include <iostream>

#define PRINTVAR(v) #v << "= " << v << " "
using namespace std;

MinuitFitterNewChi2::MinuitFitterNewChi2(int nbins) :
		fNBins(nbins) {
}

MinuitFitterNewChi2::~MinuitFitterNewChi2() {
}

void MinuitFitterNewChi2::fit() {
	mnfree(0);
	MinuitFitter::fit();
}

double MinuitFitterNewChi2::minFunction(double n, double a) {
	double chi2 = 0.;
	double M_i, D_i, m_i;
	double a_i, b_i, g_i;
	double sigma2;

	for (int i = 0; i <= fNBins; ++i) {
		M_i = 0;
		a_i = 0;
		b_i = 0;
		g_i = 0;
		FitMCSample::bContent_t b = fMCSamples->getBinContent(i);
		M_i += b.dNew;
		a_i += b.dAlpha;
		b_i += b.dBeta;
		g_i += b.dGamma;

		D_i = 0;
		D_i += fDataSamples->getBinContent(i).dSig;

		m_i = function(n, a, a_i, b_i, g_i);

		sigma2 = D_i * D_i / M_i + D_i;
		if (D_i != 0)
			chi2 += pow(D_i - m_i, 2) / sigma2;
	}

	return chi2;
}

double MinuitFitterNewChi2::function(double G, double a, double a_i, double b_i,
		double g_i) {
	return G * (a_i + 2. * a * b_i + a * a * g_i);
}
