/*
 * Fitter.cpp
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#include "../Fitter/Fitter.h"

#include <TString.h>
#include <iostream>

#include "MinuitFitterNewChi2.h"
#include "MinuitFitterNewROOT.h"

using namespace std;

Fitter::Fitter(){
}

Fitter::~Fitter() {

}

void Fitter::PrepareHistos(vector<int> allColors, vector<int> dataColors) {
	//Color histos
	InputFitDrawer drawer;
	for (unsigned int i = 0; i < fMCSamples.size(); i++) {
		vector<int> colors(allColors[i * 3], allColors[i * 3 + 2]);
		fMCSamples[i]->setPlotStyle(colors);
		fMCSamples[i]->populateStack(&drawer);
	}

	for (unsigned int i = 0; i < fDataSamples.size(); i++) {
		vector<int> colors(dataColors[i * 3], dataColors[i * 3 + 2]);
		fDataSamples[i]->setPlotStyle(colors);
		fDataSamples[i]->populateStack(&drawer);
	}

	drawer.draw();
}

void Fitter::fit(bool, bool useROOT) {
	MinuitFitter *minuit;
	if (useROOT) {
		minuit = new MinuitFitterNewROOT(fNBins);
		static_cast<MinuitFitterNewROOT*>(minuit)->init(fBinning);
		minuit->SetName("NewROOT");
	} else {
		minuit = new MinuitFitterNewChi2(fNBins);
		minuit->SetName("NewChi2");
	}
	minuit->setSamples(static_cast<FitMCSample*>(fFinalMCSample), static_cast<FitDataSample*>(fFinalDataSample));
	minuit->fit();

	minuit->printResult();
	minuit->drawResult(fMCSamples, fNBins, fBinning);
}

double Fitter::getNormalization(double a) {
	double G = 0;
	for (auto sample : fMCSamples)
		G += sample->getFFIntegral(a);
	double NSig = 0;
	for (auto sample : fDataSamples)
		NSig += sample->getFFIntegral(a);

	return NSig / G;
}
