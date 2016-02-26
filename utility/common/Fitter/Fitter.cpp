/*
 * Fitter.cpp
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#include "Fitter.h"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "../Drawer/Drawer.h"
#include "../Interface/Sample.h"
#include "../Samples/FitDataSample.h"
#include "../Samples/FitMCSample.h"
#include "MinuitFitter.h"
#include "MinuitFitterNewChi2.h"
#include "MinuitFitterNewROOT.h"

using namespace std;

Fitter::Fitter() {
}

Fitter::~Fitter() {

}

void Fitter::PrepareHistos(vector<int> allColors, vector<int> dataColors) {
	//Color histos
	for (unsigned int i = 0; i < fMCSamples.size(); i++) {
		vector<int> colors(allColors.begin() + (i * 3),
				allColors.begin() + (i * 3 + 3));
		fMCSamples[i]->setPlotStyle(colors);
	}

	for (unsigned int i = 0; i < fDataSamples.size(); i++) {
		vector<int> colors(dataColors.begin() + (i * 3),
				dataColors.begin() + (i * 3 + 3));
		fDataSamples[i]->setPlotStyle(colors);
	}

	//Drawer::drawFitPreparation(fMCSamples, fDataSamples, "Fitter");
}

void Fitter::fit(bool, bool useROOT, double maxLoss, int start, int end) {
	MinuitFitter* minuit;

	if(start==-1) start=0;
	if(end==-1) end = fFinalDataSample->getNSubSample()-1;

	for (int i = 0; i < fFinalDataSample->getNSubSample(); ++i) {
		if (useROOT) {
			minuit = new MinuitFitterNewROOT(fNBins);
			static_cast<MinuitFitterNewROOT*>(minuit)->init(fBinning);
			minuit->setName("NewROOT");
		} else {
			minuit = new MinuitFitterNewChi2(fNBins);
			minuit->setName("NewChi2");
		}
		minuit->setSamples(
				static_cast<FitMCSample*>(fFinalMCSample->getSubSample(i)),
				static_cast<FitDataSample*>(fFinalDataSample->getSubSample(i)));
		minuit->fit();
		fFitters.push_back(minuit);
		if(useROOT)
			static_cast<MinuitFitterNewROOT*>(minuit)->clean();
	}

	int dflt = fMCSamples[0]->getMainSubSample();
	cout << dflt << endl;
	vector<int> usedIndex = findWithEnoughStat(fFinalDataSample, dflt, start, end, maxLoss);
	int max = findMaxSample(fFitters, usedIndex);
	int min = findMinSample(fFitters, usedIndex);
	double average = computeAverage(fFitters, usedIndex);

	int color;
	for (unsigned int i = 0; i < fFitters.size(); ++i) {
		if (i == dflt)
			color = 1;
		else if (i == max)
			color = 2;
		else if (i == min)
			color = 3;
		else if (std::find(usedIndex.begin(), usedIndex.end(), i)
				== usedIndex.end())
			color = 4;
		else
			color = 0;
		fFitters[i]->printResult(color);
	}

	if (fFitters.size() > 1)
		Drawer::drawFitScan(fFitters, fMCSamples, fDataSamples, fFinalMCSample,
				fFinalDataSample, usedIndex);
	//Drawer::drawFitResult(vFitter, fMCSamples, fDataSamples, fFinalMCSample,
	//		fFinalDataSample, fDataSamples[0]->getMainSubSample());

	double dfltFF = fFitters[dflt]->getFormFactor();
	cout << average << " " << dfltFF << endl;
	cout << "FF difference" << endl;
	cout << "\t" << std::setw(15) << "Average: " << std::fixed
			<< std::setprecision(3) << (average - dfltFF) * 100 << endl;
	cout << "\t" << std::setw(15) << "Delta Average: " << std::fixed
			<< std::setprecision(3)
			<< (average - fFitters[max]->getFormFactor()) * 100 << endl;
	cout << "\t" << std::setw(15) << "Delta Max: " << std::fixed
			<< std::setprecision(3)
			<< (fFitters[max]->getFormFactor() - dfltFF) * 100 << endl;
	cout << "\t" << std::setw(15) << "Delta Min: " << std::fixed
			<< std::setprecision(3)
			<< (fFitters[min]->getFormFactor() - dfltFF) * 100 << endl;

}

double Fitter::getNormalization(double a) {
	double G = 0;
	for (auto sample : fMCSamples)
		G += sample->getSubSample(0)->getFFIntegral(a);
	double NSig = 0;
	for (auto sample : fDataSamples)
		NSig += sample->getSubSample(0)->getFFIntegral(a);

	return NSig / G;
}

int Fitter::findMaxSample(std::vector<MinuitFitter*> vFit, vector<int> use) {
	double max = -99;
	int maxIndex;
	for (auto index : use) {
		if (vFit[index]->getFormFactor() > max) {
			max = vFit[index]->getFormFactor();
			maxIndex = index;
		}
	}
	return maxIndex;
}

int Fitter::findMinSample(std::vector<MinuitFitter*> vFit, vector<int> use) {
	double min = 99;
	int minIndex;
	for (auto index : use) {
		if (vFit[index]->getFormFactor() < min) {
			min = vFit[index]->getFormFactor();
			minIndex = index;
		}
	}
	return minIndex;
}

double Fitter::computeAverage(std::vector<MinuitFitter*> vFit,
		vector<int> use) {
	double average = 0;
	for (auto index : use) {
		average += vFit[index]->getFormFactor();
	}
	return average / use.size();
}

vector<int> Fitter::findWithEnoughStat(Sample* finalDataSample,
		int dfltSample, int start, int end, double maxLoss) {
	vector<int> r;
	int dfltNEvents = finalDataSample->getSubSample(dfltSample)->getSelSize();
	int limitNEvents = (1-maxLoss) * dfltNEvents;
	for (unsigned int i = start; i <= end; ++i) {
		if (finalDataSample->getSubSample(i)->getSelSize() > limitNEvents)
			r.push_back(i);
	}
	return r;

}
