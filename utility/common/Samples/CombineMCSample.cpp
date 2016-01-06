/*
 * CombineMCSample.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#include "CombineMCSample.h"
#include <TFile.h>
#include <iostream>

using namespace std;

CombineMCSample::CombineMCSample() {
	// TODO Auto-generated constructor stub

}

CombineMCSample::~CombineMCSample() {
	// TODO Auto-generated destructor stub
}

CombineMCSample::CombineMCSample(int index, ConfigFile* cfg) :
		CombineSample(index, cfg) {
}

void CombineMCSample::fillHisto(ROOTPhysicsEvent* evt, ROOTRawEvent* rawEvt,
		ROOTCorrectedEvent* corrEvent, ROOTMCEvent* mcEvent, NGeom* rootGeom,
		ROOTBurst* rootBurst) {

	double weight = fWeights->applyWeights(rootBurst->nrun) * corrEvent->weight;
	CombineSample::fillHisto(evt, rawEvt, corrEvent, mcEvent, rootGeom,
			rootBurst, weight);
}

void CombineMCSample::doGet(TFile* inputFD, TFile* tempFD) {
	fitStruct totFit;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fFitBrch);

	initFitStruct(totFit);
	sumTreeFitStruct(fFitBrch, t, totFit, 1);

	doGetHisto(inputFD, tempFD);

	scale();
}

void CombineMCSample::scaleToData(bContent totalMC, double nData) {
	//Scale MC to Data
	cout << "totalMC " << totalMC.val << " nData" << nData << endl;
	double factor = nData/ totalMC.val;

	cout << d1[0]->Integral() << endl;
	for(auto plot : d1)
		plot->Scale(factor);
	for(auto plot : dMap)
		plot->Scale(factor);
	cout << d1[0]->Integral() << endl;
}

CombineMCSample::bContent CombineMCSample::getIntegrals() {
	bContent b;
	b.val = d1[0]->Integral();
	return b;
}
