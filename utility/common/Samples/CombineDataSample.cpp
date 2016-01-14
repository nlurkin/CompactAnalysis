/*
 * CombineDataSample.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#include "CombineDataSample.h"
#include <TFile.h>
#include <TTree.h>

using namespace std;

CombineDataSample::CombineDataSample() :
	fFactor(1)
{
	// TODO Auto-generated constructor stub

}

CombineDataSample::~CombineDataSample() {
	// TODO Auto-generated destructor stub
}

void CombineDataSample::fillHisto(ROOTPhysicsEvent* evt, ROOTRawEvent* rawEvt,
		ROOTCorrectedEvent* corrEvent, ROOTMCEvent* mcEvent, NGeom* rootGeom,
		ROOTBurst* rootBurst, const RunWeights *) {

	CombineSample::fillHisto(evt, rawEvt, corrEvent, mcEvent, rootGeom,
			rootBurst, 1.);
}

void CombineDataSample::doGet(TDirectory* inputFD, TFile* tempFD) {
	fitStruct totFit;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fFitBrch);

	initFitStruct(totFit);
	sumTreeFitStruct(fFitBrch, t, totFit, fFactor);

	doGetHisto(inputFD, tempFD);
}

void CombineDataSample::setPlotStyle(std::vector<int>) {
	for(auto plot : d1)
		plot->SetLineColor(kRed);

	for(auto plot : dMap)
		plot->SetLineColor(kRed);
}
