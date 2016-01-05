/*
 * CombineDataSample.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#include "CombineDataSample.h"
#include "../Drawer/CombineDrawer.h"
#include <TFile.h>

CombineDataSample::CombineDataSample() {
	// TODO Auto-generated constructor stub

}

CombineDataSample::~CombineDataSample() {
	// TODO Auto-generated destructor stub
}

CombineDataSample::CombineDataSample(int index, ConfigFile* cfg) :
		CombineSample(index, cfg) {
}

void CombineDataSample::fillHisto(ROOTPhysicsEvent* evt, ROOTRawEvent* rawEvt,
		ROOTCorrectedEvent* corrEvent, ROOTMCEvent* mcEvent, NGeom* rootGeom,
		ROOTBurst* rootBurst) {

	CombineSample::fillHisto(evt, rawEvt, corrEvent, mcEvent, rootGeom,
			rootBurst, 1.);
}

void CombineDataSample::doGet(TFile* inputFD, TFile* tempFD) {
	fitStruct totFit;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fFitBrch);

	initFitStruct(totFit);
	sumTreeFitStruct(fFitBrch, t, totFit, fFactor);

	doGetHisto(inputFD, tempFD);
}

void CombineDataSample::populateStack(HistoDrawer *drawer) {
	CombineDrawer *myDrawer = static_cast<CombineDrawer*>(drawer);

	for (unsigned int i = 0; i < d1.size(); ++i) {
		myDrawer->addHistoData(i, (TH1D*) (d1[i]->Clone()), fLegend);
	}
}
