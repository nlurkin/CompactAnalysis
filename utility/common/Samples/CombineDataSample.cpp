/*
 * CombineDataSample.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#include "CombineDataSample.h"
#include "../Drawer/CombineDrawer.h"
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

void CombineDataSample::doGet(TFile* inputFD, TFile* tempFD) {
	fitStruct totFit;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fFitBrch);

	initFitStruct(totFit);
	sumTreeFitStruct(fFitBrch, t, totFit, fFactor);

	doGetHisto(inputFD, tempFD);
}

void CombineDataSample::populateStack(HistoDrawer *drawer, string legend) {
	CombineDrawer *myDrawer = static_cast<CombineDrawer*>(drawer);

	myDrawer->addLegendData(d1[0], legend);
	for (unsigned int i = 0; i < d1.size(); ++i) {
		myDrawer->addHistoData(i, d1[i]);
	}
}

void CombineDataSample::setPlotStyle(std::vector<int>) {
	for(auto plot : d1)
		plot->SetLineColor(kRed);

	for(auto plot : dMap)
		plot->SetLineColor(kRed);
}
