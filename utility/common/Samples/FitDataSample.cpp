/*
 * FitDataSample.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#include "FitDataSample.h"
#include "exportClasses.h"
#include "ScanCuts.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include "Functions.h"

using namespace std;
FitDataSample::FitDataSample() :
		dSig(nullptr) {

}

FitDataSample::~FitDataSample() {
	// TODO Auto-generated destructor stub
}

void FitDataSample::processEvent(ROOTPhysicsEvent *eventBrch,
		ROOTBurst *, ROOTRawEvent *rawBrch,
		ROOTCorrectedEvent *corrBrch, ROOTFileHeader *,
		ROOTMCEvent *mcEvent, NGeom *geomBrch, std::vector<bool> *cutsPass,
		const ConfigFile *, const RunWeights *) {

	//Read event and fill histo
	double x, xTrue = -1;
	double weight, bweight = 1.;

	double a = fTestA;
	bool passNormal;

	passNormal = true;
	if (cutsPass) {
		if (!cutsPass->at(fScanID)) {
			fFitBrch.selEvents--;
			return;
		}
	}

	if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch,
			fFitBrch))
		return;

	x = eventBrch->x;
	if (mcEvent)
		xTrue = mcEvent->xTrue;
	weight = 1.;	//+2.*a*x+a*a*x*x;
	if (a != 0)
		bweight = (1. + 2. * a * xTrue + a * a * xTrue * xTrue)
				/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
	if (passNormal)
		dSig->Fill(x, weight * bweight);

}

void FitDataSample::doGet(TFile* inputFD, TFile* tempFD) {
	fitStruct fitBrch;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fitBrch);

	initFitStruct(fFitBrch);
	sumTreeFitStruct(fitBrch, t, fFitBrch, fFactor);

	//Create histo
	TH1D* xxx = (TH1D*) inputFD->Get("sig");
	tempFD->cd();
	dSig = ((TH1D*) xxx->Clone());
	dSig->SetName(TString::Format("sig_%i", fIndex));

	//sig->Add(dSig->at(index), factor);
}

void FitDataSample::doWrite() {
	dSig->Write();
}

void FitDataSample::doSetName() {
	//Nothing to do for data
}

void FitDataSample::initHisto(int nbins, double* bins, const ConfigFile *cfg) {
	if (cfg->isWithEqualBins())
		dSig = new TH1D("sig", "signal sample", nbins - 1, bins);
	else
		dSig = new TH1D("sig", "signal sample", NBINS, 0, MAXBIN);
}

void FitDataSample::scale() {
}

void FitDataSample::setPlotStyle(std::vector<int>) {
	dSig->SetLineColor(kRed);
}

void FitDataSample::populateStack(HistoDrawer *drawer, string legend) {
	InputFitDrawer* myDrawer = static_cast<InputFitDrawer*>(drawer);
	myDrawer->fSig->Add((TH1D*) dSig->Clone());
	myDrawer->fLegSig->AddEntry(dSig, legend.c_str());
}

FitDataSample::bContent FitDataSample::getBinContent(int bin) {
	bContent b;
	b.dSig = dSig->GetBinContent(bin);
	return b;
}

FitDataSample* FitDataSample::Add(const FitDataSample* other) {
	SubSample::Add((SubSample*) other);

	dSig->Add(other->dSig, other->fFactor);
	return this;
}
void FitDataSample::populateFit(HistoDrawer *drawer, double, double, string legend) {
	FitResultDrawer *myDrawer = static_cast<FitResultDrawer*>(drawer);

	TH1D *dSig_c = (TH1D*) dSig->Clone(
			TString::Format("dSig_c%s%i", drawer->getTitle().c_str(), fIndex));
	myDrawer->fLegFitSig->AddEntry(dSig_c, legend.c_str());
	myDrawer->fFitSig->Add(dSig_c);
}

double FitDataSample::getFFIntegral(double) {
	return dSig->Integral();
}
