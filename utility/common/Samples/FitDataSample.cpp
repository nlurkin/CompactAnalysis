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
		dSig(nullptr), fTestA(0), fFactor(1) {

}

FitDataSample::FitDataSample(int index, ConfigFile *cfg) :
		Sample(index, cfg), dSig(nullptr), fTestA(0), fFactor(1) {

}

FitDataSample::~FitDataSample() {
	// TODO Auto-generated destructor stub
}

void FitDataSample::doFill(TFile* inputFD, TFile* tempFD) {
	//Input
	int scanID;
	ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
	ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
	ROOTRawEvent *rawBrch = new ROOTRawEvent();
	NGeom *geomBrch = new NGeom();
	ROOTMCEvent *mcEvent = 0;
	vector<bool> *cutsPass = 0;
	ScanCuts *cutsLists = 0;

	TTree *t = (TTree*) inputFD->Get("event");
	TTree *tc = (TTree*) inputFD->Get("cutsDefinition");
	TTree *th = (TTree*) inputFD->Get("header");
	if (t->GetListOfBranches()->Contains("mc"))
		mcEvent = new ROOTMCEvent();
	if (t->GetListOfBranches()->Contains("cutsResult")) {
		cutsPass = new vector<bool>;
		cutsLists = new ScanCuts();
	}

	t->SetBranchAddress("pi0dEvent", &eventBrch);
	t->SetBranchAddress("corrEvent", &corrBrch);
	t->SetBranchAddress("rawEvent", &rawBrch);
	th->SetBranchAddress("geom", &geomBrch);
	if (mcEvent)
		t->SetBranchAddress("mc", &mcEvent);
	if (cutsPass) {
		t->SetBranchAddress("cutsResult", &cutsPass);
		tc->SetBranchAddress("lists", &cutsLists);
		tc->GetEntry(0);
		if (fCfg->getScanId() == -1)
			scanID = cutsLists->getDefaultIndex();
		else
			scanID = fCfg->getScanId();
		tc->GetEntry(scanID);
		cutsLists->Cuts::print();
	}

	th->GetEntry(0);
	tempFD->cd();

	// Set Number of events
	int nevt = t->GetEntries();
	int processedEvents = 0;
//	NSig = nevt;

	fFitBrch.selEvents += nevt;

	//Read event and fill histo
	int i = 0;
	double x, xTrue = -1;
	double weight, bweight = 1.;

	double a = fTestA;
	bool passNormal;

	cout << "Filling data " << fTestA << " with " << nevt << " events" << endl;
	for (i = 0; i < nevt; i++) {
		if (i % 10000 == 0)
			cout << setprecision(2) << i * 100. / (double) nevt << "% " << i
					<< "/" << nevt << "\r";
		cout.flush();
		t->GetEntry(i);
		passNormal = true;
		if (cutsPass) {
			if (!cutsPass->at(scanID)) {
				fFitBrch.selEvents--;
				continue;
			}
		}

		if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch, fFitBrch))
			continue;

		x = eventBrch->x;
		if (mcEvent)
			xTrue = mcEvent->xTrue;
		weight = 1.;	//+2.*a*x+a*a*x*x;
		if (a != 0)
			bweight = (1. + 2. * a * xTrue + a * a * xTrue * xTrue)
					/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
		if (passNormal)
			dSig->Fill(x, weight * bweight);
		processedEvents++;
	}

	cout << endl << fFitBrch.selEvents << endl;
//	return processedEvents;
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

void FitDataSample::initHisto(int nbins, double* bins) {
	if (fCfg->isWithEqualBins())
		dSig = new TH1D("sig", "signal sample", nbins - 1, bins);
	else
		dSig = new TH1D("sig", "signal sample", NBINS, 0, MAXBIN);
}

void FitDataSample::scale() {
}

void FitDataSample::setPlotStyle(std::vector<int>) {
	dSig->SetLineColor(kRed);
}

void FitDataSample::populateStack(InputFitDrawer& drawer) {
	drawer.fSig->Add((TH1D*) dSig->Clone());
	drawer.fLegSig->AddEntry(dSig, fLegend.c_str());
}

FitDataSample::bContent FitDataSample::getBinContent(int bin) {
	bContent b;
	b.dSig = dSig->GetBinContent(bin);
	return  b;
}

FitDataSample& operator +=(FitDataSample& first, const FitDataSample* other) {
	operator +=((Sample&) first, (Sample*) other);

	first.dSig->Add(other->dSig, other->fFactor);
	return first;
}
void FitDataSample::populateFit(FitResultDrawer &drawer, double, double) {
	TH1D *dSig_c = (TH1D*) dSig->Clone(
			TString::Format("dSig_c%s%i", drawer.getTitle().c_str(), fIndex));
	drawer.fLegFitSig->AddEntry(dSig_c, fLegend.c_str());
	drawer.fFitSig->Add(dSig_c);
}

