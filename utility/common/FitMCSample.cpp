/*
 * FitSample.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#include "exportClasses.h"
#include "ScanCuts.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "FitMCSample.h"
#include <iomanip>
#include "Functions.h"

using namespace std;

FitMCSample::FitMCSample(int index, ConfigFile &cfg) :
	Sample(index, cfg)
{

}

FitMCSample::~FitMCSample() {
}

void FitMCSample::doFill(TFile* inputFD, TFile* tempFD) {
	//Get the TTree
	//Input
	int scanID;
	ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
	ROOTBurst *burstBrch = new ROOTBurst();
	ROOTRawEvent *rawBrch = new ROOTRawEvent();
	ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
	ROOTFileHeader *headerBrch = new ROOTFileHeader();
	ROOTMCEvent *mcEvent = 0;
	NGeom *geomBrch = new NGeom();
	vector<bool> *cutsPass = 0;
	ScanCuts *cutsLists = 0;

	TTree *t = (TTree*) inputFD->Get("event");
	TTree *th = (TTree*) inputFD->Get("header");
	TTree *tc = (TTree*) inputFD->Get("cutsDefinition");
	if (t->GetListOfBranches()->Contains("mc"))
		mcEvent = new ROOTMCEvent();
	if (t->GetListOfBranches()->Contains("cutsResult")) {
		cutsPass = new vector<bool>;
		cutsLists = new ScanCuts();
	}

	t->SetBranchAddress("pi0dEvent", &eventBrch);
	t->SetBranchAddress("rawBurst", &burstBrch);
	t->SetBranchAddress("rawEvent", &rawBrch);
	t->SetBranchAddress("corrEvent", &corrBrch);
	th->SetBranchAddress("header", &headerBrch);
	th->SetBranchAddress("geom", &geomBrch);
	if (mcEvent)
		t->SetBranchAddress("mc", &mcEvent);
	if (cutsPass) {
		t->SetBranchAddress("cutsResult", &cutsPass);
		tc->SetBranchAddress("lists", &cutsLists);
		tc->GetEntry(0);
		if (fCfg.getScanId()== -1)
			scanID = cutsLists->getDefaultIndex();
		else
			scanID = fCfg.getScanId();
		tc->GetEntry(scanID);
		cutsLists->Cuts::print();
	}

	cout << "mc=" << mcEvent << endl;

	tempFD->cd();
	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = 0;
	for (int i = 0; i < th->GetEntries(); i++) {
		th->GetEntry(i);
		totalChanEvents += headerBrch->NProcessedEvents;
	}
	int processedEvents = 0;

	int divider = 5;

	//Read events and fill histo
	int i = 0;
	double x, xTrue = -1;
	double bweight = 1.;
	double weight, aweight;
	int mod;

	fFitBrch.totEvents += totalChanEvents;
	fFitBrch.selEvents += nevt;

	bool passNormal;
	cout << "Filling " << nevt << endl;
	for (; i < nevt; ++i) {
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
				passNormal = false;
			}
		}
		if (!fCfg.testUseRun(burstBrch->nrun, burstBrch->period))
			continue;

		if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch))
			continue;
		weight = fWeights->applyWeights(burstBrch->nrun) * corrBrch->weight;

		x = eventBrch->x;
		if (mcEvent)
			xTrue = mcEvent->xTrue;
		bweight = 1.
				/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
		mod = rawBrch->timeStamp % divider;

		aweight = 1.;

		if (passNormal) {
			dNew->Fill(x, bweight * aweight * weight);
			dAlpha->Fill(x, 1 / pow(1 + 0.032 * xTrue, 2.));
			dBeta->Fill(x, xTrue / pow(1 + 0.032 * xTrue, 2.));
			dGamma->Fill(x, pow(xTrue / (1 + 0.032 * xTrue), 2.));
		}
		if (mod == 0 || mod == 1 || mod == 2) {
			//TODO to check
			aweight = 1.;
			if (passNormal) {
				fFitBrch.n1++;
				d1->Fill(x, bweight * aweight * weight);
			}
		} else if (mod == 3) {
			//TODO to check
			aweight = xTrue;
			if (passNormal) {
				fFitBrch.nx++;
				d2->Fill(x, bweight * aweight * weight);
			}
		} else if (mod == 4) {
			//TODO to check
			aweight = xTrue * xTrue;
			if (passNormal) {
				fFitBrch.nxx++;
				d3->Fill(x, bweight * aweight * weight);
			}
		}
		processedEvents++;
	}

	cout << endl << fFitBrch.selEvents << endl;
}

void FitMCSample::doWrite() {
	d1->Write();
	d2->Write();
	d3->Write();
	dNew->Write();
	dAlpha->Write();
	dBeta->Write();
	dGamma->Write();
}

void FitMCSample::doSetName() {
//	d1->SetName(TString::Format("d1_%i", fIndex));
//	d2->SetName(TString::Format("d2_%i", fIndex));
//	d3->SetName(TString::Format("d3_%i", fIndex));
//	dNew->SetName(TString::Format("dNew_%i", fIndex));
//	dAlpha->SetName(TString::Format("dAlpha_%i", fIndex));
//	dBeta->SetName(TString::Format("dBeta_%i", fIndex));
//	dGamma->SetName(TString::Format("dGamma_%i", fIndex));
}

void FitMCSample::initHisto(int nbins, double* bins){
	if (fCfg.isWithEqualBins()) {
		d1 = new TH1D("d1", "sample 1", nbins - 1, bins);
		d1->Sumw2();
		d2 = new TH1D("d2", "sample x", nbins - 1, bins);
		d2->Sumw2();
		d3 = new TH1D("d3", "sample x^{2}", nbins - 1, bins);
		d3->Sumw2();
		dNew = new TH1D("dNew", "MC", nbins - 1, bins);
		dNew->Sumw2();
		dAlpha = new TH1D("dAlpha", "MC", nbins - 1, bins);
		dAlpha->Sumw2();
		dBeta = new TH1D("dBeta", "MC", nbins - 1, bins);
		dBeta->Sumw2();
		dGamma = new TH1D("dGamma", "MC", nbins - 1, bins);
		dGamma->Sumw2();
	}
	else{
		d1 = new TH1D("d1", "sample 1", NBINS, 0, MAXBIN);
		d1->Sumw2();
		d2 = new TH1D("d2", "sample x", NBINS, 0, MAXBIN);
		d2->Sumw2();
		d3 = new TH1D("d3", "sample x^{2}", NBINS, 0, MAXBIN);
		d3->Sumw2();
		dNew = new TH1D("dNew", "MC", NBINS, 0, MAXBIN);
		dNew->Sumw2();
		dAlpha = new TH1D("dAlpha", "MC", NBINS, 0, MAXBIN);
		dAlpha->Sumw2();
		dBeta = new TH1D("dBeta", "MC", NBINS, 0, MAXBIN);
		dBeta->Sumw2();
		dGamma = new TH1D("dGamma", "MC", NBINS, 0, MAXBIN);
		dGamma->Sumw2();
	}
}
