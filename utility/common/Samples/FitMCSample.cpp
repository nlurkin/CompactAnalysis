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
#include <TStyle.h>

using namespace std;

FitMCSample::FitMCSample() :
		d1(nullptr), d2(nullptr), d3(nullptr), dNew(nullptr), dAlpha(nullptr), dBeta(
				nullptr), dGamma(nullptr) {

}

FitMCSample::FitMCSample(int index, ConfigFile *cfg) :
		Sample(index, cfg), d1(nullptr), d2(nullptr), d3(nullptr), dNew(
				nullptr), dAlpha(nullptr), dBeta(nullptr), dGamma(nullptr) {

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
		if (fCfg->getScanId() == -1)
			scanID = cutsLists->getDefaultIndex();
		else
			scanID = fCfg->getScanId();
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
		if (!fCfg->testUseRun(burstBrch->nrun, burstBrch->period))
			continue;

		if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch, fFitBrch))
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

void FitMCSample::doGet(TFile* inputFD, TFile* tempFD) {
	fitStruct fitBrch;
	TTree *t = (TTree*) inputFD->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fitBrch);

	initFitStruct(fFitBrch);
	sumTreeFitStruct(fitBrch, t, fFitBrch, 1);

	//Create histo
	TH1D* xxx1 = (TH1D*) inputFD->Get("d1");
	TH1D* xxx2 = (TH1D*) inputFD->Get("d2");
	TH1D* xxx3 = (TH1D*) inputFD->Get("d3");
	TH1D* xxx4 = (TH1D*) inputFD->Get("dNew");
	TH1D* xxxA = (TH1D*) inputFD->Get("dAlpha");
	TH1D* xxxB = (TH1D*) inputFD->Get("dBeta");
	TH1D* xxxG = (TH1D*) inputFD->Get("dGamma");

	tempFD->cd();
	d1 = ((TH1D*) xxx1->Clone());
	d2 = ((TH1D*) xxx2->Clone());
	d3 = ((TH1D*) xxx3->Clone());
	dNew = ((TH1D*) xxx4->Clone());
	dAlpha = ((TH1D*) xxxA->Clone());
	dBeta = ((TH1D*) xxxB->Clone());
	dGamma = ((TH1D*) xxxG->Clone());
	d1->SetName(TString::Format("d1_%i", fIndex));
	d2->SetName(TString::Format("d2_%i", fIndex));
	d3->SetName(TString::Format("d3_%i", fIndex));
	dNew->SetName(TString::Format("dNew_%i", fIndex));
	dAlpha->SetName(TString::Format("dAlpha_%i", fIndex));
	dBeta->SetName(TString::Format("dBeta_%i", fIndex));
	dGamma->SetName(TString::Format("dGamma_%i", fIndex));

	scale();
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

void FitMCSample::initHisto(int nbins, double* bins) {
	if (fCfg->isWithEqualBins()) {
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
	} else {
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

void FitMCSample::scale() {
	// Rescale histo
	double totN = fFitBrch.n1 + fFitBrch.nx + fFitBrch.nxx;
	double selRatio1 = (double) fFitBrch.n1 / totN;
	double selRatiox = (double) fFitBrch.nx / totN;
	double selRatioxx = (double) fFitBrch.nxx / totN;

	cout << "NBefore scaling" << dAlpha->Integral() << endl;
	Sample::scale(d1, selRatio1);
	Sample::scale(d2, selRatiox);
	Sample::scale(d3, selRatioxx);

	Sample::scale(dNew, 1.);
	Sample::scale(dAlpha, 1.);
	Sample::scale(dBeta, 1.);
	Sample::scale(dGamma, 1.);
	cout << "Nafter scaling" << dAlpha->Integral() << endl;
}

void FitMCSample::scaleToData(bContent totalMC, double nData) {
	//Scale MC to Data
	double totalMCNorm = totalMC.d1;
	double totalMCNew = totalMC.dNew;
	double totalGreek = totalMC.dAlpha;

	double factor = ((double) (nData)) / totalMCNorm;
	double factorNew = ((double) (nData)) / totalMCNew;
	double factorGreek = ((double) (nData)) / totalGreek;

	d1->Scale(factor);
	d2->Scale(factor);
	d3->Scale(factor);
	dNew->Scale(factorNew);
	dAlpha->Scale(factorGreek);
	dBeta->Scale(factorGreek);
	dGamma->Scale(factorGreek);
}

void FitMCSample::setPlotStyle(vector<int> color) {
	d1->SetFillColor(gStyle->GetColorPalette(color[0]));
	d2->SetFillColor(gStyle->GetColorPalette(color[1]));
	d3->SetFillColor(gStyle->GetColorPalette(color[2]));
	dNew->SetFillColor(gStyle->GetColorPalette(color[0]));
	dAlpha->SetFillColor(gStyle->GetColorPalette(color[0]));
	dBeta->SetFillColor(gStyle->GetColorPalette(color[1]));
	dGamma->SetFillColor(gStyle->GetColorPalette(color[2]));
}

void FitMCSample::populateStack(HistoDrawer *drawer) {
	InputFitDrawer *myDrawer = static_cast<InputFitDrawer*>(drawer);

	myDrawer->fStd1->Add((TH1D*) d1->Clone());
	myDrawer->fStdx->Add((TH1D*) d2->Clone());
	myDrawer->fStdxx->Add((TH1D*) d3->Clone());
	myDrawer->fStdNew->Add((TH1D*) dNew->Clone());
	myDrawer->fStdAlpha->Add((TH1D*) dAlpha->Clone());
	myDrawer->fStdBeta->Add((TH1D*) dBeta->Clone());
	myDrawer->fStdGamma->Add((TH1D*) dGamma->Clone());
	myDrawer->fLeg1->AddEntry(d1, fLegend.c_str());
	myDrawer->fLegx->AddEntry(d2, fLegend.c_str());
	myDrawer->fLegxx->AddEntry(d3, fLegend.c_str());
	myDrawer->fLegNew->AddEntry(dNew, fLegend.c_str());
	myDrawer->fLegAlpha->AddEntry(dAlpha, fLegend.c_str());
	myDrawer->fLegBeta->AddEntry(dBeta, fLegend.c_str());
	myDrawer->fLegGamma->AddEntry(dGamma, fLegend.c_str());
}

FitMCSample::bContent FitMCSample::getBinContent(int bin) {
	bContent b;
	b.d1 = d1->GetBinContent(bin);
	b.d2 = d2->GetBinContent(bin);
	b.d3 = d3->GetBinContent(bin);
	b.dNew = dNew->GetBinContent(bin);
	b.dAlpha = dAlpha->GetBinContent(bin);
	b.dBeta = dBeta->GetBinContent(bin);
	b.dGamma = dGamma->GetBinContent(bin);

	return b;
}

void FitMCSample::populateFit(HistoDrawer *drawer, double norm, double a) {
	FitResultDrawer *myDrawer = static_cast<FitResultDrawer*>(drawer);

	TH1D *dAlpha_c = (TH1D*) dAlpha->Clone(
			TString::Format("dAlpha_c%s%i", drawer->getTitle().c_str(), fIndex));
	myDrawer->fLegFit->AddEntry(dAlpha_c,
			TString::Format("%s #alpha", fLegend.c_str()));
	dAlpha_c->Scale(norm);
	myDrawer->fFit->Add(dAlpha_c);
	TH1D *dBeta_c = (TH1D*) dBeta->Clone(
			TString::Format("dBeta_c%s%i", drawer->getTitle().c_str(), fIndex));
	myDrawer->fLegFit->AddEntry(dBeta_c,
			TString::Format("%s #beta", fLegend.c_str()));
	dBeta_c->Scale(norm * 2. * a);
	myDrawer->fFit->Add(dBeta_c);
	TH1D *dGamma_c = (TH1D*) dGamma->Clone(
			TString::Format("dGamma_c%s%i", drawer->getTitle().c_str(), fIndex));
	myDrawer->fLegFit->AddEntry(dGamma_c,
			TString::Format("%s #gamma", fLegend.c_str()));
	dGamma_c->Scale(norm * a * a);
	myDrawer->fFit->Add(dGamma_c);
}

FitMCSample* FitMCSample::Add(const FitMCSample* other) {
	Sample::Add((Sample*) other);
	d1->Add(other->d1, 1.);
	d2->Add(other->d2, 1.);
	d3->Add(other->d3, 1.);
	dNew->Add(other->dNew, 1.);

	dAlpha->Add(other->dAlpha, 1.);
	dBeta->Add(other->dBeta, 1.);
	dGamma->Add(other->dGamma, 1.);

	return this;
}

FitMCSample::bContent FitMCSample::getIntegrals() {
	bContent b;
	b.d1 = d1->Integral();
	b.d2 = d2->Integral();
	b.d3 = d3->Integral();
	b.dNew = dNew->Integral();
	b.dAlpha = dAlpha->Integral();
	b.dBeta = dBeta->Integral();
	b.dGamma = dGamma->Integral();

	return b;
}

double FitMCSample::getFFIntegral(double a) {
	return d1->Integral() * 1.0 + d2->Integral() * a * 2.
					+ d3->Integral() * a * a;
}
