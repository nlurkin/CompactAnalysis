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
#include "funLib.h"

extern double Mpi0;

using namespace std;
FitDataSample::FitDataSample() :
		dSig(nullptr), dXDiff(nullptr) {
}

FitDataSample::~FitDataSample() {
	// TODO Auto-generated destructor stub
}

void FitDataSample::processEvent(ROOTPhysicsEvent *eventBrch,
		ROOTBurst *burstBrch, ROOTRawEvent *rawBrch,
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

	if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch,burstBrch,
			fFitBrch))
		return;

	x = eventBrch->x;

	//double mpi0 = eventBrch->pi0.P.M();
	//double factor = Mpi0/mpi0;

	//TLorentzVector ep = eventBrch->ep.P,em = eventBrch->em.P, gamma = eventBrch->gamma.P;
	//ep *= factor;
	//em *= factor;
	//gamma *= factor;

	//double nmpi0 = (ep+em+gamma).M();
	double newx = pow((eventBrch->ep.P*1.01+eventBrch->em.P*1.01).M()/Mpi0, 2.);
	dXDiff->Fill(x-newx);
	//x = newx;
	if (mcEvent)
		xTrue = mcEvent->xTrue;
	weight = 1.;	//+2.*a*x+a*a*x*x;
//	if (a != 0)
//		bweight = (1. + 2. * a * xTrue + a * a * xTrue * xTrue)
//				/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
	if (passNormal)
		dSig->Fill(x, weight * bweight);

}

void FitDataSample::doGet(TDirectory* inputFD, TFile* tempFD) {
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
	dXDiff->Write();
}

void FitDataSample::doSetName() {
	//Nothing to do for data
}

void FitDataSample::initHisto(int nbins, double* bins, const ConfigFile *cfg) {
	if (cfg->isWithEqualBins())
		dSig = new TH1D("sig", "signal sample", nbins - 1, bins);
	else
		dSig = new TH1D("sig", "signal sample", NBINS, 0, MAXBIN);
	dXDiff = new TH1D("diffx", "diffx", 1000, -1, 1);
}

void FitDataSample::scale() {
}

void FitDataSample::setPlotStyle(std::vector<int>) {
	dSig->SetLineColor(kRed);
}

FitDataSample::bContent FitDataSample::getBinContent(int bin) {
	bContent b;
	b.dSig = dSig->GetBinContent(bin);
	return b;
}

SubSample* FitDataSample::Add(const SubSample* other) {
	SubSample::Add(other);

	const FitDataSample *myOther = static_cast<const FitDataSample*>(other);
	dSig->Add(myOther->dSig, myOther->fFactor);
	return this;
}

double FitDataSample::getFFIntegral(double) {
	return dSig->Integral();
}
