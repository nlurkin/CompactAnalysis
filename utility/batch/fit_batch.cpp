#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#define PRINTVAR(v) #v << "= " << v << " "
#define __CINT_NICO__ 1

#include "pi0DalitzHeader.h"
#include <TStyle.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMinuit.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TGraph.h>
#include "exportClasses.h"
#include "ScanCuts.h"
#include <TBrowser.h>
using namespace std;

/*************************
 * Structs
 *************************/
typedef struct fitResult_t {
	double norm;
	double normErr;
	double formFactor;
	double formFactorErr;
} fitResult;

/*************************
 * Globals
 *************************/
#define BINS 10000000
#define MAX 1
#define MAXEVENTS 0
vector<TH1D*> *d1, *d2, *d3, *dSig, *dNew, *dAlpha, *dBeta, *dGamma;
static TH1D *sig; //, *other;
vector<TH1D*> *cut_d1, *cut_d2, *cut_d3, *cut_dSig, *cut_dNew, *cut_dAlpha, *cut_dBeta, *cut_dGamma;
static TH1D *cut_sig; //, *other;
double Mpi0 = 0.1349766;
static double bins[BINS];
//static double equiBins[BINS];
static int nbins;
static double NSig, cut_NSig;

int nmc[2];
int cut_nmc[2];

namespace Fit{void minFct(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);}

/***************************
 * Mandatory from header
 * Opening/Closing files
 ***************************/
void initNewChannel() {
}

void initNewOutput(TFile **fdo, TString fileName) {
	*fdo = TFile::Open(fileName, "RECREATE");
	fitTree = new TTree("fitStruct", "fitStruct tree");
	fitTree->Branch("fitStruct", &fitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	fitBrch.n1 = 0;
	fitBrch.nx = 0;
	fitBrch.nxx = 0;
	fitBrch.selEvents = 0;
	fitBrch.totEvents = 0;

	cut_fitTree = new TTree("cut_fitStruct", "fitStruct tree");
	cut_fitTree->Branch("cut_fitStruct", &cut_fitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	cut_fitBrch.n1 = 0;
	cut_fitBrch.nx = 0;
	cut_fitBrch.nxx = 0;
	cut_fitBrch.selEvents = 0;
	cut_fitBrch.totEvents = 0;

}

void closeMCOutput(TFile *fdo, int index) {
	if (index != -1) {
		fdo->cd();
		fitTree->Fill();
		fitTree->Write();
		cut_fitTree->Fill();
		cut_fitTree->Write();
		d1->at(index)->Write();
		d2->at(index)->Write();
		d3->at(index)->Write();
		dNew->at(index)->Write();
		dAlpha->at(index)->Write();
		dBeta->at(index)->Write();
		dGamma->at(index)->Write();

		cut_d1->at(index)->Write();
		cut_d2->at(index)->Write();
		cut_d3->at(index)->Write();
		cut_dNew->at(index)->Write();
		cut_dAlpha->at(index)->Write();
		cut_dBeta->at(index)->Write();
		cut_dGamma->at(index)->Write();

		fdo->Close();

		d1->at(index)->SetName(TString::Format("d1_%i", index));
		d2->at(index)->SetName(TString::Format("d2_%i", index));
		d3->at(index)->SetName(TString::Format("d3_%i", index));
		dNew->at(index)->SetName(TString::Format("dNew_%i", index));
		dAlpha->at(index)->SetName(TString::Format("dAlpha_%i", index));
		dBeta->at(index)->SetName(TString::Format("dBeta_%i", index));
		dGamma->at(index)->SetName(TString::Format("dGamma_%i", index));

		cut_d1->at(index)->SetName(TString::Format("cut_d1_%i", index));
		cut_d2->at(index)->SetName(TString::Format("cut_d2_%i", index));
		cut_d3->at(index)->SetName(TString::Format("cut_d3_%i", index));
		cut_dNew->at(index)->SetName(TString::Format("cut_dNew_%i", index));
		cut_dAlpha->at(index)->SetName(TString::Format("cut_dAlpha_%i", index));
		cut_dBeta->at(index)->SetName(TString::Format("cut_dBeta_%i", index));
		cut_dGamma->at(index)->SetName(TString::Format("cut_dGamma_%i", index));
	}
}

void closeDataOutput(TFile *fdo, int index) {
	fdo->cd();
	fitTree->Fill();
	fitTree->Write();
	cut_fitTree->Fill();
	cut_fitTree->Write();
	dSig->at(0)->Write();
	cut_dSig->at(0)->Write();

	fdo->Close();
}

/*************************
 * Utility
 *************************/
void scaleMC(fitStruct N, int index, double br) {
	cout << "Rescaling" << endl;

	// Rescale histo
	double totN = N.n1 + N.nx + N.nxx;
	double selRatio1 = (double) N.n1 / totN;
	double selRatiox = (double) N.nx / totN;
	double selRatioxx = (double) N.nxx / totN;

	scale(d1->at(index), selRatio1, N.totEvents, br);
	scale(d2->at(index), selRatiox, N.totEvents, br);
	scale(d3->at(index), selRatioxx, N.totEvents, br);

	scale(dNew->at(index), 1., N.totEvents, br);
	scale(dAlpha->at(index), 1., N.totEvents, br);
	scale(dBeta->at(index), 1., N.totEvents, br);
	scale(dGamma->at(index), 1., N.totEvents, br);
}

void scaleMC_cut(fitStruct N, int index, double br) {
	cout << "Rescaling" << endl;

	// Rescale histo
	double totN = N.n1 + N.nx + N.nxx;
	double selRatio1 = (double) N.n1 / totN;
	double selRatiox = (double) N.nx / totN;
	double selRatioxx = (double) N.nxx / totN;

	scale(cut_d1->at(index), selRatio1, N.totEvents, br);
	scale(cut_d2->at(index), selRatiox, N.totEvents, br);
	scale(cut_d3->at(index), selRatioxx, N.totEvents, br);

	scale(cut_dNew->at(index), 1., N.totEvents, br);
	scale(cut_dAlpha->at(index), 1., N.totEvents, br);
	scale(cut_dBeta->at(index), 1., N.totEvents, br);
	scale(cut_dGamma->at(index), 1., N.totEvents, br);
}

double getNormalization(double a) {
	double G = 0;
	for (int j = 0; j < inputMCNbr; ++j) {
		G += d1->at(j)->Integral() * 1.0 + d2->at(j)->Integral() * a * 2.
				+ d3->at(j)->Integral() * a * a;
	}
	G = ((double) NSig) / G;
	return G;
}

void rebin(int binNumber = 0) {
	if (binNumber == 0) {
		//sig = (TH1D*)sig->Rebin(20, "sig_reb", 0);
		//sig->GetXaxis()->SetRange(sig->GetNbinsX()/5, sig->GetNbinsX());
		double max = sig->GetMaximum();
		//sig->GetXaxis()->SetRange();

		cout << max << endl;
		if(max<100) max = 100;
		double s;
		double sum = 0;
		double oldSum = 0;
		int j = 0;
		double var = 0.01;

		//for(int i=0; i<sig->GetNbinsX()/5; i++){
			//bins[j] = sig->GetBinLowEdge(i+1);
			//++j;
		//}
		//for(int i=sig->GetNbinsX()/5; i<=sig->GetNbinsX(); ++i){
		for (int i = 0; i <= sig->GetNbinsX(); ++i) {
			s = sig->GetBinContent(i);
			oldSum = sum;
			sum += s;
			if ((fabs(sum - max) < var * max)) {
				bins[j] = sig->GetBinLowEdge(i + 1);
				++j;
				sum = 0;
				oldSum = 0;
			} else if (sum > (max * (1. + var))) {
				if (fabs(oldSum - max) < fabs(sum - max)) {
					bins[j] = sig->GetBinLowEdge(i);
					++j;
					sum = s;
					oldSum = 0;
				} else {
					bins[j] = sig->GetBinLowEdge(i + 1);
					++j;
					sum = 0;
					oldSum = 0;
				}
			}
		}
		bins[j] = sig->GetBinLowEdge(sig->GetNbinsX() + 1);
		nbins = j;
	} else {
		int ratio = ceil(sig->GetNbinsX() / (double)binNumber);
		for (int i = 0; i <= binNumber; ++i) {
			bins[i] = sig->GetBinLowEdge(i * ratio + 1);
		}
		nbins = binNumber;
	}


	sig = (TH1D*) sig->Rebin(nbins, "sig_reb", bins);
	cut_sig = (TH1D*) cut_sig->Rebin(nbins, "cut_sig_reb", bins);
	for (int i = 0; i < inputMCNbr; ++i) {
		cut_d1->at(i) = (TH1D*) cut_d1->at(i)->Rebin(nbins,
				"cut_" + TString(cut_d1->at(i)->GetName()) + "_reb", bins);
		cut_d2->at(i) = (TH1D*) cut_d2->at(i)->Rebin(nbins,
				"cut_" + TString(cut_d2->at(i)->GetName()) + "_reb", bins);
		cut_d3->at(i) = (TH1D*) cut_d3->at(i)->Rebin(nbins,
				"cut_" + TString(cut_d3->at(i)->GetName()) + "_reb", bins);
		cut_dNew->at(i) = (TH1D*) cut_dNew->at(i)->Rebin(nbins,
				"cut_" + TString(cut_dNew->at(i)->GetName()) + "_reb", bins);
		cut_dAlpha->at(i) = (TH1D*) cut_dAlpha->at(i)->Rebin(nbins,
				"cut_" + TString(cut_dAlpha->at(i)->GetName()) + "_reb", bins);
		cut_dBeta->at(i) = (TH1D*) cut_dBeta->at(i)->Rebin(nbins,
				"cut_" + TString(cut_dBeta->at(i)->GetName()) + "_reb", bins);
		cut_dGamma->at(i) = (TH1D*) cut_dGamma->at(i)->Rebin(nbins,
				"cut_" + TString(cut_dGamma->at(i)->GetName()) + "_reb", bins);
	}

}

double chi2pValue(double chi2, int ndof) {
	double step = 0.0001;
	double currChi2 = 0;
	double integral = 0;
	//cout << currChi2 << " " << chi2 << endl;
	while (currChi2 < chi2) {
		cout << currChi2 << " " << chi2 << endl;
		integral += (1 - TMath::Prob(currChi2, ndof)) * step;
		currChi2 += step;
	}

	return 1 - integral;
}

void chi2Profile(fitResult r, TString name) {
	double chi2Min, chi2Max;
	double step;
	double chi2 = r.formFactor;
	chi2Min = chi2 * 0.8;
	chi2Max = chi2 * 1.2;
	step = (chi2Max - chi2Min) / 100.;

	TGraph* profile = new TGraph();
	profile->SetTitle("Chi2 Profile");
	profile->SetName("chi2Profile" + name);

	int npar = 3;
	double params[3] = { r.norm, 0, 1 };
	double chi2Value;
	for (int i = 0; i < 100; i++) {
		params[1] = chi2Min + i * step;
		Fit::minFct(npar, NULL, chi2Value, params, 0);
		profile->SetPoint(i, chi2Min + i * step, chi2Value);
	}

	new TCanvas("chi2Profile" + name, "chi2Profile");
	profile->Draw("A*");
}
/************************
 * Fitting
 ************************/

namespace Fit{
	double fun(double G, double a, double b1, double b2, double b3) {
		return G * (b1 + 2. * a * b2 + a * a * b3);
	}

	double funNew(double G, double a, double a_i, double b_i, double g_i) {
		return G * (a_i + 2. * a * b_i + a * a * g_i);
	}

	void minFct(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag) {
		double chi2 = 0.;
		double b1, b2, b3, s;
		double sigma, sig1 = 0, sig2 = 0, sig3 = 0, totSigma;
		double G;
		bool rootMethod = false;

		TH1D *comp, *sigComp;
		//if(par[2]==0) G = getNormalization(par[1]);
		if (par[2] != 0) {
			rootMethod = true;
			comp = new TH1D("comp", "comp", sig->GetNbinsX(), bins);
			sigComp = new TH1D("sigcomp", "sigcomp", sig->GetNbinsX(), bins);
		}
		G = par[0];

		for (int i = 0; i <= sig->GetNbinsX(); ++i) {
			//if (sig->GetBinLowEdge(i + 1) < 0.1) continue;
			b1 = 0;
			b2 = 0;
			b3 = 0;
			for (int j = 0; j < inputMCNbr; ++j) {
				b1 += d1->at(j)->GetBinContent(i);
				b2 += d2->at(j)->GetBinContent(i);
				b3 += d3->at(j)->GetBinContent(i);
				sig1 = d1->at(j)->GetBinError(i);
				sig2 = d2->at(j)->GetBinError(i);
				sig3 = d3->at(j)->GetBinError(i);
			}
			s = sig->GetBinContent(i);
			sigma = sig->GetBinError(i);

			totSigma = sigma * sigma + sig1 * sig1 + sig2 * sig2 + sig3 * sig3;
			//totSigma = fun(G, par[1], b1,b2,b3);

			if (rootMethod){
				comp->Fill(sig->GetBinCenter(i), fun(1, par[1], b1, b2, b3));
				sigComp->SetBinContent(i, s);
			}
			else if (totSigma != 0){
				chi2 += pow((s - fun(G, par[1], b1, b2, b3)), 2.) / totSigma;
			}
		}

		if (rootMethod){
			f = sigComp->Chi2Test(comp, "UW CHI2");
			delete comp;
			delete sigComp;
		}
		else
			f = chi2;
	}

	void minFctNew(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
			Int_t flag) {
		double chi2 = 0.;
		double M_i, D_i, m_i;
		double N = par[0];
		double a = par[1];
		double a_i, b_i, g_i;
		double sigma2;
		bool rootMethod = false;

		TH1D *comp, *sigComp;
		if (par[2] != 0) {
			rootMethod = true;
			comp = new TH1D("comp", "comp", sig->GetNbinsX(), bins);
			sigComp = new TH1D("sigcomp", "sigcomp", sig->GetNbinsX(), bins);
		}
		for (int i = 0; i <= sig->GetNbinsX(); ++i) {
			//if(sig->GetBinLowEdge(i+1)<0.1) continue;
			M_i = 0;
			a_i = 0;
			b_i = 0;
			g_i = 0;
			for (int j = 0; j < inputMCNbr; ++j) {
				M_i += dNew->at(j)->GetBinContent(i);
				a_i += dAlpha->at(j)->GetBinContent(i);
				b_i += dBeta->at(j)->GetBinContent(i);
				g_i += dGamma->at(j)->GetBinContent(i);
			}
			D_i = sig->GetBinContent(i);

			m_i = funNew(N, a, a_i, b_i, g_i);

			//cout << D_i << " " << m_i << " " << sigma2 << endl;
			sigma2 = D_i * D_i / M_i + D_i;
			//totSigma = sigma*sigma + sigMC*sigMC;
			if (rootMethod){
				comp->Fill(sig->GetBinCenter(i), m_i);
				sigComp->SetBinContent(i, D_i);
			}
			else if (D_i != 0)
				chi2 += pow(D_i - m_i, 2) / sigma2;
		}

		if (rootMethod){
			f = sigComp->Chi2Test(comp, "UW CHI2 P");
			delete comp;
			delete sigComp;
		}
		else
			f = chi2;
	}
	void minFctNew_cut(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
			Int_t flag) {
		double chi2 = 0.;
		double M_i, D_i, m_i;
		double N = par[0];
		double a = par[1];
		double a_i, b_i, g_i;
		double sigma2;
		bool rootMethod = false;

		TH1D *comp, *sigComp;
		if (par[2] != 0) {
			rootMethod = true;
			comp = new TH1D("cut_comp", "comp", sig->GetNbinsX(), bins);
			sigComp = new TH1D("cut_sigcomp", "sigcomp", sig->GetNbinsX(), bins);
		}
		for (int i = 0; i <= sig->GetNbinsX(); ++i) {
			//if(sig->GetBinLowEdge(i+1)<0.1) continue;
			M_i = 0;
			a_i = 0;
			b_i = 0;
			g_i = 0;

			for (int j = 0; j < inputMCNbr; ++j) {
				M_i += cut_dNew->at(j)->GetBinContent(i);
				a_i += cut_dAlpha->at(j)->GetBinContent(i);
				b_i += cut_dBeta->at(j)->GetBinContent(i);
				g_i += cut_dGamma->at(j)->GetBinContent(i);
			}
			D_i = cut_sig->GetBinContent(i);

			m_i = funNew(N, a, a_i, b_i, g_i);

			//cout << D_i << " " << m_i << " " << sigma2 << endl;
			sigma2 = D_i * D_i / M_i + D_i;
			//totSigma = sigma*sigma + sigMC*sigMC;
			if (rootMethod){
				comp->Fill(cut_sig->GetBinCenter(i), m_i);
				sigComp->SetBinContent(i, D_i);
			}
			else if (D_i != 0)
				chi2 += pow(D_i - m_i, 2) / sigma2;
		}

		if (rootMethod){
			f = sigComp->Chi2Test(comp, "UW CHI2 P");
			delete comp;
			delete sigComp;
		}
		else
			f = chi2;
	}
}

/***************************
 * Input
 ****************************/
namespace Input {

	void getInputMCFill(TFile *fd, TFile *fdout, double br, unsigned int index) {
		//Get the TTree
		//Input
		int defaultScan = -1;
		ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
		ROOTBurst *burstBrch = new ROOTBurst();
		//ROOTRawEvent *rawBrch = new ROOTRawEvent();
		ROOTRawEvent *rawBrch = xxx;
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTFileHeader *headerBrch = new ROOTFileHeader();
		ROOTMCEvent *mcEvent = 0;
		NGeom *geomBrch = new NGeom();
		vector<bool> *cutsPass = 0;
		ScanCuts *cutsLists = 0;

		TTree *t = (TTree*) fd->Get("event");
		TTree *th = (TTree*)fd->Get("header");
		TTree *tc = (TTree*)fd->Get("cutsDefinition");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();
		if(t->GetListOfBranches()->Contains("cutsResult")){
			cutsPass = new vector<bool>;
			cutsLists = new ScanCuts();
		}

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("rawBurst", &burstBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		th->SetBranchAddress("header", &headerBrch);
		th->SetBranchAddress("geom", &geomBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);
		if(cutsPass){
			t->SetBranchAddress("cutsResult", &cutsPass);
			tc->SetBranchAddress("lists", &cutsLists);
			tc->GetEntry(0);
			if(scanID==-1) scanID = cutsLists->getDefaultIndex();
			defaultScan = cutsLists->getDefaultIndex();
			tc->GetEntry(scanID);
			cutsLists->Cuts::print();
		}

		cout << "mc=" << mcEvent << endl;

		tempFD->cd();
		//Set event nb
		int nevt = t->GetEntries();
		int totalChanEvents = 0;
		for(int i=0; i<th->GetEntries(); i++){
			th->GetEntry(i);
			totalChanEvents += headerBrch->NProcessedEvents;
		}
		int processedEvents = 0;
		nevt = (MAXEVENTS>0) ? min(MAXEVENTS, nevt) : nevt;

		//int npart;
		int divider;
		//npart = nevt/5;
		divider = 5;

		//int n1 = npart*3;
		//int nx = npart;
		//int nxx = npart;

		//TODO to check
		//fitBrch.n1 = n1;
		//fitBrch.nx = nx;
		//fitBrch.nxx = nxx;

		//Create histo
		//int index = d1->size();
		if (index == d1->size()) {
			if(!withEqualBins){
				TH1D* xxx1 = new TH1D("d1", "sample 1", BINS, 0, MAX);
				xxx1->Sumw2();
				d1->push_back(xxx1);
				TH1D* xxx2 = new TH1D("d2", "sample x", BINS, 0, MAX);
				xxx2->Sumw2();
				d2->push_back(xxx2);
				TH1D* xxx3 = new TH1D("d3", "sample x^{2}", BINS, 0, MAX);
				xxx3->Sumw2();
				d3->push_back(xxx3);
				TH1D* xxxNew = new TH1D("dNew", "MC", BINS, 0, MAX);
				xxxNew->Sumw2();
				dNew->push_back(xxxNew);
				TH1D* xxxAlpha = new TH1D("dAlpha", "MC", BINS, 0, MAX);
				xxxAlpha->Sumw2();
				dAlpha->push_back(xxxAlpha);
				TH1D* xxxBeta = new TH1D("dBeta", "MC", BINS, 0, MAX);
				xxxBeta->Sumw2();
				dBeta->push_back(xxxBeta);
				TH1D* xxxGamma = new TH1D("dGamma", "MC", BINS, 0, MAX);
				xxxGamma->Sumw2();
				dGamma->push_back(xxxGamma);

				xxx1 = new TH1D("cut_d1", "sample 1", BINS, 0, MAX);
				xxx1->Sumw2();
				cut_d1->push_back(xxx1);
				xxx2 = new TH1D("cut_d2", "sample x", BINS, 0, MAX);
				xxx2->Sumw2();
				cut_d2->push_back(xxx2);
				xxx3 = new TH1D("cut_d3", "sample x^{2}", BINS, 0, MAX);
				xxx3->Sumw2();
				cut_d3->push_back(xxx3);
				xxxNew = new TH1D("cut_dNew", "MC", BINS, 0, MAX);
				xxxNew->Sumw2();
				cut_dNew->push_back(xxxNew);
				xxxAlpha = new TH1D("cut_dAlpha", "MC", BINS, 0, MAX);
				xxxAlpha->Sumw2();
				cut_dAlpha->push_back(xxxAlpha);
				xxxBeta = new TH1D("cut_dBeta", "MC", BINS, 0, MAX);
				xxxBeta->Sumw2();
				cut_dBeta->push_back(xxxBeta);
				xxxGamma = new TH1D("cut_dGamma", "MC", BINS, 0, MAX);
				xxxGamma->Sumw2();
				cut_dGamma->push_back(xxxGamma);
			}
			else{
				TH1D* xxx1 = new TH1D("d1", "sample 1", nbins-1, bins);
				xxx1->Sumw2();
				d1->push_back(xxx1);
				TH1D* xxx2 = new TH1D("d2", "sample x", nbins-1, bins);
				xxx2->Sumw2();
				d2->push_back(xxx2);
				TH1D* xxx3 = new TH1D("d3", "sample x^{2}", nbins-1, bins);
				xxx3->Sumw2();
				d3->push_back(xxx3);
				TH1D* xxxNew = new TH1D("dNew", "MC", nbins-1, bins);
				xxxNew->Sumw2();
				dNew->push_back(xxxNew);
				TH1D* xxxAlpha = new TH1D("dAlpha", "MC", nbins-1, bins);
				xxxAlpha->Sumw2();
				dAlpha->push_back(xxxAlpha);
				TH1D* xxxBeta = new TH1D("dBeta", "MC", nbins-1, bins);
				xxxBeta->Sumw2();
				dBeta->push_back(xxxBeta);
				TH1D* xxxGamma = new TH1D("dGamma", "MC", nbins-1, bins);
				xxxGamma->Sumw2();
				dGamma->push_back(xxxGamma);

				xxx1 = new TH1D("cut_d1", "sample 1", nbins-1, bins);
				xxx1->Sumw2();
				cut_d1->push_back(xxx1);
				xxx2 = new TH1D("cut_d2", "sample x", nbins-1, bins);
				xxx2->Sumw2();
				cut_d2->push_back(xxx2);
				xxx3 = new TH1D("cut_d3", "sample x^{2}", nbins-1, bins);
				xxx3->Sumw2();
				cut_d3->push_back(xxx3);
				xxxNew = new TH1D("cut_dNew", "MC", nbins-1, bins);
				xxxNew->Sumw2();
				cut_dNew->push_back(xxxNew);
				xxxAlpha = new TH1D("cut_dAlpha", "MC", nbins-1, bins);
				xxxAlpha->Sumw2();
				cut_dAlpha->push_back(xxxAlpha);
				xxxBeta = new TH1D("cut_dBeta", "MC", nbins-1, bins);
				xxxBeta->Sumw2();
				cut_dBeta->push_back(xxxBeta);
				xxxGamma = new TH1D("cut_dGamma", "MC", nbins-1, bins);
				xxxGamma->Sumw2();
				cut_dGamma->push_back(xxxGamma);
			}
		}

		//Read events and fill histo
		int i = 0;
		double x, xTrue=-1;
		double bweight = 1.;
		double weight, aweight;
		int mod;

		fitBrch.totEvents += totalChanEvents;
		fitBrch.selEvents += nevt;

		cut_fitBrch.totEvents += totalChanEvents;
		cut_fitBrch.selEvents += nevt;

		bool passNormal, passCut;
		cout << "Filling " << nevt << endl;
		for (; i < nevt; ++i) {
			if(i % 10000 == 0) cout << setprecision(2) << i*100./(double)nevt << "% " << i << "/" << nevt << "\r";
			cout.flush();
			t->GetEntry(i);
			passNormal = true;
			passCut = false;
			if(cutsPass){
//				cout << "ScanID:" << scanID << " defaultScan:" << defaultScan << endl;
				if(!cutsPass->at(scanID)){
					fitBrch.selEvents--;
					passNormal = false;
//					cout << "Pass default:" << cutsPass->at(defaultScan) << endl;
					if(cutsPass->at(defaultScan))
						passCut = true;
					else
						cut_fitBrch.selEvents--;
				}
				else
					cut_fitBrch.selEvents--;
			}
			if (!runIncluded(burstBrch->nrun, burstBrch->period))
				continue;


			if(!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch))
				continue;
			weight = applyWeights(burstBrch->nrun) * corrBrch->weight;

			x = eventBrch->x;
			if(mcEvent) xTrue = mcEvent->xTrue;
			bweight = 1.
					/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
			mod = rawBrch->timeStamp % divider;

			aweight = 1.;

//			cout << "PassNormal:" << passNormal << " passCut:" << passCut << endl;
			if(passNormal){
				dNew->at(index)->Fill(x, bweight * aweight * weight);
				dAlpha->at(index)->Fill(x, 1 / pow(1 + 0.032 * xTrue, 2.));
				dBeta->at(index)->Fill(x, xTrue / pow(1 + 0.032 * xTrue, 2.));
				dGamma->at(index)->Fill(x, pow(xTrue / (1 + 0.032 * xTrue), 2.));
			}
			if(passCut){
				cut_dNew->at(index)->Fill(x, bweight * aweight * weight);
				cut_dAlpha->at(index)->Fill(x, 1 / pow(1 + 0.032 * xTrue, 2.));
				cut_dBeta->at(index)->Fill(x, xTrue / pow(1 + 0.032 * xTrue, 2.));
				cut_dGamma->at(index)->Fill(x, pow(xTrue / (1 + 0.032 * xTrue), 2.));
			}
			if (mod == 0 || mod == 1 || mod == 2) {
				//TODO to check
				aweight = 1.;
				if(passNormal){
					fitBrch.n1++;
					d1->at(index)->Fill(x, bweight * aweight * weight);
				}
				if(passCut){
					cut_fitBrch.n1++;
					cut_d1->at(index)->Fill(x, bweight * aweight * weight);
				}
			} else if (mod == 3) {
				//TODO to check
				aweight = xTrue;
				if(passNormal){
					fitBrch.nx++;
					d2->at(index)->Fill(x, bweight * aweight * weight);
				}
				if(passCut){
					cut_fitBrch.nx++;
					cut_d2->at(index)->Fill(x, bweight * aweight * weight);
				}
			} else if (mod == 4) {
				//TODO to check
				aweight = xTrue * xTrue;
				if(passNormal){
					fitBrch.nxx++;
					d3->at(index)->Fill(x, bweight * aweight * weight);
				}
				if(passCut){
					cut_fitBrch.nxx++;
					cut_d3->at(index)->Fill(x, bweight * aweight * weight);
				}
			}
			processedEvents++;
		}

		//if(processedEvents != t->GetEntries()){
		//	double ratio = (double)processedEvents/(double)nevt;
		//	totalChanEvents = totalChanEvents*ratio;
		//}

		cout << endl << fitBrch.selEvents << endl;
	}

	int getInputDataFill(TFile *fd, TFile* fdout) {
		//Input
		int defaultScan = -1;
		ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTRawEvent *rawBrch = xxx;
		NGeom *geomBrch = new NGeom();
		ROOTMCEvent *mcEvent = 0;
		vector<bool> *cutsPass = 0;
		ScanCuts *cutsLists = 0;

		TTree *t = (TTree*) fd->Get("event");
		TTree *tc = (TTree*)fd->Get("cutsDefinition");
		TTree *th = (TTree*)fd->Get("header");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();
		if(t->GetListOfBranches()->Contains("cutsResult")){
			cutsPass = new vector<bool>;
			cutsLists = new ScanCuts();
		}

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		th->SetBranchAddress("geom", &geomBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);
		if(cutsPass){
			t->SetBranchAddress("cutsResult", &cutsPass);
			tc->SetBranchAddress("lists", &cutsLists);
			tc->GetEntry(0);
			if(scanID==-1) scanID = cutsLists->getDefaultIndex();
			defaultScan = cutsLists->getDefaultIndex();
			tc->GetEntry(scanID);
			cutsLists->Cuts::print();
		}

		th->GetEntry(0);
		tempFD->cd();

		// Set Number of events
		int nevt = t->GetEntries();
		nevt = (MAXEVENTS>0) ? min(MAXEVENTS, nevt) : nevt;
		int processedEvents = 0;
		NSig = nevt;
		cut_NSig = nevt;

		fitBrch.selEvents += nevt;
		cut_fitBrch.selEvents += nevt;

		//Create histo

		//int index = dSig->size();
		int index = 0;
		cout << "dsigsize " << dSig->size() << endl;
		if (dSig->size() == 0) {
			cout << "WithEqualBins " << withEqualBins << " " << true << " " << false << endl;
			if(!withEqualBins){
				TH1D* xxx = new TH1D("sig", "signal sample", BINS, 0, MAX);
				dSig->push_back(xxx);

				xxx = new TH1D("cut_sig", "signal sample", BINS, 0, MAX);
				cut_dSig->push_back(xxx);
			}
			else{
				TH1D* xxx = new TH1D("sig", "signal sample", nbins-1, bins);
				dSig->push_back(xxx);

				xxx = new TH1D("cut_sig", "signal sample", nbins-1, bins);
				cut_dSig->push_back(xxx);
			}
		}

		//Read event and fill histo
		int i = 0;
		double x, xTrue=-1;
		double weight, bweight = 1.;

		double a = testA;
		bool passNormal, passCut;

		cout << "Filling data " << testA << " with " << nevt << " events" << endl;
		for (i = 0; i < nevt; i++) {
			if(i % 10000 == 0) cout << setprecision(2) << i*100./(double)nevt << "% " << i << "/" << nevt << "\r";
			cout.flush();
			t->GetEntry(i);
			passNormal = true;
			passCut = false;
			if(cutsPass){
				if(!cutsPass->at(scanID)){
					fitBrch.selEvents--;
					passNormal = false;
					if(cutsPass->at(defaultScan))
						passCut = true;
					else
						cut_fitBrch.selEvents--;
				}
				else
					cut_fitBrch.selEvents--;
			}

			if(!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch)) continue;

			x = eventBrch->x;
			if(mcEvent) xTrue = mcEvent->xTrue;
			weight = 1.;	//+2.*a*x+a*a*x*x;
			if (a != 0)
				bweight = (1. + 2. * a * xTrue + a * a * xTrue * xTrue)
						/ (1. + 2. * 0.032 * xTrue + 0.032 * 0.032 * xTrue * xTrue);
			if(passNormal)
				dSig->at(index)->Fill(x, weight * bweight);
			if(passCut)
				cut_dSig->at(index)->Fill(x, weight * bweight);
			processedEvents++;
		}

		cout << endl << fitBrch.selEvents << endl;
		return processedEvents;
	}

	void getInputMCGet(TFile *fd, double br, unsigned int index) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		fitStruct cut_fitBrch, cut_totFit;
		TTree *cut_t = (TTree*) fd->Get("cut_fitStruct");
		cut_t->SetBranchAddress("cut_fitStruct", &cut_fitBrch);

		initFitStruct(totFit);
		initFitStruct(cut_totFit);
		sumTreeFitStruct(fitBrch, t, totFit, 1);
		sumTreeFitStruct(cut_fitBrch, cut_t, cut_totFit, 1);

		nmc[index] = totFit.selEvents;
		cut_nmc[index] = cut_totFit.selEvents;

		//Set event nb
		//int nevt = totFit.selEvents;
		//int totalChanEvents = totFit.totEvents;

		//Create histo
		index = d1->size();
		TH1D* xxx1 = (TH1D*) fd->Get("d1");
		TH1D* xxx2 = (TH1D*) fd->Get("d2");
		TH1D* xxx3 = (TH1D*) fd->Get("d3");
		TH1D* xxx4 = (TH1D*) fd->Get("dNew");
		TH1D* xxxA = (TH1D*) fd->Get("dAlpha");
		TH1D* xxxB = (TH1D*) fd->Get("dBeta");
		TH1D* xxxG = (TH1D*) fd->Get("dGamma");

		TH1D* cut_xxx1 = (TH1D*) fd->Get("cut_d1");
		TH1D* cut_xxx2 = (TH1D*) fd->Get("cut_d2");
		TH1D* cut_xxx3 = (TH1D*) fd->Get("cut_d3");
		TH1D* cut_xxx4 = (TH1D*) fd->Get("cut_dNew");
		TH1D* cut_xxxA = (TH1D*) fd->Get("cut_dAlpha");
		TH1D* cut_xxxB = (TH1D*) fd->Get("cut_dBeta");
		TH1D* cut_xxxG = (TH1D*) fd->Get("cut_dGamma");

		tempFD->cd();
		d1->push_back((TH1D*) xxx1->Clone());
		d2->push_back((TH1D*) xxx2->Clone());
		d3->push_back((TH1D*) xxx3->Clone());
		dNew->push_back((TH1D*) xxx4->Clone());
		dAlpha->push_back((TH1D*) xxxA->Clone());
		dBeta->push_back((TH1D*) xxxB->Clone());
		dGamma->push_back((TH1D*) xxxG->Clone());
		d1->at(index)->SetName(TString::Format("d1_%i", index));
		d2->at(index)->SetName(TString::Format("d2_%i", index));
		d3->at(index)->SetName(TString::Format("d3_%i", index));
		dNew->at(index)->SetName(TString::Format("dNew_%i", index));

		dAlpha->at(index)->SetName(TString::Format("dAlpha_%i", index));
		dBeta->at(index)->SetName(TString::Format("dBeta_%i", index));
		dGamma->at(index)->SetName(TString::Format("dGamma_%i", index));

		cut_d1->push_back((TH1D*) cut_xxx1->Clone());
		cut_d2->push_back((TH1D*) cut_xxx2->Clone());
		cut_d3->push_back((TH1D*) cut_xxx3->Clone());
		cut_dNew->push_back((TH1D*) cut_xxx4->Clone());
		cut_dAlpha->push_back((TH1D*) cut_xxxA->Clone());
		cut_dBeta->push_back((TH1D*) cut_xxxB->Clone());
		cut_dGamma->push_back((TH1D*) cut_xxxG->Clone());
		cut_d1->at(index)->SetName(TString::Format("cut_d1_%i", index));
		cut_d2->at(index)->SetName(TString::Format("cut_d2_%i", index));
		cut_d3->at(index)->SetName(TString::Format("cut_d3_%i", index));
		cut_dNew->at(index)->SetName(TString::Format("cut_dNew_%i", index));

		cut_dAlpha->at(index)->SetName(TString::Format("cut_dAlpha_%i", index));
		cut_dBeta->at(index)->SetName(TString::Format("cut_dBeta_%i", index));
		cut_dGamma->at(index)->SetName(TString::Format("cut_dGamma_%i", index));

		cut_dNew->at(index)->SaveAs("toto.png");
		scaleMC(totFit, index, br);
		scaleMC_cut(cut_totFit, index, br);
	}

	int getInputDataGet(TFile *fd, double factor) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		fitStruct cut_fitBrch, cut_totFit;
		TTree *cut_t = (TTree*) fd->Get("cut_fitStruct");
		cut_t->SetBranchAddress("cut_fitStruct", &cut_fitBrch);

		cout << fitBrch.selEvents << endl;
		initFitStruct(totFit);
		initFitStruct(cut_totFit);
		sumTreeFitStruct(fitBrch, t, totFit, 1);
		sumTreeFitStruct(cut_fitBrch, cut_t, cut_totFit, 1);

		cout << fitBrch.selEvents << endl;
		cout << totFit.selEvents << endl;
		//Set event nb
		NSig += factor*totFit.selEvents;
		cut_NSig = factor*cut_totFit.selEvents;

		//Create histo
		int index = dSig->size();
		TH1D* xxx = (TH1D*) fd->Get("sig");
		tempFD->cd();
		dSig->push_back((TH1D*) xxx->Clone());
		dSig->at(index)->SetName(TString::Format("sig_%i", index));

		sig->Add(dSig->at(index), factor);

		xxx = (TH1D*) fd->Get("cut_sig");
		tempFD->cd();
		cut_dSig->push_back((TH1D*) xxx->Clone());
		cut_dSig->at(index)->SetName(TString::Format("cut_sig_%i", index));

		cut_sig->Add(cut_dSig->at(index), factor);

		return NSig;
	}
}

/************************
 * Histograms building and drawing
 ************************/
namespace Display{
	TH1D* buildRatio(double G, double a, TString proc) {
		TH1D* sum = new TH1D("sum", "sum", nbins-1, bins);
		TH1D* r = new TH1D(TString::Format("ratio%s", proc.Data()),
				TString::Format("ratio%s", proc.Data()), nbins-1, bins);

		for (int i = 0; i < inputMCNbr; ++i) {
			sum->Add(d1->at(i), 1.);
			//sum->Add(d2->at(i), G * 2. * a);
			//sum->Add(d3->at(i), G * a * a);
		}
		sum->SetFillColor(8);

		sum->Sumw2();
		r->Sumw2();
		r->Divide(sig, sum, 1, 1, "B");
		delete sum;
		return r;
	}

	TH1D* buildRatioNew(double G, double a, TString proc) {
		TH1D* sum = new TH1D("sum", "sum", nbins-1, bins);
		TH1D* r = new TH1D(TString::Format("ratio%s", proc.Data()),
				TString::Format("ratio%s", proc.Data()), nbins-1, bins);

		for (int i = 0; i < inputMCNbr; ++i) {
			sum->Add(dAlpha->at(i), 1.);
			//sum->Add(dBeta->at(i), G * 2. * a);
			//sum->Add(dGamma->at(i), G * a * a);
		}
		sum->SetFillColor(8);

		sum->Sumw2();
		r->Sumw2();
		r->Divide(sig, sum, 1, 1, "B");
		delete sum;
		return r;
	}

	void drawFitResult(const fitResult& result, TString proc) {
		THStack *stack = new THStack("stack" + proc, "Fit Procedure " + proc);
		TCanvas *c2 = new TCanvas("fit" + proc, "Procedure " + proc);
		TLegend *leg = new TLegend(.78, 0.6, 0.98, 0.94);
		//TPaveText *fitR = new TPaveText(0.78, 0.43, 0.98, 0.59, "NDC BR");
		TPaveText *fitR = new TPaveText(0.47, 0.78, 0.77, 0.94, "NDC BR");
		fitR->AddText("Fit result");
		fitR->AddLine(0., 0.7, 1., 0.7);
		fitR->AddText(Form("G = %f #pm %f", result.norm, result.normErr));
		fitR->AddText(Form("FF = %f #pm %f", result.formFactor, result.formFactorErr));
		fitR->SetTextAlign(12);

		double gWeight = result.norm;
		double a = result.formFactor;

		for (int i = 0; i < inputMCNbr; ++i) {
			TH1D *d3_c = (TH1D*) d3->at(i)->Clone(
					TString::Format("d3_c%s%i", proc.Data(), i));
			leg->AddEntry(d3_c,
					TString::Format("%s FF=x^{2}", mcLegendTitle[i].Data()));
			d3_c->Scale(gWeight * a * a);
			stack->Add(d3_c);
		}
		for (int i = 0; i < inputMCNbr; ++i) {
			TH1D *d2_c = (TH1D*) d2->at(i)->Clone(
					TString::Format("d2_c%s%i", proc.Data(), i));
			leg->AddEntry(d2_c,
					TString::Format("%s FF=x", mcLegendTitle[i].Data()));
			d2_c->Scale(gWeight * 2. * a);
			stack->Add(d2_c);
		}
		for (int i = 0; i < inputMCNbr; ++i) {
			TH1D *d1_c = (TH1D*) d1->at(i)->Clone(
					TString::Format("d1_c%s%i", proc.Data(), i));
			leg->AddEntry(d1_c,
					TString::Format("%s FF=1", mcLegendTitle[i].Data()));
			d1_c->Scale(gWeight);
			stack->Add(d1_c);
		}

		TH1D* ratio = buildRatio(gWeight, a, proc);

		TF1 *f = new TF1("f", "(1+[0]*2.0*x+[0]*[0]*x*x)", 0, 1);
		f->SetLineColor(kRed);
		//f->SetParameter(0, result.norm);
		f->SetParameter(0, result.formFactor);
		c2->Divide(1, 2);
		c2->cd(1);
		sig->DrawClone("S E P");
		sig->GetXaxis()->SetTitle("x");
		stack->Draw("HIST SAME");
		sig->DrawClone("SAMES E P");
		leg->Draw();
		c2->cd(1)->SetLogy(true);
		c2->cd(1)->SetGrid();
		fitR->Draw();
		c2->cd(2);
		//ratio->Fit("f", "R");
		ratio->SetStats(false);
		ratio->Draw("E P");
		f->Draw("LSAME");
		c2->cd(2)->SetLogx(true);
		c2->cd(2)->SetGrid();
	}

	void drawFitNew(const fitResult& result, TString proc) {
		THStack *stack = new THStack("stack" + proc, "Fit Procedure " + proc);
		TCanvas *c2 = new TCanvas("fit" + proc, "Procedure " + proc);
		TLegend *leg = new TLegend(.78, 0.6, 0.98, 0.94);
		//TPaveText *fitR = new TPaveText(0.78, 0.43, 0.98, 0.59, "NDC BR");
		TPaveText *fitR = new TPaveText(0.47, 0.78, 0.77, 0.94, "NDC BR");
		fitR->AddText("Fit result");
		fitR->AddLine(0., 0.7, 1., 0.7);
		fitR->AddText(Form("G = %f #pm %f", result.norm, result.normErr));
		fitR->AddText(Form("FF = %f #pm %f", result.formFactor, result.formFactorErr));
		fitR->SetTextAlign(12);

		double gWeight = result.norm;
		double a = result.formFactor;

		for (int i = 0; i < inputMCNbr; ++i) {
			TH1D *dAlpha_c = (TH1D*) dAlpha->at(i)->Clone(TString::Format("dAlpha_c%s%i", proc.Data(), i));
			leg->AddEntry(dAlpha_c,	TString::Format("%s #alpha", mcLegendTitle[i].Data()));
			dAlpha_c->Scale(gWeight);
			stack->Add(dAlpha_c);
			TH1D *dBeta_c = (TH1D*) dBeta->at(i)->Clone( TString::Format("dBeta_c%s%i", proc.Data(), i));
			leg->AddEntry(dBeta_c, TString::Format("%s #beta", mcLegendTitle[i].Data()));
			dBeta_c->Scale(gWeight * 2. * a);
			stack->Add(dBeta_c);
			TH1D *dGamma_c = (TH1D*) dGamma->at(i)->Clone( TString::Format("dGamma_c%s%i", proc.Data(), i));
			leg->AddEntry(dGamma_c, TString::Format("%s #gamma", mcLegendTitle[i].Data()));
			dGamma_c->Scale(gWeight * a * a);
			stack->Add(dGamma_c);
		}

		TH1D* ratio = buildRatioNew(gWeight, a, proc);
		TF1 *f = new TF1("f", "[0]*(1+[1]*2.0*x+[1]*[1]*x*x)", 0., 1.);
		f->SetParameter(0, result.norm);
		f->SetParameter(1, result.formFactor);
		f->SetLineColor(kRed);

		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(3);
		pad1->SetGrid();
		pad1->Draw();
		pad1->cd();
		sig->DrawClone("S E P");
		stack->Draw("HIST SAME");
		sig->DrawClone("SAMES E P");
		leg->Draw();
		fitR->Draw();
		pad1->SetLogx(true);
		pad1->SetLogy(true);
		c2->cd();
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0);
		pad2->SetBottomMargin(0.2);
		pad2->SetGrid(); // vertical grid
		pad2->Draw();
		pad2->cd();
		//ratio->Fit("f", "R");
		//ratio->SetStats(0);      // No statistics on lower plot
		ratio->Draw("ep");
		f->Draw("LSAME");
		ratio->SetMarkerColor(kRed);
		pad2->SetLogx(true);

		// Y axis mc plot settings
		//stack->GetYaxis()->SetTitleSize(20);
		//stack->GetYaxis()->SetTitleFont(43);
		//stack->GetYaxis()->SetTitleOffset(1.55);

		// Ratio plot (ratio) settings
		ratio->SetTitle(""); // Remove the ratio title

		// Y axis ratio plot settings
		ratio->GetYaxis()->SetTitle("Data/MC");
		ratio->GetYaxis()->SetNdivisions(505);
		ratio->GetYaxis()->SetTitleSize(20);
		ratio->GetYaxis()->SetTitleFont(43);
		ratio->GetYaxis()->SetTitleOffset(1);
		ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		ratio->GetYaxis()->SetLabelSize(15);

		// X axis ratio plot settings
		ratio->GetXaxis()->SetTitleSize(20);
		ratio->GetXaxis()->SetTitleFont(43);
		ratio->GetXaxis()->SetTitleOffset(4.);
		ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		ratio->GetXaxis()->SetLabelSize(15);
	}

	void prepareInputHisto() {
		THStack *std1 = new THStack("std1", "Stack FF1");
		THStack *stdx = new THStack("stdx", "Stack FFx");
		THStack *stdxx = new THStack("stdxx", "Stack FFxx");
		THStack *stdNew = new THStack("stdNew", "Stack FFNew");
		THStack *stdAlpha = new THStack("stdAlpha", "Stack Alpha");
		THStack *stdBeta = new THStack("stdBeta", "Stack Beta");
		THStack *stdGamma = new THStack("stdGamma", "Stack Gamma");
		TLegend *leg1 = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legx = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legxx = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legNew = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legAlpha = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legBeta = new TLegend(.75, 0.6, 0.98, 0.82);
		TLegend *legGamma = new TLegend(.75, 0.6, 0.98, 0.82);

		//Color histos
		for (int i = 0; i < inputMCNbr; ++i) {
			d1->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i * 3]));
			d2->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i * 3 + 1]));
			d3->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i * 3 + 2]));
			dNew->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i * 3]));
			dAlpha->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i * 3]));
			dBeta->at(i)->SetFillColor(
					gStyle->GetColorPalette(mcColors[i * 3 + 1]));
			dGamma->at(i)->SetFillColor(
					gStyle->GetColorPalette(mcColors[i * 3 + 2]));

			std1->Add((TH1D*) d1->at(i)->Clone());
			stdx->Add((TH1D*) d2->at(i)->Clone());
			stdxx->Add((TH1D*) d3->at(i)->Clone());
			stdNew->Add((TH1D*) dNew->at(i)->Clone());
			stdAlpha->Add((TH1D*) dAlpha->at(i)->Clone());
			stdBeta->Add((TH1D*) dBeta->at(i)->Clone());
			stdGamma->Add((TH1D*) dGamma->at(i)->Clone());
			leg1->AddEntry(d1->at(i), mcLegendTitle[i]);
			legx->AddEntry(d2->at(i), mcLegendTitle[i]);
			legxx->AddEntry(d3->at(i), mcLegendTitle[i]);
			legNew->AddEntry(dNew->at(i), mcLegendTitle[i]);
			legAlpha->AddEntry(dAlpha->at(i), mcLegendTitle[i]);
			legBeta->AddEntry(dBeta->at(i), mcLegendTitle[i]);
			legGamma->AddEntry(dGamma->at(i), mcLegendTitle[i]);
		}
		sig->SetLineColor(kRed);

		//Plot input histo
		TCanvas *c1 = new TCanvas("cInput", "Inputs", 1600, 800);
		c1->Divide(2, 2);
		c1->cd(4);
		sig->DrawClone();
		c1->cd(1);
		std1->Draw("HIST");
		leg1->Draw();
		c1->cd(2);
		stdx->Draw("HIST");
		legx->Draw();
		c1->cd(3);
		stdxx->Draw("HIST");
		legxx->Draw();

		new TCanvas("cInputsNew", "Inputs New", 1600, 800);
		sig->DrawClone();
		stdNew->Draw("HIST");
		legNew->Draw();

		TCanvas *c3 = new TCanvas("cABG", "Alpha, Beta, gamma", 1600, 800);
		c3->Divide(2, 2);
		c3->cd(1);
		stdAlpha->Draw("HIST");
		legAlpha->Draw();
		c3->cd(2);
		stdBeta->Draw("HIST");
		legBeta->Draw();
		c3->cd(3);
		stdGamma->Draw("HIST");
		legGamma->Draw();

	}
}

/************************
 * Fit procedure
 ************************/

double fitProcedure(fitResult& result,
		void (*minimFct)(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
				Int_t flag), bool useROOT) {
	//Initialize MINUIT

	TMinuit minuit(2);
	int flag;
	double args[1];
	int fixParam = useROOT ? 1 : 0;
	minuit.SetFCN(minimFct);

	args[0] = 1;
	minuit.mnexcm("SET PRINTOUT", args, 1, flag);
	args[0] = 1;
	minuit.mnexcm("SET ERROR", args, 1, flag);
	args[0] = 1;
	minuit.mnexcm("SET STRATEGY", args, 1, flag);
	minuit.mnparm(0, "G", 1., 100, 0, 0, flag);
	minuit.mnparm(1, "a", 0.05, 0.001, 0, 0, flag);
	minuit.mnparm(2, "fix", fixParam, 0, 0, 0, flag);
	minuit.FixParameter(2);
	if (useROOT)
		minuit.FixParameter(0);

	// Chi2: 1.
	// -logl: 0.5
	minuit.SetErrorDef(0.5);

	args[0] = 100000;
	minuit.mnexcm("MIGRAD", args, 1, flag);

	//Get MINUIT results
	double gWeight, gWeightErr;
	double a, aErr;
	double fctMin;
	double edm, errdef;
	int nvpar, nparx, icstat;

	minuit.GetParameter(0, gWeight, gWeightErr);
	minuit.GetParameter(1, a, aErr);
	minuit.mnstat(fctMin, edm, errdef, nvpar, nparx, icstat);

	cout << fctMin << " " << edm << " " << errdef << " " << nvpar << " "
			<< nparx << " " << icstat << endl;
	result.norm = gWeight;
	//result.norm = getNormalization(a);
	result.normErr = gWeightErr;
	result.formFactor = a;
	result.formFactorErr = aErr;
	return fctMin;
}

/************************
 * Mains
 ************************/

void fit_batch(TString inFile) {
	struct stat buffer;
	srand(time(NULL));
	gStyle->SetOptFit(1);

	tempFileName = ".tempComb";
	tempFileName += rand() % 99999;
	tempFileName += ".root";
	while (stat(tempFileName.Data(), &buffer) == 0) {
		tempFileName = ".tempComb";
		tempFileName += rand() % 99999;
		tempFileName += ".root";
	}
	tempFD = TFile::Open(tempFileName, "RECREATE");

	d1 = new vector<TH1D*>;
	d2 = new vector<TH1D*>;
	d3 = new vector<TH1D*>;
	dNew = new vector<TH1D*>;
	dSig = new vector<TH1D*>;
	dAlpha = new vector<TH1D*>;
	dBeta = new vector<TH1D*>;
	dGamma = new vector<TH1D*>;

	cut_d1 = new vector<TH1D*>;
	cut_d2 = new vector<TH1D*>;
	cut_d3 = new vector<TH1D*>;
	cut_dNew = new vector<TH1D*>;
	cut_dSig = new vector<TH1D*>;
	cut_dAlpha = new vector<TH1D*>;
	cut_dBeta = new vector<TH1D*>;
	cut_dGamma = new vector<TH1D*>;

	//sig = new TH1D("sig", "signal sample", BINS,0,MAX);

	//Get Input
	loadWeights("/afs/cern.ch/user/n/nlurkin/Compact/pi0dalitz_weights.dat");

	if (!readConfig(inFile))
		return;

	if(withEqualBins) loadBins(bins, nbins);

	readFilesFill();

	tempFD->Close();
	remove(tempFileName.Data());

	cout << "--Done--" << endl;
}

void fit_show(TString inFile) {
	//int nbins = 100;
	srand(time(NULL));
	gStyle->SetOptFit(1);

	tempFileName = ".tempComb";
	tempFileName += rand() % 99999;
	tempFileName += ".root";

	tempFD = TFile::Open(tempFileName, "RECREATE");

	d1 = new vector<TH1D*>;
	d2 = new vector<TH1D*>;
	d3 = new vector<TH1D*>;
	dNew = new vector<TH1D*>;
	dSig = new vector<TH1D*>;
	dAlpha = new vector<TH1D*>;
	dBeta = new vector<TH1D*>;
	dGamma = new vector<TH1D*>;

	cut_d1 = new vector<TH1D*>;
	cut_d2 = new vector<TH1D*>;
	cut_d3 = new vector<TH1D*>;
	cut_dNew = new vector<TH1D*>;
	cut_dSig = new vector<TH1D*>;
	cut_dAlpha = new vector<TH1D*>;
	cut_dBeta = new vector<TH1D*>;
	cut_dGamma = new vector<TH1D*>;

	//Get Input
	readConfig(inFile);

	if(!withEqualBins){
		sig = new TH1D("sig", "signal sample", BINS, 0, MAX);
		cut_sig = new TH1D("cut_sig", "signal sample", BINS, 0, MAX);
	}
	else{
		loadBins(bins, nbins);
		sig = new TH1D("sig", "signal sample", nbins-1, bins);
		cut_sig = new TH1D("cut_sig", "signal sample", nbins-1, bins);
	}
	cout << "Initial bins: " << nbins << endl;
	readFilesGet();

	//rebin(125);

	//Scale MC to Data
	double totalMC = 0;
	double totalMCNew = 0;
	double totalGreek = 0;
	for (int i = 0; i < inputMCNbr; ++i) {
		//TODO to check
		totalMC += d1->at(i)->Integral();// + d2->at(i)->Integral() + d3->at(i)->Integral();
		//totalMC += d1->at(i)->Integral() + d2->at(i)->Integral() + d3->at(i)->Integral();
		totalMCNew += dNew->at(i)->Integral();
		totalGreek += dAlpha->at(i)->Integral();
	}

	double factor = ((double) (NSig)) / totalMC;
	double factorNew = ((double) (NSig)) / totalMCNew;
	double factorGreek = ((double) (NSig)) / totalGreek;
	//double factor = nsig/(d1->at(0)->Integral() + d2->at(0)->Integral() + d3->at(0)->Integral());
	for (int i = 0; i < inputMCNbr; ++i) {
		d1->at(i)->Scale(factor);
		d2->at(i)->Scale(factor);
		d3->at(i)->Scale(factor);
		dNew->at(i)->Scale(factorNew);
		dAlpha->at(i)->Scale(factorGreek);
		dBeta->at(i)->Scale(factorGreek);
		dGamma->at(i)->Scale(factorGreek);
	}

	//Scale MC to Data
	totalMC = 0;
	totalMCNew = 0;
	totalGreek = 0;
	for (int i = 0; i < inputMCNbr; ++i) {
		//TODO to check
		totalMC += cut_d1->at(i)->Integral();// + d2->at(i)->Integral() + d3->at(i)->Integral();
		//totalMC += d1->at(i)->Integral() + d2->at(i)->Integral() + d3->at(i)->Integral();
		totalMCNew += cut_dNew->at(i)->Integral();
		totalGreek += cut_dAlpha->at(i)->Integral();
	}

	factor = ((double) (cut_NSig)) / totalMC;
	factorNew = ((double) (cut_NSig)) / totalMCNew;
	factorGreek = ((double) (cut_NSig)) / totalGreek;
	//double factor = nsig/(d1->at(0)->Integral() + d2->at(0)->Integral() + d3->at(0)->Integral());
	for (int i = 0; i < inputMCNbr; ++i) {
		cout << PRINTVAR(factorNew) << endl;
		cut_d1->at(i)->Scale(factor);
		cut_d2->at(i)->Scale(factor);
		cut_d3->at(i)->Scale(factor);
		cut_dNew->at(i)->Scale(factorNew);
		cut_dAlpha->at(i)->Scale(factorGreek);
		cut_dBeta->at(i)->Scale(factorGreek);
		cut_dGamma->at(i)->Scale(factorGreek);
	}

	Display::prepareInputHisto();

	double chi2pv;
	double chi21, chi2ROOT, chi2New, chi2NewROOT;
	double cut_chi2pv;
	double cut_chi21, cut_chi2ROOT, cut_chi2New, cut_chi2NewROOT;

	//Fit
	fitResult result1, cut_result1;
	chi21 = fitProcedure(result1, Fit::minFct, false);
	double chi2Prob1 = TMath::Prob(chi21, 50-2);
//	cut_chi21 = fitProcedure(cut_result1, Fit::minFct, false);
//	double cut_chi2Prob1 = TMath::Prob(cut_chi21, 50-2);
	//double chi2pv = chi2pValue(chi2, nbins-2);

	//chi2Profile(result1, "1");


	//Fit
	fitResult resultROOT, cut_resultROOT;
	chi2ROOT = fitProcedure(resultROOT, Fit::minFct, true);
	double chi2ProbROOT = TMath::Prob(chi2ROOT, 50-1);
//	cut_chi2ROOT = fitProcedure(cut_resultROOT, Fit::minFct, true);
//	double cut_chi2ProbROOT = TMath::Prob(cut_chi2ROOT, 50-1);
	//chi2pv = chi2pValue(chi2, 138);

	//chi2Profile(resultROOT, "ROOT");


	//Fit
	fitResult resultNew, cut_resultNew;
	chi2New = fitProcedure(resultNew, Fit::minFctNew, false);
	double chi2ProbNew = TMath::Prob(chi2New, 50-2);
//	cut_chi2New = fitProcedure(cut_resultNew, Fit::minFctNew_cut, false);
//	double cut_chi2ProbNew = TMath::Prob(cut_chi2New, 50-2);
	//chi2pv = chi2pValue(chi2, 138);

	//chi2Profile(resultNew, "New");

	//Fit
	fitResult resultNewROOT, cut_resultNewROOT;
	chi2NewROOT = fitProcedure(resultNewROOT, Fit::minFctNew, true);
	double chi2ProbNewROOT = TMath::Prob(chi2NewROOT, 50-1);
	//cut_chi2NewROOT = fitProcedure(cut_resultNewROOT, Fit::minFctNew_cut, true);
	//double cut_chi2ProbNewROOT = TMath::Prob(cut_chi2NewROOT, 50-1);
	//chi2pv = chi2pValue(chi2, 138);

	//chi2Profile(resultNew, "New");

	Display::drawFitResult(result1, "1");
	Display::drawFitResult(resultROOT, "ROOT");
	Display::drawFitNew(resultNew, "New");
	Display::drawFitNew(resultNewROOT, "NewROOT");

	cout << "######## Procedure 1 result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << result1.norm << "+-" << result1.normErr << endl;
	cout << "Slope a : " << result1.formFactor << "+-" << result1.formFactorErr << endl;
	cout << "Chi2 : " << chi21 << " prob : " << chi2Prob1 << " p-value : " << chi2pv << endl;

	cout << "######## Procedure ROOT result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << resultROOT.norm << "+-" << resultROOT.normErr << endl;
	cout << "Slope a : " << resultROOT.formFactor << "+-" << resultROOT.formFactorErr << endl;
	cout << "Chi2 : " << chi2ROOT << " prob : " << chi2ProbROOT << " p-value : " << chi2pv << endl;

	cout << "######## Procedure New result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << resultNew.norm << "+-" << resultNew.normErr << endl;
	cout << "Slope a : " << resultNew.formFactor << "+-" << resultNew.formFactorErr << endl;
	cout << "Chi2 : " << chi2New << " prob : " << chi2ProbNew << " p-value : " << chi2pv << endl;

	cout << "######## Procedure New ROOT result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << resultNewROOT.norm << "+-" << resultNewROOT.normErr << endl;
	cout << "Slope a : " << resultNewROOT.formFactor << "+-" << resultNewROOT.formFactorErr << endl;
	cout << "Chi2 : " << chi2NewROOT << " prob : " << chi2ProbNewROOT << " p-value : " << chi2pv << endl;

//	cout << "######## Procedure New ROOT result (cut)#########" << endl << "-------------------------------------" << endl;
//	cout << "Global normalization : " << cut_resultNewROOT.norm << "+-" << cut_resultNewROOT.normErr << endl;
//	cout << "Slope a : " << cut_resultNewROOT.formFactor << "+-" << cut_resultNewROOT.formFactorErr << endl;
//	cout << "Chi2 : " << cut_chi2NewROOT << " prob : " << cut_chi2ProbNewROOT << " p-value : " << cut_chi2pv << endl;

	cout << endl << endl << "RESULTLINE:";
	cout << result1.norm << ";" << result1.normErr << ";" << result1.formFactor << ";" << result1.formFactorErr << ";" << chi21 << ";" << chi2Prob1 << ";";
	cout << resultROOT.norm << ";" << resultROOT.normErr << ";" << resultROOT.formFactor << ";" << resultROOT.formFactorErr << ";" << chi2ROOT << ";" << chi2ProbROOT << ";";
	cout << resultNew.norm << ";" << resultNew.normErr << ";" << resultNew.formFactor << ";" << resultNew.formFactorErr << ";" << chi2New << ";" << chi2ProbNew << ";";
	cout << resultNewROOT.norm << ";" << resultNewROOT.normErr << ";" << resultNewROOT.formFactor << ";" << resultNewROOT.formFactorErr << ";" << chi2NewROOT << ";" << chi2ProbNewROOT << ";";
	cout << NSig << ";" << nmc[0] << ";" << nmc[1] << endl;

	//cout << cut_result1.norm << ";" << cut_result1.normErr << ";" << cut_result1.formFactor << ";" << cut_result1.formFactorErr << ";" << cut_chi21 << ";" << cut_chi2Prob1 << ";";
	//cout << cut_resultROOT.norm << ";" << cut_resultROOT.normErr << ";" << cut_resultROOT.formFactor << ";" << cut_resultROOT.formFactorErr << ";" << cut_chi2ROOT << ";" << cut_chi2ProbROOT << ";";
	//cout << cut_resultNew.norm << ";" << cut_resultNew.normErr << ";" << cut_resultNew.formFactor << ";" << cut_resultNew.formFactorErr << ";" << cut_chi2New << ";" << cut_chi2ProbNew << ";";
	//cout << cut_resultNewROOT.norm << ";" << cut_resultNewROOT.normErr << ";" << cut_resultNewROOT.formFactor << ";" << cut_resultNewROOT.formFactorErr << ";" << cut_chi2NewROOT << ";" << cut_chi2ProbNewROOT << ";";
	//cout << cut_NSig << ";" << nmc[0] << ";" << nmc[1] << endl;


	//tempFD->Close();
	//remove(tempFileName.Data());

	cout << "--Done--" << endl;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "Missing parameter" << endl;
		return -1;
	}

	setStyle();
	gStyle->SetOptStat("11");

	signal(SIGTERM, sighandler);
	signal(SIGINT, sighandler);
	signal(SIGABRT, sighandler);

	TString config(argv[1]);
	if (argc == 2)
		fit_batch(config);
	else {
		theApp = new TApplication("combine", &argc, argv);
		fit_show(config);
		theApp->Run();
	}
}
