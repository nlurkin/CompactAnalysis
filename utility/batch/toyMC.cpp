#define __CINT_NICO__ 1

#include <TStyle.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include "pi0DalitzHeader.h"
#include "../userinc/mystructs.h"
#include <math.h>
#include <TMinuit.h>
#include <TMath.h>
using namespace std;

#define BINS 10000000
#define MAX 1
static double bins[BINS];
static int nbins;
static double NSig;
vector<TH1D*> *d1, *d2, *d3, *dSig, *dNew, *dAlpha, *dBeta, *dGamma;
static TH1D *sig, *toySig, *modelDistrib; //, *other;
TRandom3 r;


/*************************
 * Structs
 *************************/
typedef struct fitResult_t {
	double norm;
	double normErr;
	double formFactor;
	double formFactorErr;
} fitResult;

/***************************
 * Mandatory from header
 * Opening/Closing files
 ***************************/
void initNewChannel() {
}

void initNewOutput(TFile **fdo, TString fileName) {
	*fdo = TFile::Open(fileName, "RECREATE");
}

void closeMCOutput(TFile *fdo, int index) {
	if (index != -1) {
		fdo->Close();
	}
}

void closeDataOutput(TFile *fdo, int index) {
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

void rebin(int binNumber = 0) {
	if (binNumber == 0) {
		//sig = (TH1D*)sig->Rebin(20, "sig_reb", 0);
		//sig->GetXaxis()->SetRange(sig->GetNbinsX()/5, sig->GetNbinsX());
		double max = sig->GetMaximum();
		//sig->GetXaxis()->SetRange();

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
		int ratio = sig->GetNbinsX() / binNumber;
		for (int i = 0; i <= binNumber; ++i) {
			bins[i] = sig->GetBinLowEdge(i * ratio + 1);
		}
		nbins = binNumber;
	}


	sig = (TH1D*) sig->Rebin(nbins, "sig_reb", bins);
	toySig = (TH1D*) toySig->Rebin(nbins, "toySig_reb", bins);

	//cout << nbins << endl;
	//for(int i=0; i<nbins; i++) cout << bins[i] << " ";
	//cout << endl;
	for (int i = 0; i < inputMCNbr; ++i) {
		d1->at(i) = (TH1D*) d1->at(i)->Rebin(nbins,
				TString(d1->at(i)->GetName()) + "_reb", bins);
		d2->at(i) = (TH1D*) d2->at(i)->Rebin(nbins,
				TString(d2->at(i)->GetName()) + "_reb", bins);
		d3->at(i) = (TH1D*) d3->at(i)->Rebin(nbins,
				TString(d3->at(i)->GetName()) + "_reb", bins);
		dNew->at(i) = (TH1D*) dNew->at(i)->Rebin(nbins,
				TString(dNew->at(i)->GetName()) + "_reb", bins);
		dAlpha->at(i) = (TH1D*) dAlpha->at(i)->Rebin(nbins,
				TString(dAlpha->at(i)->GetName()) + "_reb", bins);
		dBeta->at(i) = (TH1D*) dBeta->at(i)->Rebin(nbins,
				TString(dBeta->at(i)->GetName()) + "_reb", bins);
		dGamma->at(i) = (TH1D*) dGamma->at(i)->Rebin(nbins,
				TString(dGamma->at(i)->GetName()) + "_reb", bins);
	}

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

		TH1D *comp;
		//if(par[2]==0) G = getNormalization(par[1]);
		if (par[2] != 0) {
			rootMethod = true;
			comp = new TH1D("comp", "comp", toySig->GetNbinsX(), bins);
		}
		G = par[0];

		for (int i = 0; i <= toySig->GetNbinsX(); ++i) {
			//if (sig->GetBinLowEdge(i + 1) < 0.005) continue;
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
			s = toySig->GetBinContent(i);
			sigma = toySig->GetBinError(i);

			totSigma = sigma * sigma + sig1 * sig1 + sig2 * sig2 + sig3 * sig3;
			//totSigma = fun(G, par[1], b1,b2,b3);

			if (rootMethod){
				comp->Fill(toySig->GetBinCenter(i), fun(1, par[1], b1, b2, b3));
			}
			else if (totSigma != 0){
				chi2 += pow((s - fun(G, par[1], b1, b2, b3)), 2.) / totSigma;
			}
		}

		if (rootMethod){
			f = toySig->Chi2Test(comp, " CHI2");
			delete comp;
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

		TH1D *comp;
		if (par[2] != 0) {
			rootMethod = true;
			comp = new TH1D("comp", "comp", toySig->GetNbinsX(), bins);
		}
		for (int i = 0; i <= toySig->GetNbinsX(); ++i) {
			//if(sig->GetBinLowEdge(i+1)<0.005) continue;
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
			D_i = toySig->GetBinContent(i);

			m_i = funNew(N, a, a_i, b_i, g_i);

			//cout << D_i << " " << m_i << " " << sigma2 << endl;
			sigma2 = D_i * D_i / M_i + D_i;
			//totSigma = sigma*sigma + sigMC*sigMC;
			if (rootMethod){
				comp->Fill(toySig->GetBinCenter(i), m_i);
			}
			else if (D_i != 0)
				chi2 += pow(D_i - m_i, 2) / sigma2;
		}

		if (rootMethod){
			f = toySig->Chi2Test(comp, " CHI2 P");
			delete comp;
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
		getInputMCGet(fd, br, index);
	}

	int getInputDataFill(TFile *fd, TFile* fdout) {
		return getInputDataGet(fd, 1);
	}

	void getInputMCGet(TFile *fd, double br, unsigned int index) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit, 1);

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

		scaleMC(totFit, index, br);
	}

	int getInputDataGet(TFile *fd, double factor) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit, factor);

		//Set event nb
		NSig = totFit.selEvents;

		//Create histo
		int index = dSig->size();
		TH1D* xxx = (TH1D*) fd->Get("sig");
		tempFD->cd();
		dSig->push_back((TH1D*) xxx->Clone());
		dSig->at(index)->SetName(TString::Format("sig_%i", index));

		cout << dSig->at(index)->GetNbinsX() << " " << sig->GetNbinsX() << endl;
		sig->Add(dSig->at(index), 1.);

		return NSig;
	}

	void readModels(){
		for(auto f : modelFiles){
			cout << f << endl;
			TFile *fd = TFile::Open(f, "READ");
			TH1D* xxxA = (TH1D*) fd->Get("dAlpha");
			modelDistrib->Add(xxxA, 1.);
			fd->Close();
		}
	}
}

/************************
 * Histograms building and drawing
 ************************/
namespace Display{
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
		c3->cd(4);

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

	args[0] = 0;
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
	minuit.SetErrorDef(1);

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
 * ToyMC Generation
 ************************/
int getBin(double value){
	int upper = nbins;
	int lower = 0;
	int currBin = upper/2;

	while(lower+1!=upper){
		if(value>=bins[currBin]) lower = currBin;
		else upper = currBin;
		currBin = (upper+lower)/2;
	}

	return lower;
}

double initPdf(double x){
	if(x<0.2) return 1.;
	else if(x<0.5) return 0.2;
	else return 0.01;
}
void genToy(TH1D *distrib, double a, int sampleSize){
	//r.SetSeed(0);
	gRandom->SetSeed(0);

	toySig->Reset();
	toySig->Sumw2();

	double x;
	for(int i=0; i<sampleSize;){
		//toySig->FillRandom(distrib, 10000);
		//break;
		//r.RndmArray(2, vals);
		//x = vals[0];
		//y = vals[1]*initPdf(vals[0]);
		//if(x<0.1) y = vals[1]*max;
		//else y = vals[1]*max/100;
		//if
		//bin = x*(nbins+1);
		//bin = getBin(x);
		//binWidth = distrib->GetBinLowEdge(bin+1)-distrib->GetBinLowEdge(bin);
		//cout << "x=" << x << " bin=" << bin << " binVal=" << bins[bin] << " y=" << y << " binContent=" << distrib->GetBinContent(bin);
		//if(vals[1]<distrib->GetBinContent(bin)/initPdf(x)){
			//toySig->Fill(x);
			//++i;
			//cout << " ==> accept " << i;
		//}
		//cout << endl;
		x = distrib->GetRandom();
		toySig->Fill(x, 1.+2*a*x+a*a*x*x);
		++i;
	}
	/*TCanvas *c = new TCanvas("ctoy" ,"ctoy");
	c->Divide(2,1);
	c->cd(1);
	distrib->Draw();
	c->cd(2);
	toySig->Draw();*/
}

void startToy(double a, double deltaA, double testNumber){
	modelDistrib->Scale(1/modelDistrib->Integral());
	double chi2pv;
	double chi21, chi2ROOT, chi2New, chi2NewROOT;

	double limit = 2;
	TH1D* fitOld = new TH1D("fitOld", "fitOld", 200, -limit, limit);
	TH1D* fitOldRoot = new TH1D("fitOldRoot", "fitOldRoot", 200, -limit, limit);
	TH1D* fitNew = new TH1D("fitNew", "fitNew", 200, -limit, limit);
	TH1D* fitNewRoot = new TH1D("fitNewRoot", "fitNewRoot", 200, -limit, limit);

	TH1D* chiOld = new TH1D("chiOld", "chiOld", 200, 20, 200);
	TH1D* chiOldRoot = new TH1D("chiOldRoot", "chiOldRoot", 200, 0, 200);
	TH1D* chiNew = new TH1D("chiNew", "chiNew", 200, 20, 200);
	TH1D* chiNewRoot = new TH1D("chiNewRoot", "chiNewRoot", 200, 0, 200);

	for(int i=0; i<testNumber; ++i){
		genToy(modelDistrib, a, sig->GetEntries());

		//Fit
		fitResult result1;
		chi21 = fitProcedure(result1, Fit::minFct, false);
		double chi2Prob1 = TMath::Prob(chi21, 50-2);

		//Fit
		fitResult resultROOT;
		chi2ROOT = fitProcedure(resultROOT, Fit::minFct, true);
		double chi2ProbROOT = TMath::Prob(chi2ROOT, 50-1);

		//Fit
		fitResult resultNew;
		chi2New = fitProcedure(resultNew, Fit::minFctNew, false);
		double chi2ProbNew = TMath::Prob(chi2New, 50-2);

		//Fit
		fitResult resultNewROOT;
		chi2NewROOT = fitProcedure(resultNewROOT, Fit::minFctNew, true);
		double chi2ProbNewROOT = TMath::Prob(chi2NewROOT, 50-1);

		fitOld->Fill((result1.formFactor-a)/deltaA);
		fitOldRoot->Fill((resultROOT.formFactor-a)/deltaA);
		fitNew->Fill((resultNew.formFactor-a)/deltaA);
		fitNewRoot->Fill((resultNewROOT.formFactor-a)/deltaA);

		chiOld->Fill(chi21);
		chiOldRoot->Fill(chi2ROOT);
		chiNew->Fill(chi2New);
		chiNewRoot->Fill(chi2NewROOT);
		//cout << "######## Procedure 1 result #########" << endl << "-------------------------------------" << endl;
		//cout << "Global normalization : " << result1.norm << "+-" << result1.normErr << endl;
		//cout << "Slope a : " << result1.formFactor << "+-" << result1.formFactorErr << endl;
		//cout << "Chi2 : " << chi21 << " prob : " << chi2Prob1 << " p-value : " << chi2pv << endl;

		//cout << "######## Procedure ROOT result #########" << endl << "-------------------------------------" << endl;
		//cout << "Global normalization : " << resultROOT.norm << "+-" << resultROOT.normErr << endl;
		//cout << "Slope a : " << resultROOT.formFactor << "+-" << resultROOT.formFactorErr << endl;
		//cout << "Chi2 : " << chi2ROOT << " prob : " << chi2ProbROOT << " p-value : " << chi2pv << endl;

		//cout << "######## Procedure New result #########" << endl << "-------------------------------------" << endl;
		//cout << "Global normalization : " << resultNew.norm << "+-" << resultNew.normErr << endl;
		//cout << "Slope a : " << resultNew.formFactor << "+-" << resultNew.formFactorErr << endl;
		//cout << "Chi2 : " << chi2New << " prob : " << chi2ProbNew << " p-value : " << chi2pv << endl;

		//cout << "######## Procedure New ROOT result #########" << endl << "-------------------------------------" << endl;
		//cout << "Global normalization : " << resultNewROOT.norm << "+-" << resultNewROOT.normErr << endl;
		//cout << "Slope a : " << resultNewROOT.formFactor << "+-" << resultNewROOT.formFactorErr << endl;
		//cout << "Chi2 : " << chi2NewROOT << " prob : " << chi2ProbNewROOT << " p-value : " << chi2pv << endl;
	}

	TCanvas *c = new TCanvas("cToyMCFit", "ToyMCFit");
	c->Divide(4,2);
	c->cd(1);
	fitOld->Draw();
	c->cd(2);
	fitOldRoot->Draw();
	c->cd(3);
	fitNew->Draw();
	c->cd(4);
	fitNewRoot->Draw();
	c->cd(5);
	chiOld->Draw();
	c->cd(6);
	chiOldRoot->Draw();
	c->cd(7);
	chiNew->Draw();
	c->cd(8);
	chiNewRoot->Draw();

}

/************************
 * Mains
 ************************/

void toyMC(TString inFile) {
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

	//Get Input
	readConfig(inFile);

	modelDistrib = new TH1D("toySource", "Toy MC Source distribution", BINS, 0, MAX);

	if(withEqualBins) loadBins(bins, nbins);

	if(!withEqualBins) sig = new TH1D("sig", "signal sample", BINS, 0, MAX);
	else sig = new TH1D("sig", "signal sample", nbins-1, bins);
	if(!withEqualBins) toySig = new TH1D("toySig", "Toy MC Source signal", BINS, 0, MAX);
	else toySig = new TH1D("toySig", "Toy MC Source signal", nbins-1, bins);

	readFilesGet();
	Input::readModels();

	cout << "Input MC Number " << inputMCNbr << endl;
	rebin(50);

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

	Display::prepareInputHisto();

	TCanvas *c = new TCanvas("toyGen", "toyGen");
	c->Divide(2,1);
	c->cd(1);
	modelDistrib->Draw();
	startToy(0.037, 0.0011, 1000);
	c->cd(2);
	toySig->Draw();
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "Missing parameter" << endl;
		return -1;
	}

	setStyle();

	signal(SIGTERM, sighandler);
	signal(SIGINT, sighandler);
	signal(SIGABRT, sighandler);

	TString config(argv[1]);
	if (argc == 2)
		toyMC(config);
	else {
		theApp = new TApplication("toyMC", &argc, argv);
		toyMC(config);
		theApp->Run();
	}
}
