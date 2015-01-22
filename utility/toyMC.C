#define __CINT_NICO__ 1

#include <TStyle.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include "pi0DalitzHeader.h"
#include "../userinc/mystructs.h"
#include <math.h>
using namespace std;

#define BINS 10000000
#define MAX 1
static double bins[BINS];
static int nbins;
static double NSig;
vector<TH1D*> *d1, *d2, *d3, *dSig, *dNew, *dAlpha, *dBeta, *dGamma;
static TH1D *sig, *toySig; //, *other;
TRandom3 r;


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


/***************************
 * Input
 ****************************/
namespace Input {
	void getInputMCFill(TFile *fd, TFile *fdout, double br, unsigned int index) {
		getInputMCGet(fd, br, index);
	}

	int getInputDataFill(TFile *fd, TFile* fdout) {
		return getInputDataGet(fd);
	}

	void getInputMCGet(TFile *fd, double br, unsigned int index) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit);

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

	int getInputDataGet(TFile *fd) {
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*) fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit);

		//Set event nb
		NSig = totFit.selEvents;

		//Create histo
		int index = dSig->size();
		TH1D* xxx = (TH1D*) fd->Get("sig");
		tempFD->cd();
		dSig->push_back((TH1D*) xxx->Clone());
		dSig->at(index)->SetName(TString::Format("sig_%i", index));

		sig->Add(dSig->at(index), 1.);

		return NSig;
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

	}
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
void genToy(TH1D *distrib, double a, int sampleSize){
	r.SetSeed(0);

	toySig->Reset();


	double vals[2];
	double x, y;
	double max = distrib->GetMaximum();
	double binWidth;
	int bin;
	cout << max << endl;
	for(int i=0; i<sampleSize;){
		r.RndmArray(2, vals);
		x = vals[0];
		y = vals[1]*max;
		//if(x<0.1) y = vals[1]*max;
		//else y = vals[1]*max/100;
		//if
		bin = x*(nbins+1);
		//bin = getBin(x);
		binWidth = distrib->GetBinLowEdge(bin+1)-distrib->GetBinLowEdge(bin);
		//cout << "x=" << x << " bin=" << bin << " binVal=" << bins[bin] << " y=" << y << " binContent=" << distrib->GetBinContent(bin);
		if(y<distrib->GetBinContent(bin)){
			toySig->Fill(distrib->GetBinLowEdge(bin));
			++i;
			//cout << " ==> accept";
		}
		//cout << endl;
	}
	TCanvas *c = new TCanvas("ctoy" ,"ctoy");
	c->Divide(2,1);
	c->cd(1);
	distrib->Draw();
	c->cd(2);
	toySig->Draw();
}

void startToy(){
	TH1D* sourceDistrib;

	if(!withEqualBins) sourceDistrib = new TH1D("toySource", "Toy MC Source distribution", BINS, 0, MAX);
	else sourceDistrib = new TH1D("toySource", "Toy MC Source distribution", nbins-1, bins);

	for(auto alpha : *dAlpha){
		sourceDistrib->Add(alpha, 1.);
	}

	genToy(sourceDistrib, 1., sig->GetEntries());
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

	if(!withEqualBins) sig = new TH1D("sig", "signal sample", BINS, 0, MAX);
	else{
		loadBins(bins, nbins);
		sig = new TH1D("sig", "signal sample", nbins-1, bins);
	}

	readFilesGet();

	//for(int i=0; i<nbins; ++i) cout << bins[i] <<  " ";
	//cout << endl;
	if(!withEqualBins) toySig = new TH1D("toySig", "Toy MC Source signal", BINS, 0, MAX);
	else toySig = new TH1D("toySig", "Toy MC Source signal", nbins-1, bins);
	startToy();

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

	//Display::prepareInputHisto();
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
