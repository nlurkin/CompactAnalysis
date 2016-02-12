#define PRINTVAR(v) #v << "= " << v << " "
#define __CINT_NICO__ 1
#include <signal.h>

#include <TStyle.h>
#include "../common/style.cpp"
#include "pi0DalitzHeader.h"
#include "../common/Fitter/Fitter.h"
#include "../common/Samples/FitMCSample.h"
#include "../common/Samples/FitDataSample.h"

using namespace std;

/*************************
 * Structs
 *************************/

/*************************
 * Globals
 *************************/
#define BINS 10000000
#define MAX 1
#define MAXEVENTS 0
vector<TH1D*> *d1, *d2, *d3, *dSig, *dNew, *dAlpha, *dBeta, *dGamma;
static TH1D *sig;
double Mpi0 = 0.1349766;
static double bins[BINS];
static int nbins;

int nmc[2];

namespace Fit {
void rebin(int binNumber = 0) {
	if (binNumber == 0) {
		double max = sig->GetMaximum();

		cout << max << endl;
		if (max < 100)
			max = 100;
		double s;
		double sum = 0;
		double oldSum = 0;
		int j = 0;
		double var = 0.01;

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
		int ratio = ceil(sig->GetNbinsX() / (double) binNumber);
		for (int i = 0; i <= binNumber; ++i) {
			bins[i] = sig->GetBinLowEdge(i * ratio + 1);
		}
		nbins = binNumber;
	}

	sig = (TH1D*) sig->Rebin(nbins, "sig_reb", bins);
}
}

/************************
 * Fitting
 ************************/

/*
 * TODO (or not) OLD
 *
 */

/*namespace Fit {
double fun(double G, double a, double b1, double b2, double b3) {
	return G * (b1 + 2. * a * b2 + a * a * b3);
}

void minFct(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
		Int_t flag) {
	double chi2 = 0.;
	double b1, b2, b3, s;
	double sigma, sig1 = 0, sig2 = 0, sig3 = 0, totSigma;
	double G;
	bool rootMethod = false;

	TH1D *comp, *sigComp;
	if (par[2] != 0) {
		rootMethod = true;
		comp = new TH1D("comp", "comp", sig->GetNbinsX(), bins);
		sigComp = new TH1D("sigcomp", "sigcomp", sig->GetNbinsX(), bins);
	}
	G = par[0];

	for (int i = 0; i <= sig->GetNbinsX(); ++i) {
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

		if (rootMethod) {
			comp->Fill(sig->GetBinCenter(i), fun(1, par[1], b1, b2, b3));
			sigComp->SetBinContent(i, s);
		} else if (totSigma != 0) {
			chi2 += pow((s - fun(G, par[1], b1, b2, b3)), 2.) / totSigma;
		}
	}

	if (rootMethod) {
		f = sigComp->Chi2Test(comp, "UW CHI2");
		delete comp;
		delete sigComp;
	} else
		f = chi2;
}
}*/

/************************
 * Histograms building and drawing
 ************************/
/*
 * TODO (or not) OLD
 *
 */

/*namespace Display {
TH1D* buildRatio(double G, double a, TString proc) {
	TH1D* sum = new TH1D("sum", "sum", nbins - 1, bins);
	TH1D* r = new TH1D(TString::Format("ratio%s", proc.Data()),
			TString::Format("ratio%s", proc.Data()), nbins - 1, bins);

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
}*/

/************************
 * Mains
 ************************/

void fit_batch() {
	gStyle->SetOptFit(1);

	if (!cfg.testAllOutputFiles())
		return;

	Fitter f;
	RunWeights weights;
	weights.loadWeights(cfg.getWeightFile());
	f.setRunWeights(&weights);
	if (cfg.isWithEqualBins()) {
		loadBins(bins, nbins);
		f.setBinning(nbins, bins);
	}

	f.prepareSamples<FitMCSample, FitDataSample>(cfg);

	f.fillSamples();

	cout << "--Done--" << endl;
}

void fit_show() {
	gStyle->SetOptFit(1);

	Fitter *f = new Fitter();
	if (cfg.isWithEqualBins()) {
		loadBins(bins, nbins);
		f->setBinning(nbins, bins);
	}

	f->prepareSamples<FitMCSample, FitDataSample>(cfg);

	f->getSamples();

	f->mergeSamples<FitMCSample, FitDataSample>();

	//rebin(125);

	f->PrepareHistos(cfg.getMcColors(), cfg.getDataColors());

	//f->fit(true, false);
	f->fit(true, true, cfg.getMaxLoss(), cfg.getStartScan(), cfg.getEndScan());
//	f->fit(false, false);
//	f->fit(false, true);

	return;

	//chi2Profile(resultNew, "New");

//	cout << "######## Procedure 1 result #########" << endl
//			<< "-------------------------------------" << endl;
//	cout << "Global normalization : " << result1.norm << "+-" << result1.normErr
//			<< endl;
//	cout << "Slope a : " << result1.formFactor << "+-" << result1.formFactorErr
//			<< endl;
//	cout << "Chi2 : " << chi21 << " prob : " << chi2Prob1 << " p-value : "
//			<< chi2pv << endl;
//
//	cout << "######## Procedure ROOT result #########" << endl
//			<< "-------------------------------------" << endl;
//	cout << "Global normalization : " << resultROOT.norm << "+-"
//			<< resultROOT.normErr << endl;
//	cout << "Slope a : " << resultROOT.formFactor << "+-"
//			<< resultROOT.formFactorErr << endl;
//	cout << "Chi2 : " << chi2ROOT << " prob : " << chi2ProbROOT << " p-value : "
//			<< chi2pv << endl;
//
//	cout << "######## Procedure New result #########" << endl
//			<< "-------------------------------------" << endl;
//	cout << "Global normalization : " << resultNew.norm << "+-"
//			<< resultNew.normErr << endl;
//	cout << "Slope a : " << resultNew.formFactor << "+-"
//			<< resultNew.formFactorErr << endl;
//	cout << "Chi2 : " << chi2New << " prob : " << chi2ProbNew << " p-value : "
//			<< chi2pv << endl;
//
//	cout << "######## Procedure New ROOT result #########" << endl
//			<< "-------------------------------------" << endl;
//	cout << "Global normalization : " << resultNewROOT.norm << "+-"
//			<< resultNewROOT.normErr << endl;
//	cout << "Slope a : " << resultNewROOT.formFactor << "+-"
//			<< resultNewROOT.formFactorErr << endl;
//	cout << "Chi2 : " << chi2NewROOT << " prob : " << chi2ProbNewROOT
//			<< " p-value : " << chi2pv << endl;

//	cout << "######## Procedure New ROOT result (cut)#########" << endl << "-------------------------------------" << endl;
//	cout << "Global normalization : " << cut_resultNewROOT.norm << "+-" << cut_resultNewROOT.normErr << endl;
//	cout << "Slope a : " << cut_resultNewROOT.formFactor << "+-" << cut_resultNewROOT.formFactorErr << endl;
//	cout << "Chi2 : " << cut_chi2NewROOT << " prob : " << cut_chi2ProbNewROOT << " p-value : " << cut_chi2pv << endl;

//	cout << endl << endl << "RESULTLINE:";
//	cout << result1.norm << ";" << result1.normErr << ";" << result1.formFactor
//			<< ";" << result1.formFactorErr << ";" << chi21 << ";" << chi2Prob1
//			<< ";";
//	cout << resultROOT.norm << ";" << resultROOT.normErr << ";"
//			<< resultROOT.formFactor << ";" << resultROOT.formFactorErr << ";"
//			<< chi2ROOT << ";" << chi2ProbROOT << ";";
//	cout << resultNew.norm << ";" << resultNew.normErr << ";"
//			<< resultNew.formFactor << ";" << resultNew.formFactorErr << ";"
//			<< chi2New << ";" << chi2ProbNew << ";";
//	cout << resultNewROOT.norm << ";" << resultNewROOT.normErr << ";"
//			<< resultNewROOT.formFactor << ";" << resultNewROOT.formFactorErr
//			<< ";" << chi2NewROOT << ";" << chi2ProbNewROOT << ";";
//	cout << NSig << ";" << nmc[0] << ";" << nmc[1] << endl;

	//cout << cut_result1.norm << ";" << cut_result1.normErr << ";" << cut_result1.formFactor << ";" << cut_result1.formFactorErr << ";" << cut_chi21 << ";" << cut_chi2Prob1 << ";";
	//cout << cut_resultROOT.norm << ";" << cut_resultROOT.normErr << ";" << cut_resultROOT.formFactor << ";" << cut_resultROOT.formFactorErr << ";" << cut_chi2ROOT << ";" << cut_chi2ProbROOT << ";";
	//cout << cut_resultNew.norm << ";" << cut_resultNew.normErr << ";" << cut_resultNew.formFactor << ";" << cut_resultNew.formFactorErr << ";" << cut_chi2New << ";" << cut_chi2ProbNew << ";";
	//cout << cut_resultNewROOT.norm << ";" << cut_resultNewROOT.normErr << ";" << cut_resultNewROOT.formFactor << ";" << cut_resultNewROOT.formFactorErr << ";" << cut_chi2NewROOT << ";" << cut_chi2ProbNewROOT << ";";
	//cout << cut_NSig << ";" << nmc[0] << ";" << nmc[1] << endl;

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

	if (!cfg.readFile(argv[1]))
		return -1;
	cfg.print();
	if (argc == 2)
		fit_batch();
	else {
		theApp = new TApplication("combine", &argc, argv);
		fit_show();
		theApp->Run();
	}
}
