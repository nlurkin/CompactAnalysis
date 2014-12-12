#include <TSystem.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TH1D.h>
#include <iostream>
#include <math.h>
#include <TGraph.h>
//#include "utility/mystructs.h"
#include <TTree.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TF1.h>
#include <TStyle.h>
#include <vector>
#include <TPaveText.h>
using namespace std;

#define BINS 10000
#define MAX 0.14
#define MAXEVENTS 0;
<<<<<<< HEAD
static int inputMCNbr = 0
=======
static int inputMCNbr = 0;
>>>>>>> ptp-push
static int inputDataNbr = 0;
static vector<TH1D*> *d1, *d2, *d3, *dSig;
static vector<int> mcColors, dataColors;
vector<TString> mcLegendTitle, dataLegendTitle;
static TH1D *sig, *other;
double Mpi0 = 0.1349766;
static double bins[BINS];
static int nbins;
static double NSig;
//static int n1, nx, nxx, nsig;
TFile *tempFD;

/************************
 * Structs
 ************************/
typedef struct fitStruct_t {
	int totEvents;
	int selEvents;
	int n1;
	int nx;
	int nxx;
} fitStruct;

typedef struct fitResult_t{
	double norm;
	double normErr;
	double formFactor;
	double formFactorErr;
} fitResult;

/************************
 * Fitting
 ************************/

double fun(double G, double a, double b1, double b2, double b3){
	return G*(b1 + 2.*a*b2 + a*a*b3);
}

double fun2(double G, double a, double b1){
	return G*(1+2.*a*b1 + a*a*b1*b1);
}

void minFct(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag){
	double chi2 = 0.;
	double b1, b2, b3, s;
	double sigma, sig1, sig2, sig3;
	double G;

	if(par[2]==0) G = getNormalization(par[1]);
	else G = par[0];

	//printIntegrals(G, par[1]);
	for(int i=0; i<=sig->GetNbinsX(); ++i){
		if(sig->GetBinLowEdge(i+1)<0.01) continue;
		b1 = 0;
		b2 = 0;
		b3 = 0;
		for(int j=0; j<inputMCNbr; ++j){
			b1 += d1->at(j)->GetBinContent(i);
			b2 += d2->at(j)->GetBinContent(i);
			b3 += d3->at(j)->GetBinContent(i);
			sig1 = d1->at(j)->GetBinError(i);
			sig2 = d2->at(j)->GetBinError(i);
			sig3 = d3->at(j)->GetBinError(i);
		}
		s = sig->GetBinContent(i);
		sigma = sig->GetBinError(i);

		if(sigma!=0) chi2 += pow((s-fun(G, par[1], b1,b2,b3)),2.)/(sigma*sigma/* + sig1*sig1 + sig2*sig2 + sig3*sig3*/);
	}

	f = chi2;
}

void minFct2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag){
	double chi2 = 0.;
	double b1;
	double sigma;
	double x;
	for(int i=0; i<=sig->GetNbinsX(); ++i){
		if(sig->GetBinLowEdge(i+1)<0.01) continue;
		if(sig->GetBinLowEdge(i+1)>0.12) continue;
		b1 = other->GetBinContent(i);
		sigma = other->GetBinError(i);
		x = other->GetBinCenter(i);

		if(sigma!=0) chi2 += pow((b1-fun2(par[0], par[1], x)),2.)/(sigma*sigma);
	}

	f = chi2;
}
/************************
 * Getting input
 ************************/
void initNewInput(TFile **fdo, fitStruct &fitBrch, TTree **fitTree, TString fileName){
	//Output
	*fdo = TFile::Open(fileName, "RECREATE");
	*fitTree = new TTree("fitStruct", "fitStruct tree");
	(*fitTree)->Branch("fitStruct", &fitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	fitBrch.n1 = 0;
	fitBrch.nx = 0;
	fitBrch.nxx = 0;
	fitBrch.selEvents = 0;
	fitBrch.totEvents = 0;
}

void closeMCFileProcess(int prevIndex, TFile *fdo, fitStruct &fitBrch, TTree *fitTree){
	if(prevIndex!=-1){
		fdo->cd();
		fitTree->Fill();
		fitTree->Write();
		d1->at(prevIndex)->Write();
		d2->at(prevIndex)->Write();
		d3->at(prevIndex)->Write();

		fdo->Close();

		d1->at(prevIndex)->SetName(TString("d1_") + (Long_t)prevIndex);
		d2->at(prevIndex)->SetName(TString("d2_") + (Long_t)prevIndex);
		d3->at(prevIndex)->SetName(TString("d3_") + (Long_t)prevIndex);
	}
}

void closeDataFileProcess(TFile *fdo, TTree *fitTree){
	fdo->cd();
	fitTree->Fill();
	fitTree->Write();
	dSig->at(0)->Write();
	fdo->Close();
}

int readConfig(TString confFile, bool test, bool fillMC, bool fillData){
	vector<double> brs;
	vector<TString> mcFileNames;
	vector<TString> mcOutputFiles;
	vector<TString> dataFileNames;
	vector<TString> dataOutputFiles;
	vector<int> mcIndexes;
	int nsig = 0;

	//read list file
	ifstream listFile(confFile.Data());

	string line;
	while(getline(listFile,line)){
		TString tline(line);
		if(tline.BeginsWith("#")) continue;
		TString key(tline(0,tline.First('=')));
		TString values(tline(tline.First('=')+1, 1000));
		TObjArray* tok = values.Tokenize(" ");
		cout << tline << endl;
		for(int i=0; i<tok->GetEntries(); ++i){
			TString entry(((TObjString*)tok->At(i))->GetString());
			if(key.CompareTo("mcfiles")==0) mcFileNames.push_back(entry);
			else if(key.CompareTo("mcout")==0) mcOutputFiles.push_back(entry);
			else if(key.CompareTo("brs")==0) brs.push_back(entry.Atof());
			else if(key.CompareTo("mccolors")==0) mcColors.push_back(entry.Atoi());
			else if(key.CompareTo("mclegends")==0) mcLegendTitle.push_back(entry);
			else if(key.CompareTo("datafiles")==0) dataFileNames.push_back(entry);
			else if(key.CompareTo("dataout")==0) dataOutputFiles.push_back(entry);
			else if(key.CompareTo("datacolors")==0) dataColors.push_back(entry.Atoi());
			else if(key.CompareTo("datalegends")==0) dataLegendTitle.push_back(entry);
			else if(key.CompareTo("mcIndex")==0) mcIndexes.push_back(entry.Atoi());
			delete tok->At(i);
		}
	}

	tempFD = TFile::Open(".tempFit.root", "RECREATE");

	vector<TFile*> fd;
	TFile* ffd;
	int prevIndex = -1;

	TTree* fitTree;
	TFile *fdo;
	fitStruct fitBrch;

	int fileNumber;

	if(fillMC) fileNumber = mcFileNames.size();
	else fileNumber = mcOutputFiles.size();

	//Getting MC
	//if(test) sig = new TH1D("sig", "signal sample", BINS,0,MAX);
	for(int i=0; i<fileNumber; ++i){
		if(fillMC){
			if(prevIndex != mcIndexes[i] && prevIndex!=-1){
				printComponentSummary(prevIndex, fitBrch, fd[prevIndex]->GetName(),brs[prevIndex]);
				closeMCFileProcess(prevIndex, fdo, fitBrch, fitTree);
				scaleMC(fitBrch, prevIndex, brs[prevIndex]);
			}
		}

		if(fillMC){
			cout << mcFileNames[i] << endl;
			ffd = TFile::Open(mcFileNames[i]);
		}
		else{
			cout << mcOutputFiles[i] << endl;
			ffd = TFile::Open(mcOutputFiles[i]);
		}
		fd.push_back(ffd);
		if(fillMC){
			if(prevIndex != mcIndexes[i]){
				prevIndex = mcIndexes[i];
				++inputMCNbr;
				initNewInput(&fdo, fitBrch, &fitTree, mcOutputFiles[prevIndex]);
			}

			getInputMCFill(ffd, fdo, fitBrch, brs[prevIndex], test, mcOutputFiles[prevIndex], prevIndex);
			ffd->Close();
		}
		else{
			getInputMCPlot(ffd, brs[i]);
			++inputMCNbr;
		}
	}

	if(fillMC){
		printComponentSummary(prevIndex, fitBrch, fd[prevIndex]->GetName(),brs[prevIndex]);
		closeMCFileProcess(prevIndex, fdo, fitBrch, fitTree);
		scaleMC(fitBrch, prevIndex, brs[prevIndex]);
	}

	if(test){
		sig->Scale(10000/sig->Integral());
		nsig = sig->Integral();
	}

	//Getting data
	if(!test){
		if(fillData){
			initNewInput(&fdo, fitBrch, &fitTree, dataOutputFiles[0]);
		}

		if(fillData) fileNumber = dataFileNames.size();
		else fileNumber = dataOutputFiles.size();

		for(int i=0; i<fileNumber; ++i){
			if(fillData){
				cout << dataFileNames[i] << endl;
				ffd = TFile::Open(dataFileNames[i]);
			}
			else{
				cout << dataOutputFiles[i] << endl;
				ffd = TFile::Open(dataOutputFiles[i]);
			}
			fd.push_back(ffd);
			++inputDataNbr;
			if(fillData){
				nsig += getInputDataFill(ffd, fdo, fitBrch);
				ffd->Close();
			}
			else nsig += getInputDataPlot(ffd);
		}

		if(fillData){
			closeDataFileProcess(fdo, fitTree);
		}
	}

	NSig = nsig;
	return nsig;
}

void getInputMCFill(TFile *fd, TFile *fdo, fitStruct &fitBrch, double br, bool test, TString outputFile, int index){

	//Get the TTree
	//Input
	pi0dEvent *eventBrch = new pi0dEvent();
	TTree *t = (TTree*)fd->Get("event");
	t->SetBranchAddress("pi0dEvent", &eventBrch);

	//Output
	//TFile *fdo = TFile::Open(outputFile, "RECREATE");
	//fitStruct fitBrch;
	//TTree* fitTree = new TTree("fitStruct", "fitStruct tree");
	//fitTree->Branch("fitStruct", &fitBrch, "totEvents/I:selEvents:n1:nx:nxx");

	tempFD->cd();
	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = ((TH1D*)fd->Get("Cuts"))->GetEntries());
	int npart;
	int divider;
	if(MAXEVENTS>0 && nevt>MAXEVENTS){
		double ratio = (double)MAXEVENTS/(double)nevt;
		nevt=MAXEVENTS;
		totalChanEvents = totalChanEvents*ratio;
	}
	if(test){
		npart = nevt/6;
		divider = 6;
	}
	else{
		npart = nevt/5;
		divider = 5;
	}
	int n1 = npart*3;
	int nx = npart;
	int nxx = npart;
	int ntest = npart;

	fitBrch.totEvents += totalChanEvents;
	fitBrch.selEvents += nevt;
	fitBrch.n1 += n1;
	fitBrch.nx += nx;
	fitBrch.nxx += nxx;

	//Create histo
	//int index = d1->size();
	if(index == d1->size()){
		TH1D* xxx1 = new TH1D("d1", "sample 1", BINS,0,MAX);
		xxx1->Sumw2();
		d1->push_back(xxx1);
		TH1D* xxx2 = new TH1D("d2", "sample x", BINS,0,MAX);
		xxx2->Sumw2();
		d2->push_back(xxx2);
		TH1D* xxx3 = new TH1D("d3", "sample x^{2}", BINS,0,MAX);
		xxx3->Sumw2();
		d3->push_back(xxx3);
	}

	//Read events and fill histo
	int i=0;
	double x;
	double bweight=1.;
	double weight;
	int mod;

	/**
	 * New
	 */
	cout << "Filling" << endl;
	for(; i<nevt; ++i){
		t->GetEntry(i);
		x = pow(eventBrch->mee/Mpi0,2.);
		bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
		mod = eventBrch->timestamp % divider;
		//cout << "Timestamp: " << eventBrch->timestamp << " % " << divider << " " << mod << endl;
		if(mod==0 || mod==1 || mod==2){
			weight = 1.;
			d1->at(index)->Fill(eventBrch->mee, bweight*weight);
		}
		else if(mod==3){
			weight = x;
			d2->at(index)->Fill(eventBrch->mee, bweight*weight);
		}
		else if(mod==4){
			weight = x*x;
			d3->at(index)->Fill(eventBrch->mee, bweight*weight);
		}
		else if (mod==5 && test) {
			double aTest = 0.5;
			//weight = 1.+2.*aTest*x+aTest*aTest*x*x;
			weight = 1.;
			//bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
			bweight = 1.;
			sig->Fill(eventBrch->mee, bweight*weight*br);
		}
	}
	/**
	 * Old
	 */
	/*cout << "Filling d1" << endl;
	for(; i<n1; i++){
		t->GetEntry(i);
		x = pow(eventBrch->mee/Mpi0,2.);
		weight = 1.;
		bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
		d1->at(index)->Fill(eventBrch->mee, bweight*weight);
	}
	cout << "Filling d2" << endl;
	for(; i<n1+nx; i++){
		t->GetEntry(i);
		x = pow(eventBrch->mee/Mpi0,2.);
		weight = x;
		bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
		d2->at(index)->Fill(eventBrch->mee, bweight*weight);
	}
	cout << "Filling d3" << endl;
	for(; i<n1+nx+nxx; i++){
		t->GetEntry(i);
		x = pow(eventBrch->mee/Mpi0,2.);
		weight = x*x;
		bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
		d3->at(index)->Fill(eventBrch->mee, bweight*weight);
	}
	if(test){
		cout << "Filling test sample" << endl;
		for(; i<n1+nx+nxx+ntest; i++){
			t->GetEntry(i);
			x = pow(eventBrch->mee/Mpi0,2.);
			double aTest = 0.5;
			//weight = 1.+2.*aTest*x+aTest*aTest*x*x;
			weight = 1.;
			//bweight = 1./(1.+2.*0.032*x+0.032*0.032*x*x);
			bweight = 1.;
			sig->Fill(eventBrch->mee, bweight*weight*br);
		}

	}*/

	//fitTree->Fill();
	//fitTree->Write();
	//xxx1->Write();
	//xxx2->Write();
	//xxx3->Write();
}

void getInputMCPlot(TFile *fd, double br){

	fitStruct fitBrch;
	TTree *t = (TTree*)fd->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fitBrch);

	t->GetEntry(0);

	//Set event nb
	int nevt = fitBrch.selEvents;
	int n1 = fitBrch.n1;
	int nx = fitBrch.nx;
	int nxx = fitBrch.nxx;
	int totalChanEvents = fitBrch.totEvents;

	//Create histo
	int index = d1->size();
	TH1D* xxx1 = (TH1D*)fd->Get("d1");
	xxx1->SetName("d1" + index);
	d1->push_back(xxx1);
	TH1D* xxx2 = (TH1D*)fd->Get("d2");
	xxx2->SetName("d2" + index);
	d2->push_back(xxx2);
	TH1D* xxx3 = (TH1D*)fd->Get("d3");
	xxx3->SetName("d3" + index);
	d3->push_back(xxx3);

	printComponentSummary(index, fitBrch, fd->GetName(), br);

	scaleMC(fitBrch, index, br);
}

int getInputDataFill(TFile *fd, TFile *fdo, fitStruct &fitBrch){
	//Input
	pi0dEvent *eventBrch = new pi0dEvent();
	TTree *t = (TTree*)fd->Get("event");
	t->SetBranchAddress("pi0dEvent", &eventBrch);

	//Output
	//TFile *fdo = TFile::Open(outFile, "RECREATE");
	//TTree* fitTree = new TTree("fitStruct", "fitStruct tree");
	//fitTree->Branch("fitStruct", &fitBrch, "totEvents/I:selEvents:n1:nx:nxx");

	tempFD->cd();

	// Set Number of events
	int nevt = t->GetEntries();
	if(MAXEVENTS>0 && nevt>MAXEVENTS) nevt=MAXEVENTS;
	int nsig = nevt;

	//Create histo
	//int index = dSig->size();
	int index = 0;
	if(dSig->size()==0){
		TH1D* xxx = new TH1D("sig", "signal sample", BINS,0,MAX);
		dSig->push_back(xxx);
	}

	//Read event and fill histo
	int i=0;
	double x;
	double bweight=1.;
	double weight;

	double a = 0.032;
	cout << "Filling data" << endl;
	for(i=0; i<nsig; i++){
		t->GetEntry(i);
		x = pow(eventBrch->mee/Mpi0,2.);
		weight = 1.;//+2.*a*x+a*a*x*x;
		sig->Fill(eventBrch->mee, weight);
		dSig->at(index)->Fill(eventBrch->mee, weight);
	}
	fitBrch.selEvents += nsig;
	return nsig;
}

int getInputDataPlot(TFile *fd){
	fitStruct fitBrch;
	TTree *t = (TTree*)fd->Get("fitStruct");
	t->SetBranchAddress("fitStruct", &fitBrch);

	t->GetEntry(0);

	//Set event nb
	int nevt = fitBrch.selEvents;
	int nsig = nevt;

	//Create histo
	int index = dSig->size();
	TH1D* xxx = (TH1D*)fd->Get("sig");
	xxx->SetName("sig" + index);
	dSig->push_back(xxx);
	sig->Add(xxx, 1.);

	return nsig;
}
/************************
 * Utility
 ************************/
void scaleMC(fitStruct N, int index, double br){
	cout << "Rescaling" << endl;

	// Rescale histo
	double selRatio1 = (double)N.selEvents/(double)N.n1;
	double selRatiox = (double)N.selEvents/(double)N.nx;
	double selRatioxx = (double)N.selEvents/(double)N.nxx;
	d1->at(index)->Scale(br/((double)N.totEvents/selRatio1));
	d2->at(index)->Scale(br/((double)N.totEvents/selRatiox));
	d3->at(index)->Scale(br/((double)N.totEvents/selRatioxx));
}

void rebin(int binNumber=0){
	if(binNumber==0){
		double max = sig->GetMaximum()*4.;
		double s;
		double sum = 0;
		double oldSum = 0;
		int j = 0;
		double var = 0.05;

		for(int i=0; i<=sig->GetNbinsX(); ++i){
			s = sig->GetBinContent(i);
			oldSum = sum;
			sum += s;
			if((fabs(sum-max)<var*max)){
				bins[j] = sig->GetBinLowEdge(i+1);
				++j;
				sum = 0;
				oldSum=0;
			}
			else if(sum>(max*(1.+var))){
				if(fabs(oldSum-max)<fabs(sum-max)){
					bins[j] = sig->GetBinLowEdge(i);
					++j;
					sum = s;
					oldSum=0;
				}
				else{
					bins[j] = sig->GetBinLowEdge(i+1);
					++j;
					sum = 0;
					oldSum=0;
				}
			}
		}
		bins[j] = sig->GetBinLowEdge(sig->GetNbinsX()+1);
		nbins = j;
	}
	else{
		int ratio = sig->GetNbinsX()/binNumber;
		for(int i=0; i<=binNumber; ++i){
			bins[i] = sig->GetBinLowEdge(i*ratio+1);
		}
		nbins = binNumber;
	}

	sig = (TH1D*)sig->Rebin(nbins, "sig_reb", bins);
	for(int i=0; i<inputMCNbr; ++i){
		d1->at(i) = (TH1D*)d1->at(i)->Rebin(nbins, TString(d1->at(i)->GetName()) + "_reb", bins);
		d2->at(i) = (TH1D*)d2->at(i)->Rebin(nbins, TString(d2->at(i)->GetName()) + "_reb", bins);
		d3->at(i) = (TH1D*)d3->at(i)->Rebin(nbins, TString(d3->at(i)->GetName()) + "_reb", bins);
	}

}

double getNormalization(double a){
	double G=0;
	for(int j=0; j<inputMCNbr; ++j){
		G += d1->at(j)->Integral()*1.0 + d2->at(j)->Integral()*a*2. + d3->at(j)->Integral()*a*a;
	}
	G = ((double)NSig)/G;
	return G;
}

void printIntegrals(double G, double a){

	double sum = 0;
	for(int j=0; j<inputMCNbr; ++j){
		sum += d1->at(j)->Integral();
		sum += 2.*a*d2->at(j)->Integral();
		sum += a*a*d3->at(j)->Integral();
	}
	cout << "Data :\t" << sig->Integral() << endl;
	cout << "\t\tBrut:" << endl;
	cout << "1 :\t" << d1->at(0)->Integral() << endl;
	cout << "x :\t" << d2->at(0)->Integral() << endl;
	cout << "xx :\t" << d3->at(0)->Integral() << endl;
	cout << "\t\tCorrection (" << a << ") :\t" << G << endl;
	cout << "1 :\t" << G*d1->at(0)->Integral() << endl;
	cout << "x :\t" << G*2.*a*d2->at(0)->Integral() << endl;
	cout << "xx :\t" << G*a*a*d3->at(0)->Integral() << endl;
	cout << "MC :\t" << G*sum << endl;
}

void printComponentSummary(int index, fitStruct fit, TString fileName, double br){
	cout << "\tChannel: " << mcLegendTitle.at(index) << endl;
	cout << "\tFile: " << fileName << endl;
	cout << "\tTotal events: " << fit.totEvents << endl;
	cout << "\tSelected events: " << fit.selEvents << endl;
	cout << "\tSample 1 events: " << fit.n1 << endl;
	cout << "\tSample x events: " << fit.nx << endl;
	cout << "\tSample x^2 events: " << fit.nxx << endl;
	cout << "\tBr: " << br << endl;
	cout << "\tScaling factor: " << (double)(br*(double)fit.selEvents)/(double)((double)fit.totEvents*(double)fit.n1) << endl;
}

/************************
 * Histograms building and drawing
 ************************/

TH1D* buildRatio(double G, double a, TString proc){
	TH1D* sum = new TH1D("sum", "sum", nbins, bins);
	TH1D* r = new TH1D("ratio" + proc, "ratio"+proc, nbins, bins);

	for(int i=0; i<inputMCNbr; ++i){
		sum->Add(d1->at(i), G);
		sum->Add(d2->at(i), G*2.*a);
		sum->Add(d3->at(i), G*a*a);
	}
	sum->SetFillColor(8);

	sum->Sumw2();
	r->Sumw2();
	r->Divide(sig, sum, 1, 1, "B");
	delete sum;
	return r;
}

void drawFitResult(const fitResult& result, TString proc){
	THStack *stack= new THStack("stack"+proc, "Fit Procedure " + proc);
	TCanvas *c2 = new TCanvas("fit" + proc, "Procedure " + proc);
	TLegend *leg = new TLegend(.78,0.4,0.98,0.74);
	TPaveText *fitR = new TPaveText(0.78, 0.23, 0.98, 0.39, "NDC BR");
	fitR->AddText("Fit result");
	fitR->AddLine(0., 0.7, 1., 0.7);
	fitR->AddText(Form("G = %f #pm %f", result.norm, result.normErr));
	fitR->AddText(Form("FF = %f #pm %f", result.formFactor, result.formFactorErr));
	fitR->SetTextAlign(12);

	double gWeight = result.norm;
	double a = result.formFactor;

	for(int i=0; i<inputMCNbr; ++i){
		TH1D *d3_c = (TH1D*)d3->at(i)->Clone(TString("d3_c")+proc + (Long_t)i);
		leg->AddEntry(d3_c, mcLegendTitle[i] + " FF=x^{2}");
		d3_c->Scale(gWeight*a*a);
		stack->Add(d3_c);
	}
	for(int i=0; i<inputMCNbr; ++i){
		TH1D *d2_c = (TH1D*)d2->at(i)->Clone(TString("d2_c")+proc + (Long_t)i);
		leg->AddEntry(d2_c, mcLegendTitle[i] + " FF=x");
		d2_c->Scale(gWeight*2.*a);
		stack->Add(d2_c);
	}
	for(int i=0; i<inputMCNbr; ++i){
		TH1D *d1_c = (TH1D*)d1->at(i)->Clone(TString("d1_c")+proc + (Long_t)i);
		leg->AddEntry(d1_c, mcLegendTitle[i] + " FF=1");
		d1_c->Scale(gWeight);
		stack->Add(d1_c);
	}

	TH1D* ratio = buildRatio(gWeight,a, proc);
	TF1 f("f", "[0]*(1+[1]*2.0*x+[1]*[1]*x*x)", 0.01, 0.14);
	f.SetLineColor(kRed);
	c2->Divide(1,2);
	c2->cd(1);
	stack->Draw("HIST");
	sig->DrawClone("SAMES E P");
	leg->Draw();
	c2->cd(1)->SetLogy(true);
	c2->cd(1)->SetGrid();
	fitR->Draw();
	c2->cd(2);
	ratio->Fit("f", "R");
	ratio->Draw("E P");
	c2->cd(2)->SetGrid();
}

void prepareInputHisto(){
	THStack *std1 = new THStack("std1", "Stack FF1");
	THStack *stdx = new THStack("stdx", "Stack FFx");
	THStack *stdxx = new THStack("stdxx", "Stack FFxx");
	TLegend *leg1 = new TLegend(.75,0.6,0.98,0.82);
	TLegend *legx = new TLegend(.75,0.6,0.98,0.82);
	TLegend *legxx = new TLegend(.75,0.6,0.98,0.82);

	//Color histos
	for(int i=0; i<inputMCNbr; ++i){
		d1->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i*3]));
		d2->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i*3+1]));
		d3->at(i)->SetFillColor(gStyle->GetColorPalette(mcColors[i*3+2]));
		std1->Add((TH1D*)d1->at(i)->Clone());
		stdx->Add((TH1D*)d2->at(i)->Clone());
		stdxx->Add((TH1D*)d3->at(i)->Clone());
		leg1->AddEntry(d1->at(i), mcLegendTitle[i]);
		legx->AddEntry(d2->at(i), mcLegendTitle[i]);
		legxx->AddEntry(d3->at(i), mcLegendTitle[i]);
	}
	sig->SetLineColor(kRed);

	//Plot input histo
	TCanvas *c1 = new TCanvas("cInput", "Inputs", 1600, 800);
	c1->Divide(2,2);
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
}

/************************
 * Fit procedure
 ************************/

void procedure1(fitResult& result){
	//Initialize MINUIT

	TMinuit minuit(2);
	int flag;
	double args[1];

	minuit.SetFCN(minFct);

	args[0] = 0;
	minuit.mnexcm("SET PRINTOUT", args,1,flag);
	args[0] = 1;
	minuit.mnexcm("SET ERROR", args,1,flag);
	args[0] = 2;
	minuit.mnexcm("SET STRATEGY", args,1,flag);
	minuit.mnparm(0, "G", 1., 100, 0, 0, flag);
	minuit.mnparm(1, "a", 0.05, 0.001, 0, 0, flag);
	minuit.mnparm(2, "fix", 1.0, 0, 0, 0, flag);
	minuit.FixParameter(2);
	//minuit.FixParameter(0);

	minuit.SetErrorDef(1);
	args[0] = 100000;
	minuit.mnexcm("MIGRAD", args,1,flag);

	//Get MINUIT results
	double gWeight, gWeightErr;
	double a, aErr;
	double ln0;
	double edm, errdef;
	int nvpar,nparx,icstat;

	minuit.GetParameter(0,gWeight, gWeightErr);
	minuit.GetParameter(1,a, aErr);
	minuit.mnstat(ln0,edm,errdef,nvpar,nparx,icstat);

	result.norm= gWeight;
	//result.norm = getNormalization(a);
	result.normErr = gWeightErr;
	result.formFactor = a;
	result.formFactorErr = aErr;
}

void procedure2(fitResult& result){
	other = buildRatio(1, 0, "2");

	//Initialize MINUIT

	TMinuit minuit(2);
	int flag;
	double args[1];

	minuit.SetFCN(minFct2);

	args[0] = 0;
	minuit.mnexcm("SET PRINTOUT", args,1,flag);
	args[0] = 1;
	minuit.mnexcm("SET ERROR", args,1,flag);
	args[0] = 2;
	minuit.mnexcm("SET STRATEGY", args,1,flag);
	minuit.mnparm(0, "G", 1., 100, 0, 0, flag);
	minuit.mnparm(1, "a", 10, 0.001, 0, 0, flag);
	//minuit.FixParameter(0);

	minuit.SetErrorDef(1);
	args[0] = 100000;
	minuit.mnexcm("MIGRAD", args,1,flag);

	//Get MINUIT results
	double gWeight, gWeightErr;
	double a, aErr;
	double ln0;
	double edm, errdef;
	int nvpar,nparx,icstat;

	minuit.GetParameter(0,gWeight, gWeightErr);
	minuit.GetParameter(1,a, aErr);
	minuit.mnstat(ln0,edm,errdef,nvpar,nparx,icstat);

	result.norm= gWeight;
	result.normErr = gWeightErr;
	result.formFactor = a;
	result.formFactorErr = aErr;

	TCanvas *c = new TCanvas("proc2Ratio", "Ratio before");
	TF1 *f = new TF1("f", "[0]*(1+[1]*2.0*x+[1]*[1]*x*x)", 0, 100);
	f->SetLineColor(kRed);
	f->SetParameter(0, gWeight);
	f->SetParameter(1, a);
	other->Draw("E P");
	f->Draw("SAME");
}

/************************
 * Main
 ************************/

int fit_test(TString inFile, bool fillMC=false, bool fillData=false){
	gSystem->Load("obj/libmystructs.so");
	gStyle->SetOptFit(1);

	d1 = new vector<TH1D*>;
	d2 = new vector<TH1D*>;
	d3 = new vector<TH1D*>;
	dSig = new vector<TH1D*>;

	sig = new TH1D("sig", "signal sample", BINS,0,MAX);

	//Get Input
	int nsig = 0;
	nsig = readConfig(inFile, false, fillMC, fillData);

	rebin(140);

	//Scale MC to Data
	double totalMC = 0;
	for(int i=0; i<inputMCNbr; ++i){
		totalMC += d1->at(i)->Integral();// + d2->at(i)->Integral() + d3->at(i)->Integral();
	}

	double factor = ((double)(nsig))/totalMC;
	//double factor = nsig/(d1->at(0)->Integral() + d2->at(0)->Integral() + d3->at(0)->Integral());
	for(int i=0; i<inputMCNbr; ++i){
		d1->at(i)->Scale(factor);
		d2->at(i)->Scale(factor);
		d3->at(i)->Scale(factor);
	}

	//printIntegrals(getNormalization(5), 5);
	prepareInputHisto();

	//Fit
	fitResult result1, result2;
	procedure1(result1);

	drawFitResult(result1, "1");

	result2.formFactor = 0.05863;
	result2.norm = 0.9856;
	drawFitResult(result2, "2");

	//procedure2(result2);

	//drawFitResult(result2, "2");

	cout << "######## Procedure 1 result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << result1.norm << "+-" << result1.normErr << endl;
	cout << "Slope a : " << result1.formFactor << "+-" << result1.formFactorErr << endl;

	/*cout << "######## Procedure 2 result #########" << endl << "-------------------------------------" << endl;
	cout << "Global normalization : " << result2.norm << "+-" << result2.normErr << endl;
	cout << "Slope a : " << result2.formFactor << "+-" << result2.formFactorErr << endl;*/
	return 0;
}
