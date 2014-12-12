#define __CINT_NICO__ 1

#include <TH1D.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TObjArray.h>
#include <TObjString.h>
#include <TTree.h>
#include <THStack.h>
#include <TSystem.h>
#include <TStyle.h>
//#include "userinc/mystructs.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
using namespace std;

#define BINS 100
#define MAX 0.14
#define MAXEVENTS 0
static int inputMCNbr = 0;
static int inputDataNbr = 0;
static vector< vector<TH1D*> > *d1;
static vector< vector<TH1D*> >*dSig;
static vector<int> mcColors, dataColors;
vector<TString> mcLegendTitle;
TString dataLegendTitle;
static double NSig;

TFile *tempFD;

Long_t iCanvas = 0;

int getInputDataFill(TFile *fd);
void getInputMCFill(TFile *fd, double br, int index);
void scaleMC(int selEvents, int totEvents, int index, double br);

/************************
 * Getting input
 ************************/
int readConfig(TString confFile){
	vector<double> brs;
	vector<TString> mcFileNames;
	vector<TString> histoList;
	vector<TString> dataFileNames;
	vector<int> mcIndexes;

	int nsig = 0;

	tempFD = TFile::Open(".tempComb.root", "RECREATE");

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
			else if(key.CompareTo("brs")==0) brs.push_back(entry.Atof());
			else if(key.CompareTo("mccolors")==0) mcColors.push_back(entry.Atoi());
			else if(key.CompareTo("mclegends")==0) mcLegendTitle.push_back(entry);
			else if(key.CompareTo("datafiles")==0) dataFileNames.push_back(entry);
			else if(key.CompareTo("datacolors")==0) dataColors.push_back(entry.Atoi());
			else if(key.CompareTo("datalegends")==0) dataLegendTitle = entry;
			else if(key.CompareTo("mcIndex")==0) mcIndexes.push_back(entry.Atoi());
			delete tok->At(i);
		}
	}

	//vector<TFile*> fd;
	TFile* ffd;
	int prevIndex = -1;

	//Getting MC
	for(int i=0; i<mcFileNames.size(); ++i){
		cout << mcFileNames[i] << endl;
		ffd = TFile::Open(mcFileNames[i]);
		//fd.push_back(ffd);
		if(prevIndex != mcIndexes[i]){
			prevIndex = mcIndexes[i];
			++inputMCNbr;
		}
		getInputMCFill(ffd, brs[prevIndex], prevIndex);
		ffd->Close();
	}

	//Getting data
	for(int i=0; i<dataFileNames.size(); ++i){
		cout << dataFileNames[i] << endl;
		ffd = TFile::Open(dataFileNames[i]);
		//fd.push_back(ffd);
		++inputDataNbr;

		nsig += getInputDataFill(ffd);
		ffd->Close();
	}

	NSig = nsig;

	cout << "Returning from config file" << endl;
	return nsig;
}

void addHisto(TString name, int index, vector<TH1D*> *v, int bins, int min, int max){
	TH1D* xxx1 = new TH1D(TString(name) + (Long_t)index, "sample 1", bins,min,max);
	xxx1->Sumw2();
	v->push_back(xxx1);
}

void addAllHisto(vector<TH1D*> *v, int index){
	tempFD->cd();
<<<<<<< HEAD
	addHisto("mK", index, v, 100, 0.4, 0.52);
	addHisto("RDCH", index, v, 100, -75, 75);
	addHisto("RLKr", index, v, 100, -75, 75);
	addHisto("Zvtx", index, v, 100, 0, 350000);
	addHisto("Pt", index, v, 100, 0, 0.1);
	addHisto("P", index, v, 100, 0, 100);
	
	//Photon
	addHisto("gEnergy", index, v, 60, 0, 60);
	addHisto("gPositionX", index, v, 120, -120, 120);
	addHisto("gPositionY", index, v, 120, -120, 120);
	addHisto("gP", index, v, 60, 0, 60);
	
	//e+/e-
	addHisto("ePMag", index, v, 60, 0, 60);
	addHisto("ePx", index, v, 60, -0.015, 0.015);
	addHisto("ePy", index, v, 60, -0.015, 0.015);
	addHisto("ePz", index, v, 110, 0, 1.1);
	addHisto("eVtxX", index, v, 60, -6, 6);
	addHisto("eVtxY", index, v, 40, -4, 4);
	addHisto("eVtxZ", index, v, 100, -2000, 8000);
	addHisto("eCDA", index, v, 100, 0, 10);
	addHisto("eEnergy", index, v, 60, 0, 60);
	addHisto("mee", index, v, 140, 0, 0.14);

	//pi+
	addHisto("pipPMag", index, v, 60, 0, 60);
	addHisto("pipPx", index, v, 60, -0.015, 0.015);
	addHisto("pipPy", index, v, 60, -0.015, 0.015);
	addHisto("pipPz", index, v, 110, 0, 1.1);
	addHisto("pipVtxX", index, v, 60, -6, 6);
	addHisto("pipVtxY", index, v, 40, -4, 4);
	addHisto("pipVtxZ", index, v, 100, -2000, 8000);
	addHisto("pipCDA", index, v, 100, 0, 10);
	addHisto("pipEnergy", index, v, 60, 0, 60);
=======
	addHisto("mee", index, v, 140, 0, 0.14);
	addHisto("mK", index, v, 100, 0.4, 0.52);
	addHisto("Pe", index, v, 100, 0, 75);
	addHisto("Ppi", index, v, 100, 0, 75);
	addHisto("RDCH", index, v, 100, -75, 75);
	addHisto("RLKr", index, v, 100, -75, 75);
	addHisto("Zvtx", index, v, 100, 0, 350000);
	addHisto("Pt", index, v, 100, 0, 5);
	addHisto("P", index, v, 100, 0, 100);
>>>>>>> ptp-push
}

void fillHistos(vector<TH1D*> *d, TObject *value, double weight=1.){
	pi0dEvent *evt = (pi0dEvent*)value;
<<<<<<< HEAD
	int i=-1;
	d->at(++i)->Fill(evt->mK, weight);
	//d->at(++i)->Fill(evt->, weight);
	//d->at(++i)->Fill(evt->mK, weight);
	//d->at(++i)->Fill(evt->zVtx, weight);
	++i;
	++i;
	++i;
	d->at(++i)->Fill(evt->pTotal.Perp(), weight);
	d->at(++i)->Fill(evt->pTotal.Mag(), weight);
	//Photon
	d->at(++i)->Fill(evt->gamma->energy, weight);
	d->at(++i)->Fill(evt->gamma->position.X(), weight);
	d->at(++i)->Fill(evt->gamma->position.Y(), weight);
	d->at(++i)->Fill(evt->pGamma.Mag(), weight);
	
	//e+/e-
	d->at(++i)->Fill(evt->em->pMag, weight);
	d->at(i)->Fill(evt->ep->pMag, weight);
	
	d->at(++i)->Fill(evt->ep->momentum.X(), weight);
	d->at(i)->Fill(evt->em->momentum.X(), weight);
	d->at(++i)->Fill(evt->ep->momentum.Y(), weight);
	d->at(i)->Fill(evt->em->momentum.Y(), weight);
	d->at(++i)->Fill(evt->ep->momentum.Z(), weight);
	d->at(i)->Fill(evt->em->momentum.Z(), weight);
	
	d->at(++i)->Fill(evt->ep->vertex.X(), weight);
	d->at(i)->Fill(evt->em->vertex.X(), weight);
	d->at(++i)->Fill(evt->ep->vertex.Y(), weight);
	d->at(i)->Fill(evt->em->vertex.Y(), weight);
	d->at(++i)->Fill(evt->ep->vertex.Z(), weight);
	d->at(i)->Fill(evt->em->vertex.Z(), weight);
	
	d->at(++i)->Fill(evt->ep->cda, weight);
	d->at(i)->Fill(evt->em->cda, weight);
	
	d->at(++i)->Fill(evt->ep->clusterEnergy, weight);
	d->at(i)->Fill(evt->em->clusterEnergy, weight);
	
	d->at(++i)->Fill(evt->mee, weight);
	
	//pi+
	d->at(++i)->Fill(evt->pip->pMag, weight);
	
	d->at(++i)->Fill(evt->pip->momentum.X(), weight);
	d->at(++i)->Fill(evt->pip->momentum.Y(), weight);
	d->at(++i)->Fill(evt->pip->momentum.Z(), weight);
	
	d->at(++i)->Fill(evt->pip->vertex.X(), weight);
	d->at(++i)->Fill(evt->pip->vertex.Y(), weight);
	d->at(++i)->Fill(evt->pip->vertex.Z(), weight);
	
	d->at(++i)->Fill(evt->pip->cda, weight);
	
	d->at(++i)->Fill(evt->pip->clusterEnergy, weight);
=======
	d->at(0)->Fill(evt->mee, weight);
	d->at(1)->Fill(evt->mK, weight);
	d->at(2)->Fill(evt->em->pMag, weight);
	d->at(2)->Fill(evt->ep->pMag, weight);
	d->at(3)->Fill(evt->pip->pMag, weight); //to enable next time
	//d->at(4)->Fill(evt->, weight);
	//d->at(5)->Fill(evt->mK, weight);
	//d->at(6)->Fill(evt->zVtx, weight);
	d->at(7)->Fill(evt->pTotal.Perp(), weight);
	d->at(8)->Fill(evt->pTotal.Mag(), weight);
>>>>>>> ptp-push
}

void getInputMCFill(TFile *fd, double br, int index){

	//Get the TTree
	//Input
	pi0dEvent *eventBrch = new pi0dEvent();
	TTree *t = (TTree*)fd->Get("event");
	t->SetBranchAddress("pi0dEvent", &eventBrch);

	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = ((TH1D*)fd->Get("Cuts"))->GetEntries();
	if(MAXEVENTS>0 && nevt>MAXEVENTS){
		double ratio = (double)MAXEVENTS/(double)nevt;
		nevt=MAXEVENTS;
		totalChanEvents = totalChanEvents*ratio;
	}

	//Create histo
	if(index == d1->size()){
		vector<TH1D*> v;
		addAllHisto(&v, index+1);
		d1->push_back(v);
	}

	//Read events and fill histo
	cout << "Filling " << nevt << endl;
	for(int i=0; i<nevt; ++i){
		t->GetEntry(i);
		fillHistos(&(d1->at(index)), eventBrch);
	}

	scaleMC(nevt, totalChanEvents, index, br);
}

int getInputDataFill(TFile *fd){
	//Input
	pi0dEvent *eventBrch = new pi0dEvent();
	TTree *t = (TTree*)fd->Get("event");
	t->SetBranchAddress("pi0dEvent", &eventBrch);

	// Set Number of events
	int nevt = t->GetEntries();
	if(MAXEVENTS>0 && nevt>MAXEVENTS) nevt=MAXEVENTS;
	int nsig = nevt;

	int index = 0;
	if(dSig->size()==0){
		vector<TH1D*> v;
		addAllHisto(&v, 0);
		dSig->push_back(v);
	}

	//Read event and fill histo
	int i=0;

	cout << "Filling data " << endl;
	for(i=0; i<nsig; i++){
		t->GetEntry(i);
		fillHistos(&(dSig->at(0)), eventBrch);
	}

	return nsig;
}

/************************
 * Utility
 ************************/
void scaleMC(int selEvents, int totEvents, int index, double br){
	cout << "Rescaling" << endl;

	// Rescale histo
	for(int i=0; i<d1->at(0).size(); ++i){
		d1->at(index).at(i)->Scale(br/((double)totEvents/(double)selEvents));
	}
}

/************************
 * Histograms building and drawing
 ************************/
void drawCanvas(TString name, THStack *stack, TH1* data, TLegend *leg){
	TCanvas *c1 = new TCanvas(TString("c") + iCanvas, name);

	stack->Draw("HIST");
	data->Draw("SAME E P");
	leg->Draw();

	c1->SetGrid();
	c1->Update();
	c1->Draw();
	++iCanvas;
}

void doPlot(int index, TString name, TString title, TLegend* leg, vector<int> colors, vector<TString> *legendTitle = NULL){
<<<<<<< HEAD
	tempFD->cd();
	
=======
>>>>>>> ptp-push
	THStack *hStack = new THStack(name, title);
	double brSum = 0;

	//Scale MC to Data
	double totalMC = 0;
	for(int i=0; i<inputMCNbr; ++i){
		totalMC += d1->at(i).at(index)->Integral();
		d1->at(i).at(index)->SetFillColor(gStyle->GetColorPalette(colors[i]));
		if(legendTitle) leg->AddEntry(d1->at(i).at(index),legendTitle->at(i).Data(),"f");
	}

	double factor = ((double)(dSig->at(0).at(index)->Integral()))/totalMC;
	for(int i=0; i<inputMCNbr; ++i){
		d1->at(i).at(index)->Scale(factor);
	}

	//Stack MC
	for(int i=0; i<d1->size();++i){
		hStack->Add(d1->at(i).at(index));
		d1->at(i).at(index)->Write();
	}

	dSig->at(0).at(index)->Write();

	//Style data
	dSig->at(0).at(index)->SetLineColor(kRed);
	if(legendTitle) leg->AddEntry(dSig->at(0).at(index),dataLegendTitle.Data(),"lep");

	drawCanvas(name, hStack, dSig->at(0).at(index), leg);

	hStack->Write();
}


/************************
 * Main
 ************************/

int combine_new(TString inFile){
	gSystem->Load("obj/libmystructs.so");
	gStyle->SetOptFit(1);

	d1 = new vector< vector<TH1D*> >;
	dSig = new vector< vector<TH1D*> >;
	TLegend *leg = new TLegend(.65,0.6,0.98,0.82);
	leg->SetHeader("The Legend Title");

	//Get Input
	int nsig = 0;
	nsig = readConfig(inFile);

	doPlot(0, "mee", "e+e- invariant mass", leg, mcColors, &mcLegendTitle);
	doPlot(1, "mK", "Kaon invariant mass", leg, mcColors);
	doPlot(2, "Pe", "e+/e- momentum", leg, mcColors);
	doPlot(3, "Ppi", "Pi+ momentum", leg, mcColors);
<<<<<<< HEAD
//	doPlot(4, "RDCH", "DCH Radius", leg, mcColors);
//	doPlot(5, "RLKr", "LKr radius", leg, mcColors);
//	doPlot(6, "Zvtx", "Z Vertex", leg, mcColors);
=======
	doPlot(4, "RDCH", "DCH Radius", leg, mcColors);
	doPlot(5, "RLKr", "LKr radius", leg, mcColors);
	doPlot(6, "Zvtx", "Z Vertex", leg, mcColors);
>>>>>>> ptp-push
	doPlot(7, "Pt", "Transverse momentum", leg, mcColors);
	doPlot(8, "P", "Total momentum", leg, mcColors);

	return 0;
}
