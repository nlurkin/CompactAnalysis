
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <TStyle.h>

using namespace std;

Long_t i = 0;
int mcInputs = 0;

#define MAXEVENTS 0;
#define BINS 10000
#define MAX 0.14
static vector<TH1D*> *d1;

void scaleHisto(TH1* h, TH1* ref, double brRatio){
	double nEvents = ref->GetEntries();
	h->Scale(brRatio);
	double factor = nEvents/h->Integral();
	h->Scale(factor);
}

void drawCanvas(TString name, THStack *stack, TH1* data, TLegend *leg){
	TCanvas *c1 = new TCanvas(TString("c") + i, name);

	stack->Draw();
	data->Draw("SAME E P");
	leg->Draw();

	c1->SetGrid();
	c1->Update();
	c1->Draw();
	++i;
}

void doPlot(TString name, TString title, vector<int> nEvents, vector<TFile *> fd, TLegend* leg, vector<int> colors, vector<double> brs, vector<TString> *legendTitle = NULL){
	cout << "Start" << endl;
	vector<TH1D*> histo;
	THStack *hStack = new THStack(name, title);
	double brSum = 0;

	cout << "Getting" << endl;
	for(int i=0; i<fd.size(); ++i){
		histo.push_back((TH1D*)fd[i]->Get(name));
		histo[i]->SetFillColor(gStyle->GetColorPalette(colors[i]));
		if(legendTitle) leg->AddEntry(histo[i],legendTitle->at(i).Data(),"f");
		brSum += brs[i];
	}

	//Normalize MC ratios
	cout << "Normalizing to MC Ratios" << endl;
	double totMCIntegral = 0;
	for(int i=0; i<fd.size()-1;++i){
		histo[i]->Scale(brs[i]/nEvents[i]);
		totMCIntegral += histo[i]->Integral();
	}

	//Normalize to data
	cout << "Normalizing to Data" << endl;
	double dataEvents = histo[mcInputs]->GetEntries();
	double factor = dataEvents/totMCIntegral;
	for(int i=0; i<fd.size()-1;++i){
		histo[i]->Scale(factor);
	}

	//Stack MC
	cout << "Stacking" << endl;
	for(int i=0; i<fd.size()-1;++i){
		hStack->Add(histo[i]);
	}

	//Style data
	cout << "Styling" << endl;
	histo[mcInputs]->SetLineColor(kRed);
	//histo[2]->SetMarkerStyle(21);
	//histo[2]->SetMarkerSize(1);
	if(legendTitle) leg->AddEntry(histo[mcInputs],legendTitle->at(mcInputs).Data(),"lep");

	cout << "Drawing" << endl;
	drawCanvas(name, hStack, histo[mcInputs], leg);
	cout << "End" << endl;
}

doPlotFromTree(){

}


void combine(TString mcList, TString data){

	vector<double> brs;// = {2.066E-1, 3.353E-2, 0};
	vector<TString> fileNames;
	vector<TString> legendTitle;
	vector<int> colors;// = {8,9,2};
	vector<int> nEvents;

	loadWeights("pi0dalitz_weights_p1.dat");
	//read list file
	ifstream listFile(mcList.Data());

	string line;
	while(getline(listFile,line)){
		TString tline(line);
		if(tline.BeginsWith("#")) continue;
		TString key(tline(0,tline.First('=')));
		TString values(tline(tline.First('=')+1, 1000));
		TObjArray* tok = values.Tokenize(" ");
		for(int i=0; i<tok->GetEntries(); ++i){
			TString entry(((TObjString*)tok->At(i))->GetString());
			if(key.CompareTo("mcfiles")==0){
				fileNames.push_back(entry);
				++mcInputs;
			}
			else if(key.CompareTo("brs")==0) brs.push_back(entry.Atof());
			else if(key.CompareTo("mccolors")==0) colors.push_back(entry.Atoi());
			else if(key.CompareTo("mclegends")==0) legendTitle.push_back(entry);
			delete tok->At(i);
		}
	}
	fileNames.push_back(data);
	brs.push_back(0);
	colors.push_back(2);
	legendTitle.push_back("Data");

	vector<TFile*> fd;
	TFile* ffd;
	for(int i=0; i<fileNames.size(); ++i){
		cout << "opening " << fileNames[i] << endl;
		ffd = TFile::Open(fileNames[i]);
		cout << "opened" << endl;
		fd.push_back(ffd);
		nEvents.push_back(((TH1D*)ffd->Get("Cuts"))->GetEntries());
		cout << "events read" << endl;
	}

	//Create Legend
	TLegend *leg = new TLegend(.65,0.6,0.98,0.82);
	leg->SetHeader("The Legend Title");

	doPlot("sel_NVertices", "Number of vertices", nEvents, fd, leg, colors, brs, &legendTitle);
	doPlot("sel_NTracks", "Number of tracks", nEvents, fd, leg, colors, brs);
	doPlot("sel_DCH1Rad", "DCH1 radius", nEvents, fd, leg, colors, brs);
	doPlot("eeMass", "ee mass", nEvents, fd, leg, colors, brs);
	doPlot("eegMass", "eeg mass", nEvents, fd, leg, colors, brs);
	doPlot("peegMass", "peeg mass", nEvents, fd, leg, colors, brs);
	doPlot("xDistrib", "x distribution", nEvents, fd, leg, colors, brs);


	pi0dEvent *eventBrch = new pi0dEvent();
	TTree *t = (TTree*)fd[0]->Get("event");
	t->SetBranchAddress("pi0dEvent", &eventBrch);

	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = ((TH1D*)fd[0]->Get("Cuts"))->GetEntries());
	int npart;
	int divider;
	if(MAXEVENTS>0 && nevt>MAXEVENTS){
		double ratio = (double)MAXEVENTS/(double)nevt;
		nevt=MAXEVENTS;
		totalChanEvents = totalChanEvents*ratio;
	}

	//Create histo
	int index = d1->size();
	TH1D* xxx1 = new TH1D("d1", "sample 1", BINS,0,MAX);
	xxx1->Sumw2();
	d1->push_back(xxx1);

	//Read events and fill histo
	int i=0;

	/**
	 * New
	 */
	cout << "Filling" << endl;
	for(; i<nevt; ++i){
		t->GetEntry(i);
		//cout << "Timestamp: " << eventBrch->timestamp << " % " << divider << " " << mod << endl;
		d1->at(index)->Fill(eventBrch->mee);
	}

	//scaleMC(fitBrch, index, br);
}
