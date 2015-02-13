#define __CINT_NICO__ 1

#include "pi0DalitzHeader.h"
#include <TSystem.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCanvas.h>
#include "../userinc/exportClasses.h"
#include <TGaxis.h>
using namespace std;

#define MAXEVENTS 0
Long_t iCanvas = 0;

static vector< vector<TH1D*> > *d1;
static vector< vector<TH1D*> > *dSig;
static vector< vector<TH2D*> > *dMap, *dSigMap;

/***************************
 * Mandatory from header
 * Opening/Closing files
 ***************************/
void initNewChannel(){
}

void initNewOutput(TFile **fdo, TString fileName){
	*fdo = TFile::Open(fileName, "RECREATE");
	fitTree = new TTree("fitStruct", "fitStruct tree");
	fitTree->Branch("fitStruct", &fitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	fitBrch.n1 = 0;
	fitBrch.nx = 0;
	fitBrch.nxx = 0;
	fitBrch.selEvents = 0;
	fitBrch.totEvents = 0;

}

void closeMCOutput(TFile *fdo, int index){
	if(index!=-1){
		fdo->cd();
		fitTree->Fill();
		fitTree->Write();

		for(unsigned int i=0; i<d1->at(index).size(); ++i){
			d1->at(index).at(i)->Write();
		}

		for(unsigned int i=0; i<dMap->at(index).size(); ++i){
			dMap->at(index).at(i)->Write();
		}

		fdo->Close();

		for(unsigned int i=0; i<d1->at(index).size(); ++i){
			d1->at(index).at(i)->SetName(TString::Format("%s%i", d1->at(index).at(i)->GetName(), index));
		}

		for(unsigned int i=0; i<dMap->at(index).size(); ++i){
			dMap->at(index).at(i)->SetName(TString::Format("%s%i", dMap->at(index).at(i)->GetName(), index));
		}
	}
}

void closeDataOutput(TFile *fdo, int index){
	fdo->cd();
	fitTree->Fill();
	fitTree->Write();

	for(unsigned int i=0; i<dSig->at(index).size(); ++i){
		dSig->at(index).at(i)->Write();
	}
	for(unsigned int i=0; i<dSigMap->at(index).size(); ++i){
		dSigMap->at(index).at(i)->Write();
	}
	fdo->Close();

	for(unsigned int i=0; i<dSig->at(index).size(); ++i){
		dSig->at(index).at(i)->SetName(TString::Format("%s%i", dSig->at(index).at(i)->GetName(), index));
	}
	for(unsigned int i=0; i<dSigMap->at(index).size(); ++i){
		dSigMap->at(index).at(i)->SetName(TString::Format("%s%i", dSigMap->at(index).at(i)->GetName(), index));
	}
}

/*******************************
 * Histogram creating/filling/fetching
 *******************************/
void addHisto(TString name, int index, vector<TH1D*> *v, int bins, double min, double max){
	TH1D* xxx1 = new TH1D(name, "sample 1", bins,min,max);
	xxx1->Sumw2();
	v->push_back(xxx1);
}

void addHisto(TString name, int index, vector<TH2D*> *v, int binsx, double minx, double maxx, int binsy, double miny, double maxy){
	TH2D* xxx1 = new TH2D(name, "sample 1", binsx,minx,maxx, binsy,miny,maxy);
	xxx1->Sumw2();
	v->push_back(xxx1);
}

void getHisto(TFile* fd, TString name, unsigned int index, vector<TH1D*> *v){
	TH1D* xxx = (TH1D*)fd->Get(name);
	if(v->size()==index){
		xxx->SetName(TString::Format("%s%i", name.Data(), index));
		tempFD->cd();
		v->push_back((TH1D*)xxx->Clone());
	}
	else{
		v->at(index)->Add(xxx, 1.);
	}
}
void getHisto(TFile* fd, TString name, unsigned int index, vector<TH2D*> *v){
	TH2D* xxx = (TH2D*)fd->Get(name);
	if(v->size()==index){
		xxx->SetName(TString::Format("%s%i", name.Data(), index));
		tempFD->cd();
		v->push_back((TH2D*)xxx->Clone());
	}
	else{
		v->at(index)->Add(xxx, 1.);
	}
}

void addAllHisto(vector<TH1D*> *v, vector<TH2D*> *vMap, int index){
	tempFD->cd();
	addHisto("mK", index, v, 100, 0.47, 0.52);
	addHisto("RDCH", index, v, 100, -75, 75);
	addHisto("RLKr", index, v, 100, -75, 75);
	addHisto("Zvtx", index, v, 100, 0, 350000);
	addHisto("Pt2", index, v, 100, 0, 0.001);
	addHisto("P", index, v, 100, 68, 80);
	addHisto("Mpi0diff", index, v, 100, -0.01, 0.01);
	
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

	addHisto("xMap", index, vMap, 1000,0,1, 1000, 0, 1);
}

void getAllHisto(TFile *fd, vector<TH1D*> *v, vector<TH2D*> *vMap){
	int i=-1;
	int iMap=-1;
	getHisto(fd, "mK", ++i, v);
	getHisto(fd, "RDCH", ++i, v);
	getHisto(fd, "RLKr", ++i, v);
	getHisto(fd, "Zvtx", ++i, v);
	getHisto(fd, "Pt2", ++i, v);
	getHisto(fd, "P", ++i, v);
	getHisto(fd, "Mpi0diff", ++i, v);
	
	//Photon
	getHisto(fd, "gEnergy", ++i, v);
	getHisto(fd, "gPositionX", ++i, v);
	getHisto(fd, "gPositionY", ++i, v);
	getHisto(fd, "gP", ++i, v);
	
	//e+/e-
	getHisto(fd, "ePMag", ++i, v);
	getHisto(fd, "ePx", ++i, v);
	getHisto(fd, "ePy", ++i, v);
	getHisto(fd, "ePz", ++i, v);
	getHisto(fd, "eVtxX", ++i, v);
	getHisto(fd, "eVtxY", ++i, v);
	getHisto(fd, "eVtxZ", ++i, v);
	getHisto(fd, "eCDA", ++i, v);
	getHisto(fd, "eEnergy", ++i, v);
	getHisto(fd, "mee", ++i, v);

	//pi+
	getHisto(fd, "pipPMag", ++i, v);
	getHisto(fd, "pipPx", ++i, v);
	getHisto(fd, "pipPy", ++i, v);
	getHisto(fd, "pipPz", ++i, v);
	getHisto(fd, "pipVtxX", ++i, v);
	getHisto(fd, "pipVtxY", ++i, v);
	getHisto(fd, "pipVtxZ", ++i, v);
	getHisto(fd, "pipCDA", ++i, v);
	getHisto(fd, "pipEnergy", ++i, v);

	getHisto(fd, "xMap", ++iMap, vMap);
}

void fillHistos(vector<TH1D*> *d, vector<TH2D*> *vMap, ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvent, ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, double weight=1.){
	//ROOTPhysicsEvent *evt = (ROOTPhysicsEvent*)value;
	int i=-1;
	int iMap=-1;

	d->at(++i)->Fill(evt->kaon.P.M(), weight);
	//d->at(++i)->Fill(evt->, weight);
	//d->at(++i)->Fill(evt->mK, weight);
	//d->at(++i)->Fill(evt->zVtx, weight);
	++i;
	++i;
	++i;
	d->at(++i)->Fill(evt->kaon.P.Perp2(corrEvent->kaonMomentum), weight);
	d->at(++i)->Fill(evt->kaon.P.Vect().Mag(), weight);
	d->at(++i)->Fill(evt->pi0.P.M()-0.1349766, weight);
	//Photon
	d->at(++i)->Fill(evt->gamma.P.E(), weight);
	d->at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.X(), weight);
	d->at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.Y(), weight);
	d->at(++i)->Fill(evt->gamma.P.Vect().Mag(), weight);
	
	//e+/e-
	d->at(++i)->Fill(evt->em.P.Vect().Mag(), weight);
	d->at(i)->Fill(evt->ep.P.Vect().Mag(), weight);
	
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().X(), weight);
	d->at(i)->Fill(evt->em.P.Vect().Unit().X(), weight);
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().Y(), weight);
	d->at(i)->Fill(evt->em.P.Vect().Unit().Y(), weight);
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().Z(), weight);
	d->at(i)->Fill(evt->em.P.Vect().Unit().Z(), weight);
	
	d->at(++i)->Fill(evt->ep.vertex.X(), weight);
	d->at(i)->Fill(evt->em.vertex.X(), weight);
	d->at(++i)->Fill(evt->ep.vertex.Y(), weight);
	d->at(i)->Fill(evt->em.vertex.Y(), weight);
	d->at(++i)->Fill(evt->ep.vertex.Z(), weight);
	d->at(i)->Fill(evt->em.vertex.Z(), weight);
	
	d->at(++i)->Fill(rawEvent->vtx[evt->ep.parentVertex].cda, weight);
	d->at(i)->Fill(rawEvent->vtx[evt->em.parentVertex].cda, weight);
	
	d->at(++i)->Fill(corrEvent->pCluster[evt->ep.parentCluster].E, weight);
	d->at(i)->Fill(corrEvent->pCluster[evt->ep.parentCluster].E, weight);
	
	d->at(++i)->Fill(evt->mee, weight);
	
	//pi+
	d->at(++i)->Fill(evt->pic.P.Vect().Mag(), weight);
	
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().X(), weight);
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().Y(), weight);
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().Z(), weight);
	
	d->at(++i)->Fill(evt->pic.vertex.X(), weight);
	d->at(++i)->Fill(evt->pic.vertex.Y(), weight);
	d->at(++i)->Fill(evt->pic.vertex.Z(), weight);
	
	d->at(++i)->Fill(rawEvent->vtx[evt->pic.parentVertex].cda, weight);
	
	d->at(++i)->Fill(corrEvent->pCluster[evt->pic.parentCluster].E, weight);

	if(mcEvent) vMap->at(++iMap)->Fill(mcEvent->xTrue, evt->x, weight);
}

/*************************
 * Utility
 *************************/
//void scaleMC(int selEvents, int totEvents, int index, double br){
void scaleMC(fitStruct N, int index, double br){
	cout << "Rescaling" << endl;
	cout << br << " " << N.selEvents << " " << N.totEvents << endl;
	// Rescale histo
	for(unsigned int i=0; i<d1->at(0).size(); ++i){
		scale(d1->at(index).at(i), 1., N.totEvents, br);
	}
	for(unsigned int i=0; i<dMap->at(0).size(); ++i){
		scale(dMap->at(index).at(i), 1., N.totEvents, br);
	}
	cout << d1->at(index).at(0)->Integral() << endl;
}

/***************************
 * Input
 ****************************/
namespace Input{
	void getInputMCFill(TFile *fd, TFile *fdout, double br, unsigned int index){

		//Get the TTree
		//Input
		ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
		ROOTBurst *burstBrch = new ROOTBurst();
		ROOTRawEvent *rawBrch = new ROOTRawEvent();
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTFileHeader *headerBrch = new ROOTFileHeader();
		ROOTMCEvent *mcEvent = 0;
		TTree *t = (TTree*)fd->Get("event");
		TTree *th = (TTree*)fd->Get("header");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("rawBurst", &burstBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		th->SetBranchAddress("header", &headerBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);

		th->GetEntry(0);
		//Set event nb
		int nevt = t->GetEntries();
		int totalChanEvents = headerBrch->NProcessedEvents;
		int processedEvents = 0;
		nevt = (MAXEVENTS>0) ? min(MAXEVENTS, nevt) : nevt;

		//Create histo
		if(index == d1->size()){
			vector<TH1D*> v;
			vector<TH2D*> vMap;
			addAllHisto(&v, &vMap, index+1);
			d1->push_back(v);
			dMap->push_back(vMap);
		}

		//Read events and fill histo
		cout << "Filling " << nevt << endl;
		double weight = 1.;
		for(int i=0; i<nevt; ++i){
			t->GetEntry(i);
			if(!runIncluded(burstBrch->nrun)) continue;
			weight = applyWeights(burstBrch->nrun);
			fillHistos(&(d1->at(index)), &(dMap->at(index)), eventBrch, rawBrch, corrBrch, mcEvent, weight);
			processedEvents++;
		}

		if(processedEvents != t->GetEntries()){
			double ratio = (double)processedEvents/(double)nevt;
			totalChanEvents = totalChanEvents*ratio;
		}
		fitBrch.totEvents += totalChanEvents;
		fitBrch.selEvents += nevt;

		//scaleMC(nevt, totalChanEvents, index, br);
	}

	int getInputDataFill(TFile *fd, TFile* fdout){
		//Input
		ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
		ROOTBurst *burstBrch = new ROOTBurst();
		ROOTRawEvent *rawBrch = new ROOTRawEvent();
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTFileHeader *headerBrch = new ROOTFileHeader();
		ROOTMCEvent *mcEvent = 0;

		TTree *t = (TTree*)fd->Get("event");
		TTree *th = (TTree*)fd->Get("header");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("rawBurst", &burstBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		th->SetBranchAddress("header", &headerBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);

		// Set Number of events
		int nevt = t->GetEntries();
		nevt = (MAXEVENTS>0) ? min(MAXEVENTS, nevt) : nevt;
		int processedEvents = 0;

		if(dSig->size()==0){
			vector<TH1D*> v;
			vector<TH2D*> vMap;
			addAllHisto(&v, &vMap, 0);
			dSig->push_back(v);
			dSigMap->push_back(vMap);
		}

		//Read event and fill histo
		int i=0;

		cout << "Filling data " << nevt << endl;
		for(i=0; i<nevt; i++){
			t->GetEntry(i);
			if(!runIncluded(burstBrch->nrun)) continue;
			fillHistos(&(dSig->at(0)), &(dSigMap->at(0)),eventBrch, rawBrch, corrBrch, mcEvent);
			processedEvents++;
		}

		fitBrch.selEvents += processedEvents;

		return processedEvents;
	}

	void getInputMCGet(TFile *fd, double br, unsigned int index){
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*)fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit);

		//Set event nb
		//int nevt = totFit.selEvents;
		//int totalChanEvents = totFit.totEvents;

		//Create histo
		//index = d1->size();
		if(index==d1->size()){
			vector<TH1D*> v;
			vector<TH2D*> vMap;
			getAllHisto(fd, &v, &vMap);
			d1->push_back(v);
			dMap->push_back(vMap);
		}
		else{
			getAllHisto(fd, &d1->at(index), &dMap->at(index));
		}

		scaleMC(totFit, index, br);
	}

	int getInputDataGet(TFile *fd){
		fitStruct fitBrch, totFit;
		TTree *t = (TTree*)fd->Get("fitStruct");
		t->SetBranchAddress("fitStruct", &fitBrch);

		initFitStruct(totFit);
		sumTreeFitStruct(fitBrch, t, totFit);

		//Set event nb
		int nsig = totFit.selEvents;

		//Create histo
		//int index = dSig->size();
		unsigned int index=0;
		if(index == dSig->size()){
			vector<TH1D*> v;
			vector<TH2D*> vMap;
			getAllHisto(fd, &v, &vMap);
			dSig->push_back(v);
			dSigMap->push_back(vMap);
		}
		else{
			getAllHisto(fd, &dSig->at(index), &dSigMap->at(index));
		}

		return nsig;
	}
}
/************************
 * Real job
 ************************/
TH1D* buildRatio(THStack* stack, TH1D* data, TString name){
	int nbins = data->GetNbinsX();
	double min = data->GetXaxis()->GetXmin();
	double max = data->GetXaxis()->GetXmax();
	TH1D* sum = new TH1D("sum", "sum", nbins, min, max);
	TH1D* r = new TH1D(TString::Format("ratio_%s", name.Data()), TString::Format("ratio_%s", name.Data()), nbins, min, max);

	for(int i=0; i<stack->GetHists()->GetEntries(); ++i){
		sum->Add((TH1D*)stack->GetHists()->At(i));
	}
	sum->SetFillColor(8);

	sum->Sumw2();
	r->Sumw2();
	r->Divide(data, sum, 1, 1, "B");
	r->SetMarkerColor(kRed);
	//r->SetMaximum(r->GetMaximum()*1.1);
	//r->SetMinimum(r->GetMinimum()*0.9);
	//r->SetMaximum(0.7);
	//r->SetMinimum(1.3);
	delete sum;
	return r;
}

void prepareRatioPlot(TCanvas *c, THStack* mc, TLegend *leg, TH1D* data, TH1D* ratio){
	//double minx = data->GetXaxis()->GetXmin();
	//double maxx = data->GetXaxis()->GetXmax();
	//double binsX = data->GetXaxis()->GetNbins();
	//double miny = data->GetMinimum();
	//double maxy = data->GetMaximum();
	//double binsY = data->GetYaxis()->GetNbins();

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(3);
	pad1->SetGrid();
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	mc->Draw("HIST");               // Draw h1
	data->Draw("SAME E P");         // Draw h2 on top of h1
	//mc->GetYaxis()->SetLabelSize(0.05);
	leg->Draw();


	c->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();
	ratio->SetStats(0);      // No statistics on lower plot
	ratio->Draw("ep");
	ratio->SetMarkerColor(kRed);
	ratio->GetYaxis()->SetRangeUser(0.5, 1.5);

	// Y axis mc plot settings
	mc->GetYaxis()->SetTitleSize(20);
	mc->GetYaxis()->SetTitleFont(43);
	mc->GetYaxis()->SetTitleOffset(1.55);

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

void drawCanvas(TString name, THStack *stack, TH1D* data, TLegend *leg){
	TCanvas *c1 = new TCanvas(TString::Format("c%li", iCanvas), name);
	TH1D* r = buildRatio(stack, data, name);

	prepareRatioPlot(c1, stack, leg, data, r);

	//c1->SetGrid();
	//c1->Update();
	//c1->Draw();
	++iCanvas;
}

void doPlot(int index, TString name, TString title, TLegend* leg, vector<int> colors, vector<TString> *legendTitle = NULL){
	tempFD->cd();
	
	THStack *hStack = new THStack(name, title);

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
	for(unsigned int i=0; i<d1->size();++i){
		hStack->Add(d1->at(i).at(index));
		d1->at(i).at(index)->Write();
	}

	dSig->at(0).at(index)->Write();

	//Style data
	dSig->at(0).at(index)->SetLineColor(kRed);
	if(legendTitle) leg->AddEntry(dSig->at(0).at(index),dataLegendTitle[0].Data(),"lep");
	
	drawCanvas(name, hStack, dSig->at(0).at(index), leg);

	hStack->Write();
}

void doPlot2(int index, TString name, TString title, TLegend* leg, vector<int> colors, vector<TString> *legendTitle = NULL){
	tempFD->cd();

	TH2D *temp = (TH2D*)dMap->at(0).at(index)->Clone(name);
	temp->SetTitle(title);
	temp->Clear();

	//Scale MC to Data
	double totalMC = 0;
	for(int i=0; i<inputMCNbr; ++i){
		totalMC += dMap->at(i).at(index)->Integral();
		dMap->at(i).at(index)->SetFillColor(gStyle->GetColorPalette(colors[i]));
		if(legendTitle) leg->AddEntry(dMap->at(i).at(index),legendTitle->at(i).Data(),"f");
	}


	double factor = ((double)(dSig->at(0).at(0)->Integral()))/totalMC;
	for(int i=0; i<inputMCNbr; ++i){
		dMap->at(i).at(index)->Scale(factor);
	}

	//Stack MC
	for(unsigned int i=0; i<dMap->size();++i){
		//hStack->Add(dMap->at(i).at(index));
		temp->Add(dMap->at(i).at(index), 1);
		dMap->at(i).at(index)->Write();
	}

	new TCanvas(TString::Format("c%li", iCanvas), name);
	temp->Draw("COLZ");
	++iCanvas;

	temp->Write();
}

/*****************************
 * Mains
 *****************************/
void combine_batch(TString inFile){
	srand(time(NULL));
	gStyle->SetOptFit(1);

	tempFileName = ".tempComb";
	tempFileName += rand() % 99999;
	tempFileName += ".root";
	tempFD = TFile::Open(tempFileName, "RECREATE");

	d1 = new vector< vector<TH1D*> >;
	dMap = new vector< vector<TH2D*> >;
	dSig = new vector< vector<TH1D*> >;
	dSigMap = new vector< vector<TH2D*> >;
	TLegend *leg = new TLegend(.65,0.6,0.98,0.82);
	leg->SetHeader("The Legend Title");

	//Get Input
	loadWeights("/afs/cern.ch/user/n/nlurkin/Compact/pi0dalitz_weights.dat");
	if(!readConfig(inFile)) return;

	readFilesFill();

	tempFD->Close();
	remove(tempFileName.Data());
}

void combine_show(TString inFile, int maxPlots){
	srand(time(NULL));
	gStyle->SetOptFit(1);

	tempFileName = ".tempComb";
	tempFileName += rand() % 99999;
	tempFileName += ".root";
	tempFD = TFile::Open(tempFileName, "RECREATE");

	d1 = new vector< vector<TH1D*> >;
	dMap = new vector< vector<TH2D*> >;
	dSig = new vector< vector<TH1D*> >;
	dSigMap = new vector< vector<TH2D*> >;
	//TLegend *leg = new TLegend(.65,0.6,0.98,0.82);
	TLegend *leg = new TLegend(.7,0.75,0.9,0.95);
	//leg->SetHeader("The Legend Title");

	//Get Input
	readConfig(inFile);

	readFilesGet();

	tempFD->cd();
	doPlot(0, "mK", "Kaon invariant mass", leg, mcColors, &mcLegendTitle);
	if(maxPlots--==0) return;
//	doPlot(1, "RDCH", "DCH Radius", leg, mcColors);
//	doPlot(2, "RLKr", "LKr radius", leg, mcColors);
//	doPlot(3, "Zvtx", "Z Vertex", leg, mcColors);
	doPlot(4, "Pt", "Squared transverse momentum", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(5, "P", "Total momentum", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(6, "Mpi0diff", "Reconstructed pi0 mass difference", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(7, "gEnergy", "Photon cluster energy", leg, mcColors);
	cout << maxPlots << endl;
	if(maxPlots--==0) return;
	doPlot(8, "gPositionX", "Photon cluster position X", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(9, "gPositionY", "Photon cluster position Y", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(10, "gP", "Photon momentum", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(11, "ePMag", "e+/e- momentum", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(12, "ePx", "e+/e- momentum X", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(13, "ePy", "e+/e- momentum Y", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(14, "ePz", "e+/e- momentum Z", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(15, "eVtxX", "e+/e- vertex X", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(16, "eVtxY", "e+/e- vertex Y", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(17, "eVtxZ", "e+/e- vertex Z", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(18, "eCDA", "e+/e- CDA", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(19, "eEnergy", "e+/e- cluster energy", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(20, "mee", "e+/e- invariant mass", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(21, "pipPMag", "Pi+ momentum", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(22, "pipPx", "Pi+ momentum X", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(23, "pipPy", "Pi+ momentum Y", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(24, "pipPz", "Pi+ momentum Z", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(25, "pipVtxX", "Vertex X", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(26, "pipVtxY", "Vertex Y", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(27, "pipVtxZ", "Vertex Z", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(28, "pipCDA", "CDA", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(29, "pipEnergy", "Pi+ cluster energy", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot2(0, "xMap", "x_reco vs. x_true", leg, mcColors);
	if(maxPlots--==0) return;
}

int main(int argc, char **argv){
	if(argc<2){
		cout << "Missing parameter" << endl;
		return -1;
	}

	setStyle();

	signal(SIGTERM, sighandler);
	signal(SIGINT, sighandler);
	signal(SIGABRT, sighandler);

	TString config(argv[1]);
	int maxPlots=-1;
	if(argc==2) combine_batch(config);
	else{
		if(argc==4) maxPlots=atoi(argv[3]);
		theApp = new TApplication("combine", &argc, argv);
		combine_show(config, maxPlots-1);
		theApp->Run();
	}
}
