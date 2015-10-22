#define __CINT_NICO__ 1

#define PRINTVAR(v) #v << "= " << v << " "

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
#include <iomanip>
#include "../userinc/funLib.h"
using namespace std;

#define MAXEVENTS 0
Long_t iCanvas = 0;

static vector< vector<TH1D*> > *d1;
static vector< vector<TH1D*> > *dSig;
static vector< vector<TH2D*> > *dMap, *dSigMap;

ROOTRawEvent *xxx = new ROOTRawEvent();
ROOTRawEvent &rawEvent = *xxx;

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

	addHisto("R_DCH1_ep", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep", index, v, 300, -150, 150);
	addHisto("R_DCH1_em", index, v, 150, 0, 150);
	addHisto("X_DCH1_em", index, v, 300, -150, 150);
	addHisto("Y_DCH1_em", index, v, 300, -150, 150);
	addHisto("R_DCH1_pip", index, v, 150, 0, 150);
	addHisto("X_DCH1_pip", index, v, 300, -150, 150);
	addHisto("Y_DCH1_pip", index, v, 300, -150, 150);
	addHisto("R_DCH1_gamma", index, v, 150, 0, 150);
	addHisto("X_DCH1_gamma", index, v, 300, -150, 150);
	addHisto("Y_DCH1_gamma", index, v, 300, -150, 150);

	addHisto("R_DCH2_ep", index, v, 150, 0, 150);
	addHisto("X_DCH2_ep", index, v, 300, -150, 150);
	addHisto("Y_DCH2_ep", index, v, 300, -150, 150);
	addHisto("R_DCH2_em", index, v, 150, 0, 150);
	addHisto("X_DCH2_em", index, v, 300, -150, 150);
	addHisto("Y_DCH2_em", index, v, 300, -150, 150);
	addHisto("R_DCH2_pip", index, v, 150, 0, 150);
	addHisto("X_DCH2_pip", index, v, 300, -150, 150);
	addHisto("Y_DCH2_pip", index, v, 300, -150, 150);
	addHisto("R_DCH2_gamma", index, v, 150, 0, 150);
	addHisto("X_DCH2_gamma", index, v, 300, -150, 150);
	addHisto("Y_DCH2_gamma", index, v, 300, -150, 150);

	addHisto("R_DCH3_ep", index, v, 150, 0, 150);
	addHisto("X_DCH3_ep", index, v, 300, -150, 150);
	addHisto("Y_DCH3_ep", index, v, 300, -150, 150);
	addHisto("R_DCH3_em", index, v, 150, 0, 150);
	addHisto("X_DCH3_em", index, v, 300, -150, 150);
	addHisto("Y_DCH3_em", index, v, 300, -150, 150);
	addHisto("R_DCH3_pip", index, v, 150, 0, 150);
	addHisto("X_DCH3_pip", index, v, 300, -150, 150);
	addHisto("Y_DCH3_pip", index, v, 300, -150, 150);
	addHisto("R_DCH3_gamma", index, v, 150, 0, 150);
	addHisto("X_DCH3_gamma", index, v, 300, -150, 150);
	addHisto("Y_DCH3_gamma", index, v, 300, -150, 150);

	addHisto("R_DCH4_ep", index, v, 150, 0, 150);
	addHisto("X_DCH4_ep", index, v, 300, -150, 150);
	addHisto("Y_DCH4_ep", index, v, 300, -150, 150);
	addHisto("R_DCH4_em", index, v, 150, 0, 150);
	addHisto("X_DCH4_em", index, v, 300, -150, 150);
	addHisto("Y_DCH4_em", index, v, 300, -150, 150);
	addHisto("R_DCH4_pip", index, v, 150, 0, 150);
	addHisto("X_DCH4_pip", index, v, 300, -150, 150);
	addHisto("Y_DCH4_pip", index, v, 300, -150, 150);
	addHisto("R_DCH4_gamma", index, v, 150, 0, 150);
	addHisto("X_DCH4_gamma", index, v, 300, -150, 150);
	addHisto("Y_DCH4_gamma", index, v, 300, -150, 150);

	addHisto("Zvtx", index, v, 100, -2000, 9000);
	addHisto("Qvtx", index, v, 10, -5, 5);
	addHisto("CDAvtx", index, v, 100, 0, 10);
	addHisto("Pt2", index, v, 100, 0, 0.001);
	addHisto("P", index, v, 100, 68, 80);

	addHisto("Mpi0", index, v, 60, 0.120, 0.150);

	//Photon
	addHisto("gEnergy", index, v, 80, 0, 80);
	addHisto("gPositionX", index, v, 300, -150, 150);
	addHisto("gPositionY", index, v, 300, -150, 150);
	addHisto("gRadius", index, v, 150, 0, 150);
	addHisto("gP", index, v, 60, 0, 60);

	//e+/e-
	addHisto("epPMag", index, v, 60, 0, 60);
	addHisto("epPx", index, v, 60, -0.015, 0.015);
	addHisto("epPy", index, v, 60, -0.015, 0.015);
	addHisto("epPz", index, v, 110, 0, 1.1);
	addHisto("epEnergy", index, v, 60, 0, 60);
	addHisto("epeop", index, v, 100, 0, 1.5);
	addHisto("epLKrX", index, v, 300, -150, 150);
	addHisto("epLKrY", index, v, 300, -150, 150);
	addHisto("epLKrR", index, v, 300, -150, 150);

	addHisto("emPMag", index, v, 60, 0, 60);
	addHisto("emPx", index, v, 60, -0.015, 0.015);
	addHisto("emPy", index, v, 60, -0.015, 0.015);
	addHisto("emPz", index, v, 110, 0, 1.1);
	addHisto("emEnergy", index, v, 60, 0, 60);
	addHisto("emeop", index, v, 100, 0, 1.5);
	addHisto("emLKrX", index, v, 300, -150, 150);
	addHisto("emLKrY", index, v, 300, -150, 150);
	addHisto("emLKrR", index, v, 300, -150, 150);

	addHisto("mee", index, v, 140, 0, 0.14);

	//pi+
	addHisto("pipPMag", index, v, 60, 0, 60);
	addHisto("pipPx", index, v, 60, -0.015, 0.015);
	addHisto("pipPy", index, v, 60, -0.015, 0.015);
	addHisto("pipPz", index, v, 110, 0, 1.1);
	addHisto("pipEnergy", index, v, 60, 0, 60);
	addHisto("pieop", index, v, 100, 0, 1.5);


	addHisto("t_epem_DCH", index, v,  150, 0, 150);
	addHisto("t_eppip_DCH", index, v, 150, 0, 150);
	addHisto("t_empip_DCH", index, v, 150, 0, 150);
	addHisto("t_epem_LKr", index, v,  400, 0, 400);
	addHisto("t_eppip_LKr", index, v, 400, 0, 400);
	addHisto("t_empip_LKr", index, v, 400, 0, 400);

	addHisto("t_gep_DCH", index, v,  150, 0, 150);
	addHisto("t_gem_DCH", index, v,  150, 0, 150);
	addHisto("t_gpip_DCH", index, v, 150, 0, 150);
	addHisto("t_gep_LKr", index, v,  400, 0, 400);
	addHisto("t_gem_LKr", index, v,  400, 0, 400);
	addHisto("t_gpip_LKr", index, v, 400, 0, 400);

	addHisto("undeft_gep_LKr", index, v,  400, 0, 400);
	addHisto("undeft_gem_LKr", index, v,  400, 0, 400);
	addHisto("undeft_gpip_LKr", index, v, 400, 0, 400);

	addHisto("L3_E_LKr_ep", index, v, 160, 0, 80);
	addHisto("L3_E_LKr_em", index, v, 160, 0, 80);
	addHisto("L3_E_LKr_gamma", index, v, 160, 0, 80);
	addHisto("L3_E_LKr", index, v, 160, 0, 80);

	cout << "Vector size: " << v->size() << endl;
	addHisto("R_DCH1_ep_0", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_0", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_0", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_0", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_1", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_1", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_1", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_1", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_2", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_2", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_2", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_2", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_3", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_3", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_3", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_3", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_4", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_4", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_4", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_4", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_5", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_5", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_5", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_5", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_6", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_6", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_6", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_6", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_7", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_7", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_7", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_7", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_8", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_8", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_8", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_8", index, v,  150, 0, 150);
	addHisto("R_DCH1_ep_9", index, v, 150, 0, 150);
	addHisto("X_DCH1_ep_9", index, v, 300, -150, 150);
	addHisto("Y_DCH1_ep_9", index, v, 300, -150, 150);
	addHisto("t_epem_DCH1_9", index, v,  150, 0, 150);
	cout << "Vector size: " << v->size() << endl;

	addHisto("xMap", 		index, vMap, 1000,0,1, 1000, 0, 1);
	addHisto("LKr_XY_ep", 	index, vMap, 600,-120,120, 600, -120, 120);
	addHisto("LKr_XY_em", 	index, vMap, 600,-120,120, 600, -120, 120);
	addHisto("LKr_XY_pip", 	index, vMap, 600,-120,120, 600, -120, 120);
	addHisto("LKr_XY_gamma",index, vMap, 600,-120,120, 600, -120, 120);
	addHisto("DCH1_XY_ep", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH1_XY_em", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH1_XY_pip", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH1_XY_gamma", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH2_XY_ep", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH2_XY_em", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH2_XY_pip", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH2_XY_gamma", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH3_XY_ep", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH3_XY_em", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH3_XY_pip", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH3_XY_gamma", index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH4_XY_ep", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH4_XY_em", 	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH4_XY_pip",	index, vMap, 500,-100,100, 500, -100, 100);
	addHisto("DCH4_XY_gamma", index, vMap, 500,-100,100, 500, -100, 100);
}

void getAllHisto(TFile *fd, vector<TH1D*> *v, vector<TH2D*> *vMap){
	int i=-1;
	int iMap=-1;
	getHisto(fd, "mK", ++i, v);

	getHisto(fd, "R_DCH1_ep", ++i, v);
	getHisto(fd, "X_DCH1_ep", ++i, v);
	getHisto(fd, "Y_DCH1_ep", ++i, v);
	getHisto(fd, "R_DCH1_em", ++i, v);
	getHisto(fd, "X_DCH1_em", ++i, v);
	getHisto(fd, "Y_DCH1_em", ++i, v);
	getHisto(fd, "R_DCH1_pip", ++i, v);
	getHisto(fd, "X_DCH1_pip", ++i, v);
	getHisto(fd, "Y_DCH1_pip", ++i, v);
	getHisto(fd, "R_DCH1_gamma", ++i, v);
	getHisto(fd, "X_DCH1_gamma", ++i, v);
	getHisto(fd, "Y_DCH1_gamma", ++i, v);

	getHisto(fd, "R_DCH2_ep", ++i, v);
	getHisto(fd, "X_DCH2_ep", ++i, v);
	getHisto(fd, "Y_DCH2_ep", ++i, v);
	getHisto(fd, "R_DCH2_em", ++i, v);
	getHisto(fd, "X_DCH2_em", ++i, v);
	getHisto(fd, "Y_DCH2_em", ++i, v);
	getHisto(fd, "R_DCH2_pip", ++i, v);
	getHisto(fd, "X_DCH2_pip", ++i, v);
	getHisto(fd, "Y_DCH2_pip", ++i, v);
	getHisto(fd, "R_DCH2_gamma", ++i, v);
	getHisto(fd, "X_DCH2_gamma", ++i, v);
	getHisto(fd, "Y_DCH2_gamma", ++i, v);

	getHisto(fd, "R_DCH3_ep", ++i, v);
	getHisto(fd, "X_DCH3_ep", ++i, v);
	getHisto(fd, "Y_DCH3_ep", ++i, v);
	getHisto(fd, "R_DCH3_em", ++i, v);
	getHisto(fd, "X_DCH3_em", ++i, v);
	getHisto(fd, "Y_DCH3_em", ++i, v);
	getHisto(fd, "R_DCH3_pip", ++i, v);
	getHisto(fd, "X_DCH3_pip", ++i, v);
	getHisto(fd, "Y_DCH3_pip", ++i, v);
	getHisto(fd, "R_DCH3_gamma", ++i, v);
	getHisto(fd, "X_DCH3_gamma", ++i, v);
	getHisto(fd, "Y_DCH3_gamma", ++i, v);

	getHisto(fd, "R_DCH4_ep", ++i, v);
	getHisto(fd, "X_DCH4_ep", ++i, v);
	getHisto(fd, "Y_DCH4_ep", ++i, v);
	getHisto(fd, "R_DCH4_em", ++i, v);
	getHisto(fd, "X_DCH4_em", ++i, v);
	getHisto(fd, "Y_DCH4_em", ++i, v);
	getHisto(fd, "R_DCH4_pip", ++i, v);
	getHisto(fd, "X_DCH4_pip", ++i, v);
	getHisto(fd, "Y_DCH4_pip", ++i, v);
	getHisto(fd, "R_DCH4_gamma", ++i, v);
	getHisto(fd, "X_DCH4_gamma", ++i, v);
	getHisto(fd, "Y_DCH4_gamma", ++i, v);

	getHisto(fd, "Zvtx", ++i, v);
	getHisto(fd, "Qvtx", ++i, v);
	getHisto(fd, "CDAvtx", ++i, v);
	getHisto(fd, "Pt2", ++i, v);
	getHisto(fd, "P", ++i, v);

	getHisto(fd, "Mpi0", ++i, v);

	//Photon
	getHisto(fd, "gEnergy", ++i, v);
	getHisto(fd, "gPositionX", ++i, v);
	getHisto(fd, "gPositionY", ++i, v);
	getHisto(fd, "gRadius", ++i, v);
	getHisto(fd, "gP", ++i, v);

	//e+/e-
	getHisto(fd, "epPMag", ++i, v);
	getHisto(fd, "epPx", ++i, v);
	getHisto(fd, "epPy", ++i, v);
	getHisto(fd, "epPz", ++i, v);
	getHisto(fd, "epEnergy", ++i, v);
	getHisto(fd, "epeop", ++i, v);
	getHisto(fd, "epLKrX", ++i, v);
	getHisto(fd, "epLKrY", ++i, v);
	getHisto(fd, "epLKrR", ++i, v);

	getHisto(fd, "emPMag", ++i, v);
	getHisto(fd, "emPx", ++i, v);
	getHisto(fd, "emPy", ++i, v);
	getHisto(fd, "emPz", ++i, v);
	getHisto(fd, "emEnergy", ++i, v);
	getHisto(fd, "emeop", ++i, v);
	getHisto(fd, "emLKrX", ++i, v);
	getHisto(fd, "emLKrY", ++i, v);
	getHisto(fd, "emLKrR", ++i, v);

	getHisto(fd, "mee", ++i, v);

	//pi+
	getHisto(fd, "pipPMag", ++i, v);
	getHisto(fd, "pipPx", ++i, v);
	getHisto(fd, "pipPy", ++i, v);
	getHisto(fd, "pipPz", ++i, v);
	getHisto(fd, "pipEnergy", ++i, v);
	getHisto(fd, "pieop", ++i, v);

	getHisto(fd, "t_epem_DCH", ++i, v);
	getHisto(fd, "t_eppip_DCH", ++i, v);
	getHisto(fd, "t_empip_DCH", ++i, v);
	getHisto(fd, "t_epem_LKr", ++i, v);
	getHisto(fd, "t_eppip_LKr", ++i, v);
	getHisto(fd, "t_empip_LKr", ++i, v);

	getHisto(fd, "t_gep_DCH", ++i, v);
	getHisto(fd, "t_gem_DCH", ++i, v);
	getHisto(fd, "t_gpip_DCH", ++i, v);
	getHisto(fd, "t_gep_LKr", ++i, v);
	getHisto(fd, "t_gem_LKr", ++i, v);
	getHisto(fd, "t_gpip_LKr", ++i, v);

	getHisto(fd, "undeft_gep_LKr", ++i, v);
	getHisto(fd, "undeft_gem_LKr", ++i, v);
	getHisto(fd, "undeft_gpip_LKr", ++i, v);

	getHisto(fd, "L3_E_LKr_ep", ++i, v);
	getHisto(fd, "L3_E_LKr_em", ++i, v);
	getHisto(fd, "L3_E_LKr_gamma", ++i, v);
	getHisto(fd, "L3_E_LKr", ++i, v);


	getHisto(fd, "R_DCH1_ep_0", ++i, v);
	getHisto(fd, "X_DCH1_ep_0", ++i, v);
	getHisto(fd, "Y_DCH1_ep_0", ++i, v);
	getHisto(fd, "t_epem_DCH1_0", ++i, v);
	getHisto(fd, "R_DCH1_ep_1", ++i, v);
	getHisto(fd, "X_DCH1_ep_1", ++i, v);
	getHisto(fd, "Y_DCH1_ep_1", ++i, v);
	getHisto(fd, "t_epem_DCH1_1", ++i, v);
	getHisto(fd, "R_DCH1_ep_2", ++i, v);
	getHisto(fd, "X_DCH1_ep_2", ++i, v);
	getHisto(fd, "Y_DCH1_ep_2", ++i, v);
	getHisto(fd, "t_epem_DCH1_2", ++i, v);
	getHisto(fd, "R_DCH1_ep_3", ++i, v);
	getHisto(fd, "X_DCH1_ep_3", ++i, v);
	getHisto(fd, "Y_DCH1_ep_3", ++i, v);
	getHisto(fd, "t_epem_DCH1_3", ++i, v);
	getHisto(fd, "R_DCH1_ep_4", ++i, v);
	getHisto(fd, "X_DCH1_ep_4", ++i, v);
	getHisto(fd, "Y_DCH1_ep_4", ++i, v);
	getHisto(fd, "t_epem_DCH1_4", ++i, v);
	getHisto(fd, "R_DCH1_ep_5", ++i, v);
	getHisto(fd, "X_DCH1_ep_5", ++i, v);
	getHisto(fd, "Y_DCH1_ep_5", ++i, v);
	getHisto(fd, "t_epem_DCH1_5", ++i, v);
	getHisto(fd, "R_DCH1_ep_6", ++i, v);
	getHisto(fd, "X_DCH1_ep_6", ++i, v);
	getHisto(fd, "Y_DCH1_ep_6", ++i, v);
	getHisto(fd, "t_epem_DCH1_6", ++i, v);
	getHisto(fd, "R_DCH1_ep_7", ++i, v);
	getHisto(fd, "X_DCH1_ep_7", ++i, v);
	getHisto(fd, "Y_DCH1_ep_7", ++i, v);
	getHisto(fd, "t_epem_DCH1_7", ++i, v);
	getHisto(fd, "R_DCH1_ep_8", ++i, v);
	getHisto(fd, "X_DCH1_ep_8", ++i, v);
	getHisto(fd, "Y_DCH1_ep_8", ++i, v);
	getHisto(fd, "t_epem_DCH1_8", ++i, v);
	getHisto(fd, "R_DCH1_ep_9", ++i, v);
	getHisto(fd, "X_DCH1_ep_9", ++i, v);
	getHisto(fd, "Y_DCH1_ep_9", ++i, v);
	getHisto(fd, "t_epem_DCH1_9", ++i, v);
	getHisto(fd, "xMap", ++iMap, vMap);
	getHisto(fd, "LKr_XY_ep", ++iMap, vMap);
	getHisto(fd, "LKr_XY_em", ++iMap, vMap);
	getHisto(fd, "LKr_XY_pip", ++iMap, vMap);
	getHisto(fd, "LKr_XY_gamma", ++iMap, vMap);
	getHisto(fd, "DCH1_XY_ep", ++iMap, vMap);
	getHisto(fd, "DCH1_XY_em", ++iMap, vMap);
	getHisto(fd, "DCH1_XY_pip", ++iMap, vMap);
	getHisto(fd, "DCH1_XY_gamma", ++iMap, vMap);
	getHisto(fd, "DCH2_XY_ep", ++iMap, vMap);
	getHisto(fd, "DCH2_XY_em", ++iMap, vMap);
	getHisto(fd, "DCH2_XY_pip", ++iMap, vMap);
	getHisto(fd, "DCH2_XY_gamma", ++iMap, vMap);
	getHisto(fd, "DCH3_XY_ep", ++iMap, vMap);
	getHisto(fd, "DCH3_XY_em", ++iMap, vMap);
	getHisto(fd, "DCH3_XY_pip", ++iMap, vMap);
	getHisto(fd, "DCH3_XY_gamma", ++iMap, vMap);
	//getHisto(fd, "DCH4_XY_ep", ++iMap, vMap);
	//getHisto(fd, "DCH4_XY_em", ++iMap, vMap);
	//getHisto(fd, "DCH4_XY_pip", ++iMap, vMap);
	//getHisto(fd, "DCH4_XY_gamma", ++iMap, vMap);
}

void fillHistos(vector<TH1D*> *d, vector<TH2D*> *vMap, ROOTPhysicsEvent *evt, ROOTRawEvent *mrawEvent, ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom, ROOTBurst *rootBurst, double weight=1.){
	//ROOTPhysicsEvent *evt = (ROOTPhysicsEvent*)value;
	int i=-1;
	int iMap=-1;
	TVector3 propPos, propPos2, propPos3;

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);

	//if(distance2D(propPos, TVector3(0,0,0))<20) return;
	//if(distance2D(propPos2, TVector3(0,0,0))<20) return;

	d->at(++i)->Fill(evt->kaon.P.M(), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[0].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[1].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[2].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[3].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);

	d->at(++i)->Fill(evt->pic.vertex.Z(), weight);
	d->at(++i)->Fill(mrawEvent->vtx[evt->pic.parentVertex].charge, weight);
	d->at(++i)->Fill(mrawEvent->vtx[evt->pic.parentVertex].cda, weight);
	d->at(++i)->Fill(evt->kaon.P.Perp2(corrEvent->kaonMomentum), weight);
	d->at(++i)->Fill(evt->kaon.P.Vect().Mag(), weight);

	d->at(++i)->Fill(evt->pi0.P.M(), weight);

	//Photon
	d->at(++i)->Fill(evt->gamma.P.E(), weight);
	d->at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.X(), weight);
	d->at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.Y(), weight);
	d->at(++i)->Fill(distance2D(corrEvent->pCluster[evt->gamma.parentCluster].position, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(evt->gamma.P.Vect().Mag(), weight);

	//e+/e-
	d->at(++i)->Fill(evt->ep.P.Vect().Mag(), weight);
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().X(), weight);
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().Y(), weight);
	d->at(++i)->Fill(evt->ep.P.Vect().Unit().Z(), weight);
	d->at(++i)->Fill(evt->ep.P.E(), weight);
	d->at(++i)->Fill(corrEvent->pTrack[evt->ep.parentTrack].E/corrEvent->pTrack[evt->ep.parentTrack].p, weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->ep.parentTrack]);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);

	d->at(++i)->Fill(evt->em.P.Vect().Mag(), weight);
	d->at(++i)->Fill(evt->em.P.Vect().Unit().X(), weight);
	d->at(++i)->Fill(evt->em.P.Vect().Unit().Y(), weight);
	d->at(++i)->Fill(evt->em.P.Vect().Unit().Z(), weight);
	d->at(++i)->Fill(evt->em.P.E(), weight);
	d->at(++i)->Fill(corrEvent->pTrack[evt->em.parentTrack].E/corrEvent->pTrack[evt->em.parentTrack].p, weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->em.parentTrack]);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);

	d->at(++i)->Fill(evt->mee, weight);

	//pi+
	d->at(++i)->Fill(evt->pic.P.Vect().Mag(), weight);
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().X(), weight);
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().Y(), weight);
	d->at(++i)->Fill(evt->pic.P.Vect().Unit().Z(), weight);
	d->at(++i)->Fill(evt->pic.P.E(), weight);
	d->at(++i)->Fill(corrEvent->pTrack[evt->pic.parentTrack].E/corrEvent->pTrack[evt->pic.parentTrack].p, weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, propPos2), weight);
	d->at(++i)->Fill(distance2D(propPos, propPos3), weight);
	d->at(++i)->Fill(distance2D(propPos2, propPos3), weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->em.parentTrack]);
	propPos3 = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, propPos2), weight);
	d->at(++i)->Fill(distance2D(propPos, propPos3), weight);
	d->at(++i)->Fill(distance2D(propPos2, propPos3), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos2, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos3, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->em.parentTrack]);
	propPos3 = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos2, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos3, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	d->at(++i)->Fill(distance2D(propPos, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos2, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);
	d->at(++i)->Fill(distance2D(propPos3, corrEvent->pCluster[evt->gamma.parentCluster].position), weight);

	double ELKr_ep = 0;
	double ELKr_em = 0;
	bool goodPBWall;
	NPhysicsTrack t_ep = corrEvent->pTrack[evt->ep.parentTrack];
	NPhysicsTrack t_em = corrEvent->pTrack[evt->em.parentTrack];

	propPos = propagateAfter(rootGeom->Lkr.z, t_ep);
	goodPBWall = true;
	if(rootBurst->pbWall && (propPos.Y()>-33.575 && propPos.Y() < -11.850)) goodPBWall = false;
	if(t_ep.lkr_acc==0 && goodPBWall && mrawEvent->track[t_ep.trackID].dDeadCell>2.) ELKr_ep = t_ep.p;

	propPos = propagateAfter(rootGeom->Lkr.z, t_em);
	goodPBWall = true;
	if(rootBurst->pbWall && (propPos.Y()>-33.575 && propPos.Y() < -11.850)) goodPBWall = false;
	if(t_em.lkr_acc==0 && goodPBWall && mrawEvent->track[t_em.trackID].dDeadCell>2.) ELKr_em = t_em.p;

	d->at(++i)->Fill(ELKr_ep, weight);
	d->at(++i)->Fill(ELKr_em, weight);
	d->at(++i)->Fill(evt->gamma.P.E(), weight);
	d->at(++i)->Fill(ELKr_ep + ELKr_em + evt->gamma.P.E(), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	i += int((propPos.Y()+150.)/30.)*4;
	d->at(++i)->Fill(distance2D(propPos, TVector3(0,0,0)), weight);
	d->at(++i)->Fill(propPos.X(), weight);
	d->at(++i)->Fill(propPos.Y(), weight);
	d->at(++i)->Fill(distance2D(propPos, propPos2), weight);

	if(mcEvent) vMap->at(++iMap)->Fill(mcEvent->xTrue, evt->x, weight);
	else ++iMap;
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->ep.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->em.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Lkr.z, corrEvent->pTrack[evt->pic.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	vMap->at(++iMap)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.X(), corrEvent->pCluster[evt->gamma.parentCluster].position.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[0].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[1].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[2].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[2].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[3].PosChamber.z, corrEvent->pTrack[evt->pic.parentTrack]);
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[3].PosChamber.z, corrEvent->pCluster[evt->gamma.parentCluster].position, evt->gamma.P.Vect());
	vMap->at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
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
		//ROOTRawEvent *rawBrch = new ROOTRawEvent();
		ROOTRawEvent *rawBrch = xxx;
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTFileHeader *headerBrch = new ROOTFileHeader();
		ROOTMCEvent *mcEvent = 0;
		TTree *t = (TTree*)fd->Get("event");
		TTree *th = (TTree*)fd->Get("header");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();
		NGeom *geomBrch = new NGeom();

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("rawBurst", &burstBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		th->SetBranchAddress("header", &headerBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);
		th->SetBranchAddress("geom", &geomBrch);

		//Set event nb
		int nevt = t->GetEntries();
		int totalChanEvents = 0;
		for(int i=0; i<th->GetEntries(); i++){
			th->GetEntry(i);
			totalChanEvents += headerBrch->NProcessedEvents;
		}
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
			if(i % 10000 == 0) cout << setprecision(2) << i*100./(double)nevt << "% " << i << "/" << nevt << "\r";
			cout.flush();
			t->GetEntry(i);
			if(!runIncluded(burstBrch->nrun)) continue;
			weight = applyWeights(burstBrch->nrun);
			fillHistos(&(d1->at(index)), &(dMap->at(index)), eventBrch, rawBrch, corrBrch, mcEvent, geomBrch, burstBrch, weight);
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
		//ROOTRawEvent *rawBrch = new ROOTRawEvent();
		ROOTRawEvent *rawBrch = xxx;
		ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
		ROOTFileHeader *headerBrch = new ROOTFileHeader();
		ROOTMCEvent *mcEvent = 0;
		NGeom *geomBrch = new NGeom();

		TTree *t = (TTree*)fd->Get("event");
		TTree *th = (TTree*)fd->Get("header");
		if(t->GetListOfBranches()->Contains("mc")) mcEvent = new ROOTMCEvent();

		t->SetBranchAddress("pi0dEvent", &eventBrch);
		t->SetBranchAddress("rawBurst", &burstBrch);
		t->SetBranchAddress("rawEvent", &rawBrch);
		t->SetBranchAddress("corrEvent", &corrBrch);
		th->SetBranchAddress("header", &headerBrch);
		if(mcEvent) t->SetBranchAddress("mc", &mcEvent);
		th->SetBranchAddress("geom", &geomBrch);

		th->GetEntry(0);

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
			if(i % 10000 == 0) cout << setprecision(2) << i*100./(double)nevt << "% " << i << "/" << nevt << "\r";
			cout.flush();
			t->GetEntry(i);
			if(!runIncluded(burstBrch->nrun)) continue;
			fillHistos(&(dSig->at(0)), &(dSigMap->at(0)),eventBrch, rawBrch, corrBrch, mcEvent, geomBrch, burstBrch);
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

TH2D* buildRatio2(TH2D* mc, TH2D* data, TString name){
	int nbinsx = data->GetNbinsX();
	double minx = data->GetXaxis()->GetXmin();
	double maxx = data->GetXaxis()->GetXmax();
	int nbinsy = data->GetNbinsY();
	double miny = data->GetYaxis()->GetXmin();
	double maxy = data->GetYaxis()->GetXmax();
	TH2D* r = new TH2D(TString::Format("ratio_%s", name.Data()), TString::Format("ratio_%s", name.Data()), nbinsx, minx, maxx, nbinsy, miny, maxy);

	r->Sumw2();
	r->Divide(data, mc, 1, 1, "B");
	r->SetMarkerColor(kRed);
	//r->SetMaximum(r->GetMaximum()*1.1);
	//r->SetMinimum(r->GetMinimum()*0.9);
	//r->SetMaximum(0.7);
	//r->SetMinimum(1.3);
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
	ratio->GetYaxis()->SetRangeUser(0.85, 1.15);

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

TCanvas* drawCanvas(TString name, THStack *stack, TH1D* data, TLegend *leg){
	TCanvas *c1 = new TCanvas(TString::Format("c%li", iCanvas), name);
	TH1D* r = buildRatio(stack, data, name);

	prepareRatioPlot(c1, stack, leg, data, r);

	//c1->SetGrid();
	//c1->Update();
	//c1->Draw();
	++iCanvas;
	return c1;
}

void doPlot(int index, TString name, TString title, TLegend* leg, vector<int> colors, vector<TString> *legendTitle = NULL){
	if(index<0) return;
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

	TCanvas *c = drawCanvas(name, hStack, dSig->at(0).at(index), leg);

	hStack->Write();
	cout << name+".png" << endl;
	c->SaveAs(name+".png");

	c->Close();
	delete c;
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

	int nbins = temp->GetYaxis()->GetNbins()/5 - 2;
	cout << nbins << endl;
	if(nbins<=0) nbins=1;
	temp = (TH2D*)temp->Rebin2D(8, nbins);
	TH2D* tempSig = (TH2D*)dSigMap->at(0).at(index)->Rebin2D(8, nbins);

	TH2D* ratio = buildRatio2(temp, tempSig, name);
	ratio->GetZaxis()->SetRangeUser(0.85, 1.15);

	TCanvas *c = new TCanvas(TString::Format("c%li", iCanvas), name);
	c->Divide(2, 2);
	c->cd(1);
	temp->Draw("COLZ");
	c->cd(2);
	tempSig->Draw("COLZ");
	c->cd(3);
	ratio->Draw("colz");
	++iCanvas;

	temp->Write();
	cout << name+".png" << endl;
	c->SaveAs(name+".png");
	//c->Close();
	//delete c;
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

void combine_show(TString inFile, int firstPlot, int maxPlots){
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


	int i=firstPlot;

	doPlot(++i, "mK", "Kaon invariant mass", leg, mcColors, &mcLegendTitle);
	if(maxPlots--==0) return;

	doPlot(++i, "R_DCH1_ep", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH1_ep", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH1_ep", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH1_em", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH1_em", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH1_em", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH1_pip", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH1_pip", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH1_pip", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH1_gamma", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH1_gamma", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH1_gamma", "DCH1 Radius", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(++i, "R_DCH2_ep", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH2_ep", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH2_ep", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH2_em", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH2_em", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH2_em", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH2_pip", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH2_pip", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH2_pip", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH2_gamma", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH2_gamma", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH2_gamma", "DCH2 Radius", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(++i, "R_DCH3_ep", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH3_ep", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH3_ep", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH3_em", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH3_em", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH3_em", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH3_pip", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH3_pip", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH3_pip", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH3_gamma", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH3_gamma", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH3_gamma", "DCH3 Radius", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(++i, "R_DCH4_ep", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH4_ep", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH4_ep", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH4_em", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH4_em", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH4_em", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH4_pip", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH4_pip", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH4_pip", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "R_DCH4_gamma", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "X_DCH4_gamma", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;
	doPlot(++i, "Y_DCH4_gamma", "DCH4 Radius", leg, mcColors);
	if(maxPlots--==0) return;

	doPlot(++i, "Zvtx", "Vertex Z", leg, mcColors);
	doPlot(++i, "Qvtx", "Vertex Charge", leg, mcColors);
	doPlot(++i, "CDAvtx", "Vertex CDA", leg, mcColors);
	doPlot(++i, "Pt2", "Square transverse momentum", leg, mcColors);
	doPlot(++i, "P", "Total momentum", leg, mcColors);

	doPlot(++i, "Mpi0", "Reconstructed Pi0 mass", leg, mcColors);

	//Photon
	doPlot(++i, "gEnergy", "Photon energy", leg, mcColors);
	doPlot(++i, "gPositionX", "Photon LKr position (X)", leg, mcColors);
	doPlot(++i, "gPositionY", "Photon LKr position (Y)", leg, mcColors);
	doPlot(++i, "gRadius", "Photon LKr radius", leg, mcColors);
	doPlot(++i, "gP", "Photon momentum", leg, mcColors);

	//e+/e-
	doPlot(++i, "epPMag", "Electron momentum", leg, mcColors);
	doPlot(++i, "epPx", "Electron momentum (X)", leg, mcColors);
	doPlot(++i, "epPy", "Electron momentum (Y)", leg, mcColors);
	doPlot(++i, "epPz", "Electron momentum (Z)", leg, mcColors);
	doPlot(++i, "epEnergy", "Electron energy", leg, mcColors);
	doPlot(++i, "epeop", "Electron E/p", leg, mcColors);
	doPlot(++i, "epLKrX", "Electron LKr position (X)", leg, mcColors);
	doPlot(++i, "epLKrY", "Electron LKr position (Y)", leg, mcColors);
	doPlot(++i, "epLKrR", "Electron LKr radius", leg, mcColors);

	doPlot(++i, "emPMag", "Electron momentum", leg, mcColors);
	doPlot(++i, "emPx", "Electron momentum (X)", leg, mcColors);
	doPlot(++i, "emPy", "Electron momentum (Y)", leg, mcColors);
	doPlot(++i, "emPz", "Electron momentum (Z)", leg, mcColors);
	doPlot(++i, "emEnergy", "Electron energy", leg, mcColors);
	doPlot(++i, "emeop", "Electron E/p", leg, mcColors);
	doPlot(++i, "emLKrX", "Electron LKr position (X)", leg, mcColors);
	doPlot(++i, "emLKrY", "Electron LKr position (Y)", leg, mcColors);
	doPlot(++i, "emLKrR", "Electron LKr radius", leg, mcColors);

	doPlot(++i, "mee", "Di-electron invariant mass", leg, mcColors);

	//pi+
	doPlot(++i, "pipPMag", "Pion momentum", leg, mcColors);
	doPlot(++i, "pipPx", "Pion momentum (X)", leg, mcColors);
	doPlot(++i, "pipPy", "Pion momentum (Y)", leg, mcColors);
	doPlot(++i, "pipPz", "Pion momentum (Z)", leg, mcColors);
	doPlot(++i, "pipEnergy", "Pion energy", leg, mcColors);
	doPlot(++i, "pieop", "Pion E/p", leg, mcColors);

	doPlot(++i, "t_epem_DCH", "Track distance DCH1", leg, mcColors);
	doPlot(++i, "t_eppip_DCH", "Track distance DCH1", leg, mcColors);
	doPlot(++i, "t_empip_DCH", "Track distance DCH1", leg, mcColors);
	doPlot(++i, "t_epem_LKr", "Track distance LKr", leg, mcColors);
	doPlot(++i, "t_eppip_LKr", "Track distance LKr", leg, mcColors);
	doPlot(++i, "t_empip_LKr", "Track distance LKr", leg, mcColors);

	doPlot(++i, "t_gep_DCH", "Track photon distance DCH1", leg, mcColors);
	doPlot(++i, "t_gem_DCH", "Track photon distance DCH1", leg, mcColors);
	doPlot(++i, "t_gpip_DCH", "Track photon distance DCH1", leg, mcColors);
	doPlot(++i, "t_gep_LKr", "Track photon distance LKr", leg, mcColors);
	doPlot(++i, "t_gem_LKr", "Track photon distance LKr", leg, mcColors);
	doPlot(++i, "t_gpip_LKr", "Track photon distance LKr", leg, mcColors);

	doPlot(++i, "undeft_gep_LKr", "Undeflected track photon distance LKr", leg, mcColors);
	doPlot(++i, "undeft_gem_LKr", "Undeflected track photon distance LKr", leg, mcColors);
	doPlot(++i, "undeft_gpip_LKr", "Undeflected track photon distance LKr", leg, mcColors);

	doPlot(++i, "L3_E_LKr_ep", "L3 Electron energy", leg, mcColors);
	doPlot(++i, "L3_E_LKr_em", "L3 Electron energy", leg, mcColors);
	doPlot(++i, "L3_E_LKr_gamma", "L3 photon energy", leg, mcColors);
	doPlot(++i, "L3_E_LKr", "L3 energy", leg, mcColors);

	doPlot(++i, "R_DCH1_ep_0", "R_DCH1_ep_0", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_0", "X_DCH1_ep_0", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_0", "Y_DCH1_ep_0", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_0", "Y_DCH1_ep_0", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_1", "R_DCH1_ep_1", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_1", "X_DCH1_ep_1", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_1", "Y_DCH1_ep_1", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_1", "Y_DCH1_ep_1", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_2", "R_DCH1_ep_2", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_2", "X_DCH1_ep_2", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_2", "Y_DCH1_ep_2", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_2", "Y_DCH1_ep_2", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_3", "R_DCH1_ep_3", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_3", "X_DCH1_ep_3", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_3", "Y_DCH1_ep_3", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_3", "Y_DCH1_ep_3", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_4", "R_DCH1_ep_4", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_4", "X_DCH1_ep_4", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_4", "Y_DCH1_ep_4", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_4", "Y_DCH1_ep_4", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_5", "R_DCH1_ep_5", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_5", "X_DCH1_ep_5", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_5", "Y_DCH1_ep_5", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_5", "Y_DCH1_ep_5", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_6", "R_DCH1_ep_6", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_6", "X_DCH1_ep_6", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_6", "Y_DCH1_ep_6", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_6", "Y_DCH1_ep_6", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_7", "R_DCH1_ep_7", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_7", "X_DCH1_ep_7", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_7", "Y_DCH1_ep_7", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_7", "Y_DCH1_ep_7", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_8", "R_DCH1_ep_8", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_8", "X_DCH1_ep_8", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_8", "Y_DCH1_ep_8", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_8", "Y_DCH1_ep_8", leg, mcColors);
	doPlot(++i, "R_DCH1_ep_9", "R_DCH1_ep_9", leg, mcColors);
	doPlot(++i, "X_DCH1_ep_9", "X_DCH1_ep_9", leg, mcColors);
	doPlot(++i, "Y_DCH1_ep_9", "Y_DCH1_ep_9", leg, mcColors);
	doPlot(++i, "t_epem_DCH1_9", "Y_DCH1_ep_9", leg, mcColors);

	int iMap = -1;
	doPlot2(++iMap, "xMap", "x_reco vs. x_true", leg, mcColors);

	doPlot2(++iMap, "LKr_XY_ep", "Electron LKr map", leg, mcColors);
	doPlot2(++iMap, "LKr_XY_em", "Electron LKr map", leg, mcColors);
	doPlot2(++iMap, "LKr_XY_pip", "Pion LKr map", leg, mcColors);
	doPlot2(++iMap, "LKr_XY_gamma", "Photon LKr map", leg, mcColors);
	doPlot2(++iMap, "DCH1_XY_ep", "Electron DCH1 map", leg, mcColors);
	doPlot2(++iMap, "DCH1_XY_em", "Electron DCH1 map", leg, mcColors);
	doPlot2(++iMap, "DCH1_XY_pip", "Pion DCH1 map", leg, mcColors);
	doPlot2(++iMap, "DCH1_XY_gamma", "Photon DCH1 map", leg, mcColors);
	doPlot2(++iMap, "DCH2_XY_ep", "Electron DCH2 map", leg, mcColors);
	doPlot2(++iMap, "DCH2_XY_em", "Electron DCH2 map", leg, mcColors);
	doPlot2(++iMap, "DCH2_XY_pip", "Pion DCH2 map", leg, mcColors);
	doPlot2(++iMap, "DCH2_XY_gamma", "Photon DCH2 map", leg, mcColors);
	doPlot2(++iMap, "DCH3_XY_ep", "Electron DCH3 map", leg, mcColors);
	doPlot2(++iMap, "DCH3_XY_em", "Electron DCH3 map", leg, mcColors);
	doPlot2(++iMap, "DCH3_XY_pip", "Pion DCH3 map", leg, mcColors);
	doPlot2(++iMap, "DCH3_XY_gamma", "Photon DCH3 map", leg, mcColors);
	/*doPlot2(++iMap, "DCH4_XY_ep", "Electron DCH4 map", leg, mcColors);
	doPlot2(++iMap, "DCH4_XY_em", "Electron DCH4 map", leg, mcColors);
	doPlot2(++iMap, "DCH4_XY_pip", "Pion DCH4 map", leg, mcColors);
	doPlot2(++iMap, "DCH4_XY_gamma", "Photon DCH4 map", leg, mcColors);*/


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
	int firstPlots=1;
	int maxPlots=-1;
	if(argc==2) combine_batch(config);
	else{
		if(argc==4) firstPlots=atoi(argv[3]);
		if(argc==5) maxPlots=atoi(argv[4]);
		theApp = new TApplication("combine", &argc, argv);
		combine_show(config, -firstPlots, maxPlots-1);
		theApp->Run();
	}
}
