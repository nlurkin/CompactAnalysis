/*
 * make_bins.C
 *
 *  Created on: Jan 13, 2015
 *      Author: ncl
 */

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <set>
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
#include "../userinc/mystructs.h"
using namespace std;

#define BINS 10000000
#define MAXEVENTS 0
static double bins[BINS];
static double NSig;
double Mpi0 = 0.1349766;

multiset<double> setVals;

/***************************
 * Mandatory from header
 * Opening/Closing files
 ***************************/
void initNewChannel() {
}

void initNewOutput(TFile **fdo, TString fileName) {
	//*fdo = TFile::Open(fileName, "RECREATE");
}

void closeMCOutput(TFile *fdo, int index) {
}

void closeDataOutput(TFile *fdo, int index) {
	//fdo->cd();
	//fdo->Close();
}

/*************************
 * Utility
 *************************/
/***************************
 * Input
 ****************************/
namespace Input {
	void getInputMCFill(TFile *fd, TFile *fdout, double br, unsigned int index) {
	}

	int getInputDataFill(TFile *fd, TFile* fdout) {
		//Input
		pi0dEvent *eventBrch = new pi0dEvent();
		TTree *t = (TTree*) fd->Get("event");
		t->SetBranchAddress("pi0dEvent", &eventBrch);

		tempFD->cd();

		// Set Number of events
		int nevt = t->GetEntries();
		if (MAXEVENTS > 0 && nevt > MAXEVENTS)
			nevt = MAXEVENTS;
		NSig = nevt;

		//Read event and fill histo
		int i = 0;
		double x;

		cout << "Filling data " << endl;
		for (i = 0; i < NSig; i++) {
			t->GetEntry(i);
			x = pow(eventBrch->mee / Mpi0, 2.);

			setVals.insert(x);
		}
		return NSig;
	}

	void getInputMCGet(TFile *fd, double br, unsigned int index) {

	}

	int getInputDataGet(TFile *fd) {
		return 0;
	}
}

/************************
 * Mains
 ************************/

void make_bins(TString inFile) {
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

	//Get Input
	readConfig(inFile);

	mcFileNames.clear();
	readFilesFill();

	int binAmount = 150;
	int nbins = ceil(setVals.size()/double(binAmount));
	multiset<double>::reverse_iterator it;
	int i = nbins;
	int j = 0;
	bins[i] = 1;
	i--;
	double currVal = -1;
	cout << "binAmount = " << binAmount << endl;
	cout << "nbins = " << nbins << endl;
	for(it = setVals.rbegin(); it != setVals.rend(); ++it){
		if(j==0 && currVal!=-1){
			cout << i << " -> " << (*it+currVal)/2. << endl;
			bins[i] = (*it+currVal)/2.;
			--i;
		}
		if(++j==binAmount){
			currVal = *it;
			j = 0;
		}
	}
	bins[0]=0;

	cout << "#Data: " << setVals.size() << endl;
	cout << "nbins: " << nbins << endl;
	ofstream fd(binsFileName);
	for(int i=0; i<=nbins; ++i){
		cout << bins[i] << endl;
		fd << bins[i] << endl;
	}
	fd.close();
	tempFD->Close();
	remove(tempFileName.Data());
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
		make_bins(config);
}
