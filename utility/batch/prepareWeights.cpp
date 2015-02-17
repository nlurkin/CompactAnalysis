#include <map>
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
#include <fstream>
#include "../userinc/exportClasses.h"
using namespace std;

map<int, double> mcRunMap, dataRunMap;
void buildRatioMap();
void readConfig(TString confFile);
void buildRunMap(TFile *fd, map<int,double> &runMap);

void readConfig(TString confFile){
	vector<TString> mcFileNames;
	vector<TString> dataFileNames;

	//read list file
	ifstream listFile(confFile.Data());

	string line;
	while(getline(listFile,line)){
		TString tline(line);
		if(tline.BeginsWith("#")) continue;
		TString key(tline(0,tline.First('=')));
		TString values(tline(tline.First('=')+1, 100000));
		TObjArray* tok = values.Tokenize(" ");
		cout << tline << endl;
		for(int i=0; i<tok->GetEntries(); ++i){
			TString entry(((TObjString*)tok->At(i))->GetString());
			cout << entry << endl;
			if(key.CompareTo("mcfiles")==0) mcFileNames.push_back(entry);
			else if(key.CompareTo("datafiles")==0) dataFileNames.push_back(entry);
			delete tok->At(i);
		}
	}

	TFile *fd;
	for(unsigned int i=0; i<mcFileNames.size(); ++i){
		fd = TFile::Open(mcFileNames[i], "READ");
		buildRunMap(fd, mcRunMap);
		fd->Close();
	}
	for(unsigned int i=0; i<dataFileNames.size(); ++i){
		fd = TFile::Open(dataFileNames[i], "READ");
		buildRunMap(fd, dataRunMap);
		fd->Close();
	}
}

void buildRunMap(TFile *fd, map<int,double> &runMap){
	//Get the TTree
	//Input
	TTree *t = (TTree*)fd->Get("event");
	ROOTBurst *burstBrch = new ROOTBurst();
	t->SetBranchAddress("rawBurst", &burstBrch);

	int nevt = t->GetEntries();
	cout << nevt << endl;
	for(int i=0; i< nevt; ++i){
		t->GetEntry(i);
		if(runMap.count(burstBrch->nrun)>0) runMap[burstBrch->nrun] += 1;
		else runMap.insert(make_pair(burstBrch->nrun, 1.));
	}
}

void buildRatioMap(){
	map<int,double>::iterator it;

	map<int,double> ratio(mcRunMap);

	for(it=dataRunMap.begin(); it!=dataRunMap.end(); ++it){
		ratio[it->first] /=it->second;
	}

	ofstream fd;
	fd.open("pi0dalitz_weights.dat");
	for(it=ratio.begin(); it!=ratio.end(); ++it){
		cout << "Run " << it->first << ": " << it->second << "(" << mcRunMap[it->first] << "/" << dataRunMap[it->first] << endl;
		fd << it->first << " " << it->second << " " << mcRunMap[it->first] << " " << dataRunMap[it->first] << endl;
	}

	fd.close();

}

void prepareWeights(TString cfgFile){
	gSystem->Load("../obj/libmystructs.so");

	readConfig(cfgFile);

	buildRatioMap();
}

int main(int argc, char** argv){
	prepareWeights(argv[1]);
}
