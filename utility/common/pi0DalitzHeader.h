#include <sys/stat.h>
#include <vector>
#include <TString.h>
#include <iostream>
#include <TObjArray.h>
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <fstream>
#include <TVector3.h>
#include <TH1.h>
#include <signal.h>
#include <stdlib.h>
#include <TApplication.h>
#include <iomanip>
#include "exportClasses.h"
#include "ScanCuts.h"
#include "../userinc/funLib.h"
using namespace std;

ROOTRawEvent *xxx = new ROOTRawEvent();
ROOTRawEvent &rawEvent = *xxx;

void setStyle();
/*************************
 * Structures
 *************************/
typedef struct fitStruct_t {
	unsigned int totEvents;
	unsigned int selEvents;
	int n1;
	int nx;
	int nxx;
} fitStruct;

void initFitStruct(fitStruct &s){
	s.n1 = 0;
	s.nx = 0;
	s.nxx = 0;
	s.selEvents = 0;
	s.totEvents = 0;
}

void sumTreeFitStruct(fitStruct &in, TTree *t, fitStruct &out){
	for(int i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		out.n1 += in.n1;
		out.nx += in.nx;
		out.nxx += in.nxx;
		out.selEvents += in.selEvents;
		out.totEvents += in.totEvents;
	}

	cout << "Total events: \t" << out.totEvents << endl;
	cout << "Sel.  events: \t" << out.selEvents << endl;
	cout << "n1    events: \t" << out.n1 << endl;
	cout << "nx    events: \t" << out.nx << endl;
	cout << "nxx   events: \t" << out.nxx << endl;
}

/*************************
 * Globals
 *************************/
int inputMCNbr = 0;
int inputDataNbr = 0;

int runStart = 0;
int runEnd = 0;

int periodStart = -1;
int periodEnd = -1;

vector<int> mcColors, dataColors;
vector<TString> mcLegendTitle, dataLegendTitle;

vector<double> brs;
vector<TString> mcFileNames;
vector<TString> dataFileNames;
vector<int> mcIndexes;
vector<TString> mcOutputFiles;
vector<TString> dataOutputFiles;
vector<TString> modelFiles;
TString binsFileName;
bool withEqualBins;

TTree* fitTree;
fitStruct fitBrch;

TString tempFileName;
TFile *tempFD;
TApplication *theApp = 0;

map<int,double> ratioMap;
double averageRatio;

double testA = 0.06;

int scanID = -1;

/*************************
 * Signal handling
 *************************/
void sighandler(int sig)
{
	cout << "Closing comme une brute" << endl;
	if(tempFD){
		tempFD->Close();
		remove(tempFileName.Data());
	}

	if(theApp) theApp->Terminate();
	exit(0);
}

/*************************
 * To be defined by user
 *************************/
void closeMCOutput(TFile *fdo, int index);
void closeDataOutput(TFile *fdo, int index);
void initNewOutput(TFile **fdo, TString fileName);

namespace Input{
	void getInputMCFill(TFile *fd, TFile *fdo, double br, unsigned int index);
	int getInputDataFill(TFile *fd, TFile *fdo);
	void getInputMCGet(TFile *fd, double br, unsigned int index);
	int getInputDataGet(TFile *fd);
}

void initNewChannel();
void scaleMC(fitStruct N, int index, double br);

/************************
 * Utility
 ************************/
void scale(TH1 *histo, double scaleFactor, double totEvents, double br){
	//histo->Scale(br/((double)totEvents/(double)selEvents)/integral);
	histo->Scale(br/(totEvents*scaleFactor));
}

bool testAllOutputs(){
	vector<TString>::iterator it;
	struct stat fStat;
	bool error = true;

	for(it=mcOutputFiles.begin(); it!= mcOutputFiles.end(); it++){
		if(stat((*it).Data(), &fStat) == 0){
			cout << "[ERROR] Output file " << *it << " already exists" << endl;
			error = false;
		}
	}

	for(it=dataOutputFiles.begin(); it!= dataOutputFiles.end(); it++){
		if(stat((*it).Data(), &fStat) == 0){
			cout << "[ERROR] Output file " << *it << " already exists" << endl;
			error = false;
		}
	}
	return error;
}

bool runIncluded(int run, int period){
	bool ret = true;
	if(runStart!=0 && run<runStart) ret = false;
	if(runEnd!=0 && run>runEnd) ret = false;
	if(periodStart!=-1 && period<periodStart) ret = false;
	if(periodEnd!=-1 && period>periodEnd) ret = false;
	return ret;
}

/************************
 * Getting input
 ************************/
bool readConfig(TString confFile){
	//read list file
	ifstream listFile(confFile.Data());

	string line;
	cout << "******** Configuration file **********" << endl;
	cout << "\t\t" << confFile << endl;
	while(getline(listFile,line)){
		TString tline(line);
		if(tline.BeginsWith("#")) continue;
		TString key(tline(0,tline.First('=')));
		TString values(tline(tline.First('=')+1, 10000));
		TObjArray* tok = values.Tokenize(" ");
		cout << tline << endl;
		for(int i=0; i<tok->GetEntries(); ++i){
			TString entry(((TObjString*)tok->At(i))->GetString());
			if(key.CompareTo("mcfiles")==0){
				mcFileNames.push_back(entry);
				/*if(entry.Contains(".root")) mcFileNames.push_back(entry);
				else{
					ifstream fd(entry);
					TString fileName;
					while(fd >> fileName){
						mcFileNames.push_back(fileName);
					}
					fd.close();
				}*/
			}
			else if(key.CompareTo("mcout")==0) mcOutputFiles.push_back(entry);
			else if(key.CompareTo("modelfiles")==0) modelFiles.push_back(entry);
			else if(key.CompareTo("brs")==0) brs.push_back(entry.Atof());
			else if(key.CompareTo("mccolors")==0) mcColors.push_back(entry.Atoi());
			else if(key.CompareTo("mclegends")==0) mcLegendTitle.push_back(entry);
			else if(key.CompareTo("datafiles")==0){
				dataFileNames.push_back(entry);
				/*if(entry.Contains(".root")) dataFileNames.push_back(entry);
				else{
					ifstream fd(entry);
					TString fileName;
					while(fd >> fileName){
						dataFileNames.push_back(fileName);
					}
					fd.close();
				}*/
			}
			else if(key.CompareTo("dataout")==0) dataOutputFiles.push_back(entry);
			else if(key.CompareTo("datacolors")==0) dataColors.push_back(entry.Atoi());
			else if(key.CompareTo("datalegends")==0) dataLegendTitle.push_back(entry);
			else if(key.CompareTo("mcIndex")==0) mcIndexes.push_back(entry.Atoi());
			else if(key.CompareTo("runstart")==0) runStart = entry.Atoi();
			else if(key.CompareTo("runend")==0) runEnd = entry.Atoi();
			else if(key.CompareTo("periodstart")==0) 	periodStart = entry.Atoi();
			else if(key.CompareTo("periodend")==0) 		periodEnd = entry.Atoi();
			else if(key.CompareTo("testA")==0) testA = entry.Atof();
			else if(key.CompareTo("binsfile")==0) binsFileName = entry;
			else if(key.CompareTo("equalbin")==0) withEqualBins = entry.CompareTo("true")==0 ? true : false;
			else if(key.CompareTo("scanid")==0) scanID = entry.Atoi();
			delete tok->At(i);
		}
	}
	cout << "*************************************" << endl;

	return testAllOutputs();
}

//Read data/MC and fill histograms. Read from multiple file and save-merge in multiple files
void readFilesFill(){
	vector<TFile*> fd;
	TFile* ffd;
	int prevIndex = -1;
	int newIndex;

	TFile *fdo;

	unsigned int fileNumber;

	fileNumber = mcFileNames.size();

	//Getting MC
	for(unsigned int i=0; i<fileNumber; ++i){
		//Do we have a new output file?
		if(i>=mcIndexes.size()) newIndex = prevIndex;
		else newIndex = mcIndexes[i];
		if(prevIndex != newIndex){
			//Need to close the previous one
			if(prevIndex!=-1){
				//scaleMC(fitBrch, prevIndex, brs[prevIndex]);
				closeMCOutput(fdo, prevIndex);
			}
			prevIndex = newIndex;
			++inputMCNbr;
			initNewOutput(&fdo, mcOutputFiles[prevIndex]);
		}

		//Open new input file
		cout << mcFileNames[i] << endl;
		vector<TString> localFiles;
		if(!mcFileNames[i].Contains(".root")){
			cout << "List file detected..." << endl;
			ifstream listFile(mcFileNames[i]);
			string buffer;
			while(getline(listFile, buffer)){
				localFiles.push_back(buffer);
			}
		}
		else{
			localFiles.push_back(mcFileNames[i]);
		}

		int iFile=0;
		double totalFiles = localFiles.size();
		for(auto files : localFiles){
			cout << "Processing file " << files << " " << setprecision(2) << std::fixed << iFile*100./totalFiles << "% " << iFile << "/" << totalFiles << endl;
			cout << "Br: " << brs[prevIndex] << endl;
			ffd = TFile::Open(files);

			//Request the TTree reading function
			Input::getInputMCFill(ffd, fdo, brs[prevIndex], prevIndex);

			//Close the input file
			ffd->Close();
			iFile++;
		}
	}

	//Close last output file
	closeMCOutput(fdo, prevIndex);

	//Getting data
	fileNumber = dataFileNames.size();

	cout << dataOutputFiles.size() << endl;
	if(fileNumber>0) initNewOutput(&fdo, dataOutputFiles[0]);

	for(unsigned int i=0; i<fileNumber; ++i){
		++inputDataNbr;
		//Open new input file
		cout << dataFileNames[i] << endl;

		vector<TString> localFiles;
		if(!dataFileNames[i].Contains(".root")){
			cout << "List file detected..." << endl;
			ifstream listFile(dataFileNames[i]);
			string buffer;
			while(getline(listFile, buffer)){
				localFiles.push_back(buffer);
			}
		}
		else{
			localFiles.push_back(dataFileNames[i]);
		}

		int iFile=0;
		double totalFiles = localFiles.size();
		for(auto files : localFiles){
			cout << "Processing file " << files << setprecision(2) << iFile*100./totalFiles << "% " << iFile << "/" << totalFiles << endl;

			ffd = TFile::Open(files);

			//Request the TTree reading function
			Input::getInputDataFill(ffd, fdo);

			//Close input file
			ffd->Close();
			iFile++;
		}
	}

	//Close last output file
	if(fileNumber>0) closeDataOutput(fdo, 0);
}

//Read data/MC and get histograms
void readFilesGet(){
	vector<TFile*> fd;
	TFile* ffd;
	int prevIndex = -1;

	int fileNumber;

	fileNumber = mcOutputFiles.size();

	//Getting MC
	for(int i=0; i<fileNumber; ++i){
		//Do we have a new channel?
		if(prevIndex != mcIndexes[i]){
			///if(prevIndex != -1){
				//scaleMC(fitBrch, prevIndex, brs[prevIndex]);
			//}

			prevIndex = mcIndexes[i];
			++inputMCNbr;
			initNewChannel();
		}

		//Open new input file
		cout << mcOutputFiles[i] << endl;
		ffd = TFile::Open(mcOutputFiles[i]);

		//Request the Histo reading function
		Input::getInputMCGet(ffd, brs[prevIndex], prevIndex);
		//cout << fitBrch.selEvents << " " << fitBrch.totEvents << endl;

		//Close the input file
		ffd->Close();
	}

	//cout << fitBrch.selEvents << " " << fitBrch.totEvents << endl;
	//scaleMC(fitBrch, prevIndex, brs[prevIndex]);
	//Close the input file
	//ffd->Close();

	//Getting data
	initNewChannel();

	fileNumber = dataOutputFiles.size();

	for(int i=0; i<fileNumber; ++i){
		++inputDataNbr;
		//Open new input file
		cout << dataOutputFiles[i] << endl;
		ffd = TFile::Open(dataOutputFiles[i]);

		//Request the Histo reading function
		Input::getInputDataGet(ffd);

		//Close input file
		ffd->Close();
	}
}

void loadWeights(TString fileName){
	cout << "Load Weights" << endl;
	ifstream fd;
	struct stat fStat;

	if(stat(fileName.Data(), &fStat) != 0){
		cout << "[ERROR] Weight file " << fileName << " does not exist" << endl;
		return;
	}
	fd.open(fileName.Data(), ios_base::in);
	int burst =-1;//, mcNumber=-1, dataNumber=-1;
	double ratio, val;
	int i;
	i = 0;
	double sumRatio = 0;
	int number = 0;

	while(fd >> val){
		if(i==0) burst = val;
		else if(i==1) ratio = val;
		//else if(i==2) mcNumber = val;
		else if(i==3){
			//dataNumber = val;
			sumRatio += ratio;
			number++;
			ratioMap.insert(std::pair<int, double>(burst, ratio));
			//cout << burst << " " << ratio << " " << mcNumber << " " << dataNumber << endl;
			i = -1;
		}
		i++;
	}

	averageRatio = sumRatio / (double)number;
	cout << "Average " << sumRatio << " " << averageRatio << endl;
	fd.close();
}

double applyWeights(int run){
	return 1.;
	map<int, double>::iterator it;

	if((it = ratioMap.find(run)) == ratioMap.end()){
		cout << "[ERROR] No weight found for run " << run << endl;
		return -1;
	}
	double ratio = ratioMap[run];
	return averageRatio / ratio;
}

void loadBins(double *bins, int& nbins){
	ifstream fd(binsFileName);

	double val;
	nbins = 0;
	while(fd >> val){
		bins[nbins] = val;
		nbins++;
	}
}

bool testAdditionalCondition(ROOTPhysicsEvent *evt, ROOTCorrectedEvent *corrEvent, NGeom *rootGeom){
	TVector3 propPos, propPos2, propPos3;

	/*if(rootBurst->period!=1){
		fitBrch.selEvents--;
		return;
	}*/
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack]);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack]);
	if(distance2D(propPos, TVector3(0,0,0))<20 && distance2D(propPos2, TVector3(0,0,0))<20 ){
		fitBrch.selEvents--;
		return false;
	}

	if(evt->x <= 0.01) {
		fitBrch.selEvents--;
		return false;
	}

	return true;
}
