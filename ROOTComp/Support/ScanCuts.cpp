/*
 * ScanCuts.cpp
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#include "ScanCuts.h"
#include <iostream>
#include <fstream>
#include <sstream>

ClassImp(ScanCuts)
ClassImp(Cuts)

ClassImp(ScanCuts::k2pi_t)
ClassImp(ScanCuts::kmu3_t)

Cuts::Cuts(){
	//Initialize with default values
	triggerMask = 0x400;
	numVertex3 = 1;
	minZVertex = -1700;
	maxZVertex = 9000;
	maxChi2Vertex = 25;
	maxExtraTracks = 0;
	maxTrackTime = 25;
	boolBadTrack = true;
	numBadTrackCombi = 0;
	numXCandidates = 1;
	boolBadECandidates = true;
	minTrackMomentum = 2.0;
	maxTrackMomentum = 74.0;
	numAddGoodCluster = 1;
	lkrAcceptance = 0;
	minGammaEnergy = 3;
	minDeadCellDist = 2;
	minGammaDCHRadius = 13;
	unDeflectedElDist = 20;
	k2pi.minTotalMomentum = 70.0;
	k2pi.maxTotalMomentum = 78.0;
	kmu3.maxTotalMomentum = 78;
	k2pi.maxPt = 0.0005;
	kmu3.minPt = 0.0005;
	kmu3.maxPt = 0.04;
	k2pi.minPi0MassDiff = 0.115;
	k2pi.maxPi0MassDiff = 0.145;
	kmu3.maxPi0MassDiff = 0.01;
	k2pi.minKaonMassDiff = 0.465;
	k2pi.maxKaonMassDiff = 0.510;
	kmu3.maxMissMassSq = 0.01;
	k2pi.pi0Mass2DiffCoarse = 0.01;
}

void Cuts::print(){
	std::cout << "Cuts definition" << std::endl;
	std::cout << "---------------" << std::endl;
	std::cout << "triggerMask\t\t--> " << triggerMask << std::endl;
	std::cout << "numVertex3\t\t--> " << numVertex3 << std::endl;
	std::cout << "minZVertex\t\t--> " << minZVertex << std::endl;
	std::cout << "maxZVertex\t\t--> " << maxZVertex << std::endl;
	std::cout << "maxChi2Vertex\t\t--> " << maxChi2Vertex << std::endl;
	std::cout << "maxExtraTracks\t\t--> " << maxExtraTracks << std::endl;
	std::cout << "maxTrackTime\t\t--> " << maxTrackTime << std::endl;
	std::cout << "boolBadTrack\t\t--> " << boolBadTrack << std::endl;
	std::cout << "numBadTrackCombi\t--> " << numBadTrackCombi << std::endl;
	std::cout << "numXCandidates\t\t--> " << numXCandidates << std::endl;
	std::cout << "boolBadECandidates\t--> " << boolBadECandidates << std::endl;
	std::cout << "minTrackMomentum\t--> " << minTrackMomentum << std::endl;
	std::cout << "maxTrackMomentum\t--> " << maxTrackMomentum << std::endl;
	std::cout << "numAddGoodCluster\t--> " << numAddGoodCluster << std::endl;
	std::cout << "lkrAcceptance\t\t--> " << lkrAcceptance << std::endl;
	std::cout << "minGammaEnergy\t\t--> " << minGammaEnergy << std::endl;
	std::cout << "minDeadCellDist\t\t--> " << minDeadCellDist << std::endl;
	std::cout << "minGammaDCHRadius\t--> " << minGammaDCHRadius << std::endl;
	std::cout << "unDeflectedElDist\t--> " << unDeflectedElDist << std::endl;
	std::cout << "  -> K2Pi specific" << std::endl;
	std::cout << "minTotalMomentum\t--> " << k2pi.minTotalMomentum << std::endl;
	std::cout << "maxTotalMomentum\t--> " << k2pi.maxTotalMomentum << std::endl;
	std::cout << "maxPt\t\t\t--> " << k2pi.maxPt << std::endl;
	std::cout << "maxPi0MassDiff\t\t--> " << k2pi.maxPi0MassDiff << std::endl;
	std::cout << "minKaonMassDiff\t\t--> " << k2pi.minKaonMassDiff << std::endl;
	std::cout << "maxKaonMassDiff\t\t--> " << k2pi.maxKaonMassDiff << std::endl;
	std::cout << "pi0Mass2DiffCoarse\t\t--> " << k2pi.pi0Mass2DiffCoarse << std::endl;
	std::cout << "  -> KMu3 specific" << std::endl;
	std::cout << "maxTotalMomentum\t--> " << kmu3.maxTotalMomentum << std::endl;
	std::cout << "minPt\t\t\t--> " << kmu3.minPt << std::endl;
	std::cout << "maxPt\t\t\t--> " << kmu3.maxPt << std::endl;
	std::cout << "maxPi0MassDiff\t\t--> " << kmu3.maxPi0MassDiff << std::endl;
	std::cout << "maxMissMassSq\t\t--> " << kmu3.maxMissMassSq << std::endl;
}

ScanCuts::ScanCuts(): defaultIndex(-1), currentList(-1){

}

ScanCuts::~ScanCuts() {
}

bool ScanCuts::addParseCutsFile(std::string fileName) {
	if(cutsLists.size()==0){
		defaultIndex = 0;
		cutsLists.push_back(Cuts());
	}
	return parseCuts(fileName);
}

void ScanCuts::print() {
	if(currentList!=-1) print(currentList);
	else if(defaultIndex!=-1) print(defaultIndex);
}

void ScanCuts::print(int index) {
	cutsLists[index].print();
}

int ScanCuts::loadNextList() {
	if(++currentList<(int)cutsLists.size()){
		return loadList(currentList);
	}
	else{
		currentList = -1;
	}
	return currentList;
}

int ScanCuts::loadList(int i) {
	Cuts::operator=(cutsLists[i]);
	currentList = i;
	return i;
}

void ScanCuts::loadDefault() {
	if(defaultIndex!=-1) loadList(defaultIndex);
}

bool ScanCuts::parseCuts(std::string fileName) {
	if(fileName.length()==0) return true;

	unsigned int id;
	std::ifstream fdCuts;
	fdCuts.open(fileName.c_str(), ifstream::in);
	std::string name, value, index;
	std::string buffer;
	if(fdCuts.is_open()){
		while(getline(fdCuts, buffer)){
			if(buffer[0]=='#') continue;
			if(buffer.length()==0) continue;

			std::stringstream ss(buffer);

			//read new format
			ss >> name >> value;
			if(defaultIndex==-1 && name.compare("defaultIndex")==0){
				defaultIndex = atoi(value.c_str());
				continue;
			}

			if(!ss.eof()) ss >> id;
			else{
				//Read old format
				id = defaultIndex;
			}

			if(id>=cutsLists.size()){
				std::cerr << "Cuts list " << id << " does not exist. Size=" << cutsLists.size() << std::endl;
				return false;
			}


			if(name.compare("triggerMask")==0){
				cutsLists[id].triggerMask = atoi(value.c_str());
			}
			else if(name.compare("numVertex3")==0){
				cutsLists[id].numVertex3 = atoi(value.c_str());
			}
			else if(name.compare("minZVertex")==0){
				cutsLists[id].minZVertex = atoi(value.c_str());
			}
			else if(name.compare("maxZVertex")==0){
				cutsLists[id].maxZVertex = atoi(value.c_str());
			}
			else if(name.compare("maxChi2Vertex")==0){
				cutsLists[id].maxChi2Vertex = atoi(value.c_str());
			}
			else if(name.compare("maxExtraTracks")==0){
				cutsLists[id].maxExtraTracks = atoi(value.c_str());
			}
			else if(name.compare("maxTrackTime")==0){
				cutsLists[id].maxTrackTime = atoi(value.c_str());
			}
			else if(name.compare("boolBadTrack")==0){
				cutsLists[id].boolBadTrack = (value.compare("true")==0) ? true : false;
			}
			else if(name.compare("numBadTrackCombi")==0){
				cutsLists[id].numBadTrackCombi = atoi(value.c_str());
			}
			else if(name.compare("numXCandidates")==0){
				cutsLists[id].numXCandidates = atoi(value.c_str());
			}
			else if(name.compare("boolBadECandidates")==0){
				cutsLists[id].boolBadECandidates = (value.compare("true")==0) ? true : false;
			}
			else if(name.compare("minTrackMomentum")==0){
				cutsLists[id].minTrackMomentum = atof(value.c_str());
			}
			else if(name.compare("maxTrackMomentum")==0){
				cutsLists[id].maxTrackMomentum = atof(value.c_str());
			}
			else if(name.compare("numAddGoodCluster")==0){
				cutsLists[id].numAddGoodCluster = atoi(value.c_str());
			}
			else if(name.compare("lkrAcceptance")==0){
				cutsLists[id].lkrAcceptance = atoi(value.c_str());
			}
			else if(name.compare("minGammaEnergy")==0){
				cutsLists[id].minGammaEnergy = atoi(value.c_str());
			}
			else if(name.compare("minDeadCellDist")==0){
				cutsLists[id].minDeadCellDist = atoi(value.c_str());
			}
			else if(name.compare("minGammaDCHRadius")==0){
				cutsLists[id].minGammaDCHRadius = atoi(value.c_str());
			}
			else if(name.compare("unDeflectedElDist")==0){
				cutsLists[id].unDeflectedElDist = atoi(value.c_str());
			}
			else if(name.compare("k2pi.minTotalMomentum")==0){
				cutsLists[id].k2pi.minTotalMomentum = atof(value.c_str());
			}
			else if(name.compare("k2pi.maxTotalMomentum")==0){
				cutsLists[id].k2pi.maxTotalMomentum = atof(value.c_str());
			}
			else if(name.compare("k2pi.maxPt")==0){
				cutsLists[id].k2pi.maxPt = atof(value.c_str());
			}
			else if(name.compare("k2pi.maxPi0MassDiff")==0){
				cutsLists[id].k2pi.maxPi0MassDiff = atof(value.c_str());
			}
			else if(name.compare("k2pi.minKaonMassDiff")==0){
				cutsLists[id].k2pi.minKaonMassDiff = atof(value.c_str());
			}
			else if(name.compare("k2pi.maxKaonMassDiff")==0){
				cutsLists[id].k2pi.maxKaonMassDiff = atof(value.c_str());
			}
			else if(name.compare("k2pi.pi0Mass2DiffCoarse")==0){
				cutsLists[id].k2pi.pi0Mass2DiffCoarse = atof(value.c_str());
			}
			else if(name.compare("kmu3.maxTotalMomentum")==0){
				cutsLists[id].kmu3.maxTotalMomentum = atoi(value.c_str());
			}
			else if(name.compare("kmu3.minPt")==0){
				cutsLists[id].kmu3.minPt = atof(value.c_str());
			}
			else if(name.compare("kmu3.maxPt")==0){
				cutsLists[id].kmu3.maxPt = atof(value.c_str());
			}
			else if(name.compare("kmu3.maxPi0MassDiff")==0){
				cutsLists[id].kmu3.maxPi0MassDiff = atof(value.c_str());
			}
			else if(name.compare("kmu3.maxMissMassSq")==0){
				cutsLists[id].kmu3.maxMissMassSq = atof(value.c_str());
			}
			else{
				std::cout << "No known " << name << " cut parameter " << value << std::endl;
			}
		}
		fdCuts.close();
		return true;
	}
	else{
		std::cout << "Unable to open cuts file" << fileName << std::endl;
		return false;
	}
}

void ScanCuts::generateLists(int number) {
	cutsLists.resize(number);
}

double ScanCuts::getDiff() const {
}
