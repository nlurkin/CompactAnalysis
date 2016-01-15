/*
 * RunWeights.cpp
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#include "RunWeights.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace std;

RunWeights::RunWeights() {
	// TODO Auto-generated constructor stub

}

RunWeights::~RunWeights() {
	// TODO Auto-generated destructor stub
}

bool RunWeights::loadWeights(string fileName){
	cout << "Load Weights" << endl;
	ifstream fd;
	struct stat fStat;

	if(stat(fileName.c_str(), &fStat) != 0){
		cout << "[ERROR] Weight file " << fileName << " does not exist" << endl;
		return false;
	}
	fd.open(fileName.c_str(), ios_base::in);
	int burst =-1;
	double ratio, val;
	int i;
	i = 0;
	double sumRatio = 0;
	int number = 0;

	while(fd >> val){
		if(i==0) burst = val;
		else if(i==1) ratio = val;
		else if(i==3){
			sumRatio += ratio;
			number++;
			fRatioMap.insert(std::pair<int, double>(burst, ratio));
			i = -1;
		}
		i++;
	}

	fAverageRatio = sumRatio / (double)number;
	cout << "Average " << sumRatio << " " << fAverageRatio << endl;
	fd.close();
	return true;
}

double RunWeights::applyWeights(int run) const{
	map<int, double>::const_iterator it;

	if((it = fRatioMap.find(run)) == fRatioMap.end()){
		cout << "[ERROR] No weight found for run " << run << endl;
		return -1;
	}
	double ratio = (*it).second;
	return fAverageRatio / ratio;
}
