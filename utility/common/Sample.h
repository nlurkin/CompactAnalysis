/*
 * Sample.h
 *
 *  Created on: Dec 21, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLE_H_
#define COMMON_SAMPLE_H_

#include <vector>
#include <string>
#include <TH1D.h>
#include "ConfigFile.h"
#include "RunWeights.h"

class TTree;
class TFile;

typedef struct fitStruct_t {
	unsigned int totEvents;
	unsigned int selEvents;
	int n1;
	int nx;
	int nxx;
} fitStruct;

class Sample {
public:
	static const int NBINS = 10000000;
	static const int MAXBIN = 1;

	Sample(int index, ConfigFile &cfg);
	virtual ~Sample();

	bool addFile(std::string);
	void fill(TFile* tempFD, int nbins, double* bins);
	void initOutput();
	void closeOutput(TFile* tempFD);

	virtual void doFill(TFile* inputFD, TFile* tempFD) = 0;
	virtual void doWrite() = 0;
	virtual void doSetName() = 0;
	virtual void initHisto(int nbins, double* bins) = 0;

	double getBr() const {
		return fBr;
	}

	void setBr(double br) {
		fBr = br;
	}

	int getSelSize() const {
		return fFitBrch.selEvents;
	}

	void setSelSize(int selSize) {
		fFitBrch.selEvents = selSize;
	}

	int getTotalSize() const {
		return fFitBrch.totEvents;
	}

	void setTotalSize(int totalSize) {
		fFitBrch.totEvents = totalSize;
	}

	void setOutputFile(const std::string& outputFile) {
		fOutputFile = outputFile;
	}

	const RunWeights* getWeights() const {
		return fWeights;
	}

	void setWeights(const RunWeights* weights) {
		fWeights = weights;
	}

protected:
	int fIndex;
	double fBr;
	std::vector<std::string> fListFiles;
	std::string fOutputFile;
	TTree* fFitTree;
	fitStruct fFitBrch;
	TFile *fOutputFD;
	ConfigFile &fCfg;
	const RunWeights *fWeights;
};

#endif /* COMMON_SAMPLE_H_ */
