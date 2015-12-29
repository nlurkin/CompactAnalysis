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
#include <TTree.h>
#include "ConfigFile.h"
#include "RunWeights.h"
#include "../Drawer/InputFitDrawer.h"
#include "../Drawer/FitResultDrawer.h"

class TTree;
class TFile;

typedef struct fitStruct_t {
	int totEvents;
	int selEvents;
	int n1;
	int nx;
	int nxx;

	fitStruct_t& operator+=(const fitStruct_t &other){
		totEvents += other.totEvents;
		selEvents += other.selEvents;
		n1 += other.n1;
		nx += other.nx;
		nxx += other.nxx;
		return *this;
	}
} fitStruct;

void initFitStruct(fitStruct &s);
void sumTreeFitStruct(fitStruct &in, TTree *t, fitStruct &out, double factor);

class Sample {
public:
	static const int NBINS = 10000000;
	static const int MAXBIN = 1;

	Sample();
	Sample(int index, ConfigFile *cfg);
	virtual ~Sample();

	bool addFile(std::string);
	void fill(TFile* tempFD, int nbins, double* bins);
	void get(TFile* tempFD);
	void initOutput();
	void closeOutput(TFile* tempFD);

	void scale(TH1 *histo, double scaleFactor);

	virtual void doFill(TFile* inputFD, TFile* tempFD) = 0;
	virtual void doGet(TFile* inputFD, TFile* tempFD) = 0;
	virtual void doWrite() = 0;
	virtual void doSetName() = 0;
	virtual void initHisto(int nbins, double* bins) = 0;
	virtual void scaleToData(double nData) = 0;

	virtual void setPlotStyle(std::vector<int> color) = 0;
	virtual void populateStack(InputFitDrawer &drawer) = 0;
	virtual void populateFit(FitResultDrawer &drawer, double norm, double a) = 0;
	virtual TH1D* getMainHisto() = 0;

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

	friend Sample& operator+=(Sample &first, const Sample* other);

	void setCfg(const ConfigFile* cfg) {
		fCfg = cfg;
	}

protected:
	int fIndex;
	double fBr;
	std::vector<std::string> fListFiles;
	std::string fOutputFile;
	TTree* fFitTree;
	fitStruct fFitBrch;
	TFile *fOutputFD;
	const ConfigFile *fCfg;
	const RunWeights *fWeights;
	std::string fLegend;
};
#endif /* COMMON_SAMPLE_H_ */
