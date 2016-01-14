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
#include "SubSample.h"
#include "DataSample.h"
#include "MCSample.h"
#include <TFile.h>

class TTree;

void inline mkdirCd(TDirectory *fd, TString name){
	if(!fd->GetKey(name)) fd->mkdir(name);
	fd->cd(name);
}

class Sample {
public:
	Sample();
	Sample(int index, ConfigFile *cfg);
	virtual ~Sample();

	bool addFile(std::string);
	void fill(TFile* tempFD, int nbins, double* bins);
	void get(TFile* tempFD);
	void initOutput();
	void closeOutput(TFile* tempFD);

	template<class SSampleType>
	void prepareNSubSamples(int N);

	void initHisto(int nbins, double* bins);
	void renameHisto() {
		for (auto ss : fSubSamples)
			ss->renameHisto();
	}
	;
	void setPlotStyle(std::vector<int> color) {
		for (auto ss : fSubSamples)
			ss->setPlotStyle(color);
	}
	;
	void populateStack(HistoDrawer *drawer) {
		for (auto ss : fSubSamples)
			ss->populateStack(drawer, fLegend);
	}
	;

	template<class SSampleType>
	std::vector<typename SSampleType::bContent> getIntegrals();
	template<class SSampleType>
	void scaleToData(const std::vector<typename SSampleType::bContent> totalMC,
			const std::vector<double> nData);

	virtual void doFill(TFile* inputFD, TFile* tempFD);
	virtual void doGet(TFile* inputFD, TFile* tempFD);

	double getBr() const {
		return fBr;
	}

	void setBr(double br) {
		fBr = br;
	}

	int getSelSize() const {
		return fSubSamples[fMainSubSample]->getSelSize();
	}

//	void setSelSize(int selSize) {
//		fFitBrch.selEvents = selSize;
//	}

	std::vector<double> getTotalSize() const {
		std::vector<double> r;
		for (auto ss : fSubSamples)
			r.push_back(ss->getTotalSize());
		return r;
	}

//	void setTotalSize(int totalSize) {
//		fFitBrch.totEvents = totalSize;
//	}

	void setOutputFile(const std::string& outputFile) {
		fOutputFile = outputFile;
	}

	const RunWeights* getWeights() const {
		return fWeights;
	}

	void setWeights(const RunWeights* weights) {
		fWeights = weights;
	}

	Sample* Add(const Sample* other);

	void setCfg(const ConfigFile* cfg) {
		fCfg = cfg;
	}

	void setLegend(const std::string& legend) {
		fLegend = legend;
	}

	void setTestA(double testA) {
		for (auto ss : fSubSamples)
			dynamic_cast<DataSample*>(ss)->setTestA(testA);
	}

	void setFactor(double factor) {
		for (auto ss : fSubSamples)
			dynamic_cast<DataSample*>(ss)->setFactor(factor);
	}

	SubSample * getSubSample(int i) {
		return fSubSamples[i];
	}
	int getNSubSample() {
		return fSubSamples.size();
	}

	const std::string& getLegend() const {
		return fLegend;
	}

	int getIndex() const {
		return fIndex;
	}

protected:
	int fIndex;
	double fBr;
	std::vector<std::string> fListFiles;
	std::string fOutputFile;
	//TTree* fFitTree;
	//fitStruct fFitBrchB;
	TFile *fOutputFD;
	const ConfigFile *fCfg;
	const RunWeights *fWeights;
	std::string fLegend;
	std::vector<SubSample*> fSubSamples;
	int fMainSubSample;
};

template<class SSampleType>
void Sample::prepareNSubSamples(int N) {
	SubSample *newSS;
	for (int i = 0; i < N; i++) {
		newSS = new SSampleType();
		newSS->setBr(fBr);
		newSS->setIndex(fIndex);
		newSS->setScanId(i);
		fSubSamples.push_back(newSS);
	}
}

template<class SSampleType>
std::vector<typename SSampleType::bContent> Sample::getIntegrals() {
	std::vector<typename SSampleType::bContent> r;
	for (auto ss : fSubSamples)
		r.push_back(static_cast<SSampleType*>(ss)->getIntegrals());
	return r;
}

template<class SSampleType>
void Sample::scaleToData(
		const std::vector<typename SSampleType::bContent> totalMC,
		const std::vector<double> nData) {
	for (unsigned int i = 0; i < fSubSamples.size(); i++) {
		static_cast<SSampleType*>(fSubSamples[i])->scaleToData(totalMC[i],
				nData[i]);
	}
}

#endif /* COMMON_SAMPLE_H_ */
