/*
 * SubSample.h
 *
 *  Created on: Jan 11, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_INTERFACE_SUBSAMPLE_H_
#define COMMON_INTERFACE_SUBSAMPLE_H_

#include "../userinc/exportClasses.h"
#include <TH1D.h>

#include "../ConfigFile.h"
#include "RunWeights.h"

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


class SubSample {
public:
	static const int NBINS = 10000000;
	static const int MAXBIN = 1;

	SubSample();
	virtual ~SubSample();

	void initNewFile(int totalChanEvents, int selEvents);

	fitStruct getFitStruct() { return fFitBrch; };
	void scale(TH1 *histo, double scaleFactor);
	virtual SubSample* Add(const SubSample* other);
	void initOutput();
	void writeTree();

	virtual void processEvent(ROOTPhysicsEvent *eventBrch, ROOTBurst *burstBrch,
			ROOTRawEvent *rawBrch, ROOTCorrectedEvent *corrBrch,
			ROOTFileHeader *headerBrch, ROOTMCEvent *mcEvent, NGeom *geomBrch,
			std::vector<bool> *cutsPass, const ConfigFile *cfg, const RunWeights *weights) = 0;
	virtual void doGet(TDirectory* inputFD, TFile* tempFD) = 0;
	virtual void doWrite() = 0;
	virtual void doSetName() = 0;
	virtual void initHisto(int nbins, double* bins, const ConfigFile *cfg) = 0;
	virtual double getFFIntegral(double a) = 0;
	virtual void renameHisto() = 0;

	virtual void setPlotStyle(std::vector<int> color) = 0;
	virtual TH1D* getMainHisto() = 0;

	void setBr(double br) {
		fBr = br;
	}

	void setIndex(int index) {
		fIndex = index;
	}

	void setScanId(int scanId) {
		fScanID = scanId;
	}

	int getSelSize() const {
		return fFitBrch.selEvents;
	}

	void setSelSize(int selSize) {
		fFitBrch.selEvents = selSize;
	}

	int getTotalSize() const {
		return fFitBrch.selEvents;
	}

	void setTotalSize(int totalSize) {
		fFitBrch.totEvents = totalSize;
	}


protected:
	fitStruct fFitBrch;
	TTree* fFitTree;
	int fScanID;
	int fIndex;
	double fBr;
};

#endif /* COMMON_INTERFACE_SUBSAMPLE_H_ */
