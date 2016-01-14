/*
 * CombineMCSample.h
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLES_COMBINESAMPLE_H_
#define COMMON_SAMPLES_COMBINESAMPLE_H_

#include "../Interface/SubSample.h"

#include <vector>
#include <TH2D.h>

#include "../userinc/exportClasses.h"

class CombineSample: public SubSample {
public:
	CombineSample();
	virtual ~CombineSample();

	virtual void processEvent(ROOTPhysicsEvent *eventBrch, ROOTBurst *burstBrch,
			ROOTRawEvent *rawBrch, ROOTCorrectedEvent *corrBrch,
			ROOTFileHeader *headerBrch, ROOTMCEvent *mcEvent, NGeom *geomBrch,
			std::vector<bool> *cutsPass, const ConfigFile *cfg, const RunWeights *weights);
	virtual void doGet(TDirectory* inputFD, TFile* tempFD) = 0;
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins, const ConfigFile *cfg);
	virtual double getFFIntegral(double a);

	virtual void setPlotStyle(std::vector<int> color);
	virtual TH1D* getMainHisto();

	void addHisto(TString name, int bins, double min, double max);
	void addHisto(TString name, int binsx, double minx, double maxx, int binsy,
			double miny, double maxy);

	void getHisto(TDirectory* fd, TFile* tempFD, TString name);
	void getHisto2(TDirectory* fd, TFile* tempFD, TString name);
	void doGetHisto(TDirectory* inputFD, TFile* tempFD);

	void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
				ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
				ROOTBurst *rootBurst, double weight);
	virtual void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
					ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
					ROOTBurst *rootBurst, const RunWeights *weights) = 0;
	void scale();
	virtual void renameHisto();
	virtual SubSample* Add(const SubSample* other);

	const std::vector<TH1D*>& getD1() const {
		return d1;
	}

	const std::vector<TH2D*>& getMap() const {
		return dMap;
	}

protected:
	std::vector<TH1D*> d1;
	std::vector<TH2D*> dMap;
};

#endif /* COMMON_SAMPLES_COMBINESAMPLE_H_ */
