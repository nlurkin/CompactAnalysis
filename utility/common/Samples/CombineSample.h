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
	virtual void doGet(TFile* inputFD, TFile* tempFD) = 0;
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins, const ConfigFile *cfg);
	virtual double getFFIntegral(double a);

	virtual void setPlotStyle(std::vector<int> color);
	virtual void populateStack(HistoDrawer *drawer, std::string legend);
	virtual void populateFit(HistoDrawer *drawer, double norm, double a, std::string legend);
	virtual TH1D* getMainHisto();

	void addHisto(TString name, int bins, double min, double max);
	void addHisto(TString name, int binsx, double minx, double maxx, int binsy,
			double miny, double maxy);

	void getHisto(TFile* fd, TFile* tempFD, TString name);
	void getHisto2(TFile* fd, TFile* tempFD, TString name);
	void doGetHisto(TFile* inputFD, TFile* tempFD);

	void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
				ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
				ROOTBurst *rootBurst, double weight);
	virtual void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
					ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
					ROOTBurst *rootBurst, const RunWeights *weights) = 0;
	void scale();
	virtual void renameHisto();
	CombineSample* Add(const CombineSample *other);

protected:
	std::vector<TH1D*> d1;
	std::vector<TH2D*> dMap;
};

#endif /* COMMON_SAMPLES_COMBINESAMPLE_H_ */
