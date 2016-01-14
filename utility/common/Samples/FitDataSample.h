/*
 * FitDataSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITDATASAMPLE_H_
#define COMMON_FITDATASAMPLE_H_

#include "../Interface/DataSample.h"
#include "../Interface/SubSample.h"

class FitDataSample: public SubSample, public DataSample{
public:
	typedef struct bContent_t{
		double dSig;
	} bContent;

	FitDataSample();
	virtual ~FitDataSample();

	virtual void processEvent(ROOTPhysicsEvent *eventBrch, ROOTBurst *burstBrch,
			ROOTRawEvent *rawBrch, ROOTCorrectedEvent *corrBrch,
			ROOTFileHeader *headerBrch, ROOTMCEvent *mcEvent, NGeom *geomBrch,
			std::vector<bool> *cutsPass, const ConfigFile *cfg, const RunWeights *weights);
	virtual void doGet(TDirectory* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins, const ConfigFile *cfg);
	virtual void scaleToData(double) {};
	virtual void setPlotStyle(std::vector<int> color);
	virtual void populateStack(HistoDrawer *drawer, std::string legend);
	virtual void populateFit(HistoDrawer *drawer, double norm, double a, std::string legend);
	virtual double getFFIntegral(double a);
	virtual void renameHisto() {};

	virtual TH1D* getMainHisto() { return dSig; }

	void scale();

	bContent getBinContent(int bin);

	virtual SubSample* Add(const SubSample* other);

	const TH1D* getSig() const {
		return dSig;
	}

private:
	TH1D *dSig;
};

#endif /* COMMON_FITDATASAMPLE_H_ */
