/*
 * CombineDataSample.h
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLES_COMBINEDATASAMPLE_H_
#define COMMON_SAMPLES_COMBINEDATASAMPLE_H_

#include "CombineSample.h"
#include "../Interface/DataSample.h"

class CombineDataSample: public CombineSample, public DataSample {
public:
	CombineDataSample();
	CombineDataSample(int index, ConfigFile *cfg);
	virtual ~CombineDataSample();

	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void populateStack(HistoDrawer *drawer);

	void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
				ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
				ROOTBurst *rootBurst);

	double getFactor() const {
		return fFactor;
	}

	void setFactor(double factor) {
		fFactor = factor;
	}

	void setTestA(double) {
	}

private:
	double fFactor;
};

#endif /* COMMON_SAMPLES_COMBINEDATASAMPLE_H_ */
