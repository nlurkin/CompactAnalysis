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
	virtual ~CombineDataSample();

	virtual void doGet(TDirectory* inputFD, TFile* tempFD);
	virtual void setPlotStyle(std::vector<int> color);

	void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
				ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
				ROOTBurst *rootBurst, const RunWeights *weigths);

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
