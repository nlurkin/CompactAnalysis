/*
 * CombineMCSample.h
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLES_COMBINEMCSAMPLE_H_
#define COMMON_SAMPLES_COMBINEMCSAMPLE_H_

#include "CombineSample.h"
#include "../Interface/MCSample.h"

class CombineMCSample: public CombineSample, public MCSample {
public:
	typedef struct bContent_t {
		double val;

		struct bContent_t &operator +=(struct bContent_t other) {
			val += other.val;
			return *this;
		}
	} bContent;

	CombineMCSample();
	CombineMCSample(int index, ConfigFile *cfg);
	virtual ~CombineMCSample();

	virtual void doGet(TFile* inputFD, TFile* tempFD);
	virtual void scaleToData(bContent totalMC, double nData);

	void fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
			ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent,
			NGeom *rootGeom, ROOTBurst *rootBurst);

	bContent getIntegrals();
};

#endif /* COMMON_SAMPLES_COMBINEMCSAMPLE_H_ */
