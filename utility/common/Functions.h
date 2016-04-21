/*
 * Functions.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FUNCTIONS_H_
#define COMMON_FUNCTIONS_H_

#include "exportClasses.h"

class TTree;

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


bool testAdditionalCondition(ROOTPhysicsEvent *evt, ROOTCorrectedEvent *corrEvent, NGeom *rootGeom, ROOTRawEvent *rawEvent, fitStruct &fitBrch);
double secondOrder_Pk(int nrun, double pkgen, int mcPartType);


#endif /* COMMON_FUNCTIONS_H_ */
