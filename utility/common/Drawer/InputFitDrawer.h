/*
 * InputFitDrawer.h
 *
 *  Created on: Dec 29, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_INPUTFITDRAWER_H_
#define COMMON_DRAWER_INPUTFITDRAWER_H_

#include "../Interface/HistoDrawer.h"

class InputFitDrawer: public HistoDrawer {
public:
	InputFitDrawer();
	virtual ~InputFitDrawer();

	virtual void draw();

public:
	THStack *fStd1;
	THStack *fStdx;
	THStack *fStdxx;
	THStack *fStdNew;
	THStack *fStdAlpha;
	THStack *fStdBeta;
	THStack *fStdGamma;
	THStack *fSig;
	TLegend *fLeg1;
	TLegend *fLegx;
	TLegend *fLegxx;
	TLegend *fLegNew;
	TLegend *fLegAlpha;
	TLegend *fLegBeta;
	TLegend *fLegGamma;
	TLegend *fLegSig;
};

#endif /* COMMON_DRAWER_INPUTFITDRAWER_H_ */
