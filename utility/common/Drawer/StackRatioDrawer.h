/*
 * DistribRatioDrawer.h
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_STACKRATIODRAWER_H_
#define COMMON_DRAWER_STACKRATIODRAWER_H_

#include "StackDrawer.h"

class StackRatioDrawer : public StackDrawer {
public:
	StackRatioDrawer();
	virtual ~StackRatioDrawer();

	void AddHisto1(TH1* h, std::string legend);
	void AddHisto2(TH1* h, std::string legend);
	void SetSecondary(TH1* h) {
		fSecondary = h;
	}

	virtual void generate(TPad* pad);
	virtual void draw();

private:
	THStack *fStack2;
	TH1 *fSecondary;
};

#endif /* COMMON_DRAWER_STACKRATIODRAWER_H_ */
