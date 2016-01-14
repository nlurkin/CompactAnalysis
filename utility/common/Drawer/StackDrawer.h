/*
 * StackDrawer.h
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_STACKDRAWER_H_
#define COMMON_DRAWER_STACKDRAWER_H_

#include <TH1.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>

class StackDrawer {
public:
	StackDrawer();
	virtual ~StackDrawer();

	void AddHisto(TH1* h, std::string legend, std::string option="");

	void setName(const std::string& name) {
		fName = name;
	}

	void setTitle(const std::string& title) {
		fTitle = title;
	}

	virtual void generate(TPad* pad);
	virtual void draw();

	virtual void save();
	virtual void free();


protected:
	THStack *fStack;
	TLegend *fLegend;
	std::string fName;
	std::string fTitle;
	TCanvas *fCanvas;
};

#endif /* COMMON_DRAWER_STACKDRAWER_H_ */
