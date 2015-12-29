/*
 * HistoDrawer.h
 *
 *  Created on: Dec 24, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_HISTODRAWER_H_
#define COMMON_HISTODRAWER_H_

#include <THStack.h>
#include <TLegend.h>

class HistoDrawer {
public:
	HistoDrawer();
	virtual ~HistoDrawer();

	virtual void draw() = 0;
	static TH1D* buildRatio(int nbins, double* binning, TH1D* mc, TH1D* data);

	void setTitle(const std::string& title) {
		fTitle = title;
	}

	const std::string& getTitle() const {
		return fTitle;
	}

public:
	std::string fTitle;
};

#endif /* COMMON_HISTODRAWER_H_ */
