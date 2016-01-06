/*
 * CombineDrawer.h
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_COMBINEDRAWER_H_
#define COMMON_DRAWER_COMBINEDRAWER_H_

#include "../Interface/HistoDrawer.h"

class CombineDrawer: public HistoDrawer {
public:
	CombineDrawer();
	virtual ~CombineDrawer();

	virtual void draw();

	void addHistoMC(unsigned int index, TH1* histo);
	void addHistoData(unsigned int index, TH1* histo);
	void addLegendMC(TH1* histo, std::string legend);
	void addLegendData(TH1* histo, std::string legend);
private:
	std::vector<THStack*> fStack, fDataStack;
	std::vector<TH1D*> fSum, fDataSum;
	TLegend *fLegend;
};

#endif /* COMMON_DRAWER_COMBINEDRAWER_H_ */
