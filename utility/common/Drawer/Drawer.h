/*
 * Drawer.h
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_DRAWER_H_
#define COMMON_DRAWER_DRAWER_H_

#include <TH1D.h>
#include "../Interface/Sample.h"
#include "../Fitter/MinuitFitter.h"

class Drawer {
public:
	Drawer();
	virtual ~Drawer();

	static TH1D* buildRatio(TH1D* mc, TH1D* data);

	static void drawFitResult(MinuitFitter * fit, std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, Sample* finalMCSample, Sample* finalDataSample);
	static void drawFitPreparation(std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, std::string title);
};

#endif /* COMMON_DRAWER_DRAWER_H_ */
