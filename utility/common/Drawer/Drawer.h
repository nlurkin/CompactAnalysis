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

	static void drawFitScan(std::vector<MinuitFitter*> fit, std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, Sample* finalMCSample, Sample* finalDataSample, std::vector<int> use);

	static void drawFitResult(std::vector<MinuitFitter*> fit, std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, Sample* finalMCSample, Sample* finalDataSample, int index);
	static void drawFitPreparation(std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, std::string title);

	static void drawCombineStack(std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples, Sample* finalMCSample, Sample* finalDataSample);
};

#endif /* COMMON_DRAWER_DRAWER_H_ */
