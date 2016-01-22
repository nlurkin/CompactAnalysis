/*
 * Fitter.h
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITTER_H_
#define COMMON_FITTER_H_

#include <sys/stat.h>
#include "../DataGetter.h"

class MinuitFitter;

class Fitter : public DataGetter{
public:
	Fitter();
	virtual ~Fitter();

	void PrepareHistos(std::vector<int> allColors, std::vector<int> dataColors);

	void fit(bool useNew, bool userROOT, double maxLoss, int start=-1, int end=-1);

	double getNormalization(double a);

	static void minFunctionStatic(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,	Int_t flag);

private:
	static int findMaxSample(std::vector<MinuitFitter*> vFit, std::vector<int> use);
	static int findMinSample(std::vector<MinuitFitter*> vFit, std::vector<int> use);
	static double computeAverage(std::vector<MinuitFitter*> vFit, std::vector<int> use);
	static std::vector<int> findWithEnoughStat(Sample* finalDataSample, int dfltSample, int start, int end, double maxLoss);

	std::vector<MinuitFitter*> fFitters;
};

#endif /* COMMON_FITTER_H_ */
