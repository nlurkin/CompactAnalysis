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

	void fit(bool useNew, bool userROOT);

	double getNormalization(double a);

	static void minFunctionStatic(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,	Int_t flag);

private:
	std::vector<MinuitFitter*> fFitters;
};

#endif /* COMMON_FITTER_H_ */
