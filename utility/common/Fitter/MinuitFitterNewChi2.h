/*
 * MinuitFitterNewChi2.h
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_MINUITFITTERNEWCHI2_H_
#define COMMON_MINUITFITTERNEWCHI2_H_

#include <TH1D.h>
#include "../Fitter/MinuitFitter.h"

class MinuitFitterNewChi2: public MinuitFitter {
public:
	MinuitFitterNewChi2(int bins);
	virtual ~MinuitFitterNewChi2();

	void fit();
	double minFunction(double n, double a);

	double function(double G, double a, double a_i, double b_i, double g_i);

private:
	int fNBins;
};

#endif /* COMMON_MINUITFITTERNEWCHI2_H_ */
