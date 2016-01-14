/*
 * MinuitFitter.h
 *
 *  Created on: Dec 28, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_MINUITFITTER_H_
#define COMMON_MINUITFITTER_H_

#include <TMinuit.h>
#include "../Samples/FitMCSample.h"
#include "../Samples/FitDataSample.h"
#include "../Interface/Sample.h"

class MinuitFitter : public TMinuit {
public:
	MinuitFitter();
	virtual ~MinuitFitter();

	void setSamples(FitMCSample* mcSamples, FitDataSample* dataSamples) {
		fMCSamples = mcSamples;
		fDataSamples = dataSamples;
	}

	virtual void fit();
	virtual double minFunction(double n, double a) = 0;

	static void minFunctionStatic(Int_t&, Double_t*, Double_t& f, Double_t* par, Int_t);

	void printResult();
	void drawResult(std::vector<Sample*> mcSamples, int nbins, double *binning);

	static double chi2pValue(double chi2, int ndof);
	void chi2Profile(TString name);

	void setName(const std::string& name) {
		fName = name;
	}

	double getFormFactor() const {
		return fFormFactor;
	}

	double getFormFactorErr() const {
		return fFormFactorErr;
	}

	double getMinimum() const {
		return fMinimum;
	}

	double getNorm() const {
		return fNorm;
	}

	double getNormErr() const {
		return fNormErr;
	}

	const std::string& getName() const {
		return fName;
	}

protected:
	double fNorm;
	double fNormErr;
	double fFormFactor;
	double fFormFactorErr;
	double fMinimum;

	FitMCSample* fMCSamples;
	FitDataSample* fDataSamples;

	std::string fName;
};

#endif /* COMMON_MINUITFITTER_H_ */
