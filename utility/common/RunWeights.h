/*
 * RunWeights.h
 *
 *  Created on: Dec 22, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_RUNWEIGHTS_H_
#define COMMON_RUNWEIGHTS_H_

#include <string>
#include <map>

class RunWeights {
public:
	RunWeights();
	virtual ~RunWeights();

	bool loadWeights(std::string fileName);
	double applyWeights(int nrun) const;

private:
	std::map<int,double> fRatioMap;
	double fAverageRatio;
};

#endif /* COMMON_RUNWEIGHTS_H_ */
