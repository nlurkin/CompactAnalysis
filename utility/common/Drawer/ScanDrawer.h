/*
 * ScanDrawer.h
 *
 *  Created on: Jan 15, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_DRAWER_SCANDRAWER_H_
#define COMMON_DRAWER_SCANDRAWER_H_

#include <vector>

class TPad;

class ScanDrawer {
public:
	ScanDrawer();
	virtual ~ScanDrawer();

	void generateResult(TPad* pad);
	void generateNSelected(TPad* pad);
	void draw();

	void addScanValue(double val, double result, double err, int ndata);
	void computeUncorrError();

	int getDefaultCutValue() const {
		return fDefaultCutValue;
	}

	void setDefaultCutValue(int defaultCutValue) {
		fDefaultCutValue = defaultCutValue;
	}

private:
	std::vector<double> fScanValues;
	std::vector<double> fResultValues;
	std::vector<double> fResultErrors;
	std::vector<double> fUncorrErrors;
	std::vector<int> fNSelected;
	int fDefaultCutValue;
};

#endif /* COMMON_DRAWER_SCANDRAWER_H_ */
