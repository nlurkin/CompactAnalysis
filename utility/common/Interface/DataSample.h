/*
 * DataSample.h
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLES_DATASAMPLE_H_
#define COMMON_SAMPLES_DATASAMPLE_H_

class DataSample{
public:
	DataSample();
	virtual ~DataSample();

	double getTestA() const {
		return fTestA;
	}

	void setTestA(double testA) {
		fTestA = testA;
	}

	double getFactor() const {
		return fFactor;
	}

	void setFactor(double factor) {
		fFactor = factor;
	}
protected:
	double fTestA;
	double fFactor;
};

#endif /* COMMON_SAMPLES_DATASAMPLE_H_ */
