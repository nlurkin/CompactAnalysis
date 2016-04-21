/*
 * MCSample.h
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_SAMPLES_MCSAMPLE_H_
#define COMMON_SAMPLES_MCSAMPLE_H_

class MCSample{
public:
	MCSample();
	virtual ~MCSample();

	void setUsePk(int usePk) {
		fUsePk = usePk;
	}

protected:
	int fUsePk;
};

#endif /* COMMON_SAMPLES_MCSAMPLE_H_ */
