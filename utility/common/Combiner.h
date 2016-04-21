/*
 * Combiner.h
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#ifndef COMMON_COMBINER_H_
#define COMMON_COMBINER_H_

#include "DataGetter.h"
#include "Samples/CombineMCSample.h"
#include "Samples/CombineDataSample.h"

class Combiner: public DataGetter {
public:
	Combiner();
	virtual ~Combiner();

	void draw(std::vector<int> allColors, std::vector<int> dataColors, int first, int last);
};

#endif /* COMMON_COMBINER_H_ */
