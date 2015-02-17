/*
 * CompactIO.cpp
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

//Local includes
#include "CompactIO.h"

CompactIO::CompactIO(): isInputList(false), doOutput(false) {
}

CompactIO::~CompactIO() {
}

bool CompactIO::openAll() {
	input.associateTrees(rawEvent, corrEvent, rootBurst, inputFileHeader,
			rootGeom, rootMC);
	if (!input.readInput(inputFileName, isInputList))
		return false;
	output.openOutput(input.getHasMC(), doOutput, rootBurst, rawEvent,
			corrEvent, rootGeom, rootMC, rootPhysics, outputFileHeader);

	return true;
}

bool CompactIO::closeAll() {
	output.close();
	return true;
}

