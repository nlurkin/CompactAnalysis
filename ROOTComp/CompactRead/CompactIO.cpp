/*
 * CompactIO.cpp
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#include "CompactIO.h"

CompactIO::CompactIO() {
	// TODO Auto-generated constructor stub

}

CompactIO::~CompactIO() {
	// TODO Auto-generated destructor stub
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
