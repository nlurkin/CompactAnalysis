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

bool CompactIO::openAll(bool doScan) {
	input.associateTrees(rawEvent, corrEvent, rootBurst, inputFileHeader,
			rootGeom, rootMC);
	if (!input.readInput(inputFileName, isInputList))
		return false;
	output.openOutput(input.getHasMC(), doOutput, doScan, rootBurst, rawEvent,
			corrEvent, rootGeom, rootMC, rootPhysics, outputFileHeader, cutsDefinition);


	if(doScan) fillCutsList();
	return true;
}

bool CompactIO::closeAll() {
	output.close();
	return true;
}

void CompactIO::fillCutsList() {
	for(int i=cutsDefinition.loadList(0); i!=-1; i=cutsDefinition.loadNextList()){
		output.fillCuts();
	}
}
