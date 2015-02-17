/*
 * CompactIO.h
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#ifndef COMPACTIO_H_
#define COMPACTIO_H_

// Local includes
#include "ROOTOutput.h"
#include "CompactImport.h"

class CompactIO {
public:
	CompactIO();
	virtual ~CompactIO();

	bool openAll();
	bool closeAll();

	const std::string& getInputFileName() const {
		return inputFileName;
	}

	void setInputFileName(const std::string& inputFileName) {
		this->inputFileName = inputFileName;
	}

	const std::string& getOutputPrefix() const {
		return outputPrefix;
	}

	void setOutputPrefix(const std::string& outputPrefix) {
		this->outputPrefix = outputPrefix;
	}

	bool isIsInputList() const {
		return isInputList;
	}

	void setIsInputList(bool isInputList) {
		this->isInputList = isInputList;
	}

	ROOTCorrectedEvent& getCorrEvent() {
		return corrEvent;
	}

	void setCorrEvent(const ROOTCorrectedEvent& corrEvent) {
		this->corrEvent = corrEvent;
	}

	bool isDoOutput() const {
		return doOutput;
	}

	void setDoOutput(bool doOutput) {
		this->doOutput = doOutput;
	}

	ROOTFileHeader& getOutputFileHeader() {
		return outputFileHeader;
	}

	void setOutputFileHeader(const ROOTFileHeader& outputFileHeader) {
		this->outputFileHeader = outputFileHeader;
	}

	ROOTRawEvent& getRawEvent() {
		return rawEvent;
	}

	void setRawEvent(const ROOTRawEvent& rawEvent) {
		this->rawEvent = rawEvent;
	}

	ROOTBurst& getRootBurst() {
		return rootBurst;
	}

	void setRootBurst(const ROOTBurst& rootBurst) {
		this->rootBurst = rootBurst;
	}

	ROOTFileHeader& getInputFileHeader() {
		return inputFileHeader;
	}

	void setInputFileHeader(const ROOTFileHeader& inputFileHeader) {
		this->inputFileHeader = inputFileHeader;
	}

	NGeom& getRootGeom() {
		return rootGeom;
	}

	void setRootGeom(const NGeom& rootGeom) {
		this->rootGeom = rootGeom;
	}

	ROOTMCEvent& getRootMc() {
		return rootMC;
	}

	void setRootMc(const ROOTMCEvent& rootMc) {
		rootMC = rootMc;
	}

	ROOTPhysicsEvent& getRootPhysics() {
		return rootPhysics;
	}

	void setRootPhysics(const ROOTPhysicsEvent& rootPhysics) {
		this->rootPhysics = rootPhysics;
	}

public:
	CompactImport input;
	ROOTOutput output;

private:
	bool isInputList;
	bool doOutput;
	std::string inputFileName;
	std::string outputPrefix;
	ROOTBurst rootBurst;
	ROOTRawEvent rawEvent;
	ROOTCorrectedEvent corrEvent;
	ROOTFileHeader inputFileHeader;
	NGeom rootGeom;
	ROOTMCEvent rootMC;
	ROOTPhysicsEvent rootPhysics;
	ROOTFileHeader outputFileHeader;
};

#endif /* COMPACTIO_H_ */
