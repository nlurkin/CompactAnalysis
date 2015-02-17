/*
 * OptionsParser.h
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#ifndef OPTIONSPARSER_H_
#define OPTIONSPARSER_H_

#include <vector>
#include "mystructs.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

class CompactIO;

class OptionsParser {
public:
	OptionsParser();
	virtual ~OptionsParser();

	bool parse(int argc, char** argv, CompactIO &io);
	void printSummary(CompactIO &io);
	bool parseFilter();

	const std::vector<eventID>& getBadEventsList() const {
		return badEventsList;
	}

	bool isExportAllEvents() const {
		return exportAllEvents;
	}

	bool isOptDebug() const {
		return optDebug;
	}

	int getOutputModulo() const {
		return outputModulo;
	}

	int getPeriodKeep() const {
		return periodKeep;
	}

	int getMaxEvents() const {
		return maxEvents;
	}

private:
	int maxEvents;
	bool optDebug;
	int outputModulo;
	int periodKeep;
	bool exportAllEvents;
	po::options_description desc;
	std::string optString;
	std::string cutsFile;
	std::string filterFile;
	std::vector<eventID> badEventsList;
};

#endif /* OPTIONSPARSER_H_ */
