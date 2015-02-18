/*
 * OptionsParser.cpp
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#include "OptionsParser.h"
#include "CompactIO.h"

OptionsParser::OptionsParser() :
		maxEvents(-1), optDebug(false), outputModulo(1), periodKeep(0), exportAllEvents(false), doScan(false), nScan(0), desc("Allowed options") {
	// Declare the supported options.
	desc.add_options()("help,h", "produce help message")("nevt,n", po::value<int>(), "max number of events")
			("file,i",po::value<std::string>(), "input file name")
			("list,l", po::value<std::string>(), "list of input files")
			("prefix,p", po::value<std::string>(), "prefix for output files")
			("debug,d",	po::value<bool>()->default_value(false), "Activate verbose debugging")
			("filter,f", po::value<std::string>(), "Filter file")
			("dooutput", po::value<bool>()->default_value(false), "Activate output text files")
			("period", po::value<int>()->default_value(0), "Keep only events from specified period")
			("mod,m", po::value<int>()->default_value(1), "Event number printing modulo")
			("cuts,c", po::value<std::string>(), "Cuts file")
			("scan", po::value<int>()->default_value(0), "Do a scan with n values")
			("eall,e", po::value<bool>()->default_value(false), "Export all events, even failed");
}

OptionsParser::~OptionsParser() {
	// TODO Auto-generated destructor stub
}

bool OptionsParser::parse(int argc, char** argv, CompactIO &io) {
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return false;
	}

	std::cout << ">>>> Received parameters:" << std::endl;
	for (const auto &it : vm) {
		std::cout << "\t" << it.first << " = ";
		auto& value = it.second.value();
		if (auto v = boost::any_cast<int>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<bool>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<std::string>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<TString>(&value))
			std::cout << *v;
		else
			std::cout << "**cast error**";
		std::cout << std::endl;
	}

	if (vm.count("nevt"))
		maxEvents = vm["nevt"].as<int>();

	/// String options
	if (vm.count("prefix"))
		io.setOutputPrefix(vm["prefix"].as<std::string>());
	if (vm.count("debug"))
		optDebug = vm["debug"].as<bool>();
	if (vm.count("period"))
		periodKeep = vm["period"].as<int>();
	if (vm.count("dooutput"))
		io.setDoOutput(vm["dooutput"].as<bool>());
	if (vm.count("mod"))
		outputModulo = vm["mod"].as<int>();
	if (vm.count("cuts"))
		cutsFile = vm["cuts"].as<std::string>();
	if (vm.count("eall"))
		exportAllEvents = vm["eall"].as<bool>();
	if (vm.count("scan")){
			nScan = vm["scan"].as<int>();
			doScan = nScan > 0;
	}
	if (vm.count("filter"))
		filterFile = vm["filter"].as<std::string>();

	/// Input (one of them mandatory)
	if (vm.count("list") and vm.count("file")) {
		std::cerr
				<< "Cannot specify both an input file and an input list at the same time"
				<< std::endl;
		return false;
	} else if (!vm.count("list") and !vm.count("file")) {
		std::cerr << "Must specify either an input file or an input list"
				<< std::endl;
		return false;
	} else {
		if (vm.count("list")) {
			io.setIsInputList(true);
			io.setInputFileName(vm["list"].as<std::string>());
		}
		if (vm.count("file"))
			io.setInputFileName(vm["file"].as<std::string>());
	}

	return true;
}

void OptionsParser::printSummary(CompactIO &io) {
	std::cout << "Output events every " << outputModulo << " events"
			<< std::endl;
	if (!cutsFile.empty())
		std::cout << "Using cuts at: " << cutsFile << std::endl;
	std::cout << "Debugging activated: " << (optDebug == true ? "Yes" : "No")
			<< std::endl;
	if (periodKeep == 0)
		std::cout << "Keeping period: All" << std::endl;
	else
		std::cout << "Keeping period: " << periodKeep << std::endl;
	if (io.isDoOutput())
		std::cout << "Text output files requested" << std::endl;
	if (exportAllEvents)
		std::cout << "Export all events requested" << std::endl;
	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cout << std::endl << std::endl;
}

bool OptionsParser::parseFilter() {
	int runNum, burstNum, timestamp;
	std::vector<eventID>::iterator it;

	if (!filterFile.empty()) {
		std::cout << ">>>> Filtering events from file " << filterFile
				<< std::endl;
		std::ifstream badEvents;
		badEvents.open(filterFile.c_str(), ifstream::in);
		if (badEvents.is_open()) {
			while (badEvents >> runNum >> burstNum >> timestamp) {
				badEventsList.push_back(eventID(runNum, burstNum, timestamp));
			}
			badEvents.close();
		} else {
			std::cerr << "Unable to open filter file" << std::endl;
			return false;
		}
		std::cout << "\t" << badEventsList.size() << " events in filter list"
				<< std::endl;
		for (auto it : badEventsList) {
			std::cout << "\t\t" << it.rnum << " " << it.bnum << " "
					<< it.timestamp << std::endl;
		}
		return true;
	}
	return false;
}
