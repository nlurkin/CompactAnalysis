/*
 * ConfigFile.cpp
 *
 *  Created on: Dec 21, 2015
 *      Author: nlurkin
 */

#include "ConfigFile.h"
#include <iostream>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <sys/stat.h>

namespace po = boost::program_options;

using namespace std;
ConfigFile::ConfigFile() :
	fDesc("Config File")
{
	// Declare the supported options.
	fDesc.add_options()
		("mcfiles", 	po::value< string>			(), 								"List of input MC files")
	    ("mcout", 		po::value< string>			(), 								"List of output MC files")
		("brs", 		po::value< string>			(), 								"List of input MC files")
		("mccolors", 	po::value< string>			(), 								"List of input MC files")
		("mclegends", 	po::value< string>			(), 								"List of input MC files")
		("mcIndex",		po::value< string>			(),				 					"List of input MC files")
		("datafiles", 	po::value< string>			(), 								"List of input MC files")
		("dataout", 	po::value< string>			(), 								"List of input MC files")
		("datacolors", 	po::value< string>			(), 								"List of input MC files")
		("datalegends", po::value< string>			(), 								"List of input MC files")
		("dataIndex",	po::value< string>			(),				 					"List of input MC files")
		("datafactor", 	po::value< string>			(), 								"List of input MC files")
		("modelfile", 	po::value< string>			(), 								"Model file for ToyMC")
		("runstart", 	po::value< int>				(&fRunStart)->default_value(0), 	"List of input MC files")
		("runend", 		po::value< int>				(&fRunEnd)->default_value(0), 		"List of input MC files")
		("periodstart", po::value< int>				(&fPeriodStart)->default_value(-1), "List of input MC files")
		("periodend", 	po::value< int>				(&fPeriodEnd)->default_value(-1), 	"List of input MC files")
		("testA", 		po::value< double>			(&fTestA)->default_value(0), 		"List of input MC files")
		("binsfile", 	po::value< string>			(&fBinsFileName), 					"List of input MC files")
		("equalbin", 	po::value< bool>			(&fWithEqualBins), 					"List of input MC files")
		("scanid", 		po::value< int>				(&fScanID)->default_value(-1),		"List of input MC files")
		("nscan", 		po::value< int>				(&fNScan)->default_value(1),	 	"List of input MC files")
		("startscan", 	po::value< int>				(&fStartScan)->default_value(-1),	"List of input MC files")
		("endscan", 	po::value< int>				(&fEndScan)->default_value(-1),	 	"List of input MC files")
		("maxloss", 	po::value< double>			(&fMaxLoss)->default_value(0.2),	"List of input MC files")
		("usepk", 		po::value< bool>			(&fUsePk)->default_value(false),	"List of input MC files")
		("weightsfile", po::value< string>			(&fWeightFile)->default_value("pi0dalits_weights.dat"),	"List of input MC files")
	;
}

ConfigFile::~ConfigFile() {
	// TODO Auto-generated destructor stub
}

vector<TString> tokenize(string s){
	TString temp(s);
	vector<TString> ret;
	TObjArray *tok = temp.Tokenize(" ");
	for(int i=0; i<tok->GetEntries(); ++i){
		TString entry(((TObjString*)tok->At(i))->GetString());
		ret.push_back(entry);
	}

	return ret;
}

bool ConfigFile::readFile(string fileName){
	po::variables_map vm;
	po::store(po::parse_config_file<char>(fileName.c_str(), fDesc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << fDesc << "\n";
	    return 1;
	}

	vector<TString> tempVec;

	if (vm.count("mcfiles")) {
		tempVec = tokenize(vm["mcfiles"].as<string>());
		for(auto it : tempVec) fMCFileNames.push_back(it.Data());
	}
	if (vm.count("mcout")) {
		tempVec = tokenize(vm["mcout"].as<string>());
		for(auto it : tempVec) fMCOutputFiles.push_back(it.Data());
	}
	if (vm.count("brs")) {
		tempVec = tokenize(vm["brs"].as<string>());
		for(auto it : tempVec) fBrs.push_back(it.Atof());
	}
	if (vm.count("mccolors")) {
		tempVec = tokenize(vm["mccolors"].as<string>());
		for(auto it : tempVec) fMCColors.push_back(it.Atoi());
	}
	if (vm.count("mclegends")) {
		tempVec = tokenize(vm["mclegends"].as<string>());
		for(auto it : tempVec) fMCLegendTitle.push_back(it.ReplaceAll("\\", "#").Data());
	}
	if (vm.count("mcIndex")) {
		tempVec = tokenize(vm["mcIndex"].as<string>());
		for(auto it : tempVec) fMCIndexes.push_back(it.Atoi());
	}

	if (vm.count("datafiles")) {
		tempVec = tokenize(vm["datafiles"].as<string>());
		for(auto it : tempVec) fDataFileNames.push_back(it.Data());
	}
	if (vm.count("dataout")) {
		tempVec = tokenize(vm["dataout"].as<string>());
		for(auto it : tempVec) fDataOutputFiles.push_back(it.Data());
	}
	if (vm.count("datacolors")) {
		tempVec = tokenize(vm["datacolors"].as<string>());
		for(auto it : tempVec) fDataColors.push_back(it.Atoi());
	}
	if (vm.count("datalegends")) {
		tempVec = tokenize(vm["datalegends"].as<string>());
		for(auto it : tempVec) fDataLegendTitle.push_back(it.ReplaceAll("\\", "#").Data());
	}
	if (vm.count("dataIndex")) {
		tempVec = tokenize(vm["dataIndex"].as<string>());
		for(auto it : tempVec) fDataIndexes.push_back(it.Atoi());
	}
	if (vm.count("datafactor")) {
		tempVec = tokenize(vm["datafactor"].as<string>());
		for(auto it : tempVec) fDataFactor.push_back(it.Atof());
	}


	if (vm.count("modelfile")) {
		tempVec = tokenize(vm["modelfile"].as<string>());
		for(auto it : tempVec) fModelFiles.push_back(it.Data());
	}

	return true;
}

void ConfigFile::print() const{
	cout << "MC Options" << endl;
	for(unsigned int i=0; i<fMCFileNames.size(); i++){
		cout << fMCFileNames[i] << " " << fMCIndexes[i] << " " << fMCOutputFiles[i] << " " << fBrs[i] << " " << fMCColors[i] << " " << fMCLegendTitle[i] << endl;
	}

	cout << "Data Options" << endl;
	for(unsigned int i=0; i<fDataFileNames.size(); i++){
		cout << fDataFileNames[i] << " " << fDataIndexes[i] << " " << fDataOutputFiles[i] << " " << fDataFactor[i] << " " << fDataColors[i] << " " << fDataLegendTitle[i] << endl;
	}

	cout << "Other options" << endl;
	cout << "Model files:";
	for(auto it : fModelFiles) cout << it << " ";
	cout << endl;
	cout << "From run " << fRunStart << " to run " << fRunEnd << endl;
	cout << "From period " << fPeriodStart << " to period " << fPeriodEnd << endl;
	cout << "TestA: " << fTestA << endl;
	cout << "Bins file: " << fBinsFileName << endl;
	cout << "With equal bins: " << fWithEqualBins << endl;
	cout << "Scan ID: " << fScanID << endl;
	cout << "# Scan: " << fNScan << endl;
	cout << "Use scan from " << fStartScan << " to " << fEndScan << endl;
}

bool ConfigFile::testAllOutputFiles() const{
	struct stat fStat;
	bool error = true;

	for(auto it: fMCOutputFiles){
		if(stat(it.c_str(), &fStat) == 0){
			cout << "[ERROR] Output file " << it << " already exists" << endl;
			error = false;
		}
	}

	for(auto it : fDataOutputFiles){
		if(stat(it.c_str(), &fStat) == 0){
			cout << "[ERROR] Output file " << it << " already exists" << endl;
			error = false;
		}
	}
	return error;
}

bool ConfigFile::testUseRun(int run, int period) const{
	bool ret = true;
	if(fRunStart!=0 && run<fRunStart) ret = false;
	if(fRunEnd!=0 && run>fRunEnd) ret = false;
	if(fPeriodStart!=-1 && period<fPeriodStart) ret = false;
	if(fPeriodEnd!=-1 && period>fPeriodEnd) ret = false;
	return ret;
}
