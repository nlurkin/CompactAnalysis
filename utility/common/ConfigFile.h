/*
 * ConfigFile.h
 *
 *  Created on: Dec 21, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_CONFIGFILE_H_
#define COMMON_CONFIGFILE_H_

#include <string>
#include <vector>
#include <boost/program_options.hpp>

class ConfigFile {
public:
	ConfigFile();
	virtual ~ConfigFile();

	bool readFile(std::string fileName);
	void print();
	bool testAllOutputFiles();
	bool testUseRun(int run, int period);

	const std::string& getBinsFileName() const {
		return fBinsFileName;
	}

	void setBinsFileName(const std::string& binsFileName) {
		fBinsFileName = binsFileName;
	}

	const std::vector<double>& getBrs() const {
		return fBrs;
	}

	void setBrs(const std::vector<double>& brs) {
		fBrs = brs;
	}

	const std::vector<int>& getDataColors() const {
		return fDataColors;
	}

	void setDataColors(const std::vector<int>& dataColors) {
		fDataColors = dataColors;
	}

	const std::vector<double>& getDataFactor() const {
		return fDataFactor;
	}

	void setDataFactor(const std::vector<double>& dataFactor) {
		fDataFactor = dataFactor;
	}

	const std::vector<std::string>& getDataFileNames() const {
		return fDataFileNames;
	}

	void setDataFileNames(const std::vector<std::string>& dataFileNames) {
		fDataFileNames = dataFileNames;
	}

	const std::vector<std::string>& getDataLegendTitle() const {
		return fDataLegendTitle;
	}

	void setDataLegendTitle(const std::vector<std::string>& dataLegendTitle) {
		fDataLegendTitle = dataLegendTitle;
	}

	const std::vector<std::string>& getDataOutputFiles() const {
		return fDataOutputFiles;
	}

	void setDataOutputFiles(const std::vector<std::string>& dataOutputFiles) {
		fDataOutputFiles = dataOutputFiles;
	}

	const std::vector<int>& getMcColors() const {
		return fMCColors;
	}

	void setMcColors(const std::vector<int>& mcColors) {
		fMCColors = mcColors;
	}

	const std::vector<std::string>& getMcFileNames() const {
		return fMCFileNames;
	}

	void setMcFileNames(const std::vector<std::string>& mcFileNames) {
		fMCFileNames = mcFileNames;
	}

	const std::vector<int>& getMcIndexes() const {
		return fMCIndexes;
	}

	void setMcIndexes(const std::vector<int>& mcIndexes) {
		fMCIndexes = mcIndexes;
	}

	const std::vector<std::string>& getMcLegendTitle() const {
		return fMCLegendTitle;
	}

	void setMcLegendTitle(const std::vector<std::string>& mcLegendTitle) {
		fMCLegendTitle = mcLegendTitle;
	}

	const std::vector<std::string>& getMcOutputFiles() const {
		return fMCOutputFiles;
	}

	void setMcOutputFiles(const std::vector<std::string>& mcOutputFiles) {
		fMCOutputFiles = mcOutputFiles;
	}

	const std::vector<std::string>& getModelFiles() const {
		return fModelFiles;
	}

	void setModelFiles(const std::vector<std::string>& modelFiles) {
		fModelFiles = modelFiles;
	}

	int getPeriodEnd() const {
		return fPeriodEnd;
	}

	void setPeriodEnd(int periodEnd) {
		fPeriodEnd = periodEnd;
	}

	int getPeriodStart() const {
		return fPeriodStart;
	}

	void setPeriodStart(int periodStart) {
		fPeriodStart = periodStart;
	}

	int getRunEnd() const {
		return fRunEnd;
	}

	void setRunEnd(int runEnd) {
		fRunEnd = runEnd;
	}

	int getRunStart() const {
		return fRunStart;
	}

	void setRunStart(int runStart) {
		fRunStart = runStart;
	}

	int getScanId() const {
		return fScanID;
	}

	void setScanId(int scanId) {
		fScanID = scanId;
	}

	double getTestA() const {
		return fTestA;
	}

	void setTestA(double testA) {
		fTestA = testA;
	}

	bool isWithEqualBins() const {
		return fWithEqualBins;
	}

	void setWithEqualBins(bool withEqualBins) {
		fWithEqualBins = withEqualBins;
	}

private:
	boost::program_options::options_description fDesc;
	int fRunStart;
	int fRunEnd;

	int fPeriodStart;
	int fPeriodEnd;

	double fTestA;
	int fScanID;

	std::vector<int> fMCColors, fDataColors;
	std::vector<std::string> fMCLegendTitle, fDataLegendTitle;

	std::vector<double> fBrs;
	std::vector<std::string> fMCFileNames;
	std::vector<std::string> fDataFileNames;
	std::vector<int> fMCIndexes;
	std::vector<double> fDataFactor;
	std::vector<std::string> fMCOutputFiles;
	std::vector<std::string> fDataOutputFiles;
	std::vector<std::string> fModelFiles;
	std::string fBinsFileName;
	bool fWithEqualBins;

};

#endif /* COMMON_CONFIGFILE_H_ */
