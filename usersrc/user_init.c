/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <sstream>
#include <iostream>
#include <map>
#include <iomanip>
#include "funLib.h"
#include "compactLib.h"
using namespace std;

int nico_ke2Init(){
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	//goodTracks.clear();
	//assocClusters.clear();
	outTree = new TTree("event", "Event");
	//outTree->Branch("goodTrack", "std::vector<CorrectedTrack>", &goodTracks, 64000, 1);
	//outTree->Branch("assocCluster", "std::vector<CorrectedCluster>", &assocClusters, 64000, 1);

	new TH1D("Cuts", "Failed Cuts", 20, 0.5, 20.5);
	new TH1D("CutsKe2", "Failed Ke2 Cuts", 20, 0.5, 20.5);

	//Histograms for common cuts
	new TH1D("timeDchOffset", "timeDchOffset", 100, 0, 100);
	new TH1D("deadCellDist", "deadCellDist", 100, 0, 15);
	new TH2D("DCHAcceptance1", "DCHAcceptance1", 400, -200, 200, 400, -200, 200);
	new TH2D("DCHAcceptance4", "DCHAcceptance4", 400, -200, 200, 400, -200, 200);
	new TH1D("RDCH1", "RDCH1", 120, 0, 120);
	new TH1D("RDCH4", "RDCH4", 140, 0, 140);
	new TH2D("excludedZone", "excludedZone", 300, -300, 300, 600, -300, 300);
	new TH2D("hotCellDistT", "hotCellDistT", 400, -200, 200, 400, -200, 200);
	new TH2D("hotCellDistC", "hotCellDistC", 400, -200, 200, 400, -200, 200);
	new TH2D("lkrAcceptance", "lkrAcceptance", 400, -200, 200, 400, -200, 200);
	new TH1D("cda", "cda", 200, 0, 20);

	new TH1D("trackMomentum", "trackMomentum", 500, 0, 100);
	new TH1D("trackQuality", "trackQuality", 200, 0, 2);

	//Histogram for ke2 Cuts
	new TH1D("mmass2", "mmass2", 400, -1, 1);
	new TH1D("clusterStatus", "clusterStatus", 10, 0, 10);
	new TH1D("trckClDist", "trckClDist", 100, 0, 10);
	new TH1D("trckClTime", "trckClTime", 300, 0, 30);
	new TH1D("eop", "eop", 100, 0, 2);

	//Final histos
	new TH1D("fmmass2", "fmmass2", 400, -1, 1);

	return 0;
}

int nico_pi0DalitzInit(){
	if(opts.count("ff")!=0) ffWeightType = atoi(opts["ff"].c_str());
	else ffWeightType = -1;

	rawEvent.clear();
	corrEvent.clear();
	outTree = new TTree("event", "Event");
	headerTree = new TTree("header", "Header");
	//outTree->Branch("goodTrack", "std::vector<CorrectedTrack*>", &goodTracks, 64000, 1);
	//outTree->Branch("assocCluster", "std::vector<CorrectedCluster*>", &assocClusters, 64000, 1);
	//outTree->Branch("cutsWord", &cutsWord, "cutsWord[19]/O");

	outTree->Branch("pi0dBurst" ,"ROOTBurst", &rootBurst);
	outTree->Branch("rawEvent" ,"ROOTRawEvent", &rawEvent);
	outTree->Branch("corrEvent" ,"ROOTCorrectedEvent", &corrEvent);
	outTree->Branch("geom" ,"NGeom", &rootGeom);
	headerTree->Branch("header" ,"ROOTFileHeader", &rootFileHeader);

	//Vertex
	new TH1D("vertexN", "vertexN", 10, 0, 10); 			///
	new TH1D("vertexZ", "vertexZ", 1300, -3000, 10000); ///
	new TH1D("vertexChi2", "vertexChi2", 1000, 0, 100); 	///
	new TH1D("vertexQ", "vertexQ", 10, -5, 5); 			///

	//Cluster
	new TH1D("clusterN", "clusterN", 10, 0, 10);					///
	new TH1D("clusterX", "clusterX", 600, -150, 150);			///
	new TH1D("clusterY", "clusterY", 600, -150, 150);			///
	new TH1D("clusterRadius", "clusterRadius", 600, 0, 150);		///
	new TH1D("clusterDPi", "clusterDPi", 500, 0, 250);			///
	new TH1D("clusterDep", "clusterDep", 500, 0, 250);			///
	new TH1D("clusterDem", "clusterDem", 500, 0, 250);			///
	new TH1D("clusterDUnep", "clusterDUnep", 500, 0, 250);		///
	new TH1D("clusterDUnem", "clusterDUnem", 500, 0, 250);		///
	new TH1D("clusterVertexTime", "clusterVertexTime", 200, 0, 200);///
	new TH1D("clusterGoodN", "clusterGoodN", 10, 0, 10);			///

	//Tracks
	new TH1D("trackN", "trackN", 10, 0, 10); 						///
	new TH1D("trackDCHTime", "trackDCHTime", 100, 0, 100);			///
	new TH1D("trackTrackTime", "trackTrackTime", 100, 0, 50); 		///
	new TH1D("trackVertexTime", "trackVertexTime", 100, 0, 100); 	///
	new TH1D("trackDCH1Radius", "trackDCH1Radius", 300, 0, 150); 	///
	new TH1D("trackDCH2Radius", "trackDCH2Radius", 300, 0, 150); 	///
	new TH1D("trackDCH4Radius", "trackDCH4Radius", 300, 0, 150); 	///
	new TH1D("trackLKrRadius", "trackLKrRadius", 400, 0, 200); 	///
	new TH1D("trackEOP", "trackEOP", 150, 0, 1.5);					///
	new TH1D("trackP", "trackP", 800, 0, 80); 						///
	new TH1D("track_ij_DCH1", "track_ij_DCH1", 400, 0, 200);		///
	new TH1D("track_ee_LKr", "track_ee_LKr", 500, 0, 250);		///
	new TH1D("track_epi_LKr", "track_epi_LKr", 500, 0, 250);		///

	//Gamma cluster
	new TH1D("gammaX", "gammaX", 600, -150, 150);				///
	new TH1D("gammaY", "gammaY", 600, -150, 150);				///
	new TH1D("gammaRadius", "gammaRadius", 300, 0, 150);			///
	new TH1D("gammaEnergy", "gammaEnergy", 800, 0, 80);				///
	new TH1D("gammaDeadCell", "gammaDeadCell", 60, 0, 30);			///
	new TH1D("gammaDCH1Radius", "gammaDCH1Radius", 300, 0, 150);	///

	//Kinematic before cuts
	new TH1D("bkinKP", "bkinKP", 100, 0, 80);			///
	new TH1D("bkinPiP", "bkinKPi", 100, 0, 80);			///
	new TH1D("bkinepP", "bkinKep", 100, 0, 80);			///
	new TH1D("bkinemP", "bkinKem", 100, 0, 80);			///
	new TH1D("bkinGammaP", "bkinKgamma", 100, 0, 80);	///
	new TH1D("bkinMee", "bkinMee", 160, 0, 0.16);			///
	new TH1D("bkinx", "bkinx", 1000, 0, 1);			///
	new TH1D("bkinMeeg", "bkinMeeg", 500, 0, 0.5);		///
	new TH1D("bkinMeegDiff", "bkinMeegDiff", 400, -0.2, 0.2);		///
	new TH1D("bkinMeegpi", "bkinMeegpi", 1000, 0, 0.6);	///
	new TH1D("bkinPTot", "bkinPTot", 300, 65, 80);		///
	new TH1D("bkinPt2", "bkinPt2", 1000, 0, 0.01);			///

	//Kinematic after cuts
	new TH1D("akinKP", "akinKP", 100, 0, 80);			///
	new TH1D("akinPiP", "akinKPi", 100, 0, 80);			///
	new TH1D("akinepP", "akinKep", 100, 0, 80);			///
	new TH1D("akinemP", "akinKem", 100, 0, 80);			///
	new TH1D("akinGammaP", "akinKgamma", 100, 0, 80);	///
	new TH1D("akinMee", "akinMee", 160, 0, 0.16);			///
	new TH1D("akinx", "akinx", 1000, 0, 1);			///
	new TH1D("akinMeeg", "akinMeeg", 500, 0, 0.5);		///
	new TH1D("akinMeegDiff", "akinMeegDiff", 400, -0.2, 0.2);		///
	new TH1D("akinMeegpi", "akinMeegpi", 1000, 0, 0.6);	///
	new TH1D("akinPTot", "akinPTot", 300, 65, 80);		///
	new TH1D("akinPt2", "akinPt2", 1000, 0, 0.01);			///



	//Selection histo
	new TH1D("Cuts", "Failed Cuts", 20, 0.5, 20.5);

	//Ana histo
	new TH1D("sel_NVertices", "Number of vertices", 10, 0, 10);
	new TH1D("sel_NTracks", "Number of tracks", 10, 0, 10);

	new TH1D("sel_DCH1Rad", "Extrapolated tracks radius on DCH1", 150, 0, 150);

	new TH1D("peegMass", "peegMass", 1000, 0.4, 0.6);
	new TH1D("eegMass", "eegMass", 1000, 0, 0.2);
	new TH1D("eeMass", "eeMass", 1000, 0, 0.2);

	new TH1D("xDistrib", "xDistrib", 1000, 0, 1);
	return 0;
}


int common_init(string filePrefix){
	int i, j, k;
	int l, m, n;
	int cpd, cell;

	int runNum, burstNum, timestamp;
	vector<eventID>::iterator it;



	string outRoot = "outfile.root";
	string outFile = "compact.txt";
	string outPass = "compactpass.txt";
	if(filePrefix.find('~')!=string::npos) filePrefix=filePrefix.replace(filePrefix.find('~'), 1, string("/afs/cern.ch/user/n/nlurkin"));
	if(filePrefix.length()>0){
		outRoot = filePrefix + ".root";
		outFile = filePrefix + ".txt";
		outPass = filePrefix + "pass.txt";
	}

	eopLoaded = false;

	CELLlength = 1.975;
	CPDlength = 8 * CELLlength;

	// Define the positions for CPDs and Cells and store them
	for (i=0; i<16; i++)
		for (j=0; j<16; j++)
		{
			k = i*16 + j;
			CPDpos_leftDownCorner[k][0] = (-1)*(7-i)*CPDlength;  // LKr RF is left-handed  --> the x sign has to be changed !!!
			CPDpos_leftDownCorner[k][1] = (7-j)*CPDlength;
			//printf ("CPD %d: position left down corner = %.2f, \t%.2f\n", k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1]);

			for (m=0; m<8; m++)
				for (n=0; n<8; n++)
				{
					l = m*8 + n;
					CELLpos_leftDownCorner[k][l][0] = CPDpos_leftDownCorner[k][0] - (7-m)*CELLlength;
					CELLpos_leftDownCorner[k][l][1] = CPDpos_leftDownCorner[k][1] + (7-n)*CELLlength;
					//printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
					//      l, k, CELLpos_leftDownCorner[k][l][0], CELLpos_leftDownCorner[k][l][1]);
				}
		}

	//Load badEvents list for debugging
	if(opts.count("filter")>0){
		cout << ">>>> Filtering events from file " << opts["filter"] << endl;
		FILE *badEvents = fopen(opts["filter"].c_str(), "r");
		if(badEvents!=NULL){
			while(fscanf(badEvents, "%i %i %i", &runNum, &burstNum, &timestamp) != EOF){
				badEventsList.push_back(eventID(runNum, burstNum, timestamp));
			}
			fclose(badEvents);
		}
		else{
			cout << "Unable to open filter file" << endl;
		}
		cout << "\t" << badEventsList.size() << " events in filter list" << endl;
		for(it=badEventsList.begin(); it!=badEventsList.end();it++){
			cout << "\t\t" << (*it).rnum << " " << (*it).bnum << " " << (*it).timestamp << endl;
		}
	}

	if(opts.count("nooutput")==0){
		noOutput = false;
		fprt=fopen(outFile.c_str(),"w");
		fprt2=fopen(outPass.c_str(),"w");
	}
	else noOutput = true;

	gFile = TFile::Open(outRoot.c_str(), "RECREATE");

	return 0;
}

int user_init() {
	map<string,string>::iterator it;

	opts = parseOptions(gString);

	cout << endl << ">>>>>>>>>>>>>>>>>>>>> Initialization" << endl;
	if(opts.count("h")!=0){
		cout << ">>>> Help " << endl;
		cout << ">>>> Syntax: param=value:param=value" << endl;
		cout << ">>>> List of parameters:" << endl;
		cout << ">>>> h: This help" << endl;
		cout << ">>>> prefix: Output file names to use (without extension)" << endl;
		cout << ">>>> can: ke2 | pi0d" << endl;
		cout << ">>>> debug: Activate debugging" << endl;
		cout << ">>>> ff: Type of form factor (0=1, 1=x, 2=x^2)" << endl;
		cout << ">>>> nooutput: Don't create output txt files" << endl;
		cout << ">>>> period: keep only events from this period" << endl;
		cout << ">>>> mod: print events index every mod events" << endl;
		cout << ">>>> cuts: specify cuts file" << endl;
		cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		exit(0);
	}
	cout << ">>>> Received parameters:" << endl;
	for(it=opts.begin(); it!=opts.end(); it++){
		cout << "\t" << it->first << " = " << it->second << endl;
	}

	common_init(opts["prefix"]);
	string chanName;

	if(strcmp(opts["can"].c_str(), "ke2")==0){
		channel=KE2;
		chanName = "KE2";
	}
	else if(strcmp(opts["can"].c_str(), "pi0d")==0){
		channel=PI0DALITZ;
		chanName = "PI0Dalitz";
	}
	else if(strcmp(opts["can"].c_str(), "none")==0){
		channel=NONE;
		chanName = "None";
	}
	else if(opts["can"].length()==0){
		channel=PI0DALITZ;
		chanName = "Pi0Dalitz";
	}

	if(opts["mod"].length()!=0) outputMod = atoi(opts["mod"].c_str());
	else outputMod = 1;

	if(opts.count("debug")!=0) optDebug = true;
	else optDebug = false;

	if(opts.count("period")!=0) periodKeep = atoi(opts["period"].c_str());
	else periodKeep = 0;

	string cutsFileName;
	if(opts.count("cuts")!=0) cutsFileName = opts["cuts"];
	else cutsFileName = "";
	parseCutsValues(cutsFileName);
	printCuts();

	cout << "Starting on channel: " << chanName << endl;
	cout << "Output events every " << outputMod << " events" << endl;
	if(cutsFileName.length()>0) cout << "Using cuts at: " << cutsFileName << endl;
	cout << "Debugging activated: " << (optDebug==true ? "Yes" : "No") << endl;
	if(periodKeep==0) cout << "Keeping period: All" << endl;
	else cout << "Keeping period: " << periodKeep << endl;
	if(noOutput) cout << "No file output requested" << endl;
	cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	cout << endl << endl;
	if(channel==KE2) nico_ke2Init();
	if(channel==PI0DALITZ) nico_pi0DalitzInit();
	/*----------- End of user C code -----------*/
	return 0;
}
