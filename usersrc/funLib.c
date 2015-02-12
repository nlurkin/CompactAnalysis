/*
 * funLib.c
 *
 *  Created on: Feb 12, 2015
 *      Author: ncl
 */
#include "funLib.h"
#include <TFile.h>
#include <cmath>

#ifdef OUTSIDECOMPACT

#else
#include "user.h"
#include "reader.h"
#endif

double Mpi0 = 0.1349766;
double Mpic = 0.139570;
double Me = 0.00051099891;

std::string printVector3(TVector3 v){
	std::ostringstream ss;
	ss.precision(7);
	ss << std::fixed;
	ss << "( " << v.X() << " , " << v.Y() << " , " << v.Z() << " )";
	return ss.str();
}

std::string printSTLvector(std::vector<TVector3> v){
	std::ostringstream ss;
	std::vector<TVector3>::iterator it;

	ss << "Vector values :" << std::endl;
	for(it=v.begin(); it!=v.end();it++){
		ss << "\t" << printVector3(*it) << std::endl;
	}
	ss << std::endl;
	return ss.str();
}

TH1D* gH(TString name){
	return (TH1D*)gFile->Get(name);
}

double distance2D(TVector3 v1, TVector3 v2){
	return sqrt(pow(v1.X()-v2.X(),2)+pow(v1.Y()-v2.Y(),2));
}

double missMass2(double m1, double m2, TVector3 p1, TVector3 p2){
	return m1*m1 + m2*m2 - 2*sqrt(m1*m1+p1.Mag2())*sqrt(m2*m2+p2.Mag2()) + 2*p1.Dot(p2);
}

double missMass2Man(double m1, double m2, float kp, float kdxdz, float kdydz, float ep, float edxdz, float edydz){
	double kEnergy = sqrt(m1*m1+kp*kp);
	double eEnergy = sqrt(m2*m2+ep*ep);
	double eMag = sqrt(edxdz*edxdz + edydz*edydz + 1);
	double kMag = sqrt(kdxdz*kdxdz + kdydz*kdydz + 1);
	double dotProduct = (kdxdz*edxdz + kdydz*edydz + 1)/(eMag*kMag);
	return m1*m1 + m2*m2 - 2*(kEnergy*eEnergy - dotProduct*kp*ep);
}
double invMass2(std::vector<double> mass, std::vector<TVector3> p){
	std::vector<double> Energy;
	double ESum = 0;
	double PSum = 0;
	double MSum = 0;

	for(int i=0; i<p.size(); i++){
		Energy.push_back(sqrt(pow(mass[i],2) + p[i].Mag2()));
	}

	for(int i=0; i<p.size()-1; i++){
		for(int j=i+1; j<p.size(); j++){
			ESum += Energy[i]*Energy[j];
			PSum += p[i].Dot(p[j]);
		}
		MSum += pow(mass[i],2);
	}
	MSum += pow(mass[p.size()-1],2);

	return MSum + 2*ESum - 2*PSum;
}

void propagateBefore(float &x, float &y, float &z, float zplane, trak t){
	x = t.bx + t.bdxdz*(zplane-Geom->DCH.bz);
	y = t.by + t.bdydz*(zplane-Geom->DCH.bz);
	z = zplane;
}
TVector3 propagateBefore(float zplane, NPhysicsTrack pt){
	NTrak t = rawEvent.track[pt.trackID];
	return t.bDetPos + t.bMomentum*((zplane-t.bDetPos.Z())/t.bMomentum.Z());
}
TVector3 propagateCorrBefore(float zplane, NPhysicsTrack pt){
	NTrak t = rawEvent.track[pt.trackID];
	return t.bDetPos + pt.momentum*((zplane-t.bDetPos.Z())/pt.momentum.Z());
}
void propagateAfter(float &x, float &y, float &z, float zplane, trak t){
	x = t.x + t.dxdz*(zplane-Geom->DCH.z);
	y = t.y + t.dydz*(zplane-Geom->DCH.z);
	z = zplane;
}
TVector3 propagateAfter(float zplane, NPhysicsTrack pt){
	NTrak t = rawEvent.track[pt.trackID];
	return t.aDetPos + t.aMomentum*((zplane-t.aDetPos.Z())/t.aMomentum.Z());
}
TVector3 propagate(float zplane, NPhysicsTrack pt){
	NTrak t = rawEvent.track[pt.trackID];
	//Ou detBpos???
	//return t->middlePos + pt.momentum*((zplane-t->middlePos.Z())/pt.momentum.Z());
	//TODO FIXME
	return TVector3();
}
TVector3 propagate(float zplane, TVector3 pos, TVector3 p){
	return pos + p*((zplane-pos.Z())/p.Z());
}


void parseCutsValues(string fileName){
	//Initialize with default values
	cutsDefinition.triggerMask = 0x400;
	cutsDefinition.numVertex3 = 1;
	cutsDefinition.minZVertex = -1800;
	cutsDefinition.maxZVertex = 9000;
	cutsDefinition.maxChi2Vertex = 25;
	cutsDefinition.maxExtraTracks = 0;
	cutsDefinition.maxTrackTime = 25;
	cutsDefinition.boolBadTrack = true;
	cutsDefinition.numBadTrackCombi = 0;
	cutsDefinition.numPiCandidates = 1;
	cutsDefinition.boolBadECandidates = true;
	cutsDefinition.minTrackMomentum = 5;
	cutsDefinition.maxTrackMomentum = 74;
	cutsDefinition.numAddGoodCluster = 1;
	cutsDefinition.lkrAcceptance = 0;
	cutsDefinition.minGammaEnergy = 3;
	cutsDefinition.minDeadCellDist = 2;
	cutsDefinition.minGammaDCHRadius = 13;
	cutsDefinition.minTotalMomentum = 70;
	cutsDefinition.maxTotalMomentum = 78;
	cutsDefinition.maxPt = 0.0005;
	cutsDefinition.maxPi0MassDiff = 0.008;
	cutsDefinition.minKaonMassDiff = 0.475;
	cutsDefinition.maxKaonMassDiff = 0.510;

	if(fileName.length()==0) return;

	FILE *fdCuts = fopen(fileName.c_str(), "r");
	char buffer[200];
	char name[200], value[200];
	if(fdCuts!=NULL){
		while(fgets(buffer, sizeof(buffer), fdCuts)){
			if(buffer[0]=='#') continue;
			sscanf(buffer, "%s %s", name, value);

			if(strcmp(name,"triggerMask")==0){
				cutsDefinition.triggerMask = atoi(value);
			}
			else if(strcmp(name,"numVertex3")==0){
				cutsDefinition.numVertex3 = atoi(value);
			}
			else if(strcmp(name,"minZVertex")==0){
				cutsDefinition.minZVertex = atoi(value);
			}
			else if(strcmp(name,"maxZVertex")==0){
				cutsDefinition.maxZVertex = atoi(value);
			}
			else if(strcmp(name,"maxChi2Vertex")==0){
				cutsDefinition.maxChi2Vertex = atoi(value);
			}
			else if(strcmp(name,"maxExtraTracks")==0){
				cutsDefinition.maxExtraTracks = atoi(value);
			}
			else if(strcmp(name,"maxTrackTime")==0){
				cutsDefinition.maxTrackTime = atoi(value);
			}
			else if(strcmp(name,"boolBadTrack")==0){
				cutsDefinition.boolBadTrack = (strcmp(value,"true")==0) ? true : false;
			}
			else if(strcmp(name,"numBadTrackCombi")==0){
				cutsDefinition.numBadTrackCombi = atoi(value);
			}
			else if(strcmp(name,"numPiCandidates")==0){
				cutsDefinition.numPiCandidates = atoi(value);
			}
			else if(strcmp(name,"boolBadECandidates")==0){
				cutsDefinition.boolBadECandidates = (strcmp(value,"true")==0) ? true : false;
			}
			else if(strcmp(name,"minTrackMomentum")==0){
				cutsDefinition.minTrackMomentum = atoi(value);
			}
			else if(strcmp(name,"maxTrackMomentum")==0){
				cutsDefinition.maxTrackMomentum = atoi(value);
			}
			else if(strcmp(name,"numAddGoodCluster")==0){
				cutsDefinition.numAddGoodCluster = atoi(value);
			}
			else if(strcmp(name,"lkrAcceptance")==0){
				cutsDefinition.lkrAcceptance = atoi(value);
			}
			else if(strcmp(name,"minGammaEnergy")==0){
				cutsDefinition.minGammaEnergy = atoi(value);
			}
			else if(strcmp(name,"minDeadCellDist")==0){
				cutsDefinition.minDeadCellDist = atoi(value);
			}
			else if(strcmp(name,"minGammaDCHRadius")==0){
				cutsDefinition.minGammaDCHRadius = atoi(value);
			}
			else if(strcmp(name,"minTotalMomentum")==0){
				cutsDefinition.minTotalMomentum = atoi(value);
			}
			else if(strcmp(name,"maxTotalMomentum")==0){
				cutsDefinition.maxTotalMomentum = atoi(value);
			}
			else if(strcmp(name,"maxPt")==0){
				cutsDefinition.maxPt = atof(value);
			}
			else if(strcmp(name,"maxPi0MassDiff")==0){
				cutsDefinition.maxPi0MassDiff = atof(value);
			}
			else if(strcmp(name,"minKaonMassDiff")==0){
				cutsDefinition.minKaonMassDiff = atof(value);
			}
			else if(strcmp(name,"maxKaonMassDiff")==0){
				cutsDefinition.maxKaonMassDiff = atof(value);
			}
			else{
				cout << "No known " << name << " cut parameter " << value << endl;
			}
		}
		fclose(fdCuts);
	}
	else{
		cout << "Unable to open cuts file" << fileName << endl;
	}
}

void printCuts(){
	cout << "Cuts definition" << endl;
	cout << "---------------" << endl;
	cout << "triggerMask\t\t--> " << cutsDefinition.triggerMask << endl;
	cout << "numVertex3\t\t--> " << cutsDefinition.numVertex3 << endl;
	cout << "minZVertex\t\t--> " << cutsDefinition.minZVertex << endl;
	cout << "maxZVertex\t\t--> " << cutsDefinition.maxZVertex << endl;
	cout << "maxChi2Vertex\t\t--> " << cutsDefinition.maxChi2Vertex << endl;
	cout << "maxExtraTracks\t\t--> " << cutsDefinition.maxExtraTracks << endl;
	cout << "maxTrackTime\t\t--> " << cutsDefinition.maxTrackTime << endl;
	cout << "boolBadTrack\t\t--> " << cutsDefinition.boolBadTrack << endl;
	cout << "numBadTrackCombi\t--> " << cutsDefinition.numBadTrackCombi << endl;
	cout << "numPiCandidates\t\t--> " << cutsDefinition.numPiCandidates << endl;
	cout << "boolBadECandidates\t--> " << cutsDefinition.boolBadECandidates << endl;
	cout << "minTrackMomentum\t--> " << cutsDefinition.minTrackMomentum << endl;
	cout << "maxTrackMomentum\t--> " << cutsDefinition.maxTrackMomentum << endl;
	cout << "numAddGoodCluster\t--> " << cutsDefinition.numAddGoodCluster << endl;
	cout << "lkrAcceptance\t\t--> " << cutsDefinition.lkrAcceptance << endl;
	cout << "minGammaEnergy\t\t--> " << cutsDefinition.minGammaEnergy << endl;
	cout << "minDeadCellDist\t\t--> " << cutsDefinition.minDeadCellDist << endl;
	cout << "minGammaDCHRadius\t--> " << cutsDefinition.minGammaDCHRadius << endl;
	cout << "minTotalMomentum\t--> " << cutsDefinition.minTotalMomentum << endl;
	cout << "maxTotalMomentum\t--> " << cutsDefinition.maxTotalMomentum << endl;
	cout << "maxPt\t\t\t--> " << cutsDefinition.maxPt << endl;
	cout << "maxPi0MassDiff\t\t--> " << cutsDefinition.maxPi0MassDiff << endl;
	cout << "minKaonMassDiff\t\t--> " << cutsDefinition.minKaonMassDiff << endl;
	cout << "maxKaonMassDiff\t\t--> " << cutsDefinition.maxKaonMassDiff << endl << endl;
}


std::map<std::string,std::string> parseOptions(std::string s){
	//string str = gString;
	std::map<std::string,std::string> opts;
	std::vector<std::string> options, keys;
	std::vector<std::string>::iterator it;

	options = tokenize(s, ':');

	for(it=options.begin(); it!=options.end(); it++){
		keys = tokenize(*it, '=');
		opts.insert(std::pair<std::string,std::string>(keys[0], keys[1]));
	}

	return opts;
}

const std::vector<std::string> tokenize(std::string s, const char delim) {
	std::vector<std::string> tokens;
	std::stringstream ss(s);
	std::string item;

	while(getline(ss, item, delim)){
		tokens.push_back(item);
	}
	return tokens;
}

bool isFilteredEvent(int nrun, int nburst, int timestamp){
	std::vector<eventID>::iterator it;

	for(it=badEventsList.begin(); it!= badEventsList.end(); it++){
		if( (nrun==(*it).rnum || (*it).rnum==0)
				&& nburst==(*it).bnum
				&& timestamp==(*it).timestamp) return true;
	}
	return false;
}


int selectOptions(std::string s){
	std::map<std::string,std::string>::iterator it;

	opts = parseOptions(s);

	std::cout << std::endl << ">>>>>>>>>>>>>>>>>>>>> Initialization" << std::endl;
	if(opts.count("h")!=0){
		std::cout << ">>>> Help " << std::endl;
		std::cout << ">>>> Syntax: param=value:param=value" << std::endl;
		std::cout << ">>>> List of parameters:" << std::endl;
		std::cout << ">>>> h: This help" << std::endl;
		std::cout << ">>>> prefix: Output file names to use (without extension)" << std::endl;
		std::cout << ">>>> can: ke2 | pi0d" << std::endl;
		std::cout << ">>>> debug: Activate debugging" << std::endl;
		std::cout << ">>>> nooutput: Don't create output txt files" << std::endl;
		std::cout << ">>>> period: keep only events from this period" << std::endl;
		std::cout << ">>>> mod: print events index every mod events" << std::endl;
		std::cout << ">>>> cuts: specify cuts file" << std::endl;
		std::cout << ">>>> eall: export all event (not only the ones passing first cuts)" << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
		exit(0);
	}
	std::cout << ">>>> Received parameters:" << std::endl;
	for(it=opts.begin(); it!=opts.end(); it++){
		std::cout << "\t" << it->first << " = " << it->second << std::endl;
	}

	common_init(opts["prefix"]);
	std::string chanName;

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

	if(opts.count("eall")!=0) exportAllEvents = true;
	else exportAllEvents = false;

	std::string cutsFileName;
	if(opts.count("cuts")!=0) cutsFileName = opts["cuts"];
	else cutsFileName = "";
	parseCutsValues(cutsFileName);
	printCuts();

	std::cout << "Starting on channel: " << chanName << std::endl;
	std::cout << "Output events every " << outputMod << " events" << std::endl;
	if(cutsFileName.length()>0) std::cout << "Using cuts at: " << cutsFileName << std::endl;
	std::cout << "Debugging activated: " << (optDebug==true ? "Yes" : "No") << std::endl;
	if(periodKeep==0) std::cout << "Keeping period: All" << std::endl;
	else std::cout << "Keeping period: " << periodKeep << std::endl;
	if(noOutput) std::cout << "No file output requested" << std::endl;
	if(exportAllEvents) std::cout << "Export all events requested" << std::endl;
	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cout << std::endl << std::endl;

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
