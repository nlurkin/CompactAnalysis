/*
 * compactLib.c
 *
 *  Created on: Feb 12, 2015
 *      Author: ncl
 */

#include "compactLib.h"
#include <string>
#include "funLib.h"
#include <TFile.h>


void loadEOPData(superBurst *sbur){
	// E/p corrections for each cell
	FILE *EopCorrfile;

	int i,j;
	int cpd, cell;

	std::string eopCorrFileName = "/afs/cern.ch/user/n/nlurkin/Compact/eopCorrFiles/eopCorrfile_p";

	if(sbur->nrun<20209){
		eopCorrFileName+= "1_v63.dat";
	}
	else if(sbur->nrun>=20209 && sbur->nrun<=20238){
		eopCorrFileName+= "2_20209-20238_v63.dat";
	}
	else if(sbur->nrun>=20254 && sbur->nrun<=20285){
		eopCorrFileName+= "2_20254-20285_v63.dat";
	}
	else if(sbur->nrun>=20286 && sbur->nrun<=20291){
		eopCorrFileName+= "3_20286-20291_v63.dat";
	}
	else if(sbur->nrun>=20296 && sbur->nrun<=20324){
		eopCorrFileName+= "3_20296-20324_v63.dat";
	}
	else if(sbur->nrun>=20332 && sbur->nrun<=20404){
		eopCorrFileName+= "4_v63.dat";
	}
	else if(sbur->nrun>=20410 && sbur->nrun<=20415){
		eopCorrFileName+= "5_20410-20415_v63.dat";
	}
	else if(sbur->nrun>=20416 && sbur->nrun<=20485){
		eopCorrFileName+= "5_20416-20485_v63.dat";
	}
	else if(sbur->nrun>=20487 && sbur->nrun<=20531){
		eopCorrFileName+= "6_v63.dat";
	}
	else if(sbur->nrun>=20611 && sbur->nrun<=20695){
		eopCorrFileName+= "7_v63.dat";
	}
	else if(sbur->nrun>=21082 && sbur->nrun<=21103){
		eopCorrFileName+= "8_21082-21103_v63.dat";
	}
	else if(sbur->nrun>=21106 && sbur->nrun<=21120){
		eopCorrFileName+= "8_21106-21120_v63.dat";
	}
	else{
		std::cout << "No eop data" << std::endl;
		user_exit();
		exit(0);
	}
	EopCorrfile = fopen (eopCorrFileName.c_str(), "r");

	if(EopCorrfile == NULL){
		std::cout << "Unable to open eop data file " << eopCorrFileName << std::endl;
		user_exit();
		exit(0);
	}
	for (i=0; i<256; i++)
		for (j=0; j<64; j++)
		{
			fscanf (EopCorrfile, "%i %i %f\n", &cpd, &cell, &EopCorr[i][j]);
		}

	fclose(EopCorrfile);
	eopLoaded = true;
}

void defineBeamCharge(superBurst *sbur){
	//P1	+-
	if(sbur->nrun>=20114 && sbur->nrun<=20203){
		rootBurst.pbWall = true;
		rootBurst.period = 1;
	}
	//P2	+-
	if(sbur->nrun>=20209 && sbur->nrun<=20285){
		rootBurst.pbWall = true;
		rootBurst.period = 2;
	}
	//P3	+-
	if(sbur->nrun>=20286 && sbur->nrun<=20303){
		rootBurst.pbWall = true;
		rootBurst.period = 3;
	}
	if(sbur->nrun>=20304 && sbur->nrun<=20324){
		rootBurst.pbWall = true;
		rootBurst.beamCharge = +1;
		rootBurst.period = 3;
	}

	//P4
	if(sbur->nrun>=20332 && sbur->nrun<=20404){
		rootBurst.pbWall = true;
		rootBurst.beamCharge = +1;
		rootBurst.period = 4;
	}

	//P4a	-/None
	if(sbur->nrun>=20385 && sbur->nrun<=20386){
		rootBurst.pbWall = true;
		rootBurst.period = 4;
	}
	if(sbur->nrun==20385){
		rootBurst.beamCharge = -1;
		rootBurst.period = 4;
	}

	//P5	+
	if(sbur->nrun>=20410 && sbur->nrun<=20485){
		rootBurst.pbWall = false;
		rootBurst.beamCharge = +1;
		rootBurst.period = 5;
	}

	//P6	-/None
	if(sbur->nrun>=20487 && sbur->nrun<=20531){
		rootBurst.pbWall = false;
		rootBurst.period = 6;
	}
	if(sbur->nrun>=20487 && sbur->nrun<=20521){
		rootBurst.beamCharge = -1;
		rootBurst.period = 6;
	}
	if(sbur->nrun>=20530 && sbur->nrun<=20531){
		rootBurst.beamCharge = -1;
		rootBurst.period = 6;
	}
	if(sbur->nrun==20525){
		rootBurst.beamCharge = -1;
		rootBurst.period = 6;
	}

	//P7	+
	if(sbur->nrun>=20611 && sbur->nrun<=20695){
		rootBurst.pbWall = false;
		rootBurst.beamCharge = +1;
		rootBurst.period = 7;
	}

	//P8a	+/-
	if(sbur->nrun>=21082 && sbur->nrun<=21085){
		rootBurst.pbWall = false;
		rootBurst.period = 8;
	}
	//P8
	if(sbur->nrun>=21088 && sbur->nrun<=21120){
		rootBurst.pbWall = false;
		rootBurst.period = 8;
		rootBurst.beamCharge = -1;
	}


	if(rootBurst.isMC){
		if(sbur->nrun>=20114 && sbur->nrun<=20133) rootBurst.alpha = 0;
		if(sbur->nrun>=20154 && sbur->nrun<=20256) rootBurst.alpha = -0.05;
		if(sbur->nrun>=20268 && sbur->nrun<=20303) rootBurst.alpha = -0.08;
		if(sbur->nrun>=20304 && sbur->nrun<=20324) rootBurst.alpha = -0.02;
		if(sbur->nrun>=20332 && sbur->nrun<=20371) rootBurst.alpha = 0.08;
		if(sbur->nrun>=20387 && sbur->nrun<=20486) rootBurst.alpha = 0.07;
		if(sbur->nrun>=20613 && sbur->nrun<=20695) rootBurst.alpha = 0.10;
	}
}

void GetCpdCellIndex(double pos_x, double pos_y, int *cpd_index, int *cell_index)
{
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                              //
	// For a given x,y coordinate at the LKr front face, this routine returns the index of the CPD  //
	// and Cell at that position (official numbering). The information is needed e.g. for the       //
	// appropriate energy recalibration per cell and the electron ID correction for Ke2 MC which    //
	// can be different for single bad CPDs.                                                        //
	//                                                                                              //
	// To reduce computing time, the definition of the corner positions of CPDs and Cells has to be //
	// done before at the beginning of the job (e.g. in user_init.c) and then stored in memory.     //
	// An example how to do it can be found on the Ke2 analysis page, section Cell-by-cell E/p      //
	// recalibration.                                                                               //
	//                                                                                              //
	// A. Winhart  -  10.11.2009                                                                    //
	//                                                                                              //
	//////////////////////////////////////////////////////////////////////////////////////////////////


	int   i, j, k;
	int   l, m, n;
	double CELLlength = 1.975;             // Cell length in cm
	double CPDlength = 8 * CELLlength;     // One CPD consists of 8x8 cells

	*cpd_index  = -1;
	*cell_index = -1;

	// Define CPD index
	for (i=0; i<16; i++)
		for (j=0; j<16; j++)
		{
			k = i*16 + j;
			if (pos_x <= CPDpos_leftDownCorner[k][0] && pos_x > (CPDpos_leftDownCorner[k][0]-CPDlength) &&
					pos_y >= CPDpos_leftDownCorner[k][1] && pos_y < (CPDpos_leftDownCorner[k][1]+CPDlength) )
			{
				//printf ("CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
				//  k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1], pos_x, pos_y);
				*cpd_index = k;
				break;
			}
		}

	// Define Cell index if CPD has been found
	if (*cpd_index >= 0)
		for (m=0; m<8; m++)
			for (n=0; n<8; n++)
			{
				l = m*8 + n;
				//printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
				//      l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1]);

				if (pos_x <= CELLpos_leftDownCorner[*cpd_index][l][0] && pos_x > (CELLpos_leftDownCorner[*cpd_index][l][0]-CELLlength) &&
						pos_y >= CELLpos_leftDownCorner[*cpd_index][l][1] && pos_y < (CELLpos_leftDownCorner[*cpd_index][l][1]+CELLlength) )
				{
					//printf ("Cell %d in CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
					//  l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1], pos_x, pos_y);
					*cell_index = l;
				}
			}
	//printf ("cpd_index = %3d \t cell_index = %3d\n", *cpd_index, *cell_index);

}
float correctedEP(superCmpEvent* sevt, trak t, float &eOverP){
	int i, j, k;
	int l, m, n;

	float x,y,z;

	propagateAfter(x,y,z,Geom->Lkr.z,t);
	// First find out to which cell is pointing the deflected track (define CPDindex and CELLindex)
	// Of course, before you had to extrapolate the track to the LKR to get the track coordinates there
	// trkatlkr[0] = x-position, trkatlkr[1] = y-position

	for (i=0; i<16; i++)
		for (j=0; j<16; j++)
		{
			k = i*16 + j;
			if (x <= CPDpos_leftDownCorner[k][0] && x > (CPDpos_leftDownCorner[k][0]-CPDlength) &&
					y >= CPDpos_leftDownCorner[k][1] && y < (CPDpos_leftDownCorner[k][1]+CPDlength) )
			{
				//printf ("CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
				//  k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1], trkatlkr[0], trkatlkr[1]);
				CPDindex = k;
				break;
			}
		}

	for (m=0; m<8; m++)
		for (n=0; n<8; n++)
		{
			l = m*8 + n;
			//printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
			//      l, CPDindex, CELLpos_leftDownCorner[CPDindex][l][0], CELLpos_leftDownCorner[CPDindex][l][1]);

			if (x <= CELLpos_leftDownCorner[CPDindex][l][0] && x > (CELLpos_leftDownCorner[CPDindex][l][0]-CELLlength) &&
					y >= CELLpos_leftDownCorner[CPDindex][l][1] && y < (CELLpos_leftDownCorner[CPDindex][l][1]+CELLlength) )
			{
				//printf ("Cell %d in CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
				//  l, CPDindex, CELLpos_leftDownCorner[CPDindex][l][0], CELLpos_leftDownCorner[CPDindex][l][1], trkatlkr[0], trkatlkr[1]);
				CELLindex = l;
			}
		}

	// Now that you know the cell hit by the track, correct the energy for the track cluster (is there is one associated to the track)
	double e_track;          // momentum and energy of the track
	int    clui_track;       // cluster index associated to the track

	clui_track = t.iClus;
	if (clui_track >= 0)
		e_track    = sevt->cluster[clui_track].energy / EopCorr[CPDindex][CELLindex];  // Ke3 E/p correction for each cell

	eOverP = EopCorr[CPDindex][CELLindex];

	return e_track;
}
float correctClusterE(superCmpEvent *sevt, NPhysicsCluster c){
	int i, j, k;
	int l, m, n;

	float x,y;

	x = c.position.X();
	y = c.position.Y();

	// First find out to which cell is pointing the deflected track (define CPDindex and CELLindex)
	// Of course, before you had to extrapolate the track to the LKR to get the track coordinates there
	// trkatlkr[0] = x-position, trkatlkr[1] = y-position

	for (i=0; i<16; i++)
		for (j=0; j<16; j++)
		{
			k = i*16 + j;
			if (x <= CPDpos_leftDownCorner[k][0] && x > (CPDpos_leftDownCorner[k][0]-CPDlength) &&
					y >= CPDpos_leftDownCorner[k][1] && y < (CPDpos_leftDownCorner[k][1]+CPDlength) )
			{
				//printf ("CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
				//  k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1], trkatlkr[0], trkatlkr[1]);
				CPDindex = k;
				break;
			}
		}

	for (m=0; m<8; m++)
		for (n=0; n<8; n++)
		{
			l = m*8 + n;
			//printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
			//      l, CPDindex, CELLpos_leftDownCorner[CPDindex][l][0], CELLpos_leftDownCorner[CPDindex][l][1]);

			if (x <= CELLpos_leftDownCorner[CPDindex][l][0] && x > (CELLpos_leftDownCorner[CPDindex][l][0]-CELLlength) &&
					y >= CELLpos_leftDownCorner[CPDindex][l][1] && y < (CELLpos_leftDownCorner[CPDindex][l][1]+CELLlength) )
			{
				//printf ("Cell %d in CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
				//  l, CPDindex, CELLpos_leftDownCorner[CPDindex][l][0], CELLpos_leftDownCorner[CPDindex][l][1], trkatlkr[0], trkatlkr[1]);
				CELLindex = l;
			}
		}

	// Now that you know the cell hit by the track, correct the energy for the track cluster (is there is one associated to the track)
	double e_cluster;          // momentum and energy of the track

	e_cluster = sevt->cluster[c.clusterID].energy / EopCorr[CPDindex][CELLindex];  // Ke3 E/p correction for each cell

	return e_cluster;
}

NPhysicsTrack correctTrack(superCmpEvent *sevt, int i){
	float eOverP;
	trak t = sevt->track[i];
	double pTrack = p_corr_ab(t.p,t.q);
	double eTrack = correctedEP(sevt, t, eOverP);

	NPhysicsTrack x;

	x.trackID = i;
	x.p = pTrack;

	if(t.iClus>=0){
		x.clusterID = t.iClus;
		if(rootBurst.isData) x.E = eTrack;
		else x.E = sevt->cluster[t.iClus].energy;
	}
	else{
		x.clusterID = -1;
		x.E = 0;
	}

	return x;
}
NPhysicsCluster correctCluster(superCmpEvent* sevt, int i){
	cluster c = sevt->cluster[i];

	NPhysicsCluster x;

	double zsh, x1, y1;

	x.clusterID = i;

	if(rootBurst.isData){
		x.position.SetZ(Geom->Lkr.z + 16.5 + 4.3*log(c.energy));
		//x->position.SetZ(Geom->Lkr.z);
		x.position.SetX((c.x + 0.136 + 0.87e-3*c.y) * (1+(x.position.Z()-Geom->Lkr.z)/10998.));
		x.position.SetY((c.y + 0.300 - 0.87e-3*c.x) * (1+(x.position.Z()-Geom->Lkr.z)/10998.));

		x.E = correctClusterE(sevt, x);
	}
	else if(rootBurst.isMC){
		x.position.SetZ(Geom->Lkr.z + 16.5 + 4.3*log(c.energy));
		x.position.SetX((c.x - 0.013) * (1+(x.position.Z()-Geom->Lkr.z)/10998.));
		x.position.SetY(c.y * (1+(x.position.Z()-Geom->Lkr.z)/10998.));

		x.E = c.energy;
	}

	return x;
}

void CreateTracks(superCmpEvent *sevt){
	for(int i=0; i<sevt->Ntrack; i++){
		NTrak x;
		x.p = sevt->track[i].p;
		x.q = sevt->track[i].q;
		x.quality = sevt->track[i].quality;
		x.chi2 = sevt->track[i].chi2;
		x.by = sevt->track[i].by;
		x.bx = sevt->track[i].bx;
		x.bdxdz = sevt->track[i].bdxdz;
		x.bdydz = sevt->track[i].bdydz;
		x.bDetPos = TVector3(sevt->track[i].bx, sevt->track[i].by, Geom->DCH.bz);
		x.aDetPos = TVector3(sevt->track[i].x, sevt->track[i].y, Geom->DCH.z);
		x.bMomentum = TVector3(x.bdxdz, x.bdydz, 1).Unit();
		x.aMomentum = TVector3(sevt->track[i].dxdz, sevt->track[i].dydz, 1).Unit();
		x.time = sevt->track[i].time;
		x.dDeadCell = sevt->track[i].dDeadCell;

		NPhysicsTrack t = correctTrack(sevt, i);

		rawEvent.track.push_back(x);

		TVector3 propPos = propagateAfter(Geom->Lkr.z, t, rawEvent);
		t.lkr_acc = LKr_acc(rootBurst.nrun, propPos.X(), propPos.Y(), 8);

		corrEvent.pTrack.push_back(t);
	}
}
void CreateClusters(superCmpEvent *sevt){
	for(int i=0; i<sevt->Ncluster; i++){
		NCluster x;
		x.position = TVector3(sevt->cluster[i].x, sevt->cluster[i].y, Geom->Lkr.z);
		x.E = sevt->cluster[i].energy;
		x.dDeadCell = sevt->cluster[i].dDeadCell;
		x.time = sevt->cluster[i].time;
		x.lkr_acc = LKr_acc(rootBurst.nrun, sevt->cluster[i].x, sevt->cluster[i].y, 8);
		x.iTrack = sevt->cluster[i].iTrack;

		NPhysicsCluster t = correctCluster(sevt, i);

		rawEvent.cluster.push_back(x);
		corrEvent.pCluster.push_back(t);
	}
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

	if(opts["nooutput"].length()!=0) noOutput = true;
	else noOutput = false;
	common_init(opts["prefix"], opts["filter"], badEventsList);
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
	parseCutsValues(cutsFileName, cutsDefinition);
	printCuts(cutsDefinition);

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


void applyEOPData(){
	int i, j, k;
	int l, m, n;
	int cpd, cell;

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
}


void openOutput(string outFile, string outPass){
	if(!noOutput){
		fprt=fopen(outFile.c_str(),"w");
		fprt2=fopen(outPass.c_str(),"w");
	}
}

int common_init(std::string filePrefix, std::string filterFile, vector<eventID> &badEventsList){
	int runNum, burstNum, timestamp;
	vector<eventID>::iterator it;

	std::string outRoot = "outfile.root";
	std::string outFile = "compact.txt";
	std::string outPass = "compactpass.txt";
	if(filePrefix.find('~')!=std::string::npos) filePrefix=filePrefix.replace(filePrefix.find('~'), 1, string("/afs/cern.ch/user/n/nlurkin"));
	if(filePrefix.length()>0){
		outRoot = filePrefix + ".root";
		outFile = filePrefix + ".txt";
		outPass = filePrefix + "pass.txt";
	}

	applyEOPData();
	openOutput(outFile, outPass);

	//Load badEvents list for debugging
	if(!filterFile.empty()){
		std::cout << ">>>> Filtering events from file " << filterFile << std::endl;
		FILE *badEvents = fopen(filterFile.c_str(), "r");
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

	gFile = TFile::Open(outRoot.c_str(), "RECREATE");

	return 0;
}
