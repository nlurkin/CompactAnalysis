/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

/*
 * TODO : Electron ID Efficiency (Ke2 MC only)
 * TODO : Pb wall
 */

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <sstream>
#include <iostream>
#include <map>
#include <sstream>
#include <iomanip>
using namespace std;

string printVector3(TVector3 v){
	ostringstream ss;
	ss.precision(7);
	ss << fixed;
	ss << "( " << v.X() << " , " << v.Y() << " , " << v.Z() << " )";
	return ss.str();
}

string printSTLvector(vector<TVector3> v){
	ostringstream ss;
	vector<TVector3>::iterator it;

	ss << "Vector values :" << endl;
	for(it=v.begin(); it!=v.end();it++){
		ss << "\t" << printVector3(*it) << endl;
	}
	ss << endl;
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
double invMass2(vector<double> mass, vector<TVector3> p){
	vector<double> Energy;
	double ESum = 0;
	double PSum = 0;
	double MSum = 0;

	for(int i=0; i<p.size(); i++){
		Energy.push_back(sqrt(pow(mass[i],2) + p[i].Mag2()));
	}

	//cout << "M " << printSTLvector<double>(mass) << endl;
	//cout << "Energy " << printSTLvector<double>(Energy) << endl;
	//cout << "P " << printSTLvector(p) << endl;

	for(int i=0; i<p.size()-1; i++){
		for(int j=i+1; j<p.size(); j++){
			ESum += Energy[i]*Energy[j];
			//cout << "E" << i << "E" << j << ": " << Energy[i]*Energy[j] << endl;
			PSum += p[i].Dot(p[j]);
			//cout << "P" << i << ".P" << j << ": " << p[i].Dot(p[j]) << endl;
		}
		MSum += pow(mass[i],2);
	}
	MSum += pow(mass[p.size()-1],2);

	//cout << "MSum " << MSum << endl;
	//cout << "ESum " << ESum << endl;
	//cout << "PSum " << PSum << endl;

	return MSum + 2*ESum - 2*PSum;
}

bool isFilteredEvent(int nrun, int nburst, int timestamp){
	vector<eventID>::iterator it;

	for(it=badEventsList.begin(); it!= badEventsList.end(); it++){
		if( (nrun==(*it).rnum || (*it).rnum==0)
				&& nburst==(*it).bnum
				&& timestamp==(*it).timestamp) return true;
	}
	return false;
}

void loadEOPData(superBurst *sbur){
	// E/p corrections for each cell
	FILE *EopCorrfile;

	int i,j;
	int cpd, cell;

	string eopCorrFileName = "/afs/cern.ch/user/n/nlurkin/Compact/eopCorrFiles/eopCorrfile_p";

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
		cout << "No eop data" << endl;
		user_exit();
		exit(0);
	}
	EopCorrfile = fopen (eopCorrFileName.c_str(), "r");

	if(EopCorrfile == NULL){
		cout << "Unable to open eop data file " << eopCorrFileName << endl;
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
	noTestCharge = false;
	beamCharge = 0;
	pbWall = false;

	//P1	+-
	if(sbur->nrun>=20114 && sbur->nrun<=20203){
		pbWall = true;
		period = 1;
	}
	//P2	+-
	if(sbur->nrun>=20209 && sbur->nrun<=20285){
		pbWall = true;
		period = 2;
	}
	//P3	+-
	if(sbur->nrun>=20286 && sbur->nrun<=20324){
		pbWall = true;
		period = 3;
	}

	//P4
	if(sbur->nrun>=20332 && sbur->nrun<=20404){
		pbWall = true;
		beamCharge = +1;
		period = 4;
	}

	//P4a	-/None
	if(sbur->nrun>=20385 && sbur->nrun<=20386){
		pbWall = true;
		period = 4;
	}
	if(sbur->nrun==20385){
		beamCharge = -1;
		period = 4;
	}

	//P5	+
	if(sbur->nrun>=20410 && sbur->nrun<=20485){
		pbWall = false;
		beamCharge = +1;
		period = 5;
	}

	//P6	-/None
	if(sbur->nrun>=20487 && sbur->nrun<=20531){
		pbWall = false;
		period = 6;
	}
	if(sbur->nrun>=20487 && sbur->nrun<=20521){
		beamCharge = -1;
		period = 6;
	}
	if(sbur->nrun>=20530 && sbur->nrun<=20531){
		beamCharge = -1;
		period = 6;
	}
	if(sbur->nrun==20525){
		beamCharge = -1;
		period = 6;
	}

	//P7	+
	if(sbur->nrun>=20611 && sbur->nrun<=20695){
		pbWall = false;
		beamCharge = +1;
		period = 7;
	}

	//P8a	+/-
	if(sbur->nrun>=21082 && sbur->nrun<=21085){
		pbWall = false;
		period = 8;
	}
	//P8
	if(sbur->nrun>=21088 && sbur->nrun<=21120){
		pbWall = false;
		period = 8;
		beamCharge = -1;
	}


	if(mcOnly){
		if(sbur->nrun>=20114 && sbur->nrun<=20133) alpha = 0;
		if(sbur->nrun>=20154 && sbur->nrun<=20256) alpha = -0.05;
		if(sbur->nrun>=20268 && sbur->nrun<=20303) alpha = -0.08;
		if(sbur->nrun>=20304 && sbur->nrun<=20324) alpha = -0.02;
		if(sbur->nrun>=20332 && sbur->nrun<=20371) alpha = 0.08;
		if(sbur->nrun>=20387 && sbur->nrun<=20486) alpha = 0.07;
		if(sbur->nrun>=20613 && sbur->nrun<=20695) alpha = 0.10;
	}
	if(beamCharge==0) noTestCharge=true;
}

void propagateBefore(float &x, float &y, float &z, float zplane, trak t){
	x = t.bx + t.bdxdz*(zplane-Geom->DCH.bz);
	y = t.by + t.bdydz*(zplane-Geom->DCH.bz);
	z = zplane;
}
TVector3 propagateBefore(float zplane, CorrectedTrack *t){
	return t->detBPos + t->unBMomentum*((zplane-t->detBPos.Z())/t->unBMomentum.Z());
}
TVector3 propagateCorrBefore(float zplane, CorrectedTrack *t){
	return t->detBPos + t->momentum*((zplane-t->detBPos.Z())/t->momentum.Z());
}
void propagateAfter(float &x, float &y, float &z, float zplane, trak t){
	x = t.x + t.dxdz*(zplane-Geom->DCH.z);
	y = t.y + t.dydz*(zplane-Geom->DCH.z);
	z = zplane;
}
TVector3 propagateAfter(float zplane, CorrectedTrack *t){
	return t->detPos + t->unMomentum*((zplane-t->detPos.Z())/t->unMomentum.Z());
}
TVector3 propagate(float zplane, CorrectedTrack *t){
	//Ou detBpos???
	return t->middlePos + t->momentum*((zplane-t->middlePos.Z())/t->momentum.Z());
}
TVector3 propagate(float zplane, TVector3 pos, TVector3 p){
	return pos + p*((zplane-pos.Z())/p.Z());
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
float correctClusterE(superCmpEvent *sevt, CorrectedCluster *c){
	int i, j, k;
	int l, m, n;

	float x,y;

	x = c->position.X();
	y = c->position.Y();

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

	e_cluster = sevt->cluster[c->clusterID].energy / EopCorr[CPDindex][CELLindex];  // Ke3 E/p correction for each cell

	return e_cluster;
}

CorrectedTrack *correctTrack(superCmpEvent *sevt, int i){
	float eOverP;
	trak t = sevt->track[i];
	double pTrack = p_corr_ab(t.p,t.q);
	double eTrack = correctedEP(sevt, t, eOverP);

	CorrectedTrack *x = new CorrectedTrack;
	x->trackID = i;
	//x.energy pas encore
	//x.momentum pas encore	-> done
	x->pMag = pTrack;
	//x.vertex pas encore -> done
	//x.cda pas encore -. done
	//x.middlePos pas encore -> done
	x->detBPos = TVector3(t.bx, t.by, Geom->DCH.bz);
	x->detPos = TVector3(t.x, t.y, Geom->DCH.z);
	//x->t = t;

	if(t.iClus>=0){
		if(dataOnly){
			x->clusterEnergy = eTrack;
		}
		else{
			x->clusterEnergy = sevt->cluster[t.iClus].energy;
		}
	}
	else{
		x->clusterEnergy = 0;
	}
	x->unMomentum = TVector3(t.dxdz, t.dydz, 1).Unit();
	//cout << "ddz " << t.dxdz << " " << t.dydz << endl;
	//cout << "unMomentum " << printVector3(TVector3(t.dxdz, t.dydz, 1)) << endl;
	//cout << "unMomentum.Unit " << printVector3(TVector3(t.dxdz, t.dydz, 1).Unit()) << endl;
	//cout << "unMomentum " << printVector3(x->unMomentum) << endl;
	x->unBMomentum = TVector3(t.bdxdz, t.bdydz, 1).Unit();
	x->unPMag = t.p;
	//x.V0 pas encore -> done
	//x.unCda pas encore -> done

	return x;
}
CorrectedCluster *correctCluster(superCmpEvent* sevt, int i){
	cluster c = sevt->cluster[i];
	CorrectedCluster *x = new CorrectedCluster;
	double zsh, x1, y1;

	x->clusterID = i;
	//x->c = c;
	x->unPosition = TVector3(c.x, c.y, Geom->Lkr.z);

	if(dataOnly){
		x->position.SetZ(Geom->Lkr.z + 16.5 + 4.3*log(c.energy));
		//x->position.SetZ(Geom->Lkr.z);
		x->position.SetX((x->unPosition.X() + 0.136 + 0.87e-3*x->unPosition.Y()) * (1+(x->position.Z()-x->unPosition.Z())/10998.));
		x->position.SetY((x->unPosition.Y() + 0.300 - 0.87e-3*x->unPosition.X()) * (1+(x->position.Z()-x->unPosition.Z())/10998.));

		x->energy = correctClusterE(sevt, x);
	}
	else if(mcOnly){
		x->position.SetZ(Geom->Lkr.z + 16.5 + 4.3*log(c.energy));
		x->position.SetX((x->unPosition.X() - 0.013) * (1+(x->position.Z()-x->unPosition.Z())/10998.));
		x->position.SetY(x->unPosition.Y() * (1+(x->position.Z()-x->unPosition.Z())/10998.));

		x->energy = c.energy;
	}
	return x;
}

vector<CorrectedTrack*> CreateTracks(superCmpEvent *sevt){
	vector<CorrectedTrack*> vtracks;
	CorrectedTrack *t;
	for(int i=0; i<sevt->Ntrack; i++){
		t = correctTrack(sevt, i);
		vtracks.push_back(t);
	}
	return vtracks;
}
vector<CorrectedCluster*> CreateClusters(superCmpEvent *sevt){
	vector<CorrectedCluster*> vClusters;
	CorrectedCluster *t;
	for(int i=0; i<sevt->Ncluster; i++){
		t = correctCluster(sevt, i);
		vClusters.push_back(t);
	}
	return vClusters;
}

const vector<string> tokenize(string s, const char delim) {
	vector<string> tokens;
	stringstream ss(s);
	string item;

	while(getline(ss, item, delim)){
		tokens.push_back(item);
	}
	return tokens;
}
map<string,string> parseOptions(){
	//string str = gString;
	map<string,string> opts;
	vector<string> options, keys;
	vector<string>::iterator it;

	options = tokenize(gString, ':');

	for(it=options.begin(); it!=options.end(); it++){
		keys = tokenize(*it, '=');
		opts.insert(pair<string,string>(keys[0], keys[1]));
	}

	return opts;
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


int nico_ke2Init(){
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	goodTracks.clear();
	assocClusters.clear();
	outTree = new TTree("event", "Event");
	outTree->Branch("goodTrack", "std::vector<CorrectedTrack>", &goodTracks, 64000, 1);
	outTree->Branch("assocCluster", "std::vector<CorrectedCluster>", &assocClusters, 64000, 1);

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

	goodTracks.clear();
	assocClusters.clear();
	outTree = new TTree("event", "Event");
	outTree->Branch("goodTrack", "std::vector<CorrectedTrack*>", &goodTracks, 64000, 1);
	//outTree->Branch("assocCluster", "std::vector<CorrectedCluster*>", &assocClusters, 64000, 1);
	//outTree->Branch("cutsWord", &cutsWord, "cutsWord[19]/O");

	outTree->Branch("pi0dEvent" ,"pi0dEvent", &fullEvent);

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

	opts = parseOptions();

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
