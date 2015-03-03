#define OUTSIDECOMPACT

/// Compact includes
#include "funLib.h"

/// Std include
#include <iomanip>
using namespace std;

// Local includes
#include "CompactIO.h"
#include "OptionsParser.h"

#include <TH2D.h>
#include <TTree.h>

CompactIO io;
OptionsParser options;

///### TTree objects
ROOTRawEvent &rawEvent = io.getRawEvent();
ROOTCorrectedEvent &corrEvent = io.getCorrEvent();
ROOTBurst &rootBurst = io.getRootBurst();
ROOTFileHeader &rootFileHeader = io.getInputFileHeader();
NGeom &rootGeom = io.getRootGeom();
ROOTMCEvent &rootMC = io.getRootMc();
ROOTPhysicsEvent &rootPhysics = io.getRootPhysics();
ROOTFileHeader &outputFileHeader = io.getOutputFileHeader();

// ### Database objects
NAbcog_params abcog_params;

TH1D meeTrue = TH1D("meeTrue", "meeTrue", 100, 0, 1);
TH1D mkTrue = TH1D("mkTrue", "mkTrue", 100, 0, 1);
TH2D meeepiTrue = TH2D("meeepiTrue", "meeepiTrue", 100, 0, 1, 100, 0, 1);
TH2D meekTrue = TH2D("meekTrue", "meekTrue", 100, 0, 1, 100, 0, 1);

TH1D meeFalse = TH1D("meeFalse", "meeFalse", 100, 0, 1);
TH1D mkFalse = TH1D("mkFalse", "mkFalse", 100, 0, 1);
TH2D meeepiFalse = TH2D("meeepiFalse", "meeepiFalse", 100, 0, 1, 100, 0, 1);
TH2D meekFalse = TH2D("meekFalse", "meekFalse", 100, 0, 1, 100, 0, 1);

TH1D meeDiffTrue = TH1D("meeDiffTrue", "meeDiffTrue", 1000, -1, 1);
TH1D mkDiffTrue = TH1D("mkDiffTrue", "mkDiffTrue", 1000, -1, 1);

TH1D meeDiffFalse = TH1D("meeDiffFalse", "meeDiffFalse", 1000, -1, 1);
TH1D mkDiffFalse = TH1D("mkDiffFalse", "mkDiffFalse", 1000, -1, 1);


TH1D nPiCandidates = TH1D("nPiCandidates", "nPiCandidates", 20, 0, 20);
TH1D nPiCandidatesNew = TH1D("nPiCandidatesNew", "nPiCandidatesNew", 20, 0, 20);

TH1D nBeamSign = TH1D("nBeamSign", "nBeamSign", 20, 0, 20);

TH2D xreco_xtrue = TH2D("xreco_xtrue", "xreco_xtrue", 1000, 0, 1, 1000, 0, 1);

struct alt_pid_res{
	//MC association
	int good=0, bad=0;
	int compareMC=0;

	//Totals
	int totEvents = 0;
	int prelimEvents = 0;

	//New
	int nNoId = 0;
	int nGoodId = 0;
	int nWrongId = 0;
	int nManyId = 0;
	int nIded = 0;
	int nPassGood = 0;
	int nPassWrong = 0;
	int nPassIded = 0;
} pid_res;

int em, ep, pic;
bool flBad;

bool associateMCTracks(){
	double limit = 2e-4;
	bool flBad = false;

	int ep1=-1, ep2=-1, ep3=-1;
	int em1=-1, em2=-1, em3=-1;
	//Find out the real tracks
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()<limit) ep1=0;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()<limit) ep2=1;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()<limit) ep3=2;

	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()<limit) em1=0;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()<limit) em2=1;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()<limit) em3=2;

	if((ep1==-1 && ep2==-1 && ep3!=-1) || (ep1==-1 && ep2!=-1 && ep3==-1) || (ep1!=-1 && ep2==-1 && ep3==-1)){
		//OK
		if(ep1!=-1) ep=ep1;
		if(ep2!=-1) ep=ep2;
		if(ep3!=-1) ep=ep3;
		pid_res.good++;
	}
	else{
		flBad = true;
		pid_res.bad++;
	}

	if((em1==-1 && em2==-1 && em3!=-1) || (em1==-1 && em2!=-1 && em3==-1) || (em1!=-1 && em2==-1 && em3==-1)){
		//OK
		if(em1!=-1) em=em1;
		if(em2!=-1) em=em2;
		if(em3!=-1) em=em3;
		pid_res.good++;
	}
	else{
		pid_res.bad++;
		flBad = true;
	}

	if(!flBad){
		if((ep==0 && em==1) || (ep==1 && em==0)) pic=2;
		else if((ep==0 && em==2) || (ep==2 && em==0)) pic=1;
		else if((ep==1 && em==2) || (ep==2 && em==1)) pic=0;
		pid_res.compareMC++;
	}
	else{
		ep = -1;
		em = -1;
		pic = -1;
	}

	return flBad;
}

int pi0d_identifyPi(int &piCandidate, bool &badElectron){
	double eop;

	int piCandidatesNb = 0;
	badElectron = false;

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		int goodTrackID = corrEvent.goodTracks[i];
		if(options.isOptDebug()) cout << "\tTrying track :\t\t" << i <<  endl;

		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		if(options.isOptDebug()) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
		if(eop<0.85){
			piCandidatesNb++;
			piCandidate = i;
		}
		else if(eop>1.15) badElectron = true;
	}

	return piCandidatesNb;
}

int pid(int &piCandidate, TLorentzVector &gamma){
	TLorentzVector tem;
	TLorentzVector t1ep, t1pi;
	TLorentzVector t2ep, t2pi;
	TLorentzVector ee1, ee2;
	TLorentzVector k1, k2;

	int goodTrack1, goodTrack2;

	int nCandidates = 0;

	int nNegative = 0;
	int vtxCharge = rawEvent.vtx[corrEvent.goodVertexID].charge;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].q==-1*vtxCharge){
		goodTrack1 = 1;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 1;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p, Me);
		nNegative++;
	}

	nBeamSign.Fill(nNegative);
	if(nNegative!=1) return 0;

	//Try e = goodTrack1, pi = goodTrack2
	//cout << endl << corrEvent.pTrack.size() << " " << corrEvent.goodTracks.size() << " " << goodTrack1 << " " << goodTrack2 << endl;
	t1ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Me);
	t1pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Mpic);
	ee1 = tem+t1ep+gamma;
	k1 = ee1+t1pi;

	//Try e = goodTrack2, pi = goodTrack1
	t2pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Mpic);
	t2ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Me);
	ee2 = tem+t2ep+gamma;
	k2 = ee2+t2pi;

	double Mk;
	double diffpi01, diffpi02;
	double diffk1, diffk2;
	double diff1, diff2;

	diffpi01 = fabs(ee1.M()-Mpi0);
	diffpi02 = fabs(ee2.M()-Mpi0);

	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1) Mk = abcog_params.mkp;
	else Mk = abcog_params.mkn;

	diffk1 = fabs(k1.M()-Mk);
	diffk2 = fabs(k2.M()-Mk);

	diffk1 = 0;
	diffk2 = 0;

	diff1 = diffpi01;// + diffk1;
	diff2 = diffpi02;// + diffk2;

	double pi0DiffLimit = io.cutsDefinition.pid_pi0Diff;
	double kDiffLimit = io.cutsDefinition.pid_kDiff;
	//take the smallest ee mass as the ep em, the other is pi
	if(ee1.M()<(Mpi0+pi0DiffLimit) && k1.M()<(Mk+kDiffLimit)){
		nCandidates++;
		piCandidate = goodTrack2;
	}
	if(ee2.M()<(Mpi0+pi0DiffLimit) && k2.M()<(Mk+kDiffLimit)){
		nCandidates++;
		piCandidate = goodTrack1;
	}

	/*if(diff1<diff2 && diff1<diff3) piCandidate = 2;
	else if(diff2<diff1 && diff2<diff3) piCandidate = 1;
	else if(diff3<diff2 && diff3<diff1) piCandidate = 0;*/

	if(!flBad){
		//Good MC association, fill the plots
		if(pic==goodTrack2){
			meeTrue.Fill(ee1.M());
			mkTrue.Fill(k1.M());
			meeepiTrue.Fill(ee1.M(), (tem+t1pi).M());
			meekTrue.Fill(ee1.M(), k1.M());

			meeDiffTrue.Fill(ee1.M()-Mpi0);
			mkDiffTrue.Fill(k1.M()-Mk);

			meeFalse.Fill(ee2.M());
			mkFalse.Fill(k2.M());
			meeepiFalse.Fill(ee2.M(), (tem+t2pi).M());
			meekFalse.Fill(ee2.M(), k2.M());

			meeDiffFalse.Fill(ee2.M()-Mpi0);
			mkDiffFalse.Fill(k2.M()-Mk);
		}
		else if(pic==goodTrack1){
			meeTrue.Fill(ee2.M());
			mkTrue.Fill(k2.M());
			meeepiTrue.Fill(ee2.M(), (tem+t2pi).M());
			meekTrue.Fill(ee2.M(), k2.M());

			meeDiffTrue.Fill(ee2.M()-Mpi0);
			mkDiffTrue.Fill(k2.M()-Mk);

			meeFalse.Fill(ee1.M());
			mkFalse.Fill(k1.M());
			meeepiFalse.Fill(ee1.M(), (tem+t1pi).M());
			meekFalse.Fill(ee1.M(), k1.M());

			meeDiffFalse.Fill(ee1.M()-Mpi0);
			mkDiffFalse.Fill(k1.M()-Mk);
		}
	}

	return nCandidates;
}

int pi0d_tracksAcceptance(){
	bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(rootGeom.Dch[0].PosChamber.x,rootGeom.Dch[0].PosChamber.y,rootGeom.Dch[0].PosChamber.z);
	TVector3 dch2(rootGeom.Dch[1].PosChamber.x,rootGeom.Dch[1].PosChamber.y,rootGeom.Dch[1].PosChamber.z);
	TVector3 dch4(rootGeom.Dch[3].PosChamber.x,rootGeom.Dch[3].PosChamber.y,rootGeom.Dch[3].PosChamber.z);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		propPos = propagateAfter(rootGeom.Lkr.z, t);
		lkrAcceptance = t.lkr_acc;
		if(options.isOptDebug()) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}

		propPos = propagateCorrBefore(rootGeom.Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		if(options.isOptDebug()) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(rootGeom.Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		if(options.isOptDebug()) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(rootGeom.Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		if(options.isOptDebug()) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;
	}

	return badTrack;
}

int pi0d_trackCombinationVeto(){
	int ntracks = 3;

	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			int trackID1 = corrEvent.goodTracks[i];
			int trackID2 = corrEvent.goodTracks[j];
			NPhysicsTrack t1 = corrEvent.pTrack[trackID1];
			NPhysicsTrack t2 = corrEvent.pTrack[trackID2];
			bad = false;
			if(options.isOptDebug()) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(options.isOptDebug()) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(rootBurst.isData){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				if(options.isOptDebug()) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(rootGeom.Lkr.z, t1);
			propPos2 = propagateAfter(rootGeom.Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(trackID1==rootPhysics.pic.parentTrack || trackID2==rootPhysics.pic.parentTrack){
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
				if(RLKr<=20) bad = true;
			}

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_goodClusters_loose(){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) conditions = 4;
	else conditions=3;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack t1 = corrEvent.pTrack[corrEvent.goodTracks[0]];
		NPhysicsTrack t2 = corrEvent.pTrack[corrEvent.goodTracks[1]];
		NPhysicsTrack t3 = corrEvent.pTrack[corrEvent.goodTracks[2]];

		// separation from pi impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, t1);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_1 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, t2);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, t3);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		/*propPos = propagate(rootGeom.Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;*/

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(options.isOptDebug()) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(options.isOptDebug()) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			rootPhysics.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

int pi0d_goodClusters_tight(){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) conditions = 6;
	else conditions=5;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack pi = corrEvent.pTrack[rootPhysics.pic.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[rootPhysics.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[rootPhysics.em.parentTrack];

		// separation from pi impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, pi);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, ep);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, em);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(options.isOptDebug()) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(options.isOptDebug()) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			rootPhysics.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

bool pi0d_failCut(int i){
	if(!options.isDoScan())corrEvent.failedCond = i;
	if(options.isOptDebug()) cout << "Event is not passing selection " << i << endl;
	if(io.isDoOutput()) io.output.f1() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << " " << i << endl;
	return false;
}
void pi0d_passSelection(){
	if(options.isOptDebug()) cout << "Event is passing selection" << endl;
	if(io.isDoOutput()) io.output.f2() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << endl;
}

bool nico_pi0DalitzSelect(){
	bool badAcceptance;

	int badCombis=0;
	int piTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int piCandNb;
	bool badElectron;

	int goodClusters;
	//int goodClusterID=-1;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	double pt;

	vector<double> vMass;
	vector<TVector3> vP;

	if(corrEvent.failedCond>=0) {
		outputFileHeader.NFailedEvents--;
		pi0d_failCut(corrEvent.failedCond);
		return false;
	}
	if(options.isOptDebug()) cout << endl;

	// 6) Track DCH time
	if(options.isOptDebug()) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(rootBurst.isData){
		if(options.isOptDebug()) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
		if(options.isOptDebug()) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
		if(options.isOptDebug()) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
	}
	else{
		if(options.isOptDebug()) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(options.isOptDebug()) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(options.isOptDebug()) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==io.cutsDefinition.boolBadTrack) {pi0d_failCut(7); return false;}

	// 8) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return false;}


	// 12) Exactly 1 good LKr cluster (loose)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_loose();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}
	TLorentzVector tempGamma;
	tempGamma.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rawEvent.vtx[corrEvent.goodVertexID].position).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);

	// 13) Photon candidate in LKr acceptance
	if(options.isOptDebug()) cout << "~~~~ Cut 13 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options.isOptDebug()) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=io.cutsDefinition.lkrAcceptance) {pi0d_failCut(13); return false;}

	// 14) E_gamma>3GeV
	if(options.isOptDebug()) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(options.isOptDebug()) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << tempGamma.E() << "\t <= 3 : rejected" << endl;
	if(tempGamma.E()<=io.cutsDefinition.minGammaEnergy) {pi0d_failCut(14); return false;}

	// 15) D_deadcell>2cm
	if(options.isOptDebug()) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(options.isOptDebug()) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t1) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	// 16) Photon DCH1 intercept >13cm
	if(options.isOptDebug()) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(rootGeom.Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, tempGamma.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options.isOptDebug()) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=io.cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16); return false;}


	//Start Kinematic cuts
	TVector3 totalP; // = t1 + t2 + t3 + tGamma
	totalP = corrEvent.pTrack[corrEvent.goodTracks[0]].momentum + corrEvent.pTrack[corrEvent.goodTracks[1]].momentum + corrEvent.pTrack[corrEvent.goodTracks[2]].momentum + tempGamma.Vect();
	//Pt
	pt = totalP.Perp2(corrEvent.kaonMomentum);

	// 11) Tracks momenta
	if(options.isOptDebug()) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_1 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[0]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_2 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[1]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_3 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[2]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[0]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[0]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[1]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[1]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[2]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[2]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}

	// 17) Total momentum 70<p<78
	if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t <70 || >78 : rejected" << endl;
	if(totalP.Mag()<io.cutsDefinition.minTotalMomentum || totalP.Mag()>io.cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return false;}

	// 18) Transverse momentum^2 < 5E-4
	if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=io.cutsDefinition.maxPt) {pi0d_failCut(18); return false;}


	// PID
	flBad = false;
	ep=-1;
	em=-1;
	pic=-1;

	bool good=false, ided=false;

	flBad = associateMCTracks();

	pid_res.prelimEvents++;


	// 9) Identify candidates
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pid(piTrack, tempGamma);
	if(piCandNb==0) pid_res.nNoId++;
	else if(piCandNb==1){
		if(!flBad){
			if(piTrack == pic){
				pid_res.nGoodId++;
				good = true;
			}
			else{
				pid_res.nWrongId++;
			}
		}
		pid_res.nIded++;
		ided = true;
	}
	else if(piCandNb>1) pid_res.nManyId++;
	nPiCandidatesNew.Fill(piCandNb);
	badElectron = false;
	//piCandNb = pi0d_identifyPi(piTrack, badElectron);
	//if(piCandNb==1 && !badElectron) pid_res.nSelectedOld++;
	//nPiCandidates.Fill(piCandNb);
	//return true;
	if(options.isOptDebug()) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=io.cutsDefinition.numPiCandidates) {pi0d_failCut(9); return false;}

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if((int)i==piTrack) continue;

		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	//Create physics event from tracks (e+, e-, pi+-)
	rootPhysics.em.parentTrack = corrEvent.goodTracks[emTrack];
	rootPhysics.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.em.parentVertex = corrEvent.goodVertexID;
	rootPhysics.ep.parentTrack = corrEvent.goodTracks[epTrack];
	rootPhysics.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.ep.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pic.parentTrack = corrEvent.goodTracks[piTrack];
	rootPhysics.pic.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pic.parentVertex = corrEvent.goodVertexID;


	rootPhysics.em.P.SetVectM(corrEvent.pTrack[rootPhysics.em.parentTrack].momentum*corrEvent.pTrack[rootPhysics.em.parentTrack].p, Me);
	rootPhysics.ep.P.SetVectM(corrEvent.pTrack[rootPhysics.ep.parentTrack].momentum*corrEvent.pTrack[rootPhysics.ep.parentTrack].p, Me);
	rootPhysics.pic.P.SetVectM(corrEvent.pTrack[rootPhysics.pic.parentTrack].momentum*corrEvent.pTrack[rootPhysics.pic.parentTrack].p, Mpic);
	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		rootPhysics.pic.pdgID = 211;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		rootPhysics.pic.pdgID = -211;
	}

	// 10) Bad electron cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options.isOptDebug()) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==io.cutsDefinition.boolBadECandidates) {pi0d_failCut(10); return false;}

	// 12) Exactly 1 good LKr cluster (tight)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_tight();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event (gamma), build reco physics (pi0, kaon)
	rootPhysics.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.gamma.parentVertex = corrEvent.goodVertexID;
	rootPhysics.gamma.P.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rootPhysics.gamma.vertex).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);
	rootPhysics.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pi0.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pi0.P = rootPhysics.em.P + rootPhysics.ep.P + rootPhysics.gamma.P;
	rootPhysics.kaon.P = rootPhysics.pic.P + rootPhysics.pi0.P;

	//Start mass cuts
	// 19) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=io.cutsDefinition.maxPi0MassDiff) {pi0d_failCut(19); return false;}

	// 20) 0.475 < M_pieeg < 0.510
	if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
	if(rootPhysics.kaon.P.M()<io.cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>io.cutsDefinition.maxKaonMassDiff) {pi0d_failCut(20); return false;}

	//pi0dalitz variables
	rootPhysics.mee = (rootPhysics.em.P + rootPhysics.ep.P).M();
	rootPhysics.x = pow(rootPhysics.mee/Mpi0, 2.);
	rootPhysics.y = 2.*(rootPhysics.em.P.E() - rootPhysics.ep.P.E())/(Mpi0*(1-rootPhysics.x));

	if(ided){
		pid_res.nPassIded++;
		if(good) pid_res.nPassGood++;
		else pid_res.nPassWrong++;
	}

	pi0d_passSelection();
	return true;
}

bool newEvent(int i, int &nevt){
	rootPhysics.clear();

	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}
	if(i==0) cout << endl;

	pid_res.totEvents++;
	++nevt;
	abcog_params = rootBurst.abcog_params;
	if(!options.isDoScan()) return nico_pi0DalitzSelect(); // Normal thing
	else{
		bool result;
		bool defaultResult = false;
		bool globalResult = false;
		io.output.resetResult();
		for(int i=io.cutsDefinition.loadList(0); i!=-1; i=io.cutsDefinition.loadNextList()){
			result = nico_pi0DalitzSelect();
			io.output.newResult(result);
			globalResult |= result;
			if(i==io.cutsDefinition.getDefaultIndex()) defaultResult = result;
		}
		corrEvent.failedCond = defaultResult;
		return globalResult;
	}
}

int main(int argc, char **argv){
	options.parse(argc, argv, io);

	if(options.isDoScan()) io.cutsDefinition.generateLists(options.getScan());
	if(!io.cutsDefinition.addParseCutsFile(options.getCutsFile())) return -1;
	io.cutsDefinition.print();
	if(!options.isDoScan()) io.cutsDefinition.loadDefault(); // No scan, load the default values

	options.printSummary(io);
	options.parseFilter();

	io.openAll(options.isDoScan());

	TTree *pid = new TTree("pid", "pid");
	pid->Branch("pid_res.good", &pid_res.good, "good/I");
	pid->Branch("pid_res.bad", &pid_res.bad, "bad/I");
	pid->Branch("pid_res.compareMC", &pid_res.compareMC, "compareMC/I");
	pid->Branch("pid_res.totEvents", &pid_res.totEvents, "totEvents/I");
	pid->Branch("pid_res.prelimEvents", &pid_res.prelimEvents, "prelimEvents/I");
	pid->Branch("pid_res.nNoId", &pid_res.nNoId, "nNoId/I");
	pid->Branch("pid_res.nGoodId", &pid_res.nGoodId, "nGoodId/I");
	pid->Branch("pid_res.nWrongId", &pid_res.nWrongId, "nWrongId/I");
	pid->Branch("pid_res.nManyId", &pid_res.nManyId, "nManyId/I");
	pid->Branch("pid_res.nIded", &pid_res.nIded, "nIded/I");
	pid->Branch("pid_res.nPassGood", &pid_res.nPassGood, "nPassGood/I");
	pid->Branch("pid_res.nPassWrong", &pid_res.nPassWrong, "nPassWrong/I");
	pid->Branch("pid_res.nPassIded", &pid_res.nPassIded, "nPassIded/I");


	int nevent = 0;

	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io.input.getNEvents() << endl;
	for(int i=io.input.firstEvent(outputFileHeader); !io.input.eof() && (options.getMaxEvents()<0 || nevent<options.getMaxEvents() ); i = io.input.nextEvent(outputFileHeader)){
		if(newEvent(i, nevent)){
			io.output.fillEvent();
			outputFileHeader.NPassedEvents++;
		}
		else{
			outputFileHeader.NFailedEvents++;
			if(options.isExportAllEvents()) io.output.fillEvent();
		}
	}
	pid->Fill();
	cout << endl << endl;

	cout << "Processed events ->\t" << outputFileHeader.NProcessedEvents << endl;
	cout << "Failed events    ->\t" << outputFileHeader.NFailedEvents << endl;
	cout << "Passed events    ->\t" << outputFileHeader.NPassedEvents << endl;

	cout << "MC Association" << endl << "--------------" << endl;
	cout << "Good = " << pid_res.good << " " << pid_res.good*100./(double)(pid_res.good+pid_res.bad) << endl;
	cout << "Bad  = " << pid_res.bad << " " << pid_res.bad*100./(double)(pid_res.good+pid_res.bad) << endl;

	cout << "Global Counters (normalized to Total Events)" << endl << "---------------" << endl;
	cout << "Total events = " << pid_res.totEvents << endl;
	cout << "Pass preliminary = " << pid_res.prelimEvents << " " << pid_res.prelimEvents*100./(double)(pid_res.totEvents) << endl;
	cout << "Identified = " << pid_res.nIded << " " << pid_res.nIded*100./(double)(pid_res.totEvents) << endl;
	cout << "Associated = " << pid_res.compareMC << " " << pid_res.compareMC*100./(double)(pid_res.totEvents) << endl;
	cout << "Pass with ID = " << pid_res.nPassIded << " " << pid_res.nPassIded*100./(double)(pid_res.totEvents) << endl;

	cout << "ID categories (normalized to Associated)" << endl << "-------------" << endl;
	cout << "No ID = " << pid_res.nNoId << " " << pid_res.nNoId*100./(double)(pid_res.compareMC) << endl;
	cout << "Good ID = " << pid_res.nGoodId << " " << pid_res.nGoodId*100./(double)(pid_res.compareMC) << endl;
	cout << "Wrong ID = " << pid_res.nWrongId << " " << pid_res.nWrongId*100./(double)(pid_res.compareMC) << endl;
	cout << "Many ID = " << pid_res.nManyId << " " << pid_res.nManyId*100./(double)(pid_res.compareMC) << endl;
	cout << "Pass with good ID = " << pid_res.nPassGood << " " << pid_res.nPassGood*100./(double)pid_res.compareMC << endl;
	cout << "Pass with wrong ID = " << pid_res.nPassWrong << " " << pid_res.nPassWrong/(double)pid_res.compareMC << endl;
	meeTrue.Write();
	mkTrue.Write();
	meeepiTrue.Write();
	meekTrue.Write();

	meeDiffTrue.Write();
	mkDiffTrue.Write();

	meeFalse.Write();
	mkFalse.Write();
	meeepiFalse.Write();
	meekFalse.Write();

	meeDiffFalse.Write();
	mkDiffFalse.Write();

	nPiCandidates.Write();
	nPiCandidatesNew.Write();
	nBeamSign.Write();


	pid->Write();


	io.closeAll();
	return 0;
}


