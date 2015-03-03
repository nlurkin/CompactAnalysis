#define OUTSIDECOMPACT

/// Compact includes
#include "funLib.h"

/// Std include
#include <iomanip>
using namespace std;

// Local includes
#include "CompactIO.h"
#include "OptionsParser.h"

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


//Selection
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

int pi0d_goodClusters(){
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

	// 9) Identify candidates
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pi0d_identifyPi(piTrack, badElectron);
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

	//Create physics event from tracks
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
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		rootPhysics.pic.pdgID = 211;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		rootPhysics.pic.pdgID = -211;
	}

	// 8) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return false;}

	// 10) Bad electron cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options.isOptDebug()) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==io.cutsDefinition.boolBadECandidates) {pi0d_failCut(10); return false;}

	// 12) Exactly 1 good LKr cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event
	rootPhysics.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.gamma.parentVertex = corrEvent.goodVertexID;
	rootPhysics.gamma.P.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rootPhysics.gamma.vertex).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);
	rootPhysics.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pi0.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pi0.P = rootPhysics.em.P + rootPhysics.ep.P + rootPhysics.gamma.P;
	rootPhysics.kaon.P = rootPhysics.pic.P + rootPhysics.pi0.P;

	// 13) Photon candidate in LKr acceptance
	if(options.isOptDebug()) cout << "~~~~ Cut 13 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options.isOptDebug()) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=io.cutsDefinition.lkrAcceptance) {pi0d_failCut(13); return false;}

	// 14) E_gamma>3GeV
	if(options.isOptDebug()) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(options.isOptDebug()) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << rootPhysics.gamma.P.E() << "\t <= 3 : rejected" << endl;
	if(rootPhysics.gamma.P.E()<=io.cutsDefinition.minGammaEnergy) {pi0d_failCut(14); return false;}

	// 15) D_deadcell>2cm
	if(options.isOptDebug()) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(options.isOptDebug()) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(pi) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.pic.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.pic.parentTrack].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(ep) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.ep.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.ep.parentTrack].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(em) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.em.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.em.parentTrack].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	// 16) Photon DCH1 intercept >13cm
	if(options.isOptDebug()) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(rootGeom.Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, rootPhysics.gamma.P.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options.isOptDebug()) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=io.cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16); return false;}

	//Start Kinematic cuts
	//Pt
	pt = rootPhysics.kaon.P.Perp2(corrEvent.kaonMomentum);

	//Mee
	rootPhysics.mee = (rootPhysics.em.P + rootPhysics.ep.P).M();
	rootPhysics.x = pow(rootPhysics.mee/Mpi0, 2.);
	rootPhysics.y = 2.*(rootPhysics.em.P.E() - rootPhysics.ep.P.E())/(Mpi0*(1-rootPhysics.x));

	// 11) Tracks momenta
	if(options.isOptDebug()) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_pi :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_e+ :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_e- :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}

	// 17) Total momentum 70<p<78
	if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_pieeg :\t\t\t" << rootPhysics.kaon.P.Vect().Mag() << "\t <70 || >78 : rejected" << endl;
	if(rootPhysics.kaon.P.Vect().Mag()<io.cutsDefinition.minTotalMomentum || rootPhysics.kaon.P.Vect().Mag()>io.cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return false;}

	// 18) Transverse momentum^2 < 5E-4
	if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=io.cutsDefinition.maxPt) {pi0d_failCut(18); return false;}

	// 19) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=io.cutsDefinition.maxPi0MassDiff) {pi0d_failCut(19); return false;}

	// 20) 0.475 < M_pieeg < 0.510
	if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
	if(rootPhysics.kaon.P.M()<io.cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>io.cutsDefinition.maxKaonMassDiff) {pi0d_failCut(20); return false;}

	pi0d_passSelection();
	return true;
}

bool newEvent(int i, int &nevt){
	rootPhysics.clear();

	if(options.getPeriodKeep()!= 0 && rootBurst.period!=options.getPeriodKeep()) return false;
	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}
	if(i==0) cout << endl;

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
	if(!options.parse(argc, argv, io)) return -1;

	if(options.isDoScan()) io.cutsDefinition.generateLists(options.getScan());
	if(!io.cutsDefinition.addParseCutsFile(options.getCutsFile())) return -1;
	io.cutsDefinition.print();
	if(!options.isDoScan()) io.cutsDefinition.loadDefault(); // No scan, load the default values

	options.printSummary(io);
	options.parseFilter();

	io.openAll(options.isDoScan());

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
	cout << endl << endl;

	cout << "Processed events ->\t" << outputFileHeader.NProcessedEvents << endl;
	cout << "Failed events    ->\t" << outputFileHeader.NFailedEvents << endl;
	cout << "Passed events    ->\t" << outputFileHeader.NPassedEvents << endl;

	io.closeAll();
	return 0;
}
