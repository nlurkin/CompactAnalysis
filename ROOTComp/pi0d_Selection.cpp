#define OUTSIDECOMPACT

/// Std include
#include <iomanip>
using namespace std;

/// Compact includes
#include "funLib.h"

// Local includes
#include "SelectionFunctions.h"
#include "pid_res.h"

struct alt_pid_res pid_res_pi, pid_res_mu;

typedef struct tempObjects_t{
	TLorentzVector tempGamma;
	ROOTPhysicsEvent piEvent;
	ROOTPhysicsEvent muEvent;
	int xPreSelected, xPreSelected_Mu;
	int nPreCandidates, nPreCandidates_Mu;
} tempObjects;

int firstCutIndex = 5;

TH1D selTrackDiff("selTrackDiff", "selTrackDiff", 10, -5, 5);

bool nico_pi0DalitzSelect_Common(tempObjects &tempObj){
	bool badAcceptance;

	int badCombis=0;

	int goodClusters;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	vector<double> vMass;
	vector<TVector3> vP;

	if(corrEvent.failedCond>=0) {
		outputFileHeader.NFailedEvents--;
		pi0d_failCut(corrEvent.failedCond);
		return false;
	}
	if(options.isOptDebug()) cout << endl;


	// 1) Track DCH time
	if(options.isOptDebug()) cout << "~~~~ Cut 1 ~~~~" << endl;
	if(rootBurst.isData){
		if(options.isOptDebug()) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>io.cutsDefinition.maxTrackTime) {pi0d_failCut(1+firstCutIndex); return false;}
		if(options.isOptDebug()) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>io.cutsDefinition.maxTrackTime) {pi0d_failCut(1+firstCutIndex); return false;}
		if(options.isOptDebug()) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>io.cutsDefinition.maxTrackTime) {pi0d_failCut(1+firstCutIndex); return false;}
	}
	else{
		if(options.isOptDebug()) cout << "\tMC: Not applicable" << endl;
	}

	// 2) Track acceptance veto
	if(options.isOptDebug()) cout << "~~~~ Cut 2 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(options.isOptDebug()) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==io.cutsDefinition.boolBadTrack) {pi0d_failCut(2+firstCutIndex); return false;}


	// 3) Track combination veto loose
	if(options.isOptDebug()) cout << "~~~~ Cut 3 ~~~~" << endl;
//	badCombis = pi0d_trackCombinationVeto(); //Wrong - depends on pi+ id
	badCombis = pi0d_trackCombinationVeto_loose();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(3+firstCutIndex); return false;}

	// 4) Michal pre identification
	tempObj.nPreCandidates = 2;
	if(options.isOptDebug()) cout << "~~~~ Cut 4 ~~~~" << endl;
	//tempObj.nPreCandidates = michal_prepid(tempObj.xPreSelected, OptionsParser::K2PI);
	//Useless for KMU3
	//tempObj.nPreCandidates_Mu = michal_prepid(tempObj.xPreSelected_Mu, OptionsParser::KMU3);
	if(options.isOptDebug()) cout << "Michal pre-id:\t\t (" << tempObj.nPreCandidates << /*" " << tempObj.nPreCandidates_Mu <<*/ ")\t == 0 : rejected" << endl;
	//if(tempObj.nPreCandidates==0 /*&& tempObj.nPreCandidates_Mu==0*/) {pi0d_failCut(4+firstCutIndex); return false;}
	//if(tempObj.nPreCandidates>1) tempObj.xPreSelected=-1;
	//if(tempObj.nPreCandidates_Mu>1) tempObj.xPreSelected_Mu=-1;

	// 5) Exactly 1 good LKr cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 5 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_loose();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(5+firstCutIndex); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}
	tempObj.tempGamma.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rawEvent.vtx[corrEvent.goodVertexID].position).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);

	if(rootBurst.pbWall){
		if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << corrEvent.pCluster[rootPhysics.gamma.parentCluster].position.Y() << " < -11.850 : reject" << endl;
		if(corrEvent.pCluster[rootPhysics.gamma.parentCluster].position.Y()>-33.575 && corrEvent.pCluster[rootPhysics.gamma.parentCluster].position.Y() < -11.850) {pi0d_failCut(5+firstCutIndex); return false;}
	}

	// 6) Photon candidate in LKr acceptance
	if(options.isOptDebug()) cout << "~~~~ Cut 6 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options.isOptDebug()) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=io.cutsDefinition.lkrAcceptance) {pi0d_failCut(6+firstCutIndex); return false;}

	// 7) E_gamma>3GeV
	if(options.isOptDebug()) cout << "~~~~ Cut 7 ~~~~" << endl;
	if(options.isOptDebug()) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << tempObj.tempGamma.E() << "\t <= 3 : rejected" << endl;
	if(tempObj.tempGamma.E()<io.cutsDefinition.minGammaEnergy) {pi0d_failCut(7+firstCutIndex); return false;}

	// 8) D_deadcell>2cm
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	if(options.isOptDebug()) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(8+firstCutIndex); return false;}

	/*if(options.isOptDebug()) cout << "d_deadcell(t1) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(8+firstCutIndex); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(8+firstCutIndex); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(8+firstCutIndex); return false;}*/

	// 9) Photon DCH1 intercept >13cm
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	propPos = propagate(rootGeom.Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, tempObj.tempGamma.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options.isOptDebug()) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=io.cutsDefinition.minGammaDCHRadius) {pi0d_failCut(9+firstCutIndex); return false;}


	//Start Kinematic cuts
	// 10) Tracks momenta
	if(options.isOptDebug()) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_1 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[0]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_2 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[1]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_3 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[2]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[0]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[0]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(10+firstCutIndex); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[1]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[1]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(10+firstCutIndex); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[2]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[2]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(10+firstCutIndex); return false;}

	return true;
}

int nico_pi0DalitzSelect_K2PI(tempObjects &tempObj, bool &good, bool &bad){
	int xTrack = -1;
	int xCandNb = 0;

	int epTrack = -1;
	int emTrack = -1;

	int goodClusters;

	TVector3 totalP; // = t1 + t2 + t3 + tGamma
	totalP =  corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p \
			+ tempObj.tempGamma.Vect();
	//Pt
	double pt = totalP.Perp2(corrEvent.kaonMomentum);
	// 11) Total momentum 70<p<78 (K2PI)
	if(options.isOptDebug()) cout << "~~~~ Cut 11 (K2PI)~~~~" << endl;
	if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t <70 || >78 : rejected" << endl;
	if(totalP.Mag()<io.cutsDefinition.k2pi.minTotalMomentum || totalP.Mag()>io.cutsDefinition.k2pi.maxTotalMomentum) return 11+firstCutIndex;

	// 12) Transverse momentum^2 < 5E-4 (K2PI)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 (K2PI) ~~~~" << endl;
	if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.k2pi.maxPt << " : rejected" << endl;
	if(pt>=io.cutsDefinition.k2pi.maxPt) return 12+firstCutIndex;

	// PID
	xTrue = 999;
	xFalse = 999;

	// Identify candidates
	xCandNb = pid(xTrack, tempObj.tempGamma, OptionsParser::K2PI);

	if(xCandNb==0){
		pid_res_pi.incNoID(!flBad);
		xMCNoID.Fill(rootMC.xTrue);
		xTruexFalseNo.Fill(xTrue, xFalse);
		xTruexMCNo.Fill(xTrue, rootMC.xTrue);
	}
	else if(xCandNb==1){
		if(!flBad){
			if(xTrack == xPart){
				pid_res_pi.incGoodId();
				good = true;
			}
			else{
				xTruexFalse.Fill(xTrue, xFalse);
				pid_res_pi.incBadId();
				bad = true;
			}
		}
		pid_res_pi.incIded(!flBad);
	}
	else if(xCandNb>1){
		pid_res_pi.incManyID(!flBad);
		xMCManyID.Fill(rootMC.xTrue);
		xTruexFalseMany.Fill(xTrue, xFalse);
		xTruexMCMany.Fill(xTrue, rootMC.xTrue);
	}
	nxCandidatesNew.Fill(xCandNb);
	// 13) Number of candidates K2PI
	if(options.isOptDebug()) cout << "~~~~ Cut 13 (K2PI)~~~~" << endl;
	if(options.isOptDebug()) cout << "Number of pi track candidates :\t" << xCandNb << "\t != 1: rejected" << endl;
	if(xCandNb!=io.cutsDefinition.numXCandidates) return 13+firstCutIndex;

	if(tempObj.xPreSelected!=-1) selTrackDiff.Fill(tempObj.xPreSelected-xTrack);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if((int)i!=xTrack){
			if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
			else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
			else{
				if(epTrack==-1) epTrack=i;
				if(emTrack==-1) emTrack=i;
			}
		}
	}

	//Create physics event from tracks (e+, e-, x+-)
	tempObj.piEvent.em.parentTrack = corrEvent.goodTracks[emTrack];
	tempObj.piEvent.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.piEvent.em.parentVertex = corrEvent.goodVertexID;
	tempObj.piEvent.ep.parentTrack = corrEvent.goodTracks[epTrack];
	tempObj.piEvent.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.piEvent.ep.parentVertex = corrEvent.goodVertexID;
	tempObj.piEvent.pic.parentTrack = corrEvent.goodTracks[xTrack];
	tempObj.piEvent.pic.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.piEvent.pic.parentVertex = corrEvent.goodVertexID;

	tempObj.piEvent.em.P.SetVectM(corrEvent.pTrack[tempObj.piEvent.em.parentTrack].momentum*corrEvent.pTrack[tempObj.piEvent.em.parentTrack].p, Me);
	tempObj.piEvent.ep.P.SetVectM(corrEvent.pTrack[tempObj.piEvent.ep.parentTrack].momentum*corrEvent.pTrack[tempObj.piEvent.ep.parentTrack].p, Me);
	tempObj.piEvent.pic.P.SetVectM(corrEvent.pTrack[tempObj.piEvent.pic.parentTrack].momentum*corrEvent.pTrack[tempObj.piEvent.pic.parentTrack].p, Mpic);
	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		tempObj.piEvent.kaon.pdgID = 321;
		tempObj.piEvent.pic.pdgID = 211;
	}
	else{
		tempObj.piEvent.kaon.pdgID = -321;
		tempObj.piEvent.pic.pdgID = -211;
	}

	//plot e/p
	double xeop = corrEvent.pTrack[tempObj.piEvent.pic.parentTrack].E/tempObj.piEvent.pic.P.Vect().Mag();
	double epeop = corrEvent.pTrack[tempObj.piEvent.ep.parentTrack].E/tempObj.piEvent.ep.P.Vect().Mag();
	double emeop = corrEvent.pTrack[tempObj.piEvent.em.parentTrack].E/tempObj.piEvent.em.P.Vect().Mag();
	eopx.Fill(xeop);
	eope.Fill(epeop);
	eope.Fill(emeop);

	if(xeop < 0.85){
		eope_goodx.Fill(epeop);
		eope_goodx.Fill(emeop);
	}
	else{
		eope_badx.Fill(epeop);
		eope_badx.Fill(emeop);
	}
	if(epeop > 0.85 && epeop < 1.15 && emeop > 0.85 && emeop < 1.15) eopx_goode.Fill(xeop);
	else eopx_bade.Fill(xeop);

	// 14) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 14 ~~~~" << endl;
	int badCombis = pi0d_trackCombinationVeto_tight(tempObj.piEvent.pic);
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) return 14+firstCutIndex;

	// 15) Exactly 1 good LKr cluster (tight)
	if(options.isOptDebug()) cout << "~~~~ Cut 15 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_tight(tempObj.piEvent.pic, tempObj.piEvent);
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) return 15+firstCutIndex;

	if(tempObj.piEvent.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event (gamma), build reco physics (pi0, kaon)
	tempObj.piEvent.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.piEvent.gamma.parentVertex = corrEvent.goodVertexID;
	tempObj.piEvent.gamma.P.SetVectM((corrEvent.pCluster[tempObj.piEvent.gamma.parentCluster].position - tempObj.piEvent.gamma.vertex).Unit()*corrEvent.pCluster[tempObj.piEvent.gamma.parentCluster].E, 0.0);
	tempObj.piEvent.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.piEvent.pi0.parentVertex = corrEvent.goodVertexID;
	tempObj.piEvent.pi0.P = tempObj.piEvent.em.P + tempObj.piEvent.ep.P + tempObj.piEvent.gamma.P;
	tempObj.piEvent.kaon.P = tempObj.piEvent.pic.P + tempObj.piEvent.pi0.P;

	//Start mass cuts
	// 16) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 16 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << tempObj.piEvent.pi0.P.M() << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(tempObj.piEvent.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(tempObj.piEvent.pi0.P.M()-Mpi0)>=io.cutsDefinition.k2pi.maxPi0MassDiff) return 16+firstCutIndex;


	// 17) 0.475 < M_pieeg < 0.510
	if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options.isOptDebug()) cout << "DB: " << MK << " mass: " << tempObj.piEvent.kaon.P.M() << endl;
	if(options.isOptDebug()) cout << "M_pieeg :\t\t" << fabs(tempObj.piEvent.kaon.P.M() - MK) << "\t >" << io.cutsDefinition.k2pi.maxKaonMassDiff << ": rejected" << endl;
	if(fabs(tempObj.piEvent.kaon.P.M() - MK) > io.cutsDefinition.k2pi.maxKaonMassDiff) return 17+firstCutIndex;

	//pi0dalitz variables
	tempObj.piEvent.mee = (tempObj.piEvent.em.P + tempObj.piEvent.ep.P).M();
	tempObj.piEvent.x = pow(tempObj.piEvent.mee/Mpi0, 2.);
	tempObj.piEvent.y = 2.*(tempObj.piEvent.em.P.E() - tempObj.piEvent.ep.P.E())/(Mpi0*(1-tempObj.piEvent.x));

	//Trigger cut
	// 18) L2: E_lkr > 14GeV
	double ELKr_ep = 0;
	double ELKr_em = 0;
	bool goodPBWall;
	TVector3 propPos;
	NPhysicsTrack t_ep = corrEvent.pTrack[tempObj.piEvent.ep.parentTrack];
	NPhysicsTrack t_em = corrEvent.pTrack[tempObj.piEvent.em.parentTrack];

	propPos = propagateAfter(rootGeom.Lkr.z, t_ep);
	goodPBWall = true;
	if(rootBurst.pbWall && (propPos.Y()>-33.575 && propPos.Y() < -11.850)) goodPBWall = false;
	if(t_ep.lkr_acc==0 && goodPBWall && rawEvent.track[t_ep.trackID].dDeadCell>2.) ELKr_ep = t_ep.p;

	propPos = propagateAfter(rootGeom.Lkr.z, t_em);
	goodPBWall = true;
	if(rootBurst.pbWall && (propPos.Y()>-33.575 && propPos.Y() < -11.850)) goodPBWall = false;
	if(t_em.lkr_acc==0 && goodPBWall && rawEvent.track[t_em.trackID].dDeadCell>2.) ELKr_em = t_em.p;

	double E_lkr = ELKr_ep + ELKr_em + tempObj.piEvent.gamma.P.E();
	if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(options.isOptDebug()) cout << "e+: lkr_acc=" << (t_ep.lkr_acc==0) << " pbwall=" << goodPBWall <<
			" dDeadCell=" << rawEvent.track[t_ep.trackID].dDeadCell <<
			" ELKr=" << ELKr_ep << endl;
	if(options.isOptDebug()) cout << "e-: lkr_acc=" << (t_em.lkr_acc==0) << " pbwall=" << goodPBWall <<
			" dDeadCell=" << rawEvent.track[t_em.trackID].dDeadCell <<
			" ELKr=" << ELKr_em << endl;
	if(options.isOptDebug()) cout << "E_LKr :\t\t\t\t" << E_lkr << "\t <14: rejected" << endl;
	if(E_lkr < 14) return 18+firstCutIndex;


	if(!pi0d_L3Trigger(t_ep) && !pi0d_L3Trigger(t_em)) return 19+firstCutIndex;

	if(tempObj.piEvent.x <= 0.01) return 20+firstCutIndex;
	return -1;
}

int nico_pi0DalitzSelect_KMU3(tempObjects &tempObj, bool &good, bool &bad){
	int xTrack = -1;
	int xCandNb = 0;

	int epTrack = -1;
	int emTrack = -1;

	int goodClusters;

	TVector3 totalP; // = t1 + t2 + t3 + tGamma
	totalP =  corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p \
			+ tempObj.tempGamma.Vect();
	//Pt
	double pt = totalP.Perp2(corrEvent.kaonMomentum);
	// 11) Total momentum p<78 (KMU3)
	if(options.isOptDebug()) cout << "~~~~ Cut 11 (KMU3)~~~~" << endl;
	if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t >78 : rejected" << endl;
	if(totalP.Mag()>io.cutsDefinition.kmu3.maxTotalMomentum) return 11+firstCutIndex;

	// 12) min value < Transverse momentum^2 < max value (KMU3)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 (KMU3) ~~~~" << endl;
	if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.kmu3.maxPt << " || <= " << io.cutsDefinition.kmu3.minPt << " : rejected" << endl;
	if(pt <= io.cutsDefinition.kmu3.minPt || pt>=io.cutsDefinition.kmu3.maxPt) return 12+firstCutIndex;

	// PID
	xTrue = 999;
	xFalse = 999;

	// Identify candidates
	xCandNb = pid(xTrack, tempObj.tempGamma, OptionsParser::KMU3);

	if(xCandNb==0){
		pid_res_mu.incNoID(!flBad);
		xMCNoID.Fill(rootMC.xTrue);
		xTruexFalseNo.Fill(xTrue, xFalse);
		xTruexMCNo.Fill(xTrue, rootMC.xTrue);
	}
	else if(xCandNb==1){
		if(!flBad){
			if(xTrack == xPart){
				pid_res_mu.incGoodId();
				good = true;
			}
			else{
				xTruexFalse.Fill(xTrue, xFalse);
				pid_res_mu.incBadId();
				bad = true;
			}
		}
		pid_res_mu.incIded(!flBad);
	}
	else if(xCandNb>1){
		pid_res_mu.incManyID(!flBad);
		xMCManyID.Fill(rootMC.xTrue);
		xTruexFalseMany.Fill(xTrue, xFalse);
		xTruexMCMany.Fill(xTrue, rootMC.xTrue);
	}
	nxCandidatesNew.Fill(xCandNb);
	// 13) Number of candidates K2PI
	if(options.isOptDebug()) cout << "~~~~ Cut 13 (KMU3)~~~~" << endl;
	if(options.isOptDebug()) cout << "Number of mu track candidates :\t" << xCandNb << "\t != 1: rejected" << endl;
	if(xCandNb!=io.cutsDefinition.numXCandidates) return 13+firstCutIndex;

	if(tempObj.xPreSelected_Mu!=-1) selTrackDiff.Fill(tempObj.xPreSelected_Mu-xTrack);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	//Create physics event from tracks (e+, e-, x+-)
	tempObj.muEvent.em.parentTrack = corrEvent.goodTracks[emTrack];
	tempObj.muEvent.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.muEvent.em.parentVertex = corrEvent.goodVertexID;
	tempObj.muEvent.ep.parentTrack = corrEvent.goodTracks[epTrack];
	tempObj.muEvent.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.muEvent.ep.parentVertex = corrEvent.goodVertexID;
	tempObj.muEvent.mu.parentTrack = corrEvent.goodTracks[xTrack];
	tempObj.muEvent.mu.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.muEvent.mu.parentVertex = corrEvent.goodVertexID;

	tempObj.muEvent.em.P.SetVectM(corrEvent.pTrack[tempObj.muEvent.em.parentTrack].momentum*corrEvent.pTrack[tempObj.muEvent.em.parentTrack].p, Me);
	tempObj.muEvent.ep.P.SetVectM(corrEvent.pTrack[tempObj.muEvent.ep.parentTrack].momentum*corrEvent.pTrack[tempObj.muEvent.ep.parentTrack].p, Me);
	tempObj.muEvent.mu.P.SetVectM(corrEvent.pTrack[tempObj.muEvent.mu.parentTrack].momentum*corrEvent.pTrack[tempObj.muEvent.mu.parentTrack].p, Mmu);

	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		tempObj.muEvent.kaon.pdgID = 321;
		tempObj.muEvent.mu.pdgID = 13;
	}
	else{
		tempObj.muEvent.kaon.pdgID = -321;
		tempObj.muEvent.mu.pdgID = -13;
	}


	//plot e/p
	double xeop = corrEvent.pTrack[tempObj.muEvent.mu.parentTrack].E/tempObj.muEvent.mu.P.Vect().Mag();
	double epeop = corrEvent.pTrack[tempObj.muEvent.ep.parentTrack].E/tempObj.muEvent.ep.P.Vect().Mag();
	double emeop = corrEvent.pTrack[tempObj.muEvent.em.parentTrack].E/tempObj.muEvent.em.P.Vect().Mag();
	eopx.Fill(xeop);
	eope.Fill(epeop);
	eope.Fill(emeop);

	if(xeop < 0.85){
		eope_goodx.Fill(epeop);
		eope_goodx.Fill(emeop);
	}
	else{
		eope_badx.Fill(epeop);
		eope_badx.Fill(emeop);
	}
	if(epeop > 0.85 && epeop < 1.15 && emeop > 0.85 && emeop < 1.15) eopx_goode.Fill(xeop);
	else eopx_bade.Fill(xeop);

	/*// 14) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 14 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto_tight();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) return 14+firstCutIndex;*/

	// 15) Exactly 1 good LKr cluster (tight)
	if(options.isOptDebug()) cout << "~~~~ Cut 15 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_tight(tempObj.muEvent.mu, tempObj.muEvent);
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) return 15+firstCutIndex;

	if(tempObj.muEvent.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event (gamma), build reco physics (pi0, kaon)
	tempObj.muEvent.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.muEvent.gamma.parentVertex = corrEvent.goodVertexID;
	tempObj.muEvent.gamma.P.SetVectM((corrEvent.pCluster[tempObj.muEvent.gamma.parentCluster].position - tempObj.muEvent.gamma.vertex).Unit()*corrEvent.pCluster[tempObj.muEvent.gamma.parentCluster].E, 0.0);
	tempObj.muEvent.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	tempObj.muEvent.pi0.parentVertex = corrEvent.goodVertexID;
	tempObj.muEvent.pi0.P = tempObj.muEvent.em.P + tempObj.muEvent.ep.P + tempObj.muEvent.gamma.P;
	tempObj.muEvent.kaon.P = tempObj.muEvent.mu.P + tempObj.muEvent.pi0.P;

	//Start mass cuts
	// 16) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 16 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << tempObj.muEvent.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(tempObj.muEvent.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(tempObj.muEvent.pi0.P.M()-Mpi0)>=io.cutsDefinition.k2pi.maxPi0MassDiff) return 16+firstCutIndex;


	// 17) Missing mass square < 0.01 MeV^2
	double mmasssq = (tempObj.muEvent.kaon.P - tempObj.muEvent.mu.P - tempObj.muEvent.pi0.P).M2();
	if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_miss^2:\t\t" << mmasssq << "\t >0.01: rejected" << endl;
	if(mmasssq > io.cutsDefinition.kmu3.maxMissMassSq) return 17+firstCutIndex;

	return -1;
}

bool nico_pi0DalitzSelect(){
	int retK2PI=0, retKMU3=0;
	tempObjects tempObj;

	if(!nico_pi0DalitzSelect_Common(tempObj)) return false;

	vector<float> eopSort = sortEOP();
	eoplowest.Fill(eopSort[0]);
	if(eopSort[0]<0.85){
		eopsecond.Fill(eopSort[1]);
		if(eopSort[1]>0.85 && eopSort[1]<1.15){
			eophighest.Fill(eopSort[2]);
		}
	}

	flBad = false;
	ep=-1;
	em=-1;
	xPart=-1;

	bool good=false, bad=false;

	flBad = associateMCTracks(pid_res_pi, &pid_res_mu);

	pid_res_pi.incPrelim();
	pid_res_mu.incPrelim();

	if(tempObj.nPreCandidates>0) retK2PI = nico_pi0DalitzSelect_K2PI(tempObj, good, bad);
	else retK2PI=4+firstCutIndex;
	if(options.getSelectionType()==OptionsParser::KMU3) retKMU3 = nico_pi0DalitzSelect_KMU3(tempObj, good, bad);

	if(retK2PI!=-1 && retKMU3!=-1) {
		int m = max(retK2PI,retKMU3);
		pid_res_pi.incNonPass(!flBad, good, bad);
		pid_res_mu.incNonPass(!flBad, good, bad);
		if(m<12+firstCutIndex) {
			pi0d_failCut(m); return false;
		}
		else {
			pi0d_failCutInc(m, !flBad, good, bad, NULL); return false;
		}
	}
	if(retK2PI==-1) pid_res_pi.incPass(!flBad, good, bad);
	else pid_res_pi.incNonPass(!flBad, good, bad);
	if(retKMU3==-1) pid_res_mu.incPass(!flBad, good, bad);
	else pid_res_mu.incNonPass(!flBad, good, bad);


	if(retK2PI==-1){
		rootPhysics = tempObj.piEvent;
		if(retKMU3==-1) rootPhysics.mu = tempObj.muEvent.mu;
	}
	else rootPhysics = tempObj.muEvent;
	pi0d_passSelection();
	return true;
}


bool newEvent(int i, int &nevt){
	//Clear the event
	rootPhysics.clear();

	// Filter events
	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}

	// Debugging info
	if(options.isOptDebug()) cout << "--------------------------------------------" << endl;
	if(options.getPeriodKeep()!= 0 && rootBurst.period!=options.getPeriodKeep()) return false;
	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(i==0) cout << endl;
	if(options.isOptDebug()) cout << endl << "--------------------------------------------" << endl;

	pid_res_pi.incTotal();
	pid_res_mu.incTotal();
	++nevt;
	abcog_params = rootBurst.abcog_params;
	if(!options.isDoScan()) return nico_pi0DalitzSelect(); // Normal thing
	else{
		// Scan cuts
		bool result;
		bool defaultResult = false;
		bool globalResult = false;

		// Reset result
		io.output.resetResult();
		pid_res_pi.ResetEntry();
		pid_res_mu.ResetEntry();

		// Loop over all cuts definition
		for(int i=io.cutsDefinition.loadList(0); i!=-1; i=io.cutsDefinition.loadNextList()){
			result = nico_pi0DalitzSelect();
			io.output.newResult(result);

			// Set as success if pass at least one cut
			globalResult |= result;
			if(i==io.cutsDefinition.getDefaultIndex()) defaultResult = result;
			pid_res_pi.NewEntry();
			pid_res_mu.NewEntry();
		}
		corrEvent.failedCond = defaultResult;
		return globalResult;
	}
}

int main(int argc, char **argv){

	// Get options

	if(!options.parse(argc, argv, io)) return 0;

	// Generate cuts
	if(options.isDoScan()) io.cutsDefinition.generateLists(options.getScan());
	if(!io.cutsDefinition.addParseCutsFile(options.getCutsFile())) return -1;
	io.cutsDefinition.print();
	if(!options.isDoScan()){
		io.cutsDefinition.loadDefault(); // No scan, load the default values
		pid_res_pi.Init(1);
		pid_res_mu.Init(1);
	}
	else{
		pid_res_pi.Init(options.getScan());
		pid_res_mu.Init(options.getScan());
	}

	options.printSummary(io);
	options.parseFilter();

	io.openAll(options.isDoScan());

	// Create the output pid_res tree
	TTree *pid_pi = new TTree("pid_pi", "pid_pi");
	//branchResTree(pid_pi, pid_res_pi);
	TTree *pid_mu = new TTree("pid_mu", "pid_mu");
	//branchResTree(pid_mu, pid_res_mu);

	// Initialise event loop
	int nevent = 0;
	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io.input.getNEvents() << endl;

	// Event loop
	for(int i=io.input.firstEvent(outputFileHeader, options.getStartEvent());
			!io.input.eof() && (options.getMaxEvents()<0 || nevent<options.getMaxEvents() );
			i = io.input.nextEvent(outputFileHeader)){
		if(newEvent(i, nevent)){
			io.output.fillEvent();
			outputFileHeader.NPassedEvents++;
		}
		else{
			outputFileHeader.NFailedEvents++;
			if(options.isExportAllEvents()) io.output.fillEvent();
		}
	}
	pid_pi->Fill();
	pid_mu->Fill();

	// Print pid result structure
	cout << endl << endl;
	//printResStruct(pid_res_pi, &pid_res_mu);

	selTrackDiff.Write();
	//pid_pi->Write();
	//pid_mu->Write();
	savePlots();
	io.closeAll();
	return 0;
}


