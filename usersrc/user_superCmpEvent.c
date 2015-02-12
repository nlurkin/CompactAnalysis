/*****************************************************************/
/* COmPACT user routine: user_superCmpEvent(superCmpEvent *sevt) */
/*                                                               */
/* User routine called everytime an event `*sevt' is             */
/* loaded. A return value of greater than zero denotes           */
/* an error condition has occured.                               */
/*                                     BH 13/2/98    RWM 20/6/97 */
/*****************************************************************/

#include <math.h>

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <iostream>
#include <iomanip>
#include "compactLib.h"
#include "funLib.h"

using namespace std;

void clearAll(){
	int i = 0;

	rawEvent.clear();
	corrEvent.clear();
	rootBurst.clear();
}

int user_superCmpEvent(superBurst *sbur,superCmpEvent *sevt) {
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	//if(!eopLoaded) loadEOPData(sbur); Done in user_superBurst
	if(periodKeep!= 0 && period!=periodKeep) return 0;
	if(iEvent==0) cout << "First event: ";
	if(iEvent % outputMod == 0) cout << iEvent << " " << sbur->nrun << " " << sbur->time << " " << sevt->timeStamp << "                 \r";
	//####### DEBUGGING
	int run = 20412;
	run = -1;
	int burst = 1188437272;
	int event = 5286224;
	//#################

	int passEvent;
	if(run!=-1){
		if(!((sbur->nrun==run) && (sbur->time==burst) && (sevt->timeStamp==event)))return 0;
	}
	if(opts.count("filter")>0){
		if(!isFilteredEvent(sbur->nrun, sbur->time, sevt->timeStamp)) return 0;
	}
	if(iEvent==0) cout << endl;

	iEvent++;

	clearAll();

	// 1) Apply all corrections
	if(dataOnly) user_lkrcalcor_SC(sbur,sevt,1);
	rootBurst = sbur;
	rawEvent = sevt;
	rootBurst.abcog_params = &abcog_params;
	rootGeom = Geom;
	rootBurst.isData = dataOnly;
	rootBurst.isMC = mcOnly;
	rootBurst.pbWall = pbWall;
	CreateTracks(sevt);
	CreateClusters(sevt);

	rootFileHeader.NProcessedEvents++;
	//Do it later: depends on the detected beam charge
	//kaonMomentum = TVector3(abcog_params.pkdxdzp, abcog_params.pkdydzp, 1.).Unit();
	//kaonP = abcog_params.pkp*(1+abcog_params.beta);



	if(channel==KE2){
		//passEvent = nico_ke2Select(sbur, sevt);
	}
	if(channel==PI0DALITZ){
		passEvent = nico_pi0DalitzSelect();
		if(passEvent==0){
			nico_pi0DalitzAna(sbur, sevt);
			rootFileHeader.NPassedEvents++;
		}
		else rootFileHeader.NFailedEvents++;

		outTree->Fill();
	}

	/*----------- End of user C code -----------*/
	return 0;
}
