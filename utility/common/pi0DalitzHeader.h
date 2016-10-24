#include <fstream>
#include <TApplication.h>
#include "ConfigFile.h"
#include "../userinc/funLib.h"
#include "Interface/Sample.h"

using namespace std;

ConfigFile cfg;
bool filterLoaded = false;
std::vector<eventID> badEventsList;

/*************************
 * Globals
 *************************/
TApplication *theApp = 0;

/*************************
 * Signal handling
 *************************/
void sighandler(int)
{
	if(theApp) theApp->Terminate();
	exit(0);
}

/************************
 * Getting input
 ************************/
void loadBins(double *bins, int& nbins){
	ifstream fd(cfg.getBinsFileName());

	double val;
	nbins = 0;
	while(fd >> val){
		bins[nbins] = val;
		nbins++;
	}
	fd.close();
}

void loadEventsFilter(){
	int run, burst, time;

	ifstream fd("/afs/cern.ch/user/n/nlurkin/Compact/filterout.dat");
	while(fd >> run >> burst >> time){
		badEventsList.push_back(eventID(run, burst, time));
	}
	fd.close();


//	fd.open("/afs/cern.ch/user/n/nlurkin/Compact/eop.dat");
//	while(fd >> run >> burst >> time){
//		badEventsList.push_back(eventID(run, burst, time));
//	}
//	fd.close();

	filterLoaded = true;
}

bool testAdditionalCondition(ROOTPhysicsEvent *evt, ROOTCorrectedEvent *corrEvent, NGeom *rootGeom, ROOTRawEvent *rawEvent, ROOTBurst *burst, fitStruct &fitBrch){
	TVector3 propPos, propPos2, propPos3;
	if(!filterLoaded) loadEventsFilter();

	if(burst->isData && isFilteredEvent(burst->nrun, burst->time, rawEvent->timeStamp, badEventsList)){
		fitBrch.selEvents--;
		return false;
	}

	//HOD efficiency
//	double eff, eff2;
//	double rndNum, rndNum2;
//	int nOk = 0;
//	TVector3 mom;
//	if(burst->isMC){
//		eff = 0.9976;
//		for(auto it : corrEvent->pTrack){
//			propPos = propagateAfter(11946.50, it, *rawEvent);
//			if (propPos.X()>-0.01 && propPos.X()<0.01) eff=0;
//			if (propPos.Y()>-0.01 && propPos.Y()<0.01) eff=0;
//
//			rndNum = (double)(rand()/(double)RAND_MAX);
//			//cout << rndNum << ">" << eff << endl;
//			if(rndNum<=eff){
//				nOk++;
//			} else{
//				propPos2 = propPos;
//				eff2 = eff;
//				rndNum2 = rndNum;
//				mom = it.momentum;
//			}
//		}
//	}
//	if(nOk==0){
//		fitBrch.selEvents--;
//		return false;
//	}

	/*if(rootBurst->period!=1){
		fitBrch.selEvents--;
		return;
	}*/
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->ep.parentTrack], *rawEvent);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z, corrEvent->pTrack[evt->em.parentTrack], *rawEvent);

	//e+ in square
	if(fabs(propPos.X())<20 && fabs(propPos.Y())<20){
		fitBrch.selEvents--;
		return false;
	}

	//e- in square
	if(fabs(propPos2.X())<20 && fabs(propPos2.Y())<20){
		fitBrch.selEvents--;
		return false;
	}

	if(rawEvent->vtx[corrEvent->goodVertexID].position.Z() < -1700){
		fitBrch.selEvents--;
		//cut_fitBrch.selEvents--;
		return false;
	}

	if(evt->x <= 0.01 || evt->x > 1) {
		fitBrch.selEvents--;
		//cut_fitBrch.selEvents--;
		return false;
	}

	return true;
}
