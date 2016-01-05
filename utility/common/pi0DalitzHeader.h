#include <fstream>
#include <TApplication.h>
#include "ConfigFile.h"
#include "../userinc/funLib.h"
#include "Interface/Sample.h"

using namespace std;

ConfigFile cfg;

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
}

bool testAdditionalCondition(ROOTPhysicsEvent *evt, ROOTCorrectedEvent *corrEvent, NGeom *rootGeom, ROOTRawEvent *rawEvent, fitStruct &fitBrch){
	TVector3 propPos, propPos2, propPos3;

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
