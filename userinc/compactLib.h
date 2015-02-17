/*
 * compactLib.h
 *
 *  Created on: Feb 12, 2015
 *      Author: ncl
 */

#ifndef COMPACTLIB_H_
#define COMPACTLIB_H_

#include "reader.h"
#include "user.h"

/*
 * compactLib.c
 *
 *  Created on: Feb 12, 2015
 *      Author: ncl
 */

void loadEOPData(superBurst *sbur);
void applyEOPData();
void openOutput(string outFile, string outPass);

void defineBeamCharge(superBurst *sbur);

void GetCpdCellIndex(double pos_x, double pos_y, int *cpd_index, int *cell_index);

float correctedEP(superCmpEvent* sevt, trak t, float &eOverP);

float correctClusterE(superCmpEvent *sevt, NPhysicsCluster c);

NPhysicsTrack correctTrack(superCmpEvent *sevt, int i);
NPhysicsCluster correctCluster(superCmpEvent* sevt, int i);

void CreateTracks(superCmpEvent *sevt);
void CreateClusters(superCmpEvent *sevt);

int selectOptions(std::string s);

#endif /* COMPACTLIB_H_ */
