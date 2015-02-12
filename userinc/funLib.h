/*
 * funLib.h
 *
 *  Created on: Feb 12, 2015
 *      Author: ncl
 */

#include <iostream>
#include <TVector3.h>
#include <TH1D.h>
#include <vector>
#include "exportClasses.h"
#include <sstream>
#include <map>

#ifndef FUNLIB_H_
#define FUNLIB_H_

//TVector related
std::string 		printVector3		(TVector3 v);
template<typename T>
std::string printSTLvector(std::vector<T> v){
	std::ostringstream ss;
	typename std::vector<T>::iterator it;

	ss << "Vector values :" << std::endl;
	for(it=v.begin(); it!=v.end();it++){
		ss << "\t" << (*it);
		ss << std::endl;
	}
	ss << std::endl;
	return ss.str();
}
std::string printSTLvector(std::vector<TVector3> v);

TH1D* gH(TString name);

double		distance2D			(TVector3 v1, TVector3 v2);
void 		propagateBefore		(float &x, float &y, float &z, float zplane, trak t);
TVector3 	propagateBefore		(float zplane, NPhysicsTrack t);
TVector3 	propagateCorrBefore	(float zplane, NPhysicsTrack t);
void 		propagateAfter		(float &x, float &y, float &z, float zplane, trak t);
TVector3 	propagateAfter		(float zplane, NPhysicsTrack t);
TVector3 	propagate			(float zplane, NPhysicsTrack t);
TVector3 	propagate			(float zplane, TVector3 pos, TVector3 p);

double 	missMass2	(double m1, double m2, TVector3 p1, TVector3 p2);
double 	missMass2Man(double m1, double m2, float kp, float kdxdz, float kdydz, float ep, float edxdz, float edydz);
double	invMass2	(std::vector<double> mass, std::vector<TVector3> p);


void parseCutsValues(std::string fileName);
void printCuts();

std::map<std::string,std::string> parseOptions(std::string s);
int selectOptions(std::string s);
const std::vector<std::string> tokenize(std::string s, const char delim);
bool isFilteredEvent(int nrun, int nburst, int timestamp);
int common_init(string filePrefix);

extern double Mpi0;
extern double Mpic;
extern double Me;

#endif /* FUNLIB_H_ */
