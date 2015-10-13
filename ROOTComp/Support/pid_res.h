/*
 * pid_res.h
 *
 *  Created on: Aug 14, 2015
 *      Author: nlurkin
 */

#ifndef PID_RES_H_
#define PID_RES_H_

#include <TTree.h>

#include "CompactIO.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

extern ROOTFileHeader &outputFileHeader;

// ### pid result struct
struct alt_pid_res{
	void incTotal(){ total.events[currentID]++;}
	void incPrelim(){ total.prelim.events[currentID]++;}
	void incNoID(bool assoc){
		if(assoc) total.prelim.associated.noID.events[currentID]++;
		else total.prelim.noID.events[currentID]++;
	}
	void incManyID(bool assoc){
		if(assoc) total.prelim.associated.manyID.events[currentID]++;
		else total.prelim.manyID.events[currentID]++;
	}
	void incAssociated(){
		total.prelim.associated.events[currentID]++;
	}
	void incIded(bool assoc){
		if(assoc) total.prelim.associated.ided.events[currentID]++;
		else total.prelim.ided.events[currentID]++;
	}
	void incGoodId(){
		total.prelim.associated.ided.good.events[currentID]++;
	}
	void incBadId(){
		total.prelim.associated.ided.bad.events[currentID]++;
	}
	void incPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.pass.events[currentID]++;
			else if(bad) total.prelim.associated.ided.bad.pass.events[currentID]++;
		}
		else{
			total.prelim.ided.pass.events[currentID]++;
		}
	}
	void incNonPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.noPass.events[currentID]++;
			else if(bad) total.prelim.associated.ided.bad.noPass.events[currentID]++;
		}
		else{
			total.prelim.ided.noPass.events[currentID]++;
		}
	}

	void NewEntry(){
		currentID++;
	}
	void ResetEntry(){
		currentID=0;
	}

	void Init(int N){
		good.resize(N, 0);
		bad.resize(N, 0);
		total.events.resize(N, 0);
		total.prelim.events.resize(N, 0);
		total.prelim.noID.events.resize(N, 0);
		total.prelim.manyID.events.resize(N, 0);
		total.prelim.associated.events.resize(N, 0);
		total.prelim.associated.noID.events.resize(N, 0);
		total.prelim.associated.manyID.events.resize(N, 0);
		total.prelim.associated.ided.events.resize(N, 0);
		total.prelim.associated.ided.good.events.resize(N, 0);
		total.prelim.associated.ided.good.pass.events.resize(N, 0);
		total.prelim.associated.ided.good.noPass.events.resize(N, 0);
		total.prelim.associated.ided.bad.events.resize(N, 0);
		total.prelim.associated.ided.bad.noPass.events.resize(N, 0);
		total.prelim.associated.ided.bad.pass.events.resize(N, 0);
		total.prelim.ided.events.resize(N, 0);
		total.prelim.ided.pass.events.resize(N, 0);
		total.prelim.ided.noPass.events.resize(N, 0);
	}

	int currentID=0;
	//MC association
	std::vector<int> good, bad;

	struct total_t{
		std::vector<int> events;

		struct prelim_t{
			std::vector<int> events;
			struct noID_t{
				std::vector<int> events;
			} noID;
			struct manyID_t{
				std::vector<int> events;
			} manyID;

			struct associated_t{
				std::vector<int> events;
				struct noID_t{
					std::vector<int> events;
				} noID;
				struct manyID_t{
					std::vector<int> events;
				} manyID;

				struct ided_t{
					std::vector<int> events;
					struct good_t{
						std::vector<int> events;
						struct pass_t{
							std::vector<int> events;
						} pass;
						struct noPass_t{
							std::vector<int> events;
						} noPass;
					} good;
					struct bad_t{
						std::vector<int> events;
						struct pass_t{
							std::vector<int> events;
						} pass;
						struct noPass_t{
							std::vector<int> events;
						} noPass;
					} bad;
				} ided;
			} associated;

			struct ided_t{
				std::vector<int> events;
				struct pass_t{
					std::vector<int> events;
				} pass;
				struct noPass_t{
					std::vector<int> events;
				}noPass;
			} ided;
		} prelim;
	} total;
};

std::string printDiv(int a, int b, int c=0, int d=0, int e=0, int f=0);
void printResStruct(struct alt_pid_res &pid_res, struct alt_pid_res *pid_res_2);
void branchResTree(TTree *pid, struct alt_pid_res &pid_res);

#endif /* PID_RES_H_ */
