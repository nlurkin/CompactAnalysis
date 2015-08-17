/*
 * pid_res.cpp
 *
 *  Created on: Aug 14, 2015
 *      Author: nlurkin
 */


#include "pid_res.h"
#include <iostream>
#include <sstream>

std::string printDiv(int a, int b, int c, int d, int e, int f){
	std::stringstream s;

	s.precision(3);
	s << std::fixed;
	s << a;
	if(b!=0) s << "\t\t" << a*100./(double)b << "%";
	if(c!=0) s << "\t\t" << a*100./(double)c << "%";
	if(d!=0) s << "\t\t" << a*100./(double)d << "%";
	if(e!=0) s << "\t\t" << a*100./(double)e << "%";
	if(f!=0) s << "\t\t" << a*100./(double)f << "%";
	return s.str();
}

void printResStruct(struct alt_pid_res &pid_res, struct alt_pid_res *pid_res_2){
	std::cout << "Processed events ->\t" << outputFileHeader.NProcessedEvents << std::endl;
	std::cout << "Failed events    ->\t" << outputFileHeader.NFailedEvents << std::endl;
	std::cout << "Passed events    ->\t" << outputFileHeader.NPassedEvents << std::endl;

	std::cout << std::endl << "MC Association (track)" << std::endl << "--------------" << std::endl;
	std::cout << "Good = " << pid_res.good[0] << " " << pid_res.good[0]*100./(double)(pid_res.good[0]+pid_res.bad[0]) << std::endl;
	std::cout << "Bad  = " << pid_res.bad[0] << " " << pid_res.bad[0]*100./(double)(pid_res.good[0]+pid_res.bad[0]) << std::endl;

	std::cout << "\t\t\t\t\t0\t\t1\t\t2\t\t3\t\t4" << std::endl;
	std::cout << "Total: \t\t\t" << printDiv(pid_res.total.events[0], pid_res.total.events[0]) << std::endl;
	std::cout << "Total: \t\t\t" << printDiv(pid_res_2->total.events[0], pid_res_2->total.events[0]) << std::endl;
	std::cout << "|-Prelim: \t\t" << printDiv(pid_res.total.prelim.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << std::endl;
	std::cout << "|-Prelim: \t\t" << printDiv(pid_res_2->total.prelim.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0]) << std::endl;
	std::cout << "  |-NoID: \t\t" << printDiv(pid_res.total.prelim.noID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << std::endl;
	std::cout << "  |-NoID: \t\t" << printDiv(pid_res_2->total.prelim.noID.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0]) << std::endl;
	std::cout << "  |-ManyID: \t\t" << printDiv(pid_res.total.prelim.manyID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << std::endl;
	std::cout << "  |-ManyID: \t\t" << printDiv(pid_res_2->total.prelim.manyID.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0]) << std::endl;
	std::cout << "  |-Ided: \t\t" << printDiv(pid_res.total.prelim.ided.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << std::endl;
	std::cout << "  |-Ided: \t\t" << printDiv(pid_res_2->total.prelim.ided.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.ided.events[0]) << std::endl;
	std::cout << "    |-Pass: \t\t" << printDiv(pid_res.total.prelim.ided.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << std::endl;
	std::cout << "    |-Pass: \t\t" << printDiv(pid_res_2->total.prelim.ided.pass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.ided.events[0]) << std::endl;
	std::cout << "    |-NoPass: \t\t" << printDiv(pid_res.total.prelim.ided.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << std::endl;
	std::cout << "    |-NoPass: \t\t" << printDiv(pid_res_2->total.prelim.ided.noPass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.ided.events[0]) << std::endl;
	std::cout << "  |-Associated: \t" << printDiv(pid_res.total.prelim.associated.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << std::endl;
	std::cout << "  |-Associated: \t" << printDiv(pid_res_2->total.prelim.associated.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0]) << std::endl;
	std::cout << "  | |-NoID: \t\t" << printDiv(pid_res.total.prelim.associated.noID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << std::endl;
	std::cout << "  | |-NoID: \t\t" << printDiv(pid_res_2->total.prelim.associated.noID.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0]) << std::endl;
	std::cout << "  | |-ManyID: \t\t" << printDiv(pid_res.total.prelim.associated.manyID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << std::endl;
	std::cout << "  | |-ManyID: \t\t" << printDiv(pid_res_2->total.prelim.associated.manyID.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0]) << std::endl;
	std::cout << "  | |-Ided: \t\t" << printDiv(pid_res.total.prelim.associated.ided.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0]) << std::endl;
	std::cout << "  | |-Ided: \t\t" << printDiv(pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0]) << std::endl;
	std::cout << "  |   |-Good: \t\t" << printDiv(pid_res.total.prelim.associated.ided.good.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   |-Good: \t\t" << printDiv(pid_res_2->total.prelim.associated.ided.good.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   | |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   | |-Pass: \t" << printDiv(pid_res_2->total.prelim.associated.ided.good.pass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   | |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   | |-NoPass: \t" << printDiv(pid_res_2->total.prelim.associated.ided.good.noPass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.good.events[0]) << std::endl;
	std::cout << "  |   |-Bad: \t\t" << printDiv(pid_res.total.prelim.associated.ided.bad.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << std::endl;
	std::cout << "  |   |-Bad: \t\t" << printDiv(pid_res_2->total.prelim.associated.ided.bad.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.bad.events[0]) << std::endl;
	std::cout << "  |     |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << std::endl;
	std::cout << "  |     |-Pass: \t" << printDiv(pid_res_2->total.prelim.associated.ided.bad.pass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.bad.events[0]) << std::endl;
	std::cout << "  |     |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << std::endl;
	std::cout << "  |     |-NoPass: \t" << printDiv(pid_res_2->total.prelim.associated.ided.bad.noPass.events[0], pid_res_2->total.events[0], pid_res_2->total.prelim.events[0], pid_res_2->total.prelim.associated.events[0], pid_res_2->total.prelim.associated.ided.events[0], pid_res_2->total.prelim.associated.ided.bad.events[0]) << std::endl;
}

void branchResTree(TTree *pid, struct alt_pid_res &pid_res){
	pid->Branch("pid_res.good", &pid_res.good);
	pid->Branch("pid_res.bad", &pid_res.bad);
	pid->Branch("pid_res.total", &pid_res.total.events);
	pid->Branch("pid_res.total.prelim", &pid_res.total.prelim.events);
	pid->Branch("pid_res.total.prelim.associated", &pid_res.total.prelim.associated.events);
	pid->Branch("pid_res.total.prelim.associated.noID", &pid_res.total.prelim.associated.noID.events);
	pid->Branch("pid_res.total.prelim.associated.manyID", &pid_res.total.prelim.associated.manyID.events);
	pid->Branch("pid_res.total.prelim.associated.ided", &pid_res.total.prelim.associated.ided.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good", &pid_res.total.prelim.associated.ided.good.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad", &pid_res.total.prelim.associated.ided.bad.events);
	pid->Branch("pid_res.total.prelim.noID", &pid_res.total.prelim.noID.events);
	pid->Branch("pid_res.total.prelim.manyID", &pid_res.total.prelim.manyID.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good.pass", &pid_res.total.prelim.associated.ided.good.pass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good.noPass", &pid_res.total.prelim.associated.ided.good.noPass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad.pass", &pid_res.total.prelim.associated.ided.bad.pass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad.noPass", &pid_res.total.prelim.associated.ided.bad.noPass.events);
	pid->Branch("pid_res.total.prelim.ided", &pid_res.total.prelim.ided.events);
	pid->Branch("pid_res.total.prelim.ided.pass", &pid_res.total.prelim.ided.pass.events);
	pid->Branch("pid_res.total.prelim.ided.noPass", &pid_res.total.prelim.ided.noPass.events);
}

