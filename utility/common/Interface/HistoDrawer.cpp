/*
 * HistoDrawer.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: nlurkin
 */

#include "HistoDrawer.h"
#include <TCanvas.h>
#include <iostream>
#include <TList.h>
using namespace std;

HistoDrawer::HistoDrawer() {

}

HistoDrawer::~HistoDrawer() {
	// TODO Auto-generated destructor stub
}

TH1D* HistoDrawer::buildRatio(int nbins, double* binning, TH1D* mc,
		TH1D* data) {
	TH1D* r;
	if (nbins != -1)
		r = new TH1D("ratio", "ratio", nbins - 1, binning);
	else
		r = new TH1D("ratio", "ratio", data->GetNbinsX(),
				data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());

	mc->SetFillColor(8);

	mc->Sumw2();
	r->Sumw2();
	r->Divide(data, mc, 1, 1, "B");
	return r;
}
