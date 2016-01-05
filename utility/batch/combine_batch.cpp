#define __CINT_NICO__ 1

#define PRINTVAR(v) #v << "= " << v << " "

#include <signal.h>
#include <TStyle.h>

#include "pi0DalitzHeader.h"
#include "../common/Combiner.h"

#include <TSystem.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCanvas.h>
#include "../userinc/exportClasses.h"
#include <TGaxis.h>
#include <iomanip>
#include <TList.h>
using namespace std;

#define MAXEVENTS 0
Long_t iCanvas = 0;

//ROOTRawEvent *xxx = new ROOTRawEvent();
//ROOTRawEvent &rawEvent = *xxx;

struct square {
	square(int y1, int y2, int x1, int x2) :
			miny(y1), maxy(y2), minx(x1), maxx(x2) {
	}
	int miny;
	int maxy;
	int minx;
	int maxx;
};

/************************
 * Real job
 ************************/
TH1D* buildRatio(THStack* stack, TH1D* data, TString name) {
	int nbins = data->GetNbinsX();
	double min = data->GetXaxis()->GetXmin();
	double max = data->GetXaxis()->GetXmax();
	TH1D* sum = new TH1D("sum", "sum", nbins, min, max);
	TH1D* r = new TH1D(TString::Format("ratio_%s", name.Data()),
			TString::Format("ratio_%s", name.Data()), nbins, min, max);

	for (int i = 0; i < stack->GetHists()->GetEntries(); ++i) {
		sum->Add((TH1D*) stack->GetHists()->At(i));
	}
	sum->SetFillColor(8);

	sum->Sumw2();
	r->Sumw2();
	r->Divide(data, sum, 1, 1, "B");
	r->SetMarkerColor(kRed);
	//r->SetMaximum(r->GetMaximum()*1.1);
	//r->SetMinimum(r->GetMinimum()*0.9);
	//r->SetMaximum(0.7);
	//r->SetMinimum(1.3);
	delete sum;
	return r;
}

TH2D* buildRatio2(TH2D* mc, TH2D* data, TString name) {
	int nbinsx = data->GetNbinsX();
	double minx = data->GetXaxis()->GetXmin();
	double maxx = data->GetXaxis()->GetXmax();
	int nbinsy = data->GetNbinsY();
	double miny = data->GetYaxis()->GetXmin();
	double maxy = data->GetYaxis()->GetXmax();
	TH2D* r = new TH2D(TString::Format("ratio_%s", name.Data()),
			TString::Format("ratio_%s", name.Data()), nbinsx, minx, maxx,
			nbinsy, miny, maxy);

	r->Sumw2();
	r->Divide(data, mc, 1, 1, "B");
	r->SetMarkerColor(kRed);
	//r->SetMaximum(r->GetMaximum()*1.1);
	//r->SetMinimum(r->GetMinimum()*0.9);
	//r->SetMaximum(0.7);
	//r->SetMinimum(1.3);
	return r;
}

void prepareRatioPlot(TCanvas *c, THStack* mc, TLegend *leg, TH1D* data,
		TH1D* ratio) {
	//double minx = data->GetXaxis()->GetXmin();
	//double maxx = data->GetXaxis()->GetXmax();
	//double binsX = data->GetXaxis()->GetNbins();
	//double miny = data->GetMinimum();
	//double maxy = data->GetMaximum();
	//double binsY = data->GetYaxis()->GetNbins();

//	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
//	pad1->SetBottomMargin(3);
//	pad1->SetGrid();
//	pad1->Draw();             // Draw the upper pad: pad1
//	pad1->cd();               // pad1 becomes the current pad
//	mc->Draw("HIST");               // Draw h1
//	data->Draw("SAME E P");         // Draw h2 on top of h1
//	//mc->GetYaxis()->SetLabelSize(0.05);
//	leg->Draw();
//
//	c->cd();
//	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
//	pad2->SetTopMargin(0);
//	pad2->SetBottomMargin(0.2);
//	pad2->SetGrid(); // vertical grid
//	pad2->Draw();
//	pad2->cd();
//	ratio->SetStats(0);      // No statistics on lower plot
//	ratio->Draw("ep");
//	ratio->SetMarkerColor(kRed);
//	ratio->GetYaxis()->SetRangeUser(0.89, 1.11);
//
//	// Y axis mc plot settings
//	mc->GetYaxis()->SetTitleSize(20);
//	mc->GetYaxis()->SetTitleFont(43);
//	mc->GetYaxis()->SetTitleOffset(1.55);
//
//	// Ratio plot (ratio) settings
//	ratio->SetTitle(""); // Remove the ratio title
//
//	// Y axis ratio plot settings
//	ratio->GetYaxis()->SetTitle("Data/MC");
//	ratio->GetYaxis()->SetNdivisions(505);
//	ratio->GetYaxis()->SetTitleSize(20);
//	ratio->GetYaxis()->SetTitleFont(43);
//	ratio->GetYaxis()->SetTitleOffset(1);
//	ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
//	ratio->GetYaxis()->SetLabelSize(15);
//
//	// X axis ratio plot settings
//	ratio->GetXaxis()->SetTitleSize(20);
//	ratio->GetXaxis()->SetTitleFont(43);
//	ratio->GetXaxis()->SetTitleOffset(4.);
//	ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
//	ratio->GetXaxis()->SetLabelSize(15);
}

TCanvas* drawCanvas(TString name, THStack *stack, TH1D* data, TLegend *leg) {
//	TCanvas *c1 = new TCanvas(TString::Format("c%li", iCanvas), name);
//	TH1D* r = buildRatio(stack, data, name);
//
//	prepareRatioPlot(c1, stack, leg, data, r);
//
//	//c1->SetGrid();
//	//c1->Update();
//	//c1->Draw();
//	++iCanvas;
//	return c1;
}

void doPlot(int index, TString name, TString title, TLegend* leg,
		vector<int> colors, vector<TString> *legendTitle = NULL) {
//	if (index < 0)
//		return;
//	tempFD->cd();
//
//	THStack *hStack = new THStack(name, title);
//
//	//Scale MC to Data
//	double totalMC = 0;
//	for (int i = 0; i < inputMCNbr; ++i) {
//		totalMC += d1->at(i).at(index)->Integral();
//		d1->at(i).at(index)->SetFillColor(gStyle->GetColorPalette(colors[i]));
//		if (legendTitle)
//			leg->AddEntry(d1->at(i).at(index), legendTitle->at(i).Data(), "f");
//	}
//
//	double factor = ((double) (dSig->at(0).at(index)->Integral())) / totalMC;
//	for (int i = 0; i < inputMCNbr; ++i) {
//		d1->at(i).at(index)->Scale(factor);
//	}
//
//	//Stack MC
//	for (unsigned int i = 0; i < d1->size(); ++i) {
//		hStack->Add(d1->at(i).at(index));
//		d1->at(i).at(index)->Write();
//	}
//
//	dSig->at(0).at(index)->Write();
//
//	//Style data
//	dSig->at(0).at(index)->SetLineColor(kRed);
//	if (legendTitle)
//		leg->AddEntry(dSig->at(0).at(index), dataLegendTitle[0].Data(), "lep");
//
//	TCanvas *c = drawCanvas(name, hStack, dSig->at(0).at(index), leg);
//
//	hStack->Write();
//	cout << name + ".png" << endl;
//	c->SaveAs(name + ".png");
//
//	c->Close();
//	delete c;
}

void doPlot2(int index, TString name, TString title, TLegend* leg,
		vector<int> colors, vector<TString> *legendTitle = NULL) {
//	tempFD->cd();
//
//	TH2D *temp = (TH2D*) dMap->at(0).at(index)->Clone(name);
//	temp->SetTitle(title);
//	temp->Clear();
//
//	//Scale MC to Data
//	double totalMC = 0;
//	for (int i = 0; i < inputMCNbr; ++i) {
//		totalMC += dMap->at(i).at(index)->Integral();
//		dMap->at(i).at(index)->SetFillColor(gStyle->GetColorPalette(colors[i]));
//		if (legendTitle)
//			leg->AddEntry(dMap->at(i).at(index), legendTitle->at(i).Data(),
//					"f");
//	}
//
//	double factor = ((double) (dSig->at(0).at(0)->Integral())) / totalMC;
//	for (int i = 0; i < inputMCNbr; ++i) {
//		dMap->at(i).at(index)->Scale(factor);
//	}
//
//	//Stack MC
//	for (unsigned int i = 0; i < dMap->size(); ++i) {
//		//hStack->Add(dMap->at(i).at(index));
//		temp->Add(dMap->at(i).at(index), 1);
//		dMap->at(i).at(index)->Write();
//	}
//
//	int nbinsy = 1;	//8;
//	int nbinsx = 1;	//temp->GetXaxis()->GetNbins()/5;
//	if (nbinsx <= 0)
//		nbinsx = 1;
//	if (nbinsy <= 0)
//		nbinsy = 1;
//	temp = (TH2D*) temp->Rebin2D(nbinsx, nbinsy);
//	TH2D* tempSig = (TH2D*) dSigMap->at(0).at(index)->Rebin2D(nbinsx, nbinsy);
//
//	TH2D* ratio = buildRatio2(temp, tempSig, name);
//	ratio->GetZaxis()->SetRangeUser(0.85, 1.15);
//
//	TCanvas *c = new TCanvas(TString::Format("c%li", iCanvas), name);
//	c->Divide(2, 2);
//	c->cd(1);
//	temp->Draw("COLZ");
//	c->cd(2);
//	tempSig->Draw("COLZ");
//	c->cd(3);
//	ratio->Draw("colz");
//	++iCanvas;
//
//	temp->Write();
//	cout << name + ".png" << endl;
//	c->SaveAs(name + ".png");
//	//c->Close();
//	//delete c;
}

/*****************************
 * Mains
 *****************************/
void combine_batch() {
	gStyle->SetOptFit(1);

	if (!cfg.testAllOutputFiles())
		return;

	Combiner combine;

	RunWeights weights;
	weights.loadWeights(
			"/afs/cern.ch/user/n/nlurkin/Compact/pi0dalitz_weights.dat");
	combine.setRunWeights(&weights);

	combine.prepareSamples<CombineMCSample, CombineDataSample>(cfg);
	combine.fillSamples();
}

void combine_show(int firstPlot, int maxPlots) {
	gStyle->SetOptFit(1);

	Combiner combine;

	RunWeights weights;
	weights.loadWeights(
			"/afs/cern.ch/user/n/nlurkin/Compact/pi0dalitz_weights.dat");
	combine.setRunWeights(&weights);

	combine.prepareSamples<CombineMCSample, CombineDataSample>(cfg);
	combine.getSamples();

	combine.mergeSamples<CombineMCSample, CombineDataSample>();
	combine.draw(cfg.getMcColors(), cfg.getDataColors());

//	TLegend *leg = new TLegend(.7, 0.75, 0.9, 0.95);
//
//	int i = firstPlot;
//
//	doPlot(++i, "mK", "Kaon invariant mass", leg, mcColors, &mcLegendTitle);
//	if (maxPlots-- == 0)
//		return;
//
//	doPlot(++i, "R_DCH1_ep", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH1_ep", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH1_ep", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH1_em", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH1_em", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH1_em", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH1_pip", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH1_pip", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH1_pip", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH1_gamma", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH1_gamma", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH1_gamma", "DCH1 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//
//	doPlot(++i, "R_DCH2_ep", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH2_ep", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH2_ep", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH2_em", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH2_em", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH2_em", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH2_pip", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH2_pip", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH2_pip", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH2_gamma", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH2_gamma", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH2_gamma", "DCH2 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//
//	doPlot(++i, "R_DCH3_ep", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH3_ep", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH3_ep", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH3_em", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH3_em", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH3_em", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH3_pip", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH3_pip", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH3_pip", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH3_gamma", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH3_gamma", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH3_gamma", "DCH3 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//
//	doPlot(++i, "R_DCH4_ep", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH4_ep", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH4_ep", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH4_em", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH4_em", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH4_em", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH4_pip", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH4_pip", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH4_pip", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "R_DCH4_gamma", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "X_DCH4_gamma", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//	doPlot(++i, "Y_DCH4_gamma", "DCH4 Radius", leg, mcColors);
//	if (maxPlots-- == 0)
//		return;
//
//	doPlot(++i, "Zvtx", "Vertex Z", leg, mcColors);
////	doPlot(++i, "Zvtx_low", "Vertex Z", leg, mcColors);
////	doPlot(++i, "Zvtx_high", "Vertex Z", leg, mcColors);
//	doPlot(++i, "Qvtx", "Vertex Charge", leg, mcColors);
//	doPlot(++i, "CDAvtx", "Vertex CDA", leg, mcColors);
//	doPlot(++i, "Pt2", "Square transverse momentum", leg, mcColors);
//	doPlot(++i, "P", "Total momentum", leg, mcColors);
//
//	doPlot(++i, "Mpi0", "Reconstructed Pi0 mass", leg, mcColors);
//
//	//Photon
//	doPlot(++i, "gEnergy", "Photon energy", leg, mcColors);
//	doPlot(++i, "gPositionX", "Photon LKr position (X)", leg, mcColors);
//	doPlot(++i, "gPositionY", "Photon LKr position (Y)", leg, mcColors);
//	doPlot(++i, "gRadius", "Photon LKr radius", leg, mcColors);
//	doPlot(++i, "gP", "Photon momentum", leg, mcColors);
//
//	//e+/e-
//	doPlot(++i, "epPMag", "Electron momentum", leg, mcColors);
//	doPlot(++i, "epPx", "Electron momentum (X)", leg, mcColors);
//	doPlot(++i, "epPy", "Electron momentum (Y)", leg, mcColors);
//	doPlot(++i, "epPz", "Electron momentum (Z)", leg, mcColors);
//	doPlot(++i, "epEnergy", "Electron energy", leg, mcColors);
//	doPlot(++i, "epeop", "Electron E/p", leg, mcColors);
//	doPlot(++i, "epLKrX", "Electron LKr position (X)", leg, mcColors);
//	doPlot(++i, "epLKrY", "Electron LKr position (Y)", leg, mcColors);
//	doPlot(++i, "epLKrR", "Electron LKr radius", leg, mcColors);
//
//	doPlot(++i, "emPMag", "Electron momentum", leg, mcColors);
//	doPlot(++i, "emPx", "Electron momentum (X)", leg, mcColors);
//	doPlot(++i, "emPy", "Electron momentum (Y)", leg, mcColors);
//	doPlot(++i, "emPz", "Electron momentum (Z)", leg, mcColors);
//	doPlot(++i, "emEnergy", "Electron energy", leg, mcColors);
//	doPlot(++i, "emeop", "Electron E/p", leg, mcColors);
//	doPlot(++i, "emLKrX", "Electron LKr position (X)", leg, mcColors);
//	doPlot(++i, "emLKrY", "Electron LKr position (Y)", leg, mcColors);
//	doPlot(++i, "emLKrR", "Electron LKr radius", leg, mcColors);
//
//	doPlot(++i, "mee", "Di-electron invariant mass", leg, mcColors);
//
//	//pi+
//	doPlot(++i, "pipPMag", "Pion momentum", leg, mcColors);
//	doPlot(++i, "pipPx", "Pion momentum (X)", leg, mcColors);
//	doPlot(++i, "pipPy", "Pion momentum (Y)", leg, mcColors);
//	doPlot(++i, "pipPz", "Pion momentum (Z)", leg, mcColors);
//	doPlot(++i, "pipEnergy", "Pion energy", leg, mcColors);
//	doPlot(++i, "pieop", "Pion E/p", leg, mcColors);
//
//	doPlot(++i, "t_epem_DCH", "Track distance DCH1", leg, mcColors);
//	doPlot(++i, "t_eppip_DCH", "Track distance DCH1", leg, mcColors);
//	doPlot(++i, "t_empip_DCH", "Track distance DCH1", leg, mcColors);
//	doPlot(++i, "t_epem_LKr", "Track distance LKr", leg, mcColors);
//	doPlot(++i, "t_eppip_LKr", "Track distance LKr", leg, mcColors);
//	doPlot(++i, "t_empip_LKr", "Track distance LKr", leg, mcColors);
//
//	doPlot(++i, "t_gep_DCH", "Track photon distance DCH1", leg, mcColors);
//	doPlot(++i, "t_gem_DCH", "Track photon distance DCH1", leg, mcColors);
//	doPlot(++i, "t_gpip_DCH", "Track photon distance DCH1", leg, mcColors);
//	doPlot(++i, "t_gep_LKr", "Track photon distance LKr", leg, mcColors);
//	doPlot(++i, "t_gem_LKr", "Track photon distance LKr", leg, mcColors);
//	doPlot(++i, "t_gpip_LKr", "Track photon distance LKr", leg, mcColors);
//
//	doPlot(++i, "undeft_gep_LKr", "Undeflected track photon distance LKr", leg,
//			mcColors);
//	doPlot(++i, "undeft_gem_LKr", "Undeflected track photon distance LKr", leg,
//			mcColors);
//	doPlot(++i, "undeft_gpip_LKr", "Undeflected track photon distance LKr", leg,
//			mcColors);
//
//	doPlot(++i, "L3_E_LKr_ep", "L3 Electron energy", leg, mcColors);
//	doPlot(++i, "L3_E_LKr_em", "L3 Electron energy", leg, mcColors);
//	doPlot(++i, "L3_E_LKr_gamma", "L3 photon energy", leg, mcColors);
//	doPlot(++i, "L3_E_LKr", "L3 energy", leg, mcColors);
//
//	doPlot(++i, "R_DCH1_ep_0", "R_DCH1_ep_0", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_0", "X_DCH1_ep_0", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_0", "Y_DCH1_ep_0", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_0", "Y_DCH1_ep_0", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_1", "R_DCH1_ep_1", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_1", "X_DCH1_ep_1", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_1", "Y_DCH1_ep_1", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_1", "Y_DCH1_ep_1", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_2", "R_DCH1_ep_2", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_2", "X_DCH1_ep_2", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_2", "Y_DCH1_ep_2", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_2", "Y_DCH1_ep_2", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_3", "R_DCH1_ep_3", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_3", "X_DCH1_ep_3", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_3", "Y_DCH1_ep_3", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_3", "Y_DCH1_ep_3", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_4", "R_DCH1_ep_4", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_4", "X_DCH1_ep_4", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_4", "Y_DCH1_ep_4", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_4", "Y_DCH1_ep_4", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_5", "R_DCH1_ep_5", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_5", "X_DCH1_ep_5", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_5", "Y_DCH1_ep_5", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_5", "Y_DCH1_ep_5", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_6", "R_DCH1_ep_6", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_6", "X_DCH1_ep_6", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_6", "Y_DCH1_ep_6", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_6", "Y_DCH1_ep_6", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_7", "R_DCH1_ep_7", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_7", "X_DCH1_ep_7", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_7", "Y_DCH1_ep_7", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_7", "Y_DCH1_ep_7", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_8", "R_DCH1_ep_8", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_8", "X_DCH1_ep_8", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_8", "Y_DCH1_ep_8", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_8", "Y_DCH1_ep_8", leg, mcColors);
//	doPlot(++i, "R_DCH1_ep_9", "R_DCH1_ep_9", leg, mcColors);
//	doPlot(++i, "X_DCH1_ep_9", "X_DCH1_ep_9", leg, mcColors);
//	doPlot(++i, "Y_DCH1_ep_9", "Y_DCH1_ep_9", leg, mcColors);
//	doPlot(++i, "t_epem_DCH1_9", "Y_DCH1_ep_9", leg, mcColors);
//
//	int iMap = -1;
//	//doPlot2(++iMap, "xMap", "x_reco vs. x_true", leg, mcColors);
//	++iMap;
//
//	doPlot2(++iMap, "LKr_XY_ep", "Electron LKr map", leg, mcColors);
//	doPlot2(++iMap, "LKr_XY_em", "Electron LKr map", leg, mcColors);
//	doPlot2(++iMap, "LKr_XY_pip", "Pion LKr map", leg, mcColors);
//	//iMap+=3;
//	++iMap;
//	//doPlot2(++iMap, "LKr_XY_gamma", "Photon LKr map", leg, mcColors);
//	doPlot2(++iMap, "DCH1_XY_ep", "Electron DCH1 map", leg, mcColors);
//	doPlot2(++iMap, "DCH1_XY_em", "Electron DCH1 map", leg, mcColors);
//	doPlot2(++iMap, "DCH1_XY_pip", "Pion DCH1 map", leg, mcColors);
//	++iMap;
//	//doPlot2(++iMap, "DCH1_XY_gamma", "Photon DCH1 map", leg, mcColors);
//	doPlot2(++iMap, "DCH2_XY_ep", "Electron DCH2 map", leg, mcColors);
//	doPlot2(++iMap, "DCH2_XY_em", "Electron DCH2 map", leg, mcColors);
//	doPlot2(++iMap, "DCH2_XY_pip", "Pion DCH2 map", leg, mcColors);
//	++iMap;
//	//doPlot2(++iMap, "DCH2_XY_gamma", "Photon DCH2 map", leg, mcColors);
//	doPlot2(++iMap, "DCH3_XY_ep", "Electron DCH3 map", leg, mcColors);
//	doPlot2(++iMap, "DCH3_XY_em", "Electron DCH3 map", leg, mcColors);
//	doPlot2(++iMap, "DCH3_XY_pip", "Pion DCH3 map", leg, mcColors);
//	++iMap;
//	//doPlot2(++iMap, "DCH3_XY_gamma", "Photon DCH3 map", leg, mcColors);
//	doPlot2(++iMap, "DCH4_XY_ep", "Electron DCH4 map", leg, mcColors);
//	doPlot2(++iMap, "DCH4_XY_em", "Electron DCH4 map", leg, mcColors);
//	doPlot2(++iMap, "DCH4_XY_pip", "Pion DCH4 map", leg, mcColors);
//	++iMap;
//	//doPlot2(++iMap, "DCH4_XY_gamma", "Photon DCH4 map", leg, mcColors);

}

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "Missing parameter" << endl;
		return -1;
	}

//	setStyle();

	signal(SIGTERM, sighandler);
	signal(SIGINT, sighandler);
	signal(SIGABRT, sighandler);

	int firstPlots = 1;
	int maxPlots = -1;

	if (!cfg.readFile(argv[1]))
		return -1;
	cfg.print();
	if (argc == 2)
		combine_batch();
	else {
		if (argc >= 4)
			firstPlots = atoi(argv[3]);
		if (argc >= 5)
			maxPlots = atoi(argv[4]);
		theApp = new TApplication("combine", &argc, argv);
		combine_show(-firstPlots, maxPlots - 1);
		theApp->Run();
	}
}
