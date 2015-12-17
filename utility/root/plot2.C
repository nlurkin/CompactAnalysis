#include <TGraphErrors.h>
#include <TColor.h>
#include <TROOT.h>
#include <TStyle.h>

#include <vector>
using namespace std;

// This is the file rootlogon.C
void setStyle(){

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

  // from ROOT plain stylei
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleColor(0);
  myStyle->SetStatColor(0);
  //myStyle->SetPalette(1);
  myStyle->SetLabelSize(0.025,"xyz"); // size of axis values

  // default canvas positioning
  myStyle->SetCanvasDefX(900);
  myStyle->SetCanvasDefY(20);
  myStyle->SetCanvasDefH(550);
  myStyle->SetCanvasDefW(540);
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadLeftMargin(0.125);
  myStyle->SetPadRightMargin(0.125);

  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);

  myStyle->SetFrameBorderMode(0);


  //Nice colors
  const Int_t NRGBs = 5;
  const Int_t NCont = 999;
  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  myStyle->SetNumberContours(NCont);

  gROOT->SetStyle("MyStyle"); //uncomment to set this style

  TCanvas *c = new TCanvas("c","Marker colors",0,0,500,200);
  c->DrawColorTable();
}


void plotpt2(){
	//setStyle();
	double pointsOld	[11] = 		{0.0365678,0.0365678,0.0343745,0.0343745,0.0331968,0.0331968,0.0322055,0.0322055,0.0321715,0.0321715,0.0320711};
	double errOld		[11] = 		{0.00793701,0.00793701,0.00687002,0.00687002,0.00660519,0.00660519,0.00657244,0.00657244,0.00656821,0.00656821,0.00656661};

	double pointsOldRoot[11] = 	{0.0384388,0.0384388,0.0356803,0.0356803,0.0342473,0.0342473,0.0332669,0.0332669,0.0332091,0.0332091,0.0331162};
	double errOldRoot	[11] = 	{0.00756655,0.00756655,0.00766112,0.00766112,0.00771578,0.00771578,0.00772124,0.00772124,0.00772294,0.00772294,0.0077232};

	double pointsNew	[11] = 		{0.040984,0.040984,0.0375896,0.0375896,0.0362318,0.0362318,0.0352833,0.0352833,0.035254,0.035254,0.0351612};
	double errNew		[11] = 		{0.00765139,0.00765139,0.00774394,0.00774394,0.00779874,0.00779874,0.00780497,0.00780497,0.00780531,0.00780531,0.00780551};

	double pointsNewRoot[11] = 	{0.0408873,0.0408873,0.0375794,0.0375794,0.0361593,0.0361593,0.0351997,0.0351997,0.0351612,0.0351612,0.0350738};
	double errNewRoot	[11] = 	{0.00757255,0.00757255,0.00766464,0.00766464,0.00772138,0.00772138,0.00772708,0.00772708,0.00772828,0.00772828,0.00772844};

	double pointsNEvents[11] = 	{424190,424190,356547,356547,335555,335555,332869,332869,332512,332512,332383};

	//Creation
	TGraphErrors *pt2Old = new TGraphErrors();
	pt2Old->SetName("pt2Old");
	//Points
	for(int i=0; i<10; ++i){
		pt2Old->SetPoint(i, 0.0001*(i+1), pointsOld[i]);
	}
	//Errors
	for(int i=0; i<10; ++i){
		pt2Old->SetPointError(i, 0, errOld[i]);
	}
	//Style
	pt2Old->SetMarkerStyle(20);
	pt2Old->SetMarkerColor(2);
	pt2Old->SetLineColor(2);
	pt2Old->SetFillStyle(0);


	//Creation
	TGraphErrors *pt2OldRoot = new TGraphErrors();
	pt2OldRoot->SetName("pt2OldRoot");
	//Points
	for(int i=0; i<10; ++i){
		pt2OldRoot->SetPoint(i, 0.0001*(i+1)+0.00002, pointsOldRoot[i]);
	}
	//Errors
	for(int i=0; i<10; ++i){
		pt2OldRoot->SetPointError(i, 0, errOldRoot[i]);
	}
	//Style
	pt2OldRoot->SetMarkerStyle(21);
	pt2OldRoot->SetMarkerColor(3);
	pt2OldRoot->SetLineColor(3);
	pt2OldRoot->SetFillStyle(0);


	//Creation
	TGraphErrors *pt2New = new TGraphErrors();
	pt2New->SetName("pt2New");
	//Points
	for(int i=0; i<10; ++i){
		pt2New->SetPoint(i, 0.0001*(i+1)+0.00004, pointsNew[i]);
	}
	//Errors
	for(int i=0; i<10; ++i){
		pt2New->SetPointError(i, 0, errNew[i]);
	}
	//Style
	pt2New->SetMarkerStyle(22);
	pt2New->SetMarkerColor(4);
	pt2New->SetLineColor(4);
	pt2New->SetFillStyle(0);


	//Creation
	TGraphErrors *pt2NewRoot = new TGraphErrors();
	pt2NewRoot->SetName("pt2NewRoot");
	//Points
	for(int i=0; i<10; ++i){
		pt2NewRoot->SetPoint(i, 0.0001*(i+1)+0.00006, pointsNewRoot[i]);
	}
	//Errors
	for(int i=0; i<10; ++i){
		pt2NewRoot->SetPointError(i, 0, errNewRoot[i]);
	}
	//Style
	pt2NewRoot->SetMarkerStyle(23);
	pt2NewRoot->SetMarkerColor(6);
	pt2NewRoot->SetLineColor(6);
	pt2NewRoot->SetFillStyle(0);

	//Creation
	TGraph *pt2Selected = new TGraph();
	pt2Selected->SetName("pt2Selected");
	//Points
	for(int i=0; i<10; ++i){
		pt2Selected->SetPoint(i, 0.0001*(i+1), pointsNEvents[i]);
	}
	//Style
	pt2Selected->SetMarkerStyle(20);
	pt2Selected->SetMarkerColor(4);
	pt2Selected->SetLineColor(4);
	pt2Selected->SetFillStyle(0);

	//Plotting
	TCanvas *c = new TCanvas("plot");
	c->SetGrid(1, 1);

	TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.95);
	leg->AddEntry(pt2Old, "Pt^{2} Old Method");
	leg->AddEntry(pt2OldRoot, "Pt^{2} Old Method with ROOT #chi^{2}");
	leg->AddEntry(pt2New, "Pt^{2} New Method");
	leg->AddEntry(pt2NewRoot, "Pt^{2} New Method with ROOT #chi^{2}");

	pt2Old->GetYaxis()->SetRangeUser(0.02, 0.07);
	pt2Old->SetTitle("Pt^{2} scan");
	pt2Old->GetXaxis()->SetTitle("Pt^{2} cut");
	pt2Old->GetYaxis()->SetTitle("FF Slope a");
	pt2Old->GetYaxis()->SetTitleOffset(1.5);
	pt2Old->Draw("AP");
	pt2OldRoot->Draw("SAME P");
	pt2New->Draw("SAME P");
	pt2NewRoot->Draw("SAME P");
	leg->Draw("LEP");

	//Plotting
	TCanvas *c = new TCanvas("Selected");
	c->SetGrid(1, 1);

	pt2Selected->SetTitle("Pt^{2} scan");
	pt2Selected->GetXaxis()->SetTitle("Pt^{2} cut");
	pt2Selected->GetYaxis()->SetTitle("Selected events");
	pt2Selected->GetYaxis()->SetTitleOffset(1.5);
	pt2Selected->Draw("APL");
}

void plotBins(){
	//setStyle();
	int nPoints = 6;
	double x[6] = {11, 25, 55, 125, 275, 1375};
	double pointsOld[6] = {0.0489317, 0.0477602, 0.0453805, 0.0435394, 0.0435394, 0.0435394};
	double errOld[6] = {0.0117164, 0.0108575, 0.0104976, 0.010366, 0.010366, 0.010366};

	double pointsOldRoot[6] = {0.0489442, 0.0477102, 0.0451809, 0.042954, 0.042954, 0.042954};
	double errOldRoot[6] = {0.0147777, 0.013616, 0.0131078, 0.0128983, 0.0128983, 0.0128983};

	double pointsNew[6] = {0.0489511, 0.0485804, 0.0460614, 0.0445144, 0.0445144, 0.0445144};
	double errNew[6] = {0.0147983, 0.0136198, 0.0130741, 0.0128599, 0.0128599, 0.0128599};

	double pointsNewRoot[6] = {0.0489394, 0.0485061, 0.0460815, 0.0442785, 0.0442785, 0.0442785};
	double errNewRoot[6] = {0.0147336, 0.0135508, 0.013022, 0.0128102, 0.0128102, 0.0128102};

	//Creation
	TGraphErrors *binsOld = new TGraphErrors();
	binsOld->SetName("binsOld");
	//Points
	for(int i=0; i<nPoints; ++i){
		binsOld->SetPoint(i, x[i], pointsOld[i]);
	}
	//Errors
	for(int i=0; i<nPoints; ++i){
		binsOld->SetPointError(i, 0, errOld[i]);
	}
	//Style
	binsOld->SetMarkerStyle(20);
	binsOld->SetMarkerColor(2);
	binsOld->SetLineColor(2);
	binsOld->SetFillStyle(0);


	//Creation
	TGraphErrors *binsOldRoot = new TGraphErrors();
	binsOldRoot->SetName("binsOldRoot");
	//Points
	for(int i=0; i<nPoints; ++i){
		binsOldRoot->SetPoint(i, x[i]+2, pointsOldRoot[i]);
	}
	//Errors
	for(int i=0; i<nPoints; ++i){
		binsOldRoot->SetPointError(i, 0, errOldRoot[i]);
	}
	//Style
	binsOldRoot->SetMarkerStyle(21);
	binsOldRoot->SetMarkerColor(3);
	binsOldRoot->SetLineColor(3);
	binsOldRoot->SetFillStyle(0);


	//Creation
	TGraphErrors *binsNew = new TGraphErrors();
	binsNew->SetName("binsNew");
	//Points
	for(int i=0; i<nPoints; ++i){
		binsNew->SetPoint(i, x[i]+4, pointsNew[i]);
	}
	//Errors
	for(int i=0; i<nPoints; ++i){
		binsNew->SetPointError(i, 0, errNew[i]);
	}
	//Style
	binsNew->SetMarkerStyle(22);
	binsNew->SetMarkerColor(4);
	binsNew->SetLineColor(4);
	binsNew->SetFillStyle(0);


	//Creation
	TGraphErrors *binsNewRoot = new TGraphErrors();
	binsNewRoot->SetName("binsNewRoot");
	//Points
	for(int i=0; i<nPoints; ++i){
		binsNewRoot->SetPoint(i, x[i]+6, pointsNewRoot[i]);
	}
	//Errors
	for(int i=0; i<nPoints; ++i){
		binsNewRoot->SetPointError(i, 0, errNewRoot[i]);
	}
	//Style
	binsNewRoot->SetMarkerStyle(23);
	binsNewRoot->SetMarkerColor(6);
	binsNewRoot->SetLineColor(6);
	binsNewRoot->SetFillStyle(0);

	//Plotting
	TCanvas *c = new TCanvas("plotBins");
	c->SetGrid(1, 1);

	TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.95);
	leg->AddEntry(binsOld, "Old Method");
	leg->AddEntry(binsOldRoot, "Old Method with ROOT #chi^{2}");
	leg->AddEntry(binsNew, "New Method");
	leg->AddEntry(binsNewRoot, "New Method with ROOT #chi^{2}");

	binsOld->GetYaxis()->SetRangeUser(0.02, 0.07);
	binsOld->SetTitle("Binning scan");
	binsOld->GetXaxis()->SetTitle("Number of bins");
	binsOld->GetYaxis()->SetTitle("FF Slope a");
	binsOld->GetYaxis()->SetTitleOffset(1.5);
	binsOld->Draw("AP");
	binsOldRoot->Draw("SAME P");
	binsNew->Draw("SAME P");
	binsNewRoot->Draw("SAME P");
	leg->Draw("LEP");
}


void plotpi0mass(){
	//setStyle();
	int NPoints = 10;
	double points[10] = {0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011};

	double pointsOld[10] = {0.0257525, 0.0299572, 0.0347466, 0.0410449, 0.0436154, 0.0455023, 0.0467107, 0.0472451, 0.048017, 0.0484429};
	double errOld[10] = {0.012046, 0.011257, 0.0109968, 0.0108846, 0.0108221, 0.0107851, 0.0107605, 0.0107444, 0.0107404, 0.0107241};

	double pointsOldRoot[10] = {0.0256596, 0.029835, 0.0345065, 0.0406962, 0.0432381, 0.0450902, 0.0462382, 0.0467681, 0.0474945, 0.0479093};
	double errOldRoot[10] = {0.0149477, 0.0139699, 0.0136587, 0.0135292, 0.0134542, 0.0134128, 0.0133885, 0.0133693, 0.0133716, 0.01335};

	double pointsNew[10] = {0.0259345, 0.0311638, 0.0359116, 0.0423887, 0.0449354, 0.0466156, 0.0477873, 0.0481367, 0.0488021, 0.0491589};
	double errNew[10] = {0.0149385, 0.0139718, 0.0136461, 0.0135038, 0.0134258, 0.0133747, 0.0133404, 0.0133222, 0.0133192, 0.0132983};

	double pointsNewRoot[10] = {0.0257485, 0.0309926, 0.035792, 0.0423119, 0.0448136, 0.0465055, 0.0476842, 0.0480334, 0.0487064, 0.0490642};
	double errNewRoot[10] = {0.0148516, 0.0138768, 0.0135627, 0.0134277, 0.0133522, 0.0133057 ,0.0132792 ,0.0132628 ,0.0132658, 0.0132472};

	double pointsNEvents[10] = {756544, 874179, 919446, 940313, 951571, 958237, 962005, 964212, 963751, 966306};

	//Creation
/*	TGraphErrors *pi0massOld = new TGraphErrors();
	pi0massOld->SetName("pi0massOld");
	//Points
	for(int i=0; i<NPoints; ++i){
		pi0massOld->SetPoint(i, points[i], pointsOld[i]);
	}
	//Errors
	for(int i=0; i<NPoints; ++i){
		pi0massOld->SetPointError(i, 0, errOld[i]);
	}
	//Style
	pi0massOld->SetMarkerStyle(20);
	pi0massOld->SetMarkerColor(2);
	pi0massOld->SetLineColor(2);
	pi0massOld->SetFillStyle(0);*/


	//Creation
	TGraphErrors *pi0massOldRoot = new TGraphErrors();
	pi0massOldRoot->SetName("pi0massOldRoot");
	//Points
	for(int i=0; i<NPoints; ++i){
		pi0massOldRoot->SetPoint(i, points[i]+0.00002, pointsOldRoot[i]);
	}
	//Errors
	for(int i=0; i<NPoints; ++i){
		pi0massOldRoot->SetPointError(i, 0, errOldRoot[i]);
	}
	//Style
	pi0massOldRoot->SetMarkerStyle(21);
	pi0massOldRoot->SetMarkerColor(3);
	pi0massOldRoot->SetLineColor(3);
	pi0massOldRoot->SetFillStyle(0);


	//Creation
	TGraphErrors *pi0massNew = new TGraphErrors();
	pi0massNew->SetName("pi0massNew");
	//Points
	for(int i=0; i<NPoints; ++i){
		pi0massNew->SetPoint(i, points[i]+0.00004, pointsNew[i]);
	}
	//Errors
	for(int i=0; i<NPoints; ++i){
		pi0massNew->SetPointError(i, 0, errNew[i]);
	}
	//Style
	pi0massNew->SetMarkerStyle(22);
	pi0massNew->SetMarkerColor(4);
	pi0massNew->SetLineColor(4);
	pi0massNew->SetFillStyle(0);

	//Creation
	TGraphErrors *pi0massNewRoot = new TGraphErrors();
	pi0massNewRoot->SetName("pi0massNewRoot");
	//Points
	for(int i=0; i<NPoints; ++i){
		pi0massNewRoot->SetPoint(i, points[i]+0.00006, pointsNewRoot[i]);
	}
	//Errors
	for(int i=0; i<NPoints; ++i){
		pi0massNewRoot->SetPointError(i, 0, errNewRoot[i]);
	}
	//Style
	pi0massNewRoot->SetMarkerStyle(23);
	pi0massNewRoot->SetMarkerColor(6);
	pi0massNewRoot->SetLineColor(6);
	pi0massNewRoot->SetFillStyle(0);

	//Creation
	TGraph *pi0massSelected = new TGraph();
	pi0massSelected->SetName("pi0massSelected");
	//Points
	for(int i=0; i<10; ++i){
		pi0massSelected->SetPoint(i, points[i], pointsNEvents[i]);
	}
	//Style
	pi0massSelected->SetMarkerStyle(20);
	pi0massSelected->SetMarkerColor(4);
	pi0massSelected->SetLineColor(4);
	pi0massSelected->SetFillStyle(0);


	//Plotting
	TCanvas *c = new TCanvas("mPi0plot");
	c->SetGrid(1, 1);

	TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.95);
//	leg->AddEntry(pi0massOld, "m_{ee#gamma} Old Method");
	leg->AddEntry(pi0massOldRoot, "m_{ee#gamma} Old Method with ROOT #chi^{2}");
	leg->AddEntry(pi0massNew, "m_{ee#gamma} New Method");
	leg->AddEntry(pi0massNewRoot, "m_{ee#gamma} New Method with ROOT #chi^{2}");

	pi0massOldRoot->GetYaxis()->SetRangeUser(0, 0.08);
	pi0massOldRoot->SetTitle("m_{ee#gamma} scan");
	pi0massOldRoot->GetXaxis()->SetTitle("m_{ee#gamma} cut");
	pi0massOldRoot->GetYaxis()->SetTitle("FF Slope a");
	pi0massOldRoot->GetYaxis()->SetTitleOffset(1.5);
//	pi0massOld->Draw("AP");
//	pi0massOldRoot->Draw("SAME P");
	pi0massOldRoot->Draw("AP");
	pi0massNew->Draw("SAME P");
	pi0massNewRoot->Draw("SAME P");
	leg->Draw("LEP");

	//Plotting
	TCanvas *c = new TCanvas("mPi0Selected");
	c->SetGrid(1, 1);

	pi0massSelected->SetTitle("m_{ee#gamma} scan");
	pi0massSelected->GetXaxis()->SetTitle("m_{ee#gamma} cut");
	pi0massSelected->GetYaxis()->SetTitle("Selected events");
	pi0massSelected->GetYaxis()->SetTitleOffset(1.5);
	pi0massSelected->Draw("APL");

}

void plot_only1(){
	int npoints = 		  9;
	double cutValue		 [9] = 	{-0.5,0
,1
,2
,3
,4
,5
,6
,7};

	double pointsNewRoot [9] = 	{1000,3.52929
,3.5285
,3.62514
,3.72007
,3.82021
,3.79847
,3.77486
,3.74425};

	double errNewRoot	 [9] = 	{0,0.525758
,0.525757
,0.530455
,0.546722
,0.571993
,0.606369
,0.651326
,0.707002};

	double pointsNEvents [9] = 	{0,1075340
,1075340
,1053280
,977979
,875727
,761218
,642452
,530199};

	double errUncorrRoot [9] = 	{0,0.070434597
,0.070442061
,0
,0.132372347
,0.213994124
,0.293770075
,0.377945831
,0.467407019};


	//Creation
	TGraphErrors *pt2NewRoot = new TGraphErrors();
	pt2NewRoot->SetName("pt2NewRoot");
	TGraphErrors *pt2NewRootUncorr = new TGraphErrors();
	pt2NewRootUncorr->SetName("pt2NewRootUncorr");
	//Points
	for(int i=0; i<npoints; ++i){
		pt2NewRoot->SetPoint(i, cutValue[i], pointsNewRoot[i]/100.);
		pt2NewRootUncorr->SetPoint(i, cutValue[i], pointsNewRoot[i]/100.);
	}
	//Errors
	for(int i=0; i<npoints; ++i){
		pt2NewRoot->SetPointError(i, 0, errNewRoot[i]/100.);
		pt2NewRootUncorr->SetPointError(i, 0, errUncorrRoot[i]/100.);
	}
	//Style
	pt2NewRoot->SetMarkerStyle(23);
	pt2NewRoot->SetMarkerColor(4);
	pt2NewRoot->SetLineColor(4);
	pt2NewRoot->SetFillStyle(0);

	pt2NewRootUncorr->SetMarkerStyle(0);
	pt2NewRootUncorr->SetMarkerColor(4);
	pt2NewRootUncorr->SetLineColor(6);
	pt2NewRootUncorr->SetFillStyle(0);

	//Creation
	TGraph *pt2Selected = new TGraph();
	pt2Selected->SetName("fffit");
	//Points
	for(int i=0; i<npoints; ++i){
		pt2Selected->SetPoint(i, cutValue[i], pointsNEvents[i]);
	}
	//Style
	pt2Selected->SetMarkerStyle(20);
	pt2Selected->SetMarkerColor(4);
	pt2Selected->SetLineColor(4);
	pt2Selected->SetFillStyle(0);

	//Plotting
	TCanvas *c = new TCanvas("plot");
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	//pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();

	pt2NewRoot->SetTitle("FF Slope fit result");
	pt2NewRoot->GetXaxis()->SetTitle("Cut value");
	pt2NewRoot->GetYaxis()->SetTitle("FF Slope a");
	pt2NewRoot->GetYaxis()->SetTitleOffset(1.5);

	pt2NewRootUncorr->SetTitle("");
	pt2NewRoot->SetTitle("");
	pt2NewRoot->SetLineWidth(2);
	pt2NewRootUncorr->SetLineWidth(2);
	pt2NewRoot->GetYaxis()->SetNdivisions(505);
	pt2NewRoot->GetYaxis()->SetTitleSize(15);
	pt2NewRoot->GetYaxis()->SetTitleFont(43);
	pt2NewRoot->GetYaxis()->SetTitleOffset(1.4);
	pt2NewRoot->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	pt2NewRoot->GetYaxis()->SetLabelSize(15);
	pt2NewRoot->GetXaxis()->SetLabelFont(43);
	pt2NewRoot->GetXaxis()->SetLabelSize(15);
	pt2NewRoot->GetYaxis()->SetRangeUser(0.025, 0.05);
	pt2NewRootUncorr->GetYaxis()->SetRangeUser(0.025, 0.05);

	pt2NewRoot->Draw("AP");
	pt2NewRootUncorr->Draw("PSAME");

//	TCanvas *c = new TCanvas("plotUncorr");
//	pt2NewRootUncorr->SetTitle("FF Slope fit result (Uncorrelated errors)");
//	pt2NewRootUncorr->GetXaxis()->SetTitle("Cut value");
//	pt2NewRootUncorr->GetYaxis()->SetTitle("FF Slope a");
//	pt2NewRootUncorr->GetYaxis()->SetTitleOffset(1.5);

	//Plotting

	TCanvas *c = new TCanvas("Selected");
	pt2Selected->SetTitle("Number of selected events");
	pt2Selected->GetXaxis()->SetTitle("Cut value");
	pt2Selected->GetYaxis()->SetTitle("Selected events");
	pt2Selected->GetYaxis()->SetTitleOffset(1.5);
	pt2Selected->Draw("APL");
}

void plot2(){
	//plotpt2();
	//plotpi0mass();
	//plotBins();
	plot_only1();
}

