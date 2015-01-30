#include <TGraphErrors.h>
#include <TColor.h>
#include <TROOT.h>
#include <TStyle.h>

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

	//Creation
	TGraphErrors *pt2Old = new TGraphErrors();
	pt2Old->SetName("pt2Old");

	//Points
	pt2Old->SetPoint(0, 0.0001, 0.024681);
	pt2Old->SetPoint(1, 0.0002, 0.033619);
	pt2Old->SetPoint(2, 0.0003, 0.035015);
	pt2Old->SetPoint(3, 0.0004, 0.036419);
	pt2Old->SetPoint(4, 0.0005, 0.035884);
	pt2Old->SetPoint(5, 0.0006, 0.035861);
	pt2Old->SetPoint(6, 0.0007, 0.035496);
	pt2Old->SetPoint(7, 0.0008, 0.035786);
	pt2Old->SetPoint(8, 0.0009, 0.035620);
	pt2Old->SetPoint(9, 0.0010, 0.035405);

	//Errors
	pt2Old->SetPointError(0, 0, 0.0109678);
	pt2Old->SetPointError(1, 0, 0.0107274);
	pt2Old->SetPointError(2, 0, 0.0106531);
	pt2Old->SetPointError(3, 0, 0.0106226);
	pt2Old->SetPointError(4, 0, 0.0106064);
	pt2Old->SetPointError(5, 0, 0.0105968);
	pt2Old->SetPointError(6, 0, 0.0105909);
	pt2Old->SetPointError(7, 0, 0.0105868);
	pt2Old->SetPointError(8, 0, 0.0105842);
	pt2Old->SetPointError(9, 0, 0.0105821);

	//Style
	pt2Old->SetMarkerStyle(20);
	pt2Old->SetMarkerColor(2);
	pt2Old->SetLineColor(2);
	pt2Old->SetFillStyle(0);


	//Creation
	TGraphErrors *pt2OldRoot = new TGraphErrors();
	pt2Old->SetName("pt2OldRoot");


	//Points
	pt2OldRoot->SetPoint(0, 0.00012, 0.0246773);
	pt2OldRoot->SetPoint(1, 0.00022, 0.0335654);
	pt2OldRoot->SetPoint(2, 0.00032, 0.0349375);
	pt2OldRoot->SetPoint(3, 0.00042, 0.0363305);
	pt2OldRoot->SetPoint(4, 0.00052, 0.0357917);
	pt2OldRoot->SetPoint(5, 0.00062, 0.0357685);
	pt2OldRoot->SetPoint(6, 0.00072, 0.0354012);
	pt2OldRoot->SetPoint(7, 0.00082, 0.0356871);
	pt2OldRoot->SetPoint(8, 0.00092, 0.0355217);
	pt2OldRoot->SetPoint(9, 0.00102, 0.0353101);

	//Errors
	pt2OldRoot->SetPointError(0, 0, 0.0137759);
	pt2OldRoot->SetPointError(1, 0, 0.0134755);
	pt2OldRoot->SetPointError(2, 0, 0.0133764);
	pt2OldRoot->SetPointError(3, 0, 0.0133344);
	pt2OldRoot->SetPointError(4, 0, 0.0133112);
	pt2OldRoot->SetPointError(5, 0, 0.0132987);
	pt2OldRoot->SetPointError(6, 0, 0.0132893);
	pt2OldRoot->SetPointError(7, 0, 0.0132858);
	pt2OldRoot->SetPointError(8, 0, 0.0132825);
	pt2OldRoot->SetPointError(9, 0, 0.0132792);

	//Style
	pt2OldRoot->SetMarkerStyle(21);
	pt2OldRoot->SetMarkerColor(3);
	pt2OldRoot->SetLineColor(3);
	pt2OldRoot->SetFillStyle(0);


	//Creation
	TGraphErrors *pt2New = new TGraphErrors();
	pt2Old->SetName("pt2New");


	//Points
	pt2New->SetPoint(0, 0.00014, 0.0263395);
	pt2New->SetPoint(1, 0.00024, 0.0349747);
	pt2New->SetPoint(2, 0.00034, 0.0364357);
	pt2New->SetPoint(3, 0.00044, 0.0377355);
	pt2New->SetPoint(4, 0.00054, 0.0371447);
	pt2New->SetPoint(5, 0.00064, 0.0370775);
	pt2New->SetPoint(6, 0.00074, 0.036707);
	pt2New->SetPoint(7, 0.00084, 0.0370406);
	pt2New->SetPoint(8, 0.00094, 0.0368565);
	pt2New->SetPoint(9, 0.00104, 0.0367007);

	//Errors
	pt2New->SetPointError(0, 0, 0.0138228);
	pt2New->SetPointError(1, 0, 0.0134967);
	pt2New->SetPointError(2, 0, 0.0133997);
	pt2New->SetPointError(3, 0, 0.0133588);
	pt2New->SetPointError(4, 0, 0.0133362);
	pt2New->SetPointError(5, 0, 0.0133233);
	pt2New->SetPointError(6, 0, 0.0133156);
	pt2New->SetPointError(7, 0, 0.013318);
	pt2New->SetPointError(8, 0, 0.0133087);
	pt2New->SetPointError(9, 0, 0.0133054);

	//Style
	pt2New->SetMarkerStyle(22);
	pt2New->SetMarkerColor(4);
	pt2New->SetLineColor(4);
	pt2New->SetFillStyle(0);

	//Creation
	TGraphErrors *pt2NewRoot = new TGraphErrors();
	pt2Old->SetName("pt2NewRoot");


	//Points
	pt2NewRoot->SetPoint(0, 0.00016, 0.0262831);
	pt2NewRoot->SetPoint(1, 0.00026, 0.0349844);
	pt2NewRoot->SetPoint(2, 0.00036, 0.0364254);
	pt2NewRoot->SetPoint(3, 0.00046, 0.0377165);
	pt2NewRoot->SetPoint(4, 0.00056, 0.0371151);
	pt2NewRoot->SetPoint(5, 0.00066, 0.037052);
	pt2NewRoot->SetPoint(6, 0.00076, 0.0366747);
	pt2NewRoot->SetPoint(7, 0.00086, 0.0370075);
	pt2NewRoot->SetPoint(8, 0.00096, 0.0368243);
	pt2NewRoot->SetPoint(9, 0.00106, 0.0366688);

	//Errors
	pt2NewRoot->SetPointError(0, 0, 0.0137436);
	pt2NewRoot->SetPointError(1, 0, 0.0134296);
	pt2NewRoot->SetPointError(2, 0, 0.0133354);
	pt2NewRoot->SetPointError(3, 0, 0.0132939);
	pt2NewRoot->SetPointError(4, 0, 0.0132721);
	pt2NewRoot->SetPointError(5, 0, 0.0132602);
	pt2NewRoot->SetPointError(6, 0, 0.0132522);
	pt2NewRoot->SetPointError(7, 0, 0.0132479);
	pt2NewRoot->SetPointError(8, 0, 0.013246);
	pt2NewRoot->SetPointError(9, 0, 0.0132423);

	//Style
	pt2NewRoot->SetMarkerStyle(23);
	pt2NewRoot->SetMarkerColor(6);
	pt2NewRoot->SetLineColor(6);
	pt2NewRoot->SetFillStyle(0);

	//Plotting
	TCanvas *c = new TCanvas("plot");
	c->SetGrid(1, 1);

	TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.95);
	leg->AddEntry(pt2Old, "Pt^{2} Old Method");
	leg->AddEntry(pt2OldRoot, "Pt^{2} Old Method with ROOT #chi^{2}");
	leg->AddEntry(pt2New, "Pt^{2} New Method");
	leg->AddEntry(pt2NewRoot, "Pt^{2} New Method with ROOT #chi^{2}");

	pt2Old->GetYaxis()->SetRangeUser(0.01, 0.06);
	pt2Old->SetTitle("Pt^{2} scan");
	pt2Old->GetXaxis()->SetTitle("Pt^{2} cut");
	pt2Old->GetYaxis()->SetTitle("FF Slope a");
	pt2Old->GetYaxis()->SetTitleOffset(1.5);
	pt2Old->Draw("AP");
	pt2OldRoot->Draw("SAME P");
	pt2New->Draw("SAME P");
	pt2NewRoot->Draw("SAME P");
	leg->Draw("LEP");
}

void plotpi0mass(){
	//setStyle();

	//Creation
	TGraphErrors *pi0massOld = new TGraphErrors();
	pi0massOld->SetName("pi0massOld");

	//Points
	pi0massOld->SetPoint(0, 0.0001, 0.024681);
	pi0massOld->SetPoint(1, 0.0002, 0.033619);
	pi0massOld->SetPoint(2, 0.0003, 0.035015);
	pi0massOld->SetPoint(3, 0.0004, 0.036419);
	pi0massOld->SetPoint(4, 0.0005, 0.035884);
	pi0massOld->SetPoint(5, 0.0006, 0.035861);
	pi0massOld->SetPoint(6, 0.0007, 0.035496);
	pi0massOld->SetPoint(7, 0.0008, 0.035786);
	pi0massOld->SetPoint(8, 0.0009, 0.035620);
	pi0massOld->SetPoint(9, 0.0010, 0.035405);

	//Errors
	pi0massOld->SetPointError(0, 0, 0.0109678);
	pi0massOld->SetPointError(1, 0, 0.0107274);
	pi0massOld->SetPointError(2, 0, 0.0106531);
	pi0massOld->SetPointError(3, 0, 0.0106226);
	pi0massOld->SetPointError(4, 0, 0.0106064);
	pi0massOld->SetPointError(5, 0, 0.0105968);
	pi0massOld->SetPointError(6, 0, 0.0105909);
	pi0massOld->SetPointError(7, 0, 0.0105868);
	pi0massOld->SetPointError(8, 0, 0.0105842);
	pi0massOld->SetPointError(9, 0, 0.0105821);

	//Style
	pi0massOld->SetMarkerStyle(20);
	pi0massOld->SetMarkerColor(2);
	pi0massOld->SetLineColor(2);
	pi0massOld->SetFillStyle(0);


	//Creation
	TGraphErrors *pi0massOldRoot = new TGraphErrors();
	pi0massOld->SetName("pi0massOldRoot");


	//Points
	pi0massOldRoot->SetPoint(0, 0.00012, 0.0246773);
	pi0massOldRoot->SetPoint(1, 0.00022, 0.0335654);
	pi0massOldRoot->SetPoint(2, 0.00032, 0.0349375);
	pi0massOldRoot->SetPoint(3, 0.00042, 0.0363305);
	pi0massOldRoot->SetPoint(4, 0.00052, 0.0357917);
	pi0massOldRoot->SetPoint(5, 0.00062, 0.0357685);
	pi0massOldRoot->SetPoint(6, 0.00072, 0.0354012);
	pi0massOldRoot->SetPoint(7, 0.00082, 0.0356871);
	pi0massOldRoot->SetPoint(8, 0.00092, 0.0355217);
	pi0massOldRoot->SetPoint(9, 0.00102, 0.0353101);

	//Errors
	pi0massOldRoot->SetPointError(0, 0, 0.0137759);
	pi0massOldRoot->SetPointError(1, 0, 0.0134755);
	pi0massOldRoot->SetPointError(2, 0, 0.0133764);
	pi0massOldRoot->SetPointError(3, 0, 0.0133344);
	pi0massOldRoot->SetPointError(4, 0, 0.0133112);
	pi0massOldRoot->SetPointError(5, 0, 0.0132987);
	pi0massOldRoot->SetPointError(6, 0, 0.0132893);
	pi0massOldRoot->SetPointError(7, 0, 0.0132858);
	pi0massOldRoot->SetPointError(8, 0, 0.0132825);
	pi0massOldRoot->SetPointError(9, 0, 0.0132792);

	//Style
	pi0massOldRoot->SetMarkerStyle(21);
	pi0massOldRoot->SetMarkerColor(3);
	pi0massOldRoot->SetLineColor(3);
	pi0massOldRoot->SetFillStyle(0);


	//Creation
	TGraphErrors *pi0massNew = new TGraphErrors();
	pi0massOld->SetName("pi0massNew");


	//Points
	pi0massNew->SetPoint(0, 0.00014, 0.0263395);
	pi0massNew->SetPoint(1, 0.00024, 0.0349747);
	pi0massNew->SetPoint(2, 0.00034, 0.0364357);
	pi0massNew->SetPoint(3, 0.00044, 0.0377355);
	pi0massNew->SetPoint(4, 0.00054, 0.0371447);
	pi0massNew->SetPoint(5, 0.00064, 0.0370775);
	pi0massNew->SetPoint(6, 0.00074, 0.036707);
	pi0massNew->SetPoint(7, 0.00084, 0.0370406);
	pi0massNew->SetPoint(8, 0.00094, 0.0368565);
	pi0massNew->SetPoint(9, 0.00104, 0.0367007);

	//Errors
	pi0massNew->SetPointError(0, 0, 0.0138228);
	pi0massNew->SetPointError(1, 0, 0.0134967);
	pi0massNew->SetPointError(2, 0, 0.0133997);
	pi0massNew->SetPointError(3, 0, 0.0133588);
	pi0massNew->SetPointError(4, 0, 0.0133362);
	pi0massNew->SetPointError(5, 0, 0.0133233);
	pi0massNew->SetPointError(6, 0, 0.0133156);
	pi0massNew->SetPointError(7, 0, 0.013318);
	pi0massNew->SetPointError(8, 0, 0.0133087);
	pi0massNew->SetPointError(9, 0, 0.0133054);

	//Style
	pi0massNew->SetMarkerStyle(22);
	pi0massNew->SetMarkerColor(4);
	pi0massNew->SetLineColor(4);
	pi0massNew->SetFillStyle(0);

	//Creation
	TGraphErrors *pi0massNewRoot = new TGraphErrors();
	pi0massOld->SetName("pi0massNewRoot");


	//Points
	pi0massNewRoot->SetPoint(0, 0.00016, 0.0262831);
	pi0massNewRoot->SetPoint(1, 0.00026, 0.0349844);
	pi0massNewRoot->SetPoint(2, 0.00036, 0.0364254);
	pi0massNewRoot->SetPoint(3, 0.00046, 0.0377165);
	pi0massNewRoot->SetPoint(4, 0.00056, 0.0371151);
	pi0massNewRoot->SetPoint(5, 0.00066, 0.037052);
	pi0massNewRoot->SetPoint(6, 0.00076, 0.0366747);
	pi0massNewRoot->SetPoint(7, 0.00086, 0.0370075);
	pi0massNewRoot->SetPoint(8, 0.00096, 0.0368243);
	pi0massNewRoot->SetPoint(9, 0.00106, 0.0366688);

	//Errors
	pi0massNewRoot->SetPointError(0, 0, 0.0137436);
	pi0massNewRoot->SetPointError(1, 0, 0.0134296);
	pi0massNewRoot->SetPointError(2, 0, 0.0133354);
	pi0massNewRoot->SetPointError(3, 0, 0.0132939);
	pi0massNewRoot->SetPointError(4, 0, 0.0132721);
	pi0massNewRoot->SetPointError(5, 0, 0.0132602);
	pi0massNewRoot->SetPointError(6, 0, 0.0132522);
	pi0massNewRoot->SetPointError(7, 0, 0.0132479);
	pi0massNewRoot->SetPointError(8, 0, 0.013246);
	pi0massNewRoot->SetPointError(9, 0, 0.0132423);

	//Style
	pi0massNewRoot->SetMarkerStyle(23);
	pi0massNewRoot->SetMarkerColor(6);
	pi0massNewRoot->SetLineColor(6);
	pi0massNewRoot->SetFillStyle(0);

	//Plotting
	TCanvas *c = new TCanvas("plot");
	c->SetGrid(1, 1);

	TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.95);
	leg->AddEntry(pi0massOld, "m_{ee\gamma} Old Method");
	leg->AddEntry(pi0massOldRoot, "m_{ee\gamma} Old Method with ROOT #chi^{2}");
	leg->AddEntry(pi0massNew, "m_{ee\gamma} New Method");
	leg->AddEntry(pi0massNewRoot, "m_{ee\gamma} New Method with ROOT #chi^{2}");

	pi0massOld->GetYaxis()->SetRangeUser(0.01, 0.06);
	pi0massOld->SetTitle("m_{ee\gamma} scan");
	pi0massOld->GetXaxis()->SetTitle("m_{ee\gamma} cut");
	pi0massOld->GetYaxis()->SetTitle("FF Slope a");
	pi0massOld->GetYaxis()->SetTitleOffset(1.5);
	pi0massOld->Draw("AP");
	pi0massOldRoot->Draw("SAME P");
	pi0massNew->Draw("SAME P");
	pi0massNewRoot->Draw("SAME P");
	leg->Draw("LEP");
}

void plot(){
	plotpt2();
}

