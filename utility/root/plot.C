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
	double pointsOld[10] = {0.0397792, 0.041743, 0.0428235, 0.0438376, 0.0435394, 0.043297, 0.0433736, 0.0435701, 0.0434316, 0.0432505};
	double errOld[10] = {0.010694, 0.0104719, 0.0104066, 0.0103798, 0.010366, 0.0103568, 0.0103507, 0.0103469, 0.0103442, 0.0103415};
	
	double pointsOldRoot[10] = {0.0394773, 0.0411424, 0.0422285, 0.0432323, 0.042954, 0.0427035, 0.0427771, 0.0429641, 0.0428285, 0.0426402};
	double errOldRoot[10] = {0.013312, 0.013045, 0.0129558, 0.0129198, 0.0128983, 0.0128855, 0.0128762, 0.0128733, 0.0128695, 0.0128659}; 

	double pointsNew[10] = {0.0408979, 0.0426648, 0.043845, 0.0448861, 0.0445144, 0.0443114, 0.0443202, 0.0445949, 0.0444315, 0.0442956};
	double errNew[10] = {0.0133229, 0.0130023, 0.0129171, 0.0128793, 0.0128599, 0.0128485, 0.0128416, 0.0128373, 0.0128355, 0.0128317};
	
	double pointsNewRoot[10] = {0.0407679, 0.0424753, 0.0436177, 0.0446595, 0.0442785, 0.0440594, 0.0440673, 0.0443324, 0.0441651, 0.0440268};
	double errNewRoot[10] = {0.013248, 0.0129565, 0.0128696, 0.0128313, 0.0128102, 0.0127985, 0.0127918, 0.0127874, 0.0127854, 0.0127822};
	
	double pointsNEvents[10] = {952820, 1005159, 1020619, 1027590, 1031026, 1033011, 1034237, 1035008, 1035359, 1035927};

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

	double pointsOld[10] = {0.0149355, 0.0173647, 0.0238414, 0.032847, 0.0362169, 0.0383385, 0.0403324, 0.0411539, 0.0412099, 0.0446113};
	double errOld[10] = {0.0118117, 0.0110085, 0.0107409, 0.0106093, 0.0105521, 0.0105102, 0.0104584, 0.0104281, 0.0104336, 0.0103517};
		
	double pointsOldRoot[10] = {0.0152272, 0.0174402, 0.0237519, 0.0325828, 0.035904, 0.0380045, 0.0399277, 0.0406624, 0.0406641, 0.0439828};
	double errOldRoot[10] = {0.0145287, 0.0135896, 0.0132953, 0.0131689, 0.0131064, 0.0130632, 0.013027, 0.0129856, 0.0130037, 0.0128912};

	double pointsNew[10] = {0.0146534, 0.0186988, 0.0251068, 0.0338209, 0.0372407, 0.0400598, 0.0415006, 0.0422051, 0.04217, 0.0452986};
	double errNew[10] = {0.014576, 0.0136149, 0.0133031, 0.013167, 0.0130938, 0.0130392, 0.0130014, 0.0129531, 0.0129684, 0.01285};
	
	double pointsNewRoot[10] = {0.0145182, 0.0184307, 0.0248997, 0.033599, 0.0370024, 0.0398518, 0.0412759, 0.0419821, 0.0418912, 0.0450634};
	double errNewRoot[10] = {0.0144619, 0.0135061, 0.013209, 0.0130864, 0.0130184, 0.0129677, 0.0129375, 0.0128975, 0.012916, 0.0128079};
	
	double pointsNEvents[10] = {792280, 916082, 963504, 985346, 997136, 1005345, 1009929, 1016579, 1013040, 1030818};

	//Creation
	TGraphErrors *pi0massOld = new TGraphErrors();
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
	pi0massOld->SetFillStyle(0);


	//Creation
	TGraphErrors *pi0massOldRoot = new TGraphErrors();
	pi0massOld->SetName("pi0massOldRoot");
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
	pi0massOld->SetName("pi0massNew");
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
	pi0massOld->SetName("pi0massNewRoot");
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
	leg->AddEntry(pi0massOld, "m_{ee#gamma} Old Method");
	leg->AddEntry(pi0massOldRoot, "m_{ee#gamma} Old Method with ROOT #chi^{2}");
	leg->AddEntry(pi0massNew, "m_{ee#gamma} New Method");
	leg->AddEntry(pi0massNewRoot, "m_{ee#gamma} New Method with ROOT #chi^{2}");

	pi0massOld->GetYaxis()->SetRangeUser(0, 0.08);
	pi0massOld->SetTitle("m_{ee#gamma} scan");
	pi0massOld->GetXaxis()->SetTitle("m_{ee#gamma} cut");
	pi0massOld->GetYaxis()->SetTitle("FF Slope a");
	pi0massOld->GetYaxis()->SetTitleOffset(1.5);
	pi0massOld->Draw("AP");
	pi0massOldRoot->Draw("SAME P");
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

void plot(){
	plotpt2();
	plotpi0mass();
	plotBins();
}

