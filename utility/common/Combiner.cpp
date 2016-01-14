/*
 * Combiner.cpp
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#include "Combiner.h"
#include "Drawer/Drawer.h"

using namespace std;

Combiner::Combiner() {
	// TODO Auto-generated constructor stub

}

Combiner::~Combiner() {
	// TODO Auto-generated destructor stub
}

void Combiner::draw(vector<int> allColors, vector<int> dataColors) {
	for (unsigned int i = 0; i < fMCSamples.size(); i++) {
		vector<int> colors(allColors.begin() + (i*3), allColors.begin() + (i*3+3));
		fMCSamples[i]->setPlotStyle(colors);
	}

	for (unsigned int i = 0; i < fDataSamples.size(); i++) {
		vector<int> colors(dataColors.begin() + (i*3), dataColors.begin() + (i*3+3));
		fDataSamples[i]->setPlotStyle(colors);
	}

	Drawer::drawCombineStack(fMCSamples, fDataSamples, fFinalMCSample, fFinalDataSample);
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
//	//	doPlot(++i, "Zvtx_low", "Vertex Z", leg, mcColors);
//	//	doPlot(++i, "Zvtx_high", "Vertex Z", leg, mcColors);
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
	//doPlot2(++iMap, "DCH4_XY_gamma", "Photon DCH4 map", leg, mcColors);
}
