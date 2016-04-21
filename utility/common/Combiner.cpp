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

void Combiner::draw(vector<int> allColors, vector<int> dataColors, int first, int last) {
	for (unsigned int i = 0; i < fMCSamples.size(); i++) {
		vector<int> colors(allColors.begin() + (i*3), allColors.begin() + (i*3+3));
		std::cout << "toto" << endl;
		for(auto it: colors) std::cout << it << " ";
		std::cout << std::endl;
		fMCSamples[i]->setPlotStyle(colors);
	}

	for (unsigned int i = 0; i < fDataSamples.size(); i++) {
		vector<int> colors(dataColors.begin() + (i*3), dataColors.begin() + (i*3+3));
		fDataSamples[i]->setPlotStyle(colors);
	}

	Drawer::drawCombineStack(fMCSamples, fDataSamples, fFinalMCSample, fFinalDataSample, first, last);
}
