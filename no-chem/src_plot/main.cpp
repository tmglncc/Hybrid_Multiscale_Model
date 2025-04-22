#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

#include "PlotFactory.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 4)
		return EXIT_FAILURE;

	if (strcmp(argv[1], "EGF") == 0)
		PlotFactory::makePlot(argv[2], argv[3], "EGF");
	else if (strcmp(argv[1], "NUT") == 0)
		PlotFactory::makePlot(argv[2], argv[3], "Oxygen");
	else if (strcmp(argv[1], "CHEM") == 0)
		PlotFactory::makePlot(argv[2], argv[3], "Chemotherapy");
	else if (strcmp(argv[1], "OUT") == 0)
		PlotFactory::makeFramePlot(argv[2], argv[3]);
	else if (strcmp(argv[1], "PHEN") == 0) {
		PlotFactory::makePhenotypeData(argv[2]);
		PlotFactory::makePhenotypePlot(argv[2], argv[3]);
	} else if (strcmp(argv[1], "STOC") == 0) {
		PlotFactory::makeStochasticData(argv[2]);
		PlotFactory::makeStochasticPlot(argv[2], argv[3]);
	}

	return EXIT_SUCCESS;
}