#ifndef PLOT_FACTORY
#define PLOT_FACTORY

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "Util.hpp"
#include "Vector.hpp"

class PlotFactory {
public:
	static void makePlot(char* inputFilename, char* outputFilename, const char* outputMode);
	static void makeFramePlot(char* inputFilename, char* outputFilename, 
		const char* bvFilename = "../../input/inside-blood-vessel1-2D.dat", double bvRadius = 9.953);
	static void makePhenotypeData(char* outputFilename, int timeMax = 800, int inputFilenamesNumber = 2);
	static void makePhenotypePlot(char* inputFilename, char* outputFilename);
	static void makeChemMeanData(char* outputFilename, int timeMax = 700, int inputFilenamesNumber = 2,
		const char* chemProtocol = "B");
	static void makeStochasticData(char* outputFilename, int seeds = 10, int timeMax = 700);
	static void makeStochasticPlot(char* inputFilename, char* outputFilename);
	static void makeStochasticChemData(char* outputFilename, int seeds = 10, int timeMax = 700);
	static void makeStochasticChemPlot(char* inputFilename, char* outputFilename);
};

#endif
