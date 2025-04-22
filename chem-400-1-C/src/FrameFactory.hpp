#ifndef FRAME_FACTORY
#define FRAME_FACTORY

#include <fstream>
#include <cmath>
#include <string>
#include <iostream>

#include "Frame.hpp"
#include "Pathway.hpp"
#include "NR3.hpp"
#include "Ran.hpp"
#include "Chemotherapy.hpp"

class FrameFactory {
public:
	static Frame* makeFrame3D(std::string fileName, Ran* ran, CellCyclePathway* cellCyclePathway = NULL, Chemotherapy* chemotherapy = NULL);
	static Frame* makeFrame2D(std::string fileName, Ran* ran, CellCyclePathway* cellCyclePathway = NULL, Chemotherapy* chemotherapy = NULL);
};

#endif /* end of include guard: FRAME_FACTORY */
