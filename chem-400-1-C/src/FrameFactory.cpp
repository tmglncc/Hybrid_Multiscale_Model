#include "FrameFactory.hpp"

Frame* FrameFactory::makeFrame3D(std::string fileName, Ran* ran, CellCyclePathway* cellCyclePathway, Chemotherapy* chemotherapy) {
	// FILE
	std::ifstream file;

	// FRAME
	Frame* frame = new Frame();

	Cell c;
	// int state;
	int numCells = 0;

	file.open(fileName.c_str());

	file >> frame->domain.x >> frame->domain.y >> frame->domain.z;
	file >> numCells >> frame->time;
	file >> frame->outCells >> frame->tumorCells;

	// Reading all cells
	for (std::size_t i = 0; i < numCells; i++) {
		file >> c.state;
		file >> c.coordinates.x >> c.coordinates.y >> c.coordinates.z;
		file >> c.nucleusRadius >> c.radius >> c.actionRadius >> c.lifetime >> c.previousState;
		file >> c.oConsumption >> c.egfSource >> c.calcification >> c.sigmaEGF >> c.sigmaO;
		file >> c.speed.x >> c.speed.y >> c.speed.z;

		if (chemotherapy != NULL)
			c.sigmaChem = std::vector<double>(chemotherapy->chemDrugs.size(), 0.0);
		if (cellCyclePathway != NULL)
			cellCyclePathway->setInitialAttributes(ran, &c.cellCycleConc[0], &c.HIF, &c.mu);

		frame->cells.push_back(c);
	}

	file.close();

	return frame;
}

Frame* FrameFactory::makeFrame2D(std::string fileName, Ran* ran, CellCyclePathway* cellCyclePathway, Chemotherapy* chemotherapy) {
	// FILE
	std::ifstream file;

	// FRAME
	Frame* frame = new Frame();

	Cell c;
	// int state;
	int numCells = 0;

	file.open(fileName.c_str());

	file >> frame->domain.x >> frame->domain.y; frame->domain.z = 0.0;
	file >> numCells >> frame->time;
	file >> frame->outCells >> frame->tumorCells;

	// Reading all cells
	for (std::size_t i = 0; i < numCells; i++) {
		file >> c.state;
		file >> c.coordinates.x >> c.coordinates.y; c.coordinates.z = 0.0;
		file >> c.nucleusRadius >> c.radius >> c.actionRadius >> c.lifetime >> c.previousState;
		file >> c.oConsumption >> c.egfSource >> c.calcification >> c.sigmaEGF >> c.sigmaO;
		file >> c.speed.x >> c.speed.y; c.speed.z = 0.0;

		if (chemotherapy != NULL)
			c.sigmaChem = std::vector<double>(chemotherapy->chemDrugs.size(), 0.0);
		if (cellCyclePathway != NULL)
			cellCyclePathway->setInitialAttributes(ran, &c.cellCycleConc[0], &c.HIF, &c.mu);

		frame->cells.push_back(c);
	}

	file.close();

	return frame;
}
