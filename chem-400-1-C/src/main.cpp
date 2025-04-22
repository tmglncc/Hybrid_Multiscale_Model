#include <iostream>
#include <cmath>

#include "Util.hpp"
#include "Cell.hpp"
#include "CellManipulation.hpp"
#include "FileFactory.hpp"
#include "Frame.hpp"
#include "FrameFactory.hpp"
#include "Mesh.hpp"
#include "ConfigHandler.hpp"

#include "Pathway.hpp"
#include "Ran.hpp"
#include "NR3.hpp"
#include "Macro.hpp"
#include "BloodVesselFactory.hpp"
#include "Chemotherapy.hpp"

int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	ConfigHandler config(argv[1], atoi(argv[2]));

	if (!config.created())
		return EXIT_FAILURE;

	Ran ran(config.input.seed);
	
	EGFRPathway egfrPathway;
	CellCyclePathway cellCyclePathway;

	Chemotherapy chemotherapy(config.continuum.chemotherapy.chemTimes, config.continuum.chemotherapy.chemDosages, config.continuum.chemotherapy.chemProtocol);
	
	/*PharmacologicModel pharmacologic(&config);

	pharmacologic.buildPKCurve();
	pharmacologic.buildPDCurve(1.0e-2);

	if (config.output.files) {
		FileFactory::makePKFile(&pharmacologic, &config);
		FileFactory::makePDFile(&pharmacologic, &config);
	}*/

	std::vector<Vector3> bloodVessel;
	Frame* frame;
	Mesh* mesh;
	Macro* macro;

	if (config.input.fileFormat == BI_D) {
		bloodVessel = BloodVesselFactory::makeBloodVessel2D(config.input.bloodVessel);
		frame = FrameFactory::makeFrame2D(config.input.initialCondition, &ran, &cellCyclePathway, &chemotherapy); // initial condition
		mesh = new Mesh2D(frame->domain, &config, &chemotherapy);
		macro = new Macro2D(mesh, frame, &config, bloodVessel, &chemotherapy);
	} else {
		frame = FrameFactory::makeFrame3D(config.input.initialCondition, &ran, &cellCyclePathway, &chemotherapy); // initial condition
		mesh = new Mesh3D(frame->domain, &config, &chemotherapy);
		macro = new Macro3D(mesh, frame, &config, bloodVessel, &chemotherapy);
	}

	if (config.output.files)
		FileFactory::makeFrameFile(frame, &config);

	macro->reaction();
	macro->differenceO();
	macro->differenceEGF();
	macro->differenceChem();
	
	double* F_hap = CellManipulation::generateF_hap(frame, mesh, &ran);

	while (frame->time < config.input.timeMax) {
		frame->time += config.parameters.delta_tt;

		if (config.input.fileFormat == BI_D)
			CellManipulation2D::updateFrame(frame, &config, mesh, &egfrPathway, &cellCyclePathway, &ran, &chemotherapy);
		else
			CellManipulation3D::updateFrame(frame, &config, mesh, &egfrPathway, &cellCyclePathway, &ran, &chemotherapy);

		if (config.output.prints)
			std::cout << frame->to_string();
		else if (config.output.justLastPrint && frame->time == config.input.timeMax)
			std::cout << frame->to_string();

		CellManipulation::calculateGrad_F_hap(frame, mesh, F_hap);

		if (config.input.fileFormat == BI_D)
			CellManipulation2D::force(frame, &config, bloodVessel);
		else
			CellManipulation3D::force(frame, &config, bloodVessel);

		if (config.output.files && frame->time%config.output.filesFrequency == 0)
			FileFactory::makeFrameFile(frame, &config);
		else if (config.output.justLastFile && frame->time == config.input.timeMax)
			FileFactory::makeFrameFile(frame, &config);

		macro->reaction();
		macro->differenceO();
		macro->differenceEGF();
		macro->differenceChem();
	}

	delete F_hap;

	delete frame;
	delete mesh;
	delete macro;

	return EXIT_SUCCESS;
}

/*int main(int argc, char *argv[]) {
	ConfigHandler config(argv[1], atoi(argv[2]));

	if (!config.created())
		return EXIT_FAILURE;

	Ran ran(config.input.seed);

	std::vector<Vector3> bloodVessel;
	Frame* frame;

	if (config.input.fileFormat == BI_D) {
		bloodVessel = BloodVesselFactory::makeBloodVessel2D(config.input.bloodVessel);
		frame = FrameFactory::makeFrame2D(config.input.initialCondition, &ran); // initial condition
	} else {
		frame = FrameFactory::makeFrame3D(config.input.initialCondition, &ran); // initial condition
	}

	CellManipulation::updateInitialFrame(frame, &config, &ran);

	if (config.output.prints)
		std::cout << frame->to_string();

	if (config.input.fileFormat == BI_D)
		CellManipulation2D::force(frame, &config, bloodVessel);
	else
		CellManipulation3D::force(frame, &config, bloodVessel);

	if (config.output.files && frame->time%config.output.filesFrequency == 0)
		FileFactory::makeFrameFile(frame, &config);

	delete frame;

	return EXIT_SUCCESS;
}*/
