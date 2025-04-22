#include "BloodVesselFactory.hpp"

std::vector<Vector3> BloodVesselFactory::makeBloodVessel2D(std::string fileName) {
	// FILE
	std::ifstream file;

	// BLOOD VESSEL
	std::vector<Vector3> bloodVessel = std::vector<Vector3>();

	Vector3 coordinates;
	int numVessels = 0;

	file.open(fileName.c_str());

	file >> numVessels;

	// Reading all vessels
	for (std::size_t i = 0; i < numVessels; i++) {
		file >> coordinates.x >> coordinates.y; coordinates.z = 0.0;

		bloodVessel.push_back(coordinates);
	}

	file.close();

	return bloodVessel;
}