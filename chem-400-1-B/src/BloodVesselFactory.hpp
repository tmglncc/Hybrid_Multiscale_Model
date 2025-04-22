#ifndef BLOOD_VESSEL_FACTORY
#define BLOOD_VESSEL_FACTORY

#include <fstream>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>

#include "Vector.hpp"

class BloodVesselFactory {
public:
	static std::vector<Vector3> makeBloodVessel2D(std::string fileName);
};

#endif /* end of include guard: BLOOD_VESSEL_FACTORY */