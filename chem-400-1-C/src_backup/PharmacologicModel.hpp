#ifndef PHARMACOLOGIC_MODEL
#define PHARMACOLOGIC_MODEL

#include <cmath>

#include "../src/ConfigHandler.hpp"

class PharmacologicModel {
private:
	ConfigHandler* config;

public:
	int pkSize;
	double* time;
	double* pkCurve;

	int pdSize;
	double* concentration;
	double* pdCurve;

	PharmacologicModel(ConfigHandler* config);
	double computeConcentration(double t);
	double computeDrugEffectFactor(double u);
	void buildPKCurve(double delta_t = 1);
	void buildPDCurve(double delta_u = 1);
};

#endif