#include "PharmacologicModel.hpp"

PharmacologicModel::PharmacologicModel(ConfigHandler* config) {
	this->config = config;
}

double PharmacologicModel::computeConcentration(double t) {
	double u = 0.0;
	double smallValue = 1.0e-8;
	double bigValue = 4.819e42;

	if (t <= this->config->pharmacologic.pk.p1)
		u = smallValue*exp(this->config->pharmacologic.pk.lambda_a*(fmod(t, this->config->pharmacologic.pk.tau)));
	else if (t <= this->config->pharmacologic.pk.p2)
		u = this->config->pharmacologic.pk.zeta;
	else if (t <= this->config->pharmacologic.pk.p3)
		u = bigValue*exp(-this->config->pharmacologic.pk.lambda_d*(fmod(t, this->config->pharmacologic.pk.tau)));

	return u;
}

double PharmacologicModel::computeDrugEffectFactor(double u) {
	double eta_drug;

	if (u < this->config->pharmacologic.pd.lowerThreshold)
		eta_drug = 0.0;
	else if (u <= this->config->pharmacologic.pd.upperThreshold)
		eta_drug = this->config->pharmacologic.pd.q_u*(u - this->config->pharmacologic.pd.lowerThreshold);
	else
		eta_drug = this->config->pharmacologic.pd.q_u*(this->config->pharmacologic.pd.upperThreshold - this->config->pharmacologic.pd.lowerThreshold);

	return eta_drug;
}

void PharmacologicModel::buildPKCurve(double delta_t) {
	this->pkSize = (int) (this->config->pharmacologic.pk.tau/delta_t);

	this->time = new double[this->pkSize];
	this->pkCurve = new double[this->pkSize];

	for (int i = 0; i < this->pkSize; i++) {
		this->time[i] = i*delta_t;
		this->pkCurve[i] = computeConcentration(this->time[i]);
	}
}

void PharmacologicModel::buildPDCurve(double delta_u) {
	this->pdSize = (int) (this->config->pharmacologic.pd.maxConcentration/delta_u);

	this->concentration = new double[this->pdSize];
	this->pdCurve = new double[this->pdSize];

	for (int i = 0; i < this->pdSize; i++) {
		this->concentration[i] = i*delta_u;
		this->pdCurve[i] = computeDrugEffectFactor(this->concentration[i]);
	}
}