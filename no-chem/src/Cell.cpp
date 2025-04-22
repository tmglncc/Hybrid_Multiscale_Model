#include "Cell.hpp"

Cell::Cell(int state, int previousState, int lifetime, 
		Vector3 coordinates, Vector3 speed, Vector3 force, 
		double nucleusRadius, double radius, double actionRadius, 
		double oConsumption, double egfSource, double chemUptake,
		double calcification, 
		double sigmaO, double sigmaEGF,
		int HIF, double mu):
state(state), previousState(previousState), lifetime(lifetime),
coordinates(coordinates), speed(speed), force(force),
nucleusRadius(nucleusRadius), radius(radius), actionRadius(actionRadius), 
oConsumption(oConsumption), egfSource(egfSource), chemUptake(chemUptake),
calcification(calcification),
sigmaO(sigmaO), sigmaEGF(sigmaEGF),
HIF(HIF), mu(mu) {
	this->grad_F_hap[0] = 0.0;
	this->grad_F_hap[1] = 0.0;

	for (size_t i = 0; i < 6; i++) {
		this->cellCycleConc[i] = 0.0;
	}
}

std::string Cell::to_string() {
	std::string out;
	switch (this->state) {
		case NECROTIC:
		out = "State = Necrotic\n";
		break;

		case QUIESCENT:
		out = "State = Quiescent (G0)\n";
		break;

		case PROLIFERATIVE2:
		out = "State = Proliferative (S-G2-M)\n";
		break;

		case HYPOXIC:
		out = "State = Hypoxic\n";
		break;

		case APOPTOTIC:
		out = "State = Apoptotic\n";
		break;

		case PROLIFERATIVE1:
		out = "State = Proliferative (G1)\n";
		break;

		case NORMOXIC:
		out = "State = Normoxic\n";
		break;

		case MIGRATE:
		out = "State = Migrate\n";
		break;

		case KILLED:
		out = "State = Killed\n";
		break;

		default:
		out = "State = " + std::to_string(this->state) + "\n";
	}

	out += "Time = " + std::to_string(this->lifetime) + "\n";
	out += "Coordinates = " + this->coordinates.to_string() + "\n";
	out += "Speed = " + this->speed.to_string() + "\n";
	out += "Force = " + this->force.to_string() + "\n";
	out += "Radius | Nucleus = " + std::to_string(this->nucleusRadius) + ", Cell = " + std::to_string(this->radius) + ", Action = " + std::to_string(this->actionRadius) + "\n";
	out += "Calcification = " + std::to_string(this->calcification) + "\n";
	
	return out;
}
