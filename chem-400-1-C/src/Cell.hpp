#ifndef CELL
#define CELL

#include <vector>

#include "Vector.hpp"
#include "Util.hpp"

class Cell {
public:
	int state, previousState, lifetime;
	Vector3 coordinates, speed, force;
	double nucleusRadius, radius, actionRadius;
	double oConsumption, egfSource, chemUptake;
	double calcification;
	double sigmaO, sigmaEGF;
	std::vector<double> sigmaChem;
	double grad_F_hap[2];

	double cellCycleConc[6];
	int HIF;
	double mu;

	Cell(int state = -1, int previousState = 0, int lifetime = 0, 
		Vector3 coordinates = Vector3(), Vector3 speed = Vector3(), Vector3 force = Vector3(), 
		double nucleusRadius = 0.0, double radius = 0.0, double actionRadius = 0.0, 
		double oConsumption = 0.0, double egfSource = 0.0, double chemUptake = 0.0,
		double calcification = 0.0, 
		double sigmaO = 0.0, double sigmaEGF = 0.0,
		int HIF = -1, double mu = 0.0);
	std::string to_string();

};

#endif /* end of include guard: CELL */
