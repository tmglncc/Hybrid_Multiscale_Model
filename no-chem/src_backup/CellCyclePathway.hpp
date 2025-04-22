#ifndef CELL_CYCLE_PATHWAY
#define CELL_CYCLE_PATHWAY

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

#include "NR3.hpp"
#include "Ran.hpp"

class CellCyclePathway {
private:
	//Reaction Number
	double V1(double* y, double k1, double k2_1, double k2_2, double p27p21, int HIF);
	double V2(double* y, double k3_1, double k3_2, double k4, double J3, double J4);
	double V3(double* y, double k5_1, double k5_2, double k6, double J5, int n);
	double V4(double* y, double k6, double k7, double k8, double J7, double J8, double Mad);
	double V5(double* y, double k9, double k10);
	double V6(double* y, double mu, double m_ast);
	double fl(int equation, double* y, int HIF, double mu);

	double k[13] = {0.12, 0.12, 4.5, 3.0, 30.0, 105.0, 0.015, 0.6, 0.3, 3.0, 1.5, 0.3, 0.06};
	double J[5] = {0.04, 0.04, 0.3, 0.001, 0.001};
	int n = 4;
	double p27p21 = 1.05;
	double Mad = 1.0;
	double m_ast = 10.0;
	
	double IC[6] = {0.0938864628820961, 0.9978010536294271, 0.9978010536294271, 0.8755303112713486, 0.5742202676032262, 0.4519495252451473};
	double mu_plus = 0.03;
	double epsilon = 0.006;
	double mu_hat[2] = {-1.0, 1.0};

	double dt, tMax;
	int eqnum;

public:
	double CycB_th = 0.1;

	CellCyclePathway(double dt = 0.001, double tMax = 1.0, int eqnum = 6);
	void runge_kutta(double* IC, int HIF, double mu);
	void setInitialAttributes(Ran* ran, double* cellCycleConc, int* HIF, double* mu);
	double generateMu(double randomNumber);
};

#endif /* end of include guard: PATHWAY */
