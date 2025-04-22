#ifndef PATHWAY
#define PATHWAY

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

#include "NR3.hpp"
#include "Ran.hpp"

class Pathway {
protected:
	virtual double fl(int equation, double t, double* y) = 0;
	void runge_kutta(double* yk);

	double dt, tMax;
	int eqnum;

public:
	Pathway(double dt, double tMax, int eqnum);
};

class EGFRPathway: public Pathway {
private:
	//Reaction Number
	double V1(double* y, double P1, double P2);
	double V2(double* y, double P1, double P2);
	double V3(double* y, double P1, double P2);
	double V4(double* y, double P1, double P2);
	double V5(double* y, double P1, double P2);
	double V6(double* y, double P1, double P2);
	double V7(double* y, double P1, double P2);
	double V8(double* y, double P1, double P2);
	double V9(double* y, double P1, double P2);
	double V10(double* y, double P1, double P2);
	double V11(double* y, double P1, double P2);
	double V12(double* y, double P1, double P2, double P3);
	double V13(double* y, double P1, double P2, double P3);
	double V14(double* y, double P1, double P2, double P3);
	double V15(double* y, double P1, double P2, double P3);
	double V16(double* y, double P1, double P2, double P3);
	double V17(double* y, double P1, double P2, double P3);
	double V18(double* y, double P1, double P2, double P3);
	double V19(double* y, double P1, double P2, double P3);
	double fl(int equation, double t, double* y);

	// double kV[19] = {10.8, 36.0, 3600.0, 1.62e6, 216.0, 3600.0, 1080.0, 3600.0, 3600.0, 770.4, 1.44e4, 1.26e4, 208.8, 1.044e4, 208.8, 3.42e4, 1080.0, 5.76e4, 972.0};
	// double kmK[19] = {216.0, 360.0, 36.0, 50.0, 720.0, 180.0, 21.6, 100.0, 108.0, 1.89e4, 64.0, 317.0, 2200.0, 317.0, 60.0, 1.46e5, 160.0, 1.46e5, 60.0};
	double kV[19] = {0.003, 0.01, 1.0, 450.0, 0.06, 1.0, 0.3, 1.0, 1.0, 0.214, 4.0, 3.5, 0.058, 2.9, 0.058, 9.5, 0.3, 16.0, 0.27};
	double kmK[19] = {0.06, 0.1, 0.01, 50.0, 0.2, 0.05, 0.006, 100.0, 0.03, 5.25, 64.0, 317.0, 2200.0, 317.0, 60.0, 1.46e5, 160.0, 1.46e5, 60.0};
	double IC[20] = {0, 80.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 100.0, 0.0, 120.0, 0.0, 0.0, 100.0, 0.0, 0.0};

public:
	double ROC_PLC_th = -1e-4;
	double ROC_ERK_th = 1.2e-5;

	EGFRPathway(double dt = 0.1, double tMax = 3600, int eqnum = 20);
	void runge_kutta(double sigmaEGF, double* ROC_PLC, double* ROC_ERK);
};

class CellCyclePathway: public Pathway {
private:
	//Reaction Number
	double V1(double* y, double k1, double k2_1, double k2_2, double p27p21, int HIF);
	double V2(double* y, double k3_1, double k3_2, double k4, double J3, double J4);
	double V3(double* y, double k5_1, double k5_2, double k6, double J5, int n);
	double V4(double* y, double k6, double k7, double k8, double J7, double J8, double Mad);
	double V5(double* y, double k9, double k10);
	double V6(double* y, double mu, double m_ast);
	double fl(int equation, double t, double* y);

	double k[13] = {0.12, 0.12, 4.5, 3.0, 30.0, 105.0, 0.015, 0.6, 0.3, 3.0, 1.5, 0.3, 0.06};
	double J[5] = {0.04, 0.04, 0.3, 0.001, 0.001};
	int n = 4;
	double p27p21 = 1.05;
	double Mad = 1.0;
	double m_ast = 10.0;
	int HIF;
	double mu;
	
	double IC[6] = {0.0938864628820961, 0.9978010536294271, 0.9978010536294271, 0.8755303112713486, 0.5742202676032262, 0.4519495252451473};
	double mu_plus = 0.03;
	double epsilon = 0.006;
	double mu_hat[2] = {-1.0, 1.0};

public:
	double CycB_th = 0.1;

	CellCyclePathway(double dt = 0.001, double tMax = 1.0, int eqnum = 6);
	void runge_kutta(double* cellCycleConc, int HIF, double mu);
	void setInitialAttributes(Ran* ran, double* cellCycleConc, int* HIF, double* mu);
	double generateMu(double randomNumber);
};

#endif /* end of include guard: PATHWAY */