#include "CellCyclePathway.hpp"

using namespace std;

double CellCyclePathway::fl(int equation, double* y, int HIF, double mu) {
	switch(equation) {
		case 0: return this->V1(y, this->k[0], this->k[1], this->k[2], this->p27p21, HIF);
		case 1: return this->V2(y, this->k[3], this->k[4], this->k[5], this->J[0], this->J[1]);
		case 2: return this->V3(y, this->k[6], this->k[7], this->k[8], this->J[2], this->n);
		case 3: return this->V4(y, this->k[8], this->k[9], this->k[10], this->J[3], this->J[4], this->Mad);
		case 4: return this->V5(y, this->k[11], this->k[12]);
		case 5: return this->V6(y, mu, this->m_ast);
		default: return 0;
	}
}

//Reaction Number: V1
double CellCyclePathway::V1(double* y, double k1, double k2_1, double k2_2, double p27p21, int HIF) {
	return k1 - (k2_1 + k2_2*y[1] + p27p21*HIF)*y[0];
}

//Reaction Number: V2
double CellCyclePathway::V2(double* y, double k3_1, double k3_2, double k4, double J3, double J4) {
	return ((k3_1 + k3_2*y[3])*(1.0 - y[1]))/(J3 + 1.0 - y[1]) - (k4*y[5]*y[0]*y[1])/(J4 + y[1]);
}

//Reaction Number: V3
double CellCyclePathway::V3(double* y, double k5_1, double k5_2, double k6, double J5, int n) {
	return k5_1 + k5_2*(pow(y[0]*y[5], n))/(pow(J5, n) + pow(y[0]*y[5], n)) - k6*y[2];
}

//Reaction Number: V4
double CellCyclePathway::V4(double* y, double k6, double k7, double k8, double J7, double J8, double Mad) {
	return (k7*y[4]*(y[2] - y[3]))/(J7 + y[2] - y[3]) - (k8*Mad*y[3])/(J8 + y[3]) - k6*y[3];
}

//Reaction Number: V5
double CellCyclePathway::V5(double* y, double k9, double k10) {
	return k9*y[5]*y[0]*(1.0 - y[4]) - k10*y[4];
}

//Reaction Number: V6
double CellCyclePathway::V6(double* y, double mu, double m_ast) {
	return mu*y[5]*(1.0 - y[5]/m_ast);
}

CellCyclePathway::CellCyclePathway(double dt, double tMax, int eqnum): dt(dt), tMax(tMax), eqnum(eqnum) {}

///Method: Runge-Kutta 4th Order
void CellCyclePathway::runge_kutta(double* IC, int HIF, double mu) {
	//Parameters: RK4;
	double** K = new double*[4];
	for (int i = 0; i < 4; i++)
		K[i] = new double[this->eqnum];

	double* yk = new double[this->eqnum]; // (double*) malloc(this->eqnum * sizeof(double));
	for (size_t i = 0; i < this->eqnum; i++) {
		yk[i] = IC[i];
	}
	double* yref = new double[this->eqnum]; // (double*) malloc(this->eqnum * sizeof(double));

	double tk = 0.0;
	int subd = this->tMax/this->dt;
	
	for(int i = 0; i < subd; i++) {
		for(int eq = 0; eq < this->eqnum; eq++)
			yref[eq] = yk[eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[0][eq] = this->dt*fl(eq, yk, HIF, mu);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + 0.5*K[0][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[1][eq] = this->dt*fl(eq, yk, HIF, mu);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + 0.5*K[1][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[2][eq] = this->dt*fl(eq, yk, HIF, mu);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + K[2][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[3][eq] = this->dt*fl(eq, yk, HIF, mu);

		//Update
		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + (K[0][eq] + 2.0*K[1][eq] + 2.0*K[2][eq] + K[3][eq])/6.0;

		tk += this->dt;
	}

	for (size_t i = 0; i < this->eqnum; i++) {
		IC[i] = yk[i];
	}

	for (int i = 0; i < 4; i++)
		delete K[i];
	delete K;
	delete yk; // free(yk);
	delete yref; // free(yref);
}

void CellCyclePathway::setInitialAttributes(Ran* ran, double* cellCycleConc, int* HIF, double* mu) {
	for (size_t i = 0; i < this->eqnum; i++) {
		cellCycleConc[i] = this->IC[i];
	}

	*HIF = 0;
	*mu = this->generateMu(ran->doub());
}

double CellCyclePathway::generateMu(double randomNumber) {
	return this->mu_plus + this->epsilon*((this->mu_hat[1] - this->mu_hat[0])*randomNumber + this->mu_hat[0]);
}
