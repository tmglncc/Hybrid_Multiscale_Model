#include "Pathway.hpp"

Pathway::Pathway(double dt, double tMax, int eqnum): dt(dt), tMax(tMax), eqnum(eqnum) {}

void Pathway::runge_kutta(double* yk) {
	//Parameters: RK4;
	double** K = new double*[4];
	for (int i = 0; i < 4; i++)
		K[i] = new double[this->eqnum];

	double* yref = new double[this->eqnum]; // (double*) malloc(this->eqnum * sizeof(double));

	double tk = 0.0;
	int subd = this->tMax/this->dt;
	
	for(int i = 0; i < subd; i++) {
		for(int eq = 0; eq < this->eqnum; eq++)
			yref[eq] = yk[eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[0][eq] = this->dt*fl(eq, tk, yk);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + 0.5*K[0][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[1][eq] = this->dt*fl(eq, tk, yk);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + 0.5*K[1][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[2][eq] = this->dt*fl(eq, tk, yk);

		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + K[2][eq];

		for(int eq = 0; eq < this->eqnum; eq++)
			K[3][eq] = this->dt*fl(eq, tk, yk);

		//Update
		for(int eq = 0; eq < this->eqnum; eq++)
			yk[eq] = yref[eq] + (K[0][eq] + 2.0*K[1][eq] + 2.0*K[2][eq] + K[3][eq])/6.0;

		tk += this->dt;
	}

	for (int i = 0; i < 4; i++)
		delete K[i];
	delete K;
	delete yref; // free(yref);
}

double EGFRPathway::fl(int equation, double t, double* y) {
	switch(equation) {
		case 0: return - this->V1(y, this->kV[0], this->kmK[0]);
		case 1: return - this->V1(y, this->kV[0], this->kmK[0]);
		case 2: return this->V1(y, this->kV[0], this->kmK[0]) - 2.0*this->V2(y, this->kV[1], this->kmK[1]);
		case 3: return this->V2(y, this->kV[1], this->kmK[1]) + this->V4(y, this->kV[3], this->kmK[3]) - this->V3(y, this->kV[2], this->kmK[2]);
		case 4: return this->V3(y, this->kV[2], this->kmK[2]) + this->V7(y, this->kV[6], this->kmK[6]) - this->V4(y, this->kV[3], this->kmK[3]) - this->V5(y, this->kV[4], this->kmK[4]);
		case 5: return this->V8(y, this->kV[7], this->kmK[7]) - this->V5(y, this->kV[4], this->kmK[4]);
		case 6: return this->V5(y, this->kV[4], this->kmK[4]) - this->V6(y, this->kV[5], this->kmK[5]);
		case 7: return this->V6(y, this->kV[5], this->kmK[5]) - this->V7(y, this->kV[6], this->kmK[6]);
		case 8: return this->V7(y, this->kV[6], this->kmK[6]) - this->V8(y, this->kV[7], this->kmK[7]) - this->V9(y, this->kV[8], this->kmK[8]) - this->V10(y, this->kV[9], this->kmK[9]);
		case 9: return this->V9(y, this->kV[8], this->kmK[8]);
		case 10: return - this->V10(y, this->kV[9], this->kmK[9]);
		case 11: return this->V10(y, this->kV[9], this->kmK[9]) - this->V11(y, this->kV[10], this->kmK[10]);
		case 12: return - this->V11(y, this->kV[10], this->kmK[10]);
		case 13: return this->V11(y, this->kV[10], this->kmK[10]) - this->V12(y, this->kV[11], this->kmK[11], this->kmK[13]) - this->V14(y, this->kV[13], this->kmK[13], this->kmK[11]);
		case 14: return this->V13(y, this->kV[12], this->kmK[12], this->kmK[14]) - this->V12(y, this->kV[11], this->kmK[11], this->kmK[13]);
		case 15: return this->V12(y, this->kV[11], this->kmK[11], this->kmK[13]) + this->V15(y, this->kV[14], this->kmK[14], this->kmK[12]) - this->V13(y, this->kV[12], this->kmK[12], this->kmK[14]) - this->V14(y, this->kV[13], this->kmK[13], this->kmK[11]);
		case 16: return this->V14(y, this->kV[13], this->kmK[13], this->kmK[11]) - this->V15(y, this->kV[14], this->kmK[14], this->kmK[12]) - this->V16(y, this->kV[15], this->kmK[15], this->kmK[17]) - this->V18(y, this->kV[17], this->kmK[17], this->kmK[15]);
		case 17: return this->V17(y, this->kV[16], this->kmK[16], this->kmK[18]) - this->V16(y, this->kV[15], this->kmK[15], this->kmK[17]);
		case 18: return this->V16(y, this->kV[15], this->kmK[15], this->kmK[17]) + this->V19(y, this->kV[18], this->kmK[18], this->kmK[16]) - this->V17(y, this->kV[16], this->kmK[16], this->kmK[18]) - this->V18(y, this->kV[17], this->kmK[17], this->kmK[15]);
		case 19: return this->V18(y, this->kV[17], this->kmK[17], this->kmK[15]) - this->V19(y, this->kV[18], this->kmK[18], this->kmK[16]);
		default: return 0;
	}
}

//Reaction Number: V1
double EGFRPathway::V1(double* y, double P1, double P2) {
	return P1*y[0]*y[1] - P2*y[2];
}

//Reaction Number: V2
double EGFRPathway::V2(double* y, double P1, double P2) {
	return P1*y[2]*y[2] - P2*y[3];
}

//Reaction Number: V3
double EGFRPathway::V3(double* y, double P1, double P2) {
	return P1*y[3] - P2*y[4];
}

//Reaction Number: V4
double EGFRPathway::V4(double* y, double P1, double P2) {
	return P1*y[4]/(P2 + y[4]);
}

//Reaction Number: V5
double EGFRPathway::V5(double* y, double P1, double P2) {
	return P1*y[4]*y[5] - P2*y[6];
}

//Reaction Number: V6
double EGFRPathway::V6(double* y, double P1, double P2) {
	return P1*y[6] - P2*y[7];
}

//Reaction Number: V7
double EGFRPathway::V7(double* y, double P1, double P2) {
	return P1*y[7] - P2*y[4]*y[8];
}

//Reaction Number: V8
double EGFRPathway::V8(double* y, double P1, double P2) {
	return P1*y[8]/(P2 + y[8]);
}

//Reaction Number: V9
double EGFRPathway::V9(double* y, double P1, double P2) {
	return P1*y[8] - P2*y[9];
}

//Reaction Number: V10
double EGFRPathway::V10(double* y, double P1, double P2) {
	return P1*y[8]*y[10] - P2*y[11];
}

//Reaction Number: V11
double EGFRPathway::V11(double* y, double P1, double P2) {
	return P1*y[11]*y[12]/(P2 + y[12]);
}

//Reaction Number: V12
double EGFRPathway::V12(double* y, double P1, double P2, double P3) {
	return P1*y[13]*y[14]/(P2*(1.0 + y[15]/P3) + y[14]);
}

//Reaction Number: V13
double EGFRPathway::V13(double* y, double P1, double P2, double P3) {
	return P1*y[15]/(P2*(1.0 + y[16]/P3) + y[15]);
}

//Reaction Number: V14
double EGFRPathway::V14(double* y, double P1, double P2, double P3) {
	return P1*y[13]*y[15]/(P2*(1.0 + y[14]/P3) + y[15]);
}

//Reaction Number: V15
double EGFRPathway::V15(double* y, double P1, double P2, double P3) {
	return P1*y[16]/(P2*(1.0 + y[15]/P3) + y[16]);
}

//Reaction Number: V16
double EGFRPathway::V16(double* y, double P1, double P2, double P3) {
	return P1*y[16]*y[17]/(P2*(1.0 + y[18]/P3) + y[17]);
}

//Reaction Number: V17
double EGFRPathway::V17(double* y, double P1, double P2, double P3) {
	return P1*y[18]/(P2*(1.0 + y[19]/P3) + y[18]);
}

//Reaction Number: V18
double EGFRPathway::V18(double* y, double P1, double P2, double P3) {
	return P1*y[16]*y[18]/(P2*(1.0 + y[17]/P3) + y[18]);
}

//Reaction Number: V19
double EGFRPathway::V19(double* y, double P1, double P2, double P3) {
	return P1*y[19]/(P2*(1.0 + y[18]/P3) + y[19]);
}

EGFRPathway::EGFRPathway(double dt, double tMax, int eqnum): Pathway(dt, tMax, eqnum) {}

///Method: Runge-Kutta 4th Order
void EGFRPathway::runge_kutta(double sigmaEGF, double* ROC_PLC, double* ROC_ERK) {
	double* yk = new double[this->eqnum]; // (double*) malloc(this->eqnum * sizeof(double));
	for (size_t i = 0; i < this->eqnum; i++) {
		yk[i] = this->IC[i];
	}
	yk[0] = sigmaEGF;
	
	Pathway::runge_kutta(yk);

	*ROC_PLC = V8(yk, 1.0, 100.0) - V5(yk, 0.06, 0.2);
	*ROC_ERK = V17(yk, 0.3, 160.0, 60.0) - V16(yk, 9.5, 1.46e5, 1.46e5);

	delete yk; // free(yk);
}

double CellCyclePathway::fl(int equation, double t, double* y) {
	switch(equation) {
		case 0: return this->V1(y, this->k[0], this->k[1], this->k[2], this->p27p21, this->HIF);
		case 1: return this->V2(y, this->k[3], this->k[4], this->k[5], this->J[0], this->J[1]);
		case 2: return this->V3(y, this->k[6], this->k[7], this->k[8], this->J[2], this->n);
		case 3: return this->V4(y, this->k[8], this->k[9], this->k[10], this->J[3], this->J[4], this->Mad);
		case 4: return this->V5(y, this->k[11], this->k[12]);
		case 5: return this->V6(y, this->mu, this->m_ast);
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

CellCyclePathway::CellCyclePathway(double dt, double tMax, int eqnum): Pathway(dt, tMax, eqnum) {}

///Method: Runge-Kutta 4th Order
void CellCyclePathway::runge_kutta(double* cellCycleConc, int HIF, double mu) {
	double* yk = new double[this->eqnum]; // (double*) malloc(this->eqnum * sizeof(double));
	for (size_t i = 0; i < this->eqnum; i++) {
		yk[i] = cellCycleConc[i];
	}

	this->HIF = HIF;
	this->mu = mu;

	Pathway::runge_kutta(yk);

	for (size_t i = 0; i < this->eqnum; i++) {
		cellCycleConc[i] = yk[i];
	}

	delete yk; // free(yk);
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
