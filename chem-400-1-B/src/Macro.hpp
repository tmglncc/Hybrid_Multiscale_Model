#ifndef MACRO
#define MACRO

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "Cell.hpp"
#include "Mesh.hpp"
#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "CRS.hpp"
#include "mgmres.hpp"
#include "FileFactory.hpp"
#include "Chemotherapy.hpp"

class Macro {
protected:
	double* oProduction;
	double* oUptake; // remove
	double* egfSource; // remove
	std::vector<double*> chemSupply;
	std::vector<double*> chemUptake;

	Mesh* mesh;
	Frame* frame;
	ConfigHandler* config;
	std::vector<Vector3> bloodVessel;
	Chemotherapy* chemotherapy;
	
public:
	Macro(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy);
	~Macro();
	virtual void reaction() = 0;
	virtual void differenceO(int itr_max = 1000, int mr = 5, double tol_abs = 1.0e-5, double tol_rel = 1.0e-5) = 0;
	virtual void differenceEGF(int itr_max = 1000, int mr = 5, double tol_abs = 1.0e-5, double tol_rel = 1.0e-5) = 0;
	virtual void differenceChem(int itr_max = 1000, int mr = 5, double tol_abs = 1.0e-5, double tol_rel = 1.0e-5) = 0;

	/*int m(Vector3 coordinates, std::vector<Vector3> bloodVessel);
	int cell(Vector3 coordinates, std::vector<Cell> cells);*/
};

class Macro3D: public Macro {
private:
	std::vector<std::vector<std::vector<double>>> oUptakeMesh, egfSourceMesh;

	void homogenize(int kBegin, int kEnd, int lBegin, int lEnd, int jBegin, int jEnd, int i, double homogenizationFactor);

public:
	Macro3D(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy);
	void reaction();
	void differenceO(int itr_max, int mr, double tol_abs, double tol_rel);
	void differenceEGF(int itr_max, int mr, double tol_abs, double tol_rel);
	void differenceChem(int itr_max, int mr, double tol_abs, double tol_rel);
};

class Macro2D: public Macro {
private:
	std::vector<std::vector<double>> oProductionMesh, oUptakeMesh, egfSourceMesh;
	std::vector<std::vector<std::vector<double>>> chemSupplyMesh, chemUptakeMesh;

	void homogenize(int lBegin, int lEnd, int jBegin, int jEnd, int i, double homogenizationFactor);

public:
	Macro2D(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy);
	void reaction();
	void differenceO(int itr_max, int mr, double tol_abs, double tol_rel);
	void differenceEGF(int itr_max, int mr, double tol_abs, double tol_rel);
	void differenceChem(int itr_max, int mr, double tol_abs, double tol_rel);
};

#endif /* end of include guard: MACRO */
