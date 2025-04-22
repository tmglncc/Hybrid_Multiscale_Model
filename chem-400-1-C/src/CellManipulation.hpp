#ifndef CELL_MANIPULATION
#define CELL_MANIPULATION

#include "NR3.hpp"
#include "Ran.hpp"

#include <vector>
#include <string>
#include <cmath>

#include "Cell.hpp"
#include "Vector.hpp"
#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "Pathway.hpp"
#include "Mesh.hpp"
#include "Ran.hpp"
#include "Chemotherapy.hpp"
#include "../src_backup/PharmacologicModel.hpp"

#include <iostream>
#include <algorithm>

class CellManipulation {
protected:
	static void changeCellState(Cell* cell, int previousState, int state, int lifetime, double oConsumption, double egfSource, double chemUptake);
	static double norma(Vector3 pos);
	static Vector3 func_var_phi(Vector3 normal, double actionRadius, int n);
	static Vector3 func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m);

public:
	static void updateInitialFrame(Frame* frame, ConfigHandler* config, Ran* ran);
	static double* generateF_hap(Frame* frame, Mesh* mesh, Ran* ran);
	static void calculateGrad_F_hap(Frame* frame, Mesh* mesh, double* F_hap);
};

class CellManipulation3D: public CellManipulation {
private:
	static void calculateCellSigmas(Cell* cell, Mesh* mesh);
	static Cell divide(Cell* cell, double rand1, double rand2);
	static Vector3 normal(Vector3 coordinates, double domainRadius);

public:
	static void updateFrame(Frame* frame, ConfigHandler* config, Mesh* mesh, EGFRPathway* egfrPathway, CellCyclePathway* cellCyclePathway, Ran* ran, Chemotherapy* chemotherapy, PharmacologicModel* pharmacologic = NULL);
	static void force(Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel);
};

class CellManipulation2D: public CellManipulation {
private:
	static void calculateCellSigmas(Cell* cell, Mesh* mesh);
	static Cell divide(Cell* cell, double rand1, CellCyclePathway* cellCyclePathway);
	static Vector3 normal(Vector3 coordinates, double domainRadius);

public:
	static void updateFrame(Frame* frame, ConfigHandler* config, Mesh* mesh, EGFRPathway* egfrPathway, CellCyclePathway* cellCyclePathway, Ran* ran, Chemotherapy* chemotherapy, PharmacologicModel* pharmacologic = NULL);
	static void force(Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel);

	// static double norma(Vector3 pos);
	// static Vector3 func_var_phi(Vector3 normal, double actionRadius, int n);
	// static Vector3 func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m);
};

#endif /* end of include guard: CELL_MANIPULATION */
