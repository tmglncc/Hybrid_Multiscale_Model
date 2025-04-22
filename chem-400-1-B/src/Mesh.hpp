#ifndef MESH
#define MESH

#include "Vector.hpp"
#include "CRS.hpp"
#include "ConfigHandler.hpp"
#include "Chemotherapy.hpp"

class Mesh {
protected:
	virtual void createCRS(Vector3 domain, Chemotherapy* chemotherapy) = 0;
	virtual void setBoundaryConditions(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy) = 0;

public:
	double hCoarse, hRefined, hRatio;
	double deltaT;
	int matrixSize, nzNum;
	Vector3i unityCoarse, unityRefined;
	std::vector<Vector3> refined;
	std::vector<Vector3> pos;

	std::vector<double> uO;
	std::vector<double> uEGF;
	std::vector<std::vector<double>> uChem;

	double sigmaO;
	double sigmaEGF;
	std::vector<double> sigmaChem;
	std::vector<double> chemDecay;

	double* bO;
	double* bEGF;
	std::vector<double*> bChem;

	/*std::vector<int> rowPtr;
	int* colInd;
	double* val;
	int* rowPtrEGF;
	int* colIndEGF;
	double* valEGF;*/
	CRS* crsO;
	CRS* crsEGF;
	std::vector<CRS*> crsChem;

	// Mesh(Vector3 domain, double oDiffusion, double egfDiffusion, double oConsumptionBorder, double hCoarse=10.0, double hRefined=2.5, double deltaT=1.0);
	void printAttributes();
	double calculateMax_uO();
	double calculateMax_uChem(int chemDrug);
};

class Mesh3D: public Mesh {
private:
	void createCRS(Vector3 domain, Chemotherapy* chemotherapy);
	void setBoundaryConditions(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy);

public:
	Mesh3D(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy);
	~Mesh3D();
};

class Mesh2D: public Mesh {
private:
	void createCRS(Vector3 domain, Chemotherapy* chemotherapy);
	void setBoundaryConditions(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy);

public:
	Mesh2D(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy);
	~Mesh2D();
};

#endif /* end of include guard: MESH */
