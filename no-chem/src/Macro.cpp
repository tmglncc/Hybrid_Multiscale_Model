#include "Macro.hpp"

Macro::Macro(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy): mesh(mesh), frame(frame), config(config), bloodVessel(bloodVessel), chemotherapy(chemotherapy) {
	this->oProduction = new double[this->mesh->matrixSize];
	this->oUptake = new double[this->mesh->matrixSize];
	this->egfSource = new double[this->mesh->matrixSize];

	this->chemSupply.resize(this->chemotherapy->chemDrugs.size());
	this->chemUptake.resize(this->chemotherapy->chemDrugs.size());
	for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
		this->chemSupply[chemDrug] = new double[this->mesh->matrixSize];
		this->chemUptake[chemDrug] = new double[this->mesh->matrixSize];
	}
}

Macro::~Macro() {
	delete this->oProduction;
	delete this->oUptake;
	delete this->egfSource;

	for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
		delete this->chemSupply[chemDrug];
		delete this->chemUptake[chemDrug];
	}
}

/*int Macro::m(Vector3 coordinates, std::vector<Vector3> bloodVessel) {
	for (int i = 0; i < bloodVessel.size(); i++) {
		if (sqrt(pow(coordinates.x - bloodVessel[i].x, 2) + pow(coordinates.y - bloodVessel[i].y, 2) + pow(coordinates.z - bloodVessel[i].z, 2)) <= this->config->bloodVessel.radius)
			return 1;
	}

	return 0;
}

int Macro::cell(Vector3 coordinates, std::vector<Cell> cells) {
	for (int i = 0; i < cells.size(); i++) {
		if (cells[i].state != NORMOXIC &&
			sqrt(pow(coordinates.x - cells[i].coordinates.x, 2) + pow(coordinates.y - cells[i].coordinates.y, 2) + pow(coordinates.z - cells[i].coordinates.z, 2)) <= cells[i].radius)
			return 1;
	}

	return 0;
}*/

// GUSTAVO: Conferir depois
Macro3D::Macro3D(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy): Macro(mesh, frame, config, bloodVessel, chemotherapy) {
	this->oUptakeMesh.resize(this->mesh->unityRefined.x);
	this->egfSourceMesh.resize(this->mesh->unityRefined.x);
	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		this->oUptakeMesh[i].resize(this->mesh->unityRefined.y);
		this->egfSourceMesh[i].resize(this->mesh->unityRefined.y);
		for (int j = 0; j < this->mesh->unityRefined.y; ++j) {
			this->oUptakeMesh[i][j].resize(this->mesh->unityRefined.z);
			this->egfSourceMesh[i][j].resize(this->mesh->unityRefined.z);
		}
	}
}

void Macro3D::differenceO(int itr_max, int mr, double tol_abs, double tol_rel) {
	// double sigma = this->config->continuum.oxygen.oDiffusion/pow(this->mesh->hCoarse,2);

	if (config->output.prints) std::cout << "Oxygen Concentration..." << std::endl;

	for (int i = 0; i < this->mesh->matrixSize; i++){
		if (((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0 &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse != this->frame->domain.x &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0 &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse!=this->frame->domain.y &&
			(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse != 0 &&
			(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse!=this->frame->domain.z)
			this->mesh->crsO->set(i, i, this->oUptake[i]+6*this->mesh->sigmaO);
	}
	// double* a = &this->mesh->uO[0];
	//Calcula U[]
	MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, &mesh->crsO->rowPtr[0], mesh->crsO->colInd, mesh->crsO->val, &this->mesh->uO[0], mesh->bO, itr_max, mr, tol_abs, tol_rel);

	if(this->config->output.nut && this->frame->time%this->config->output.filesFrequency == 0)
		FileFactory::makeFile(this->frame, this->config, this->mesh, this->chemotherapy, NUT);
}

void Macro3D::differenceEGF(int itr_max, int mr, double tol_abs, double tol_rel) {
	double time = 0.0, timeMax = 1.0;

	// double * B = new double[this->mesh->matrixSize];

	if (config->output.prints) std::cout << "EGF Concentration..." << std::endl;

	//Condição inicial
	this->mesh->uEGF = std::vector<double>(this->mesh->matrixSize, this->config->continuum.egf.egfSourceBorder);
	//Testa se tem EGF
	for(int i=0;i<this->mesh->matrixSize;i++){
		if (this->egfSource[i] != 0){ time = 0; break;}
		else{time = timeMax;}
	}

	while (time < timeMax){
		for(int i = 0; i < this->mesh->matrixSize; i++){
			this->mesh->bEGF[i] = this->mesh->uEGF[i] + mesh->deltaT*this->egfSource[i];
			if (((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse == 0 ||
				((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse == this->frame->domain.x ||
				((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse == 0 ||
				((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse == this->frame->domain.y ||
				(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse == 0 ||
				(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse == this->frame->domain.z)
				this->mesh->bEGF[i]=this->config->continuum.egf.egfSourceBorder;
		}

		MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, mesh->crsEGF->rowPtr, mesh->crsEGF->colInd, mesh->crsEGF->val, &this->mesh->uEGF[0], this->mesh->bEGF, itr_max, mr, tol_abs, tol_rel);

		time += mesh->deltaT;
	}

	if(this->config->output.egf && this->frame->time%this->config->output.filesFrequency == 0)
		FileFactory::makeFile(this->frame, this->config, this->mesh, this->chemotherapy, EGF);
}

void Macro3D::differenceChem(int itr_max, int mr, double tol_abs, double tol_rel) {

}

void Macro3D::reaction() {
	if (config->output.prints) std::cout << "Reaction/Source term ..." << std::endl;

	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		for (int j = 0; j < this->mesh->unityRefined.y; ++j) {
			for (int k = 0; k < this->mesh->unityRefined.z; ++k) {
				this->oUptakeMesh[i][j][k] = 0.0;
				this->egfSourceMesh[i][j][k] = 0.0;
			}
		}
	}

	//Calcula termo de captação
	for(int l=0;l<this->frame->cells.size();l++){
		for(int k = ( (this->frame->cells[l].coordinates.z - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;k < (this->frame->cells[l].coordinates.z + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;k++){
			for(int i = ( (this->frame->cells[l].coordinates.y - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;i< (this->frame->cells[l].coordinates.y + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;i++){
				for(int j = ( (this->frame->cells[l].coordinates.x - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;j< (this->frame->cells[l].coordinates.x + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;j++){
					if (k >= this->mesh->unityRefined.z || i >= this->mesh->unityRefined.y || j >= this->mesh->unityRefined.x || k < 0 || i < 0 || j < 0) continue;
					if ((this->frame->cells[l].state != NORMOXIC) && (sqrt(pow(this->mesh->refined[j].x-this->frame->cells[l].coordinates.x,2) + pow(this->mesh->refined[i].y-this->frame->cells[l].coordinates.y,2) + pow(this->mesh->refined[k].z-this->frame->cells[l].coordinates.z,2)) <= this->frame->cells[l].radius)){
						this->oUptakeMesh[k][i][j] += this->frame->cells[l].oConsumption; // EXECUTAR OS DOIS
						this->egfSourceMesh[k][i][j] += this->frame->cells[l].egfSource;
					}
				}
			}
		}
	}

	//Calcula o termo reação na malha do nutriente
	double x, y, z;
	for(int i = 0; i < this->mesh->matrixSize; i++){
		x = ((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		y = ((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		z = (i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse;
		//cout << "Nó = " << i << " PONTO = (" << x << ", " << y << ", " << z << ") ";

		//Nós internos
		if (x != 0.0 && x != this->frame->domain.x && y != 0.0 && y != this->frame->domain.y && z != 0.0 && z != this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, pow(this->mesh->hRatio+1,3));
			//cout << "Nós internos!" << endl;
			continue;
		}

		//Nós (---,---,0) ■ Face 1
		if (x != 0.0 && x != this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 1!" << endl;
			continue;
		}

		//Nós (0,---,---) ■ Face 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0 && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 2!" << endl;
			continue;
		}

		//Nós (---,0,---) ■ Face 3
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 3!" << endl;
			continue;
		}

		//Nós (this->frame->domain.x,---,---) ■ Face 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 4!" << endl;
			continue;
		}

		//Nós (---,this->frame->domain.y,---) ■ Face 5
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 5!" << endl;
			continue;
		}

		//Nós (---,---,this->frame->domain.z) ■ Face 6
		if (x != 0.0 && x != this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(this->mesh->hRatio+1,2)*((this->mesh->hRatio/2)+1)));
			//cout << "Face 6!" << endl;
			continue;
		}

		//Nós (---,0,0) -- Aresta 1
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 1!" << endl;
			continue;
		}

		//Nós (0,---,0) -- Aresta 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 2!" << endl;
			continue;
		}

		//Nós (---,this->frame->domain.y,0) -- Aresta 3
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 3!" << endl;
			continue;
		}

		//Nós (this->frame->domain.x,---,0) -- Aresta 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 4!" << endl;
			continue;
		}

		//Nós (---,0,this->frame->domain.z) -- Aresta 5
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 6!" << endl;
			continue;
		}

		//Nós (---,this->frame->domain.y,this->frame->domain.z) -- Aresta 7
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 7!" << endl;
			continue;
		}

		//Nós (this->frame->domain.x,---,this->frame->domain.z) -- Aresta 8
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 8!" << endl;
			continue;
		}

		//Nós (0,0,---) -- Aresta 9
		if (x == 0.0 && y == 0.0 && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 9!" << endl;
			continue;
		}

		//Nós (0,this->frame->domain.y,---) -- Aresta 10
		if (x == 0.0 && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 10!" << endl;
			continue;
		}

		//Nós (this->frame->domain.x,0,---) -- Aresta 11
		if (x == this->frame->domain.x && y == 0.0 && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 11!" << endl;
			continue;
		}

		//Nós (this->frame->domain.x,this->frame->domain.y,---) -- Aresta 12
		if (x == this->frame->domain.x && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (pow(((this->mesh->hRatio/2))+1,2)*(this->mesh->hRatio+1)));
			//cout << "Aresta 12!" << endl;
			continue;
		}

		//Nó (0,0,0) * Ponto 1
		if (x == 0.0 && y == 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 1!" << endl;
			continue;
		}

		//Nó (this->frame->domain.x,0,0) * Ponto 2
		if (x == this->frame->domain.x && y == 0.0 && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 2!" << endl;
			continue;
		}

		//Nó (0,this->frame->domain.y,0) * Ponto 3
		if (x == 0.0 && y == this->frame->domain.y && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 3!" << endl;
			continue;
		}

		//Nó (this->frame->domain.x,this->frame->domain.y,0) * Ponto 4
		if (x == this->frame->domain.x && y == this->frame->domain.y && z == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 5!" << endl;
			continue;
		}

		//Nó (this->frame->domain.x,0,this->frame->domain.z) * Ponto 6
		if (x == this->frame->domain.x && y == 0.0 && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 6!" << endl;
			continue;
		}

		//Nó (0,this->frame->domain.y,this->frame->domain.z) * Ponto 7
		if (x == 0.0 && y == this->frame->domain.y && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 7!" << endl;
			continue;
		}

		//Nó (this->frame->domain.x,this->frame->domain.y,this->frame->domain.z) * Ponto 8
		if (x == this->frame->domain.x && y == this->frame->domain.y && z == this->frame->domain.z) {
			this->homogenize((z/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (z/this->mesh->hCoarse)*this->mesh->hRatio,
				(y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,3));
			//cout << "Ponto 8!" << endl;
			continue;
		}
	}
}

void Macro3D::homogenize(int kBegin, int kEnd, int lBegin, int lEnd, int jBegin, int jEnd, int i, double homogenizationFactor) {
	double oUptakeTemp = 0.0;
	double egfSourceTemp = 0.0;

	for (int k = kBegin; k <= kEnd; k++) {
		for (int l = lBegin; l <= lEnd; l++) {
			for (int j = jBegin; j <= jEnd; j++) {
				if (this->oUptakeMesh[k][l][j] == 0)
					this->oUptakeMesh[k][l][j] = this->config->continuum.oxygen.oConsumptionBg;
				oUptakeTemp += this->oUptakeMesh[k][l][j];

				if (this->egfSourceMesh[k][l][j] == 0)
					this->egfSourceMesh[k][l][j] = this->config->continuum.egf.egfSourceBg;
				egfSourceTemp += this->egfSourceMesh[k][l][j];
			}
		}
	}
	
	this->oUptake[i] = oUptakeTemp/homogenizationFactor;
	this->egfSource[i] = egfSourceTemp/homogenizationFactor;
}

Macro2D::Macro2D(Mesh* mesh, Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel, Chemotherapy* chemotherapy): Macro(mesh, frame, config, bloodVessel, chemotherapy) {
	this->oProductionMesh.resize(this->mesh->unityRefined.x);
	this->oUptakeMesh.resize(this->mesh->unityRefined.x);
	this->egfSourceMesh.resize(this->mesh->unityRefined.x);
	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		this->oProductionMesh[i].resize(this->mesh->unityRefined.y);
		this->oUptakeMesh[i].resize(this->mesh->unityRefined.y);
		this->egfSourceMesh[i].resize(this->mesh->unityRefined.y);
	}

	this->chemSupplyMesh.resize(this->chemotherapy->chemDrugs.size());
	this->chemUptakeMesh.resize(this->chemotherapy->chemDrugs.size());
	for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
		this->chemSupplyMesh[chemDrug].resize(this->mesh->unityRefined.x);
		this->chemUptakeMesh[chemDrug].resize(this->mesh->unityRefined.x);
		for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
			this->chemSupplyMesh[chemDrug][i].resize(this->mesh->unityRefined.y);
			this->chemUptakeMesh[chemDrug][i].resize(this->mesh->unityRefined.y);
		}
	}
}

void Macro2D::differenceO(int itr_max, int mr, double tol_abs, double tol_rel) {
	if (config->output.prints)
		std::cout << "Oxygen Concentration..." << std::endl;

	// double sigma = this->config->continuum.oxygen.oDiffusion/pow(this->mesh->hCoarse, 2); // Rocha, 2018
	// double sigma = this->mesh->deltaT*this->config->continuum.oxygen.oDiffusion/pow(this->mesh->hCoarse, 2); // Powathil, 2012

	// Rocha, 2018
	/*for (int i = 0; i < this->mesh->matrixSize; i++) {
		// GUSTAVO:
		// this->mesh->pos[i].y == (i/this->mesh->unityCoarse.x)*this->mesh->hCoarse
		// this->mesh->pos[i].x == (i%this->mesh->unityCoarse.x)*this->mesh->hCoarse
		if (this->mesh->pos[i].x != 0 && this->mesh->pos[i].x != this->frame->domain.x &&
			this->mesh->pos[i].y != 0 && this->mesh->pos[i].y != this->frame->domain.y)
			this->mesh->crsO->set(i, i, this->oUptake[i] + 4.0*this->mesh->sigmaO);
	}

	// double* a = &this->mesh->uO[0];
	//Calcula U[]
	MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, &this->mesh->crsO->rowPtr[0], this->mesh->crsO->colInd, this->mesh->crsO->val, &this->mesh->uO[0], this->mesh->bO, itr_max, mr, tol_abs, tol_rel);*/

	// Powathil, 2012
	double time = 0.0, timeMax = 1.0;

	// double* B = new double[this->mesh->matrixSize];

	//Condição inicial
	// this->mesh->uO = std::vector<double>(this->mesh->matrixSize, 0.0);

	while (time < timeMax) {
		for (int i = 0; i < this->mesh->matrixSize; i++) {
			this->mesh->crsO->set(i, i, 1.0 + 4.0*this->mesh->sigmaO + this->mesh->deltaT*this->oUptake[i]);

			this->mesh->bO[i] = this->mesh->deltaT*this->oProduction[i] + this->mesh->uO[i];
			if (this->mesh->pos[i].x == 0 || this->mesh->pos[i].y == 0)
				this->mesh->bO[i] += -2.0*this->mesh->sigmaO*this->mesh->hCoarse*this->config->continuum.oxygen.oConsumptionBorder;
			else if (this->mesh->pos[i].x == this->frame->domain.x || this->mesh->pos[i].y == this->frame->domain.y)
				this->mesh->bO[i] += 2.0*this->mesh->sigmaO*this->mesh->hCoarse*this->config->continuum.oxygen.oConsumptionBorder;
		}

		MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, &this->mesh->crsO->rowPtr[0], this->mesh->crsO->colInd, this->mesh->crsO->val, &this->mesh->uO[0], this->mesh->bO, itr_max, mr, tol_abs, tol_rel);

		time += this->mesh->deltaT;
	}

	// Powathil, 2012 (Quasi-stationary)
	/*// double* B = new double[this->mesh->matrixSize];

	for (int i = 0; i < this->mesh->matrixSize; i++) {
		this->mesh->crsO->set(i, i, 4.0*this->mesh->sigmaO + this->oUptake[i]);

		this->mesh->bO[i] = this->oProduction[i];
		if (this->mesh->pos[i].x == 0 || this->mesh->pos[i].y == 0)
			this->mesh->bO[i] += -2.0*this->mesh->sigmaO*this->mesh->hCoarse*this->config->continuum.oxygen.oConsumptionBorder;
		else if (this->mesh->pos[i].x == this->frame->domain.x || this->mesh->pos[i].y == this->frame->domain.y)
			this->mesh->bO[i] += 2.0*this->mesh->sigmaO*this->mesh->hCoarse*this->config->continuum.oxygen.oConsumptionBorder;
	}

	MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, &this->mesh->crsO->rowPtr[0], this->mesh->crsO->colInd, this->mesh->crsO->val, &this->mesh->uO[0], this->mesh->bO, itr_max, mr, tol_abs, tol_rel);*/

	if (this->config->output.nut && this->frame->time%this->config->output.filesFrequency == 0)
		FileFactory::makeFile(this->frame, this->config, this->mesh, this->chemotherapy, NUT);
}

void Macro2D::differenceEGF(int itr_max, int mr, double tol_abs, double tol_rel) {
	if (config->output.prints)
		std::cout << "EGF Concentration..." << std::endl;

	double time = 0.0, timeMax = 1.0;

	// double* B = new double[this->mesh->matrixSize];

	//Condição inicial
	this->mesh->uEGF = std::vector<double>(this->mesh->matrixSize, this->config->continuum.egf.egfSourceBorder);

	//Testa se tem EGF
	for (int i = 0; i < this->mesh->matrixSize; i++) {
		if (this->egfSource[i] != 0) {
			time = 0.0;
			break;
		} else
			time = timeMax;
	}

	while (time < timeMax) {
		for (int i = 0; i < this->mesh->matrixSize; i++) {
			this->mesh->bEGF[i] = this->mesh->uEGF[i] + this->mesh->deltaT*this->egfSource[i];
			if (this->mesh->pos[i].x == 0 || this->mesh->pos[i].x == this->frame->domain.x ||
				this->mesh->pos[i].y == 0 || this->mesh->pos[i].y == this->frame->domain.y)
				this->mesh->bEGF[i] = this->config->continuum.egf.egfSourceBorder;
		}

		MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, this->mesh->crsEGF->rowPtr, this->mesh->crsEGF->colInd, this->mesh->crsEGF->val, &this->mesh->uEGF[0], this->mesh->bEGF, itr_max, mr, tol_abs, tol_rel);

		time += this->mesh->deltaT;
	}

	if (this->config->output.egf && this->frame->time%this->config->output.filesFrequency == 0)
		FileFactory::makeFile(this->frame, this->config, this->mesh, this->chemotherapy, EGF);
}

void Macro2D::differenceChem(int itr_max, int mr, double tol_abs, double tol_rel) {
	if (config->output.prints)
		std::cout << "Chemotherapy Concentration..." << std::endl;

	// Powathil, 2012
	for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
		double time = 0.0, timeMax = 1.0;

		for (int i = 0; i < this->mesh->matrixSize; i++) {
			if (this->mesh->uChem[chemDrug][i] != 0 || this->chemSupply[chemDrug][i] != 0) {
				time = 0.0;
				break;
			} else
				time = timeMax;
		}

		while (time < timeMax) {
			for (int i = 0; i < this->mesh->matrixSize; i++) {
				this->mesh->crsChem[chemDrug]->set(i, i, 1.0 + 4.0*this->mesh->sigmaChem[chemDrug] + this->mesh->deltaT*this->chemUptake[chemDrug][i] + this->mesh->deltaT*this->mesh->chemDecay[chemDrug]);

				this->mesh->bChem[chemDrug][i] = this->mesh->deltaT*this->chemSupply[chemDrug][i] + this->mesh->uChem[chemDrug][i];
				if (this->mesh->pos[i].x == 0 || this->mesh->pos[i].y == 0)
					this->mesh->bChem[chemDrug][i] += -2.0*this->mesh->sigmaChem[chemDrug]*this->mesh->hCoarse*this->config->continuum.chemotherapy.chemUptakeBorder;
				else if (this->mesh->pos[i].x == this->frame->domain.x || this->mesh->pos[i].y == this->frame->domain.y)
					this->mesh->bChem[chemDrug][i] += 2.0*this->mesh->sigmaChem[chemDrug]*this->mesh->hCoarse*this->config->continuum.chemotherapy.chemUptakeBorder;
			}

			MGMRES::pmgmres_ilu_cr(this->mesh->matrixSize, this->mesh->nzNum, &this->mesh->crsChem[chemDrug]->rowPtr[0], this->mesh->crsChem[chemDrug]->colInd, this->mesh->crsChem[chemDrug]->val, &this->mesh->uChem[chemDrug][0], this->mesh->bChem[chemDrug], itr_max, mr, tol_abs, tol_rel);

			time += this->mesh->deltaT;
		}

		if (this->config->output.chem && this->frame->time%this->config->output.filesFrequency == 0)
			FileFactory::makeFile(this->frame, this->config, this->mesh, this->chemotherapy, CHEM, chemDrug);
	}
}

void Macro2D::reaction() {
	if (config->output.prints)
		std::cout << "Reaction/Source term ..." << std::endl;

	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		for (int j = 0; j < this->mesh->unityRefined.y; ++j) {
			this->oProductionMesh[i][j] = 0.0;
			this->oUptakeMesh[i][j] = 0.0;
			this->egfSourceMesh[i][j] = 0.0;

			for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
				this->chemSupplyMesh[chemDrug][i][j] = 0.0;
				this->chemUptakeMesh[chemDrug][i][j] = 0.0;
			}
		}
	}

	//Powathil, 2012
	bool applyInjection = this->chemotherapy->applyInjection(this->frame->time);
	for (int l = 0; l < this->bloodVessel.size(); l++) {
		for (int i = (this->bloodVessel[l].y - this->config->bloodVessel.actionRadius)/this->mesh->hRefined; i < (this->bloodVessel[l].y + this->config->bloodVessel.actionRadius)/this->mesh->hRefined + 1; i++) {
			for (int j = (this->bloodVessel[l].x - this->config->bloodVessel.actionRadius)/this->mesh->hRefined; j < (this->bloodVessel[l].x + this->config->bloodVessel.actionRadius)/this->mesh->hRefined + 1; j++) {
				if (i >= this->mesh->unityRefined.y || j >= this->mesh->unityRefined.x || i < 0 || j < 0)
					continue;

				if (sqrt(pow(this->mesh->refined[j].x - this->bloodVessel[l].x, 2) + pow(this->mesh->refined[i].y - this->bloodVessel[l].y, 2)) <= this->config->bloodVessel.radius) {
					this->oProductionMesh[i][j] += this->config->bloodVessel.oProduction;

					for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
						if (applyInjection && this->config->continuum.chemotherapy.chemProtocol[this->chemotherapy->nextInjection - 1].compare(this->chemotherapy->chemDrugs[chemDrug]) == 0)
							this->chemSupplyMesh[chemDrug][i][j] += this->config->continuum.chemotherapy.chemDosages[this->chemotherapy->nextInjection - 1]*this->mesh->matrixSize/this->bloodVessel.size();
					}
				}
			}
		}
	}

	//Calcula termo de captação
	for (int l = 0; l < this->frame->cells.size(); l++) {
		if (this->frame->cells[l].state == NORMOXIC)
			continue;

		for (int i = (this->frame->cells[l].coordinates.y - this->frame->cells[l].actionRadius)/this->mesh->hRefined; i < (this->frame->cells[l].coordinates.y + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1; i++) {
			for (int j = (this->frame->cells[l].coordinates.x - this->frame->cells[l].actionRadius)/this->mesh->hRefined; j < (this->frame->cells[l].coordinates.x + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1; j++) {
				if (i >= this->mesh->unityRefined.y || j >= this->mesh->unityRefined.x || i < 0 || j < 0)
					continue;

				if (sqrt(pow(this->mesh->refined[j].x - this->frame->cells[l].coordinates.x, 2) + pow(this->mesh->refined[i].y - this->frame->cells[l].coordinates.y, 2)) <= this->frame->cells[l].radius) {
					this->oUptakeMesh[i][j] += this->frame->cells[l].oConsumption; // EXECUTAR OS DOIS
					this->egfSourceMesh[i][j] += this->frame->cells[l].egfSource;
					for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++)
						this->chemUptakeMesh[chemDrug][i][j] += this->frame->cells[l].chemUptake;
				}
			}
		}
	}

	//Calcula o termo reação na malha do nutriente
	double x, y;
	for(int i = 0; i < this->mesh->matrixSize; i++) {
		x = this->mesh->pos[i].x;
		y = this->mesh->pos[i].y;
		//std::cout << "Nó = " << i << " PONTO = (" << x << ", " << y << ") ";

		//Nós internos
		if (x != 0.0 && x != this->frame->domain.x && y != 0.0 && y != this->frame->domain.y) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, pow(this->mesh->hRatio+1,2));
			//std::cout << "Nós internos!" << std::endl;
			continue;
		}

		//Nós (---,0) -- Aresta 1
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (((this->mesh->hRatio/2)+1)*(this->mesh->hRatio+1)));
			//std::cout << "Aresta 1!" << std::endl;
			continue;
		}

		//Nós (0,---) -- Aresta 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, (((this->mesh->hRatio/2)+1)*(this->mesh->hRatio+1)));
			//std::cout << "Aresta 2!" << std::endl;
			continue;
		}

		//Nós (---,this->frame->domain.y) -- Aresta 3
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				i, (((this->mesh->hRatio/2)+1)*(this->mesh->hRatio+1)));
			//std::cout << "Aresta 3!" << std::endl;
			continue;
		}

		//Nós (this->frame->domain.x,---) -- Aresta 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio + (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, (((this->mesh->hRatio/2)+1)*(this->mesh->hRatio+1)));
			//std::cout << "Aresta 4!" << std::endl;
			continue;
		}

		//Nó (0,0) * Ponto 1
		if (x == 0.0 && y == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				0, (this->mesh->hRatio/2),
				i, pow(((this->mesh->hRatio/2))+1,2));
			//std::cout << "Ponto 1!" << std::endl;
			continue;
		}

		//Nó (this->frame->domain.x,0) * Ponto 2
		if (x == this->frame->domain.x && y == 0.0) {
			this->homogenize(0, (this->mesh->hRatio/2),
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,2));
			//std::cout << "Ponto 2!" << std::endl;
			continue;
		}

		//Nó (0,this->frame->domain.y) * Ponto 3
		if (x == 0.0 && y == this->frame->domain.y) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				0, (this->mesh->hRatio/2),
				i, pow(((this->mesh->hRatio/2))+1,2));
			//std::cout << "Ponto 3!" << std::endl;
			continue;
		}

		//Nó (this->frame->domain.x,this->frame->domain.y) * Ponto 4
		if (x == this->frame->domain.x && y == this->frame->domain.y) {
			this->homogenize((y/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (y/this->mesh->hCoarse)*this->mesh->hRatio,
				(x/this->mesh->hCoarse)*this->mesh->hRatio - (this->mesh->hRatio/2), (x/this->mesh->hCoarse)*this->mesh->hRatio,
				i, pow(((this->mesh->hRatio/2))+1,2));
			//std::cout << "Ponto 4!" << std::endl;
			continue;
		}
	}
}

void Macro2D::homogenize(int lBegin, int lEnd, int jBegin, int jEnd, int i, double homogenizationFactor) {
	double oProductionTemp = 0.0;
	double oUptakeTemp = 0.0;
	double egfSourceTemp = 0.0;
	std::vector<double> chemSupplyTemp = std::vector<double>(this->chemotherapy->chemDrugs.size(), 0.0);
	std::vector<double> chemUptakeTemp = std::vector<double>(this->chemotherapy->chemDrugs.size(), 0.0);

	for (int l = lBegin; l <= lEnd; l++) {
		for (int j = jBegin; j <= jEnd; j++) {
			oProductionTemp += this->oProductionMesh[l][j];

			if (this->oUptakeMesh[l][j] == 0)
				this->oUptakeMesh[l][j] = this->config->continuum.oxygen.oConsumptionBg;
			oUptakeTemp += this->oUptakeMesh[l][j];

			if (this->egfSourceMesh[l][j] == 0)
				this->egfSourceMesh[l][j] = this->config->continuum.egf.egfSourceBg;
			egfSourceTemp += this->egfSourceMesh[l][j];

			for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
				chemSupplyTemp[chemDrug] += this->chemSupplyMesh[chemDrug][l][j];

				if (this->chemUptakeMesh[chemDrug][l][j] == 0)
					this->chemUptakeMesh[chemDrug][l][j] = this->config->continuum.chemotherapy.chemUptakeBg;
				chemUptakeTemp[chemDrug] += this->chemUptakeMesh[chemDrug][l][j];
			}
		}
	}

	this->oProduction[i] = oProductionTemp/homogenizationFactor;
	this->oUptake[i] = oUptakeTemp/homogenizationFactor;
	this->egfSource[i] = egfSourceTemp/homogenizationFactor;

	for (int chemDrug = 0; chemDrug < this->chemotherapy->chemDrugs.size(); chemDrug++) {
		this->chemSupply[chemDrug][i] = chemSupplyTemp[chemDrug]/homogenizationFactor;
		this->chemUptake[chemDrug][i] = chemUptakeTemp[chemDrug]/homogenizationFactor;
	}
}
