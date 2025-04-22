#include "Mesh.hpp"

using namespace std;

void Mesh::printAttributes() {
	cout << "hCoarse = " << this->hCoarse << endl;
	cout << "hRefined = " << this->hRefined << endl;
	cout << "hRatio = " << this->hRatio << endl;
	cout << "deltaT = " << this->deltaT << endl;
	cout << "matrixSize = " << this->matrixSize << endl;
	cout << "nzNum = " << this->nzNum << endl;
	cout << "unityCoarse = " << this->unityCoarse.to_string() << endl;
	cout << "unityRefined = " << this->unityRefined.to_string() << endl;

	cout << "refined = ";
	for (vector<Vector3>::iterator it = this->refined.begin(); it != this->refined.end(); it++)
		cout << (*it).to_string() << "\t";
	cout << endl;

	cout << "pos = ";
	for (vector<Vector3>::iterator it = this->pos.begin(); it != this->pos.end(); it++)
		cout << (*it).to_string() << "\t";
	cout << endl;

	cout << "uO = ";
	for (vector<double>::iterator it = this->uO.begin(); it != this->uO.end(); it++)
		cout << *it << "\t";
	cout << endl;

	cout << "uEGF = ";
	for (vector<double>::iterator it = this->uEGF.begin(); it != this->uEGF.end(); it++)
		cout << *it << "\t";
	cout << endl;

	for (int i = 0; i < this->uChem.size(); i++) {
		cout << "uChem[" << i << "] = ";
		for (vector<double>::iterator it = this->uChem[i].begin(); it != this->uChem[i].end(); it++)
			cout << *it << "\t";
		cout << endl;
	}

	cout << "sigmaO = " << this->sigmaO << endl;
	cout << "sigmaEGF = " << this->sigmaEGF << endl;
	for (int i = 0; i < this->sigmaChem.size(); i++)
		cout << "sigmaChem[" << i << "] = " << this->sigmaChem[i] << endl;
	for (int i = 0; i < this->chemDecay.size(); i++)
		cout << "chemDecay[" << i << "] = " << this->chemDecay[i] << endl;

	cout << "bO = ";
	for (int i = 0; i < this->matrixSize; i++)
		cout << this->bO[i] << "\t";
	cout << endl;

	cout << "bEGF = ";
	for (int i = 0; i < this->matrixSize; i++)
		cout << this->bEGF[i] << "\t";
	cout << endl;

	for (int i = 0; i < this->bChem.size(); i++) {
		cout << "bChem[" << i << "] = ";
		for (int j = 0; j < this->matrixSize; j++)
			cout << this->bChem[i][j] << "\t";
		cout << endl;
	}

	cout << "crsO = " << endl;
	this->crsO->printAsMatrix();

	cout << "crsEGF = " << endl;
	this->crsEGF->printAsMatrix();

	for (int i = 0; i < this->crsChem.size(); i++) {
		cout << "crsChem[" << i << "] = " << endl;
		this->crsChem[i]->printAsMatrix();
	}
}

double Mesh::calculateMax_uO() {
	double max_uO = this->uO[0];
	for (int i = 1; i < this->matrixSize; i++) {
		if (this->uO[i] > max_uO)
			max_uO = this->uO[i];
	}

	return (max_uO == 0.0) ? 1.0 : max_uO;
}

double Mesh::calculateMax_uChem(int chemDrug) {
	double max_uChem = this->uChem[chemDrug][0];
	for (int i = 1; i < this->matrixSize; i++) {
		if (this->uChem[chemDrug][i] > max_uChem)
			max_uChem = this->uChem[chemDrug][i];
	}

	return (max_uChem < 1.0e-10) ? 1.0 : max_uChem;
}

// GUSTAVO: Conferir depois
Mesh3D::Mesh3D(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy) {
	this->hCoarse = config->continuum.hCoarse;
	this->hRefined = config->continuum.hRefined;

	this->unityCoarse = Vector3i((domain/this->hCoarse)+1);
	this->unityRefined = Vector3i((domain/this->hRefined)+1);
	this->matrixSize = this->unityCoarse.z * this->unityCoarse.y * this->unityCoarse.z;
	this->uO = std::vector<double>(this->matrixSize, 0.0);
	this->uEGF = std::vector<double>(this->matrixSize, 0.0);
	this->nzNum = 8*4 + 4*(this->unityCoarse.x-2)*5 + 4*(this->unityCoarse.y-2)*5 + 4*(this->unityCoarse.z-2)*5 + 2*(this->unityCoarse.y-2)*(this->unityCoarse.z-2)*6 + 2*(this->unityCoarse.x-2)*(this->unityCoarse.z-2)*6 + 2*(this->unityCoarse.x-2)*(this->unityCoarse.y-2)*6 + (this->unityCoarse.x-2)*(this->unityCoarse.y-2)*(this->unityCoarse.z-2)*7;
	this->hRatio = this->hCoarse/this->hRefined;
	this->deltaT = config->continuum.deltaT;

	this->sigmaO = config->continuum.oxygen.oDiffusion/pow(this->hCoarse, 2);
	this->sigmaEGF = this->deltaT*config->continuum.egf.egfDiffusion/pow(this->hCoarse, 2);

	this->bO = new double[this->matrixSize];
	this->bEGF = new double[this->matrixSize];
	
	this->pos = std::vector<Vector3>(this->matrixSize, Vector3());

	for (size_t i = 0; i < std::max(this->unityRefined.x, std::max(this->unityRefined.y, this->unityRefined.z)); i++) {
		this->refined.push_back(Vector3(i*this->hRefined, i*this->hRefined, i*this->hRefined));
	}

	this->createCRS(domain, chemotherapy);
	this->setBoundaryConditions(domain, config, chemotherapy);
}

Mesh3D::~Mesh3D() {
	delete this->bO;
	delete this->bEGF;
	
	delete this->crsO;
	delete this->crsEGF;
}

void Mesh3D::createCRS(Vector3 domain, Chemotherapy* chemotherapy) {
	/*this->rowPtr = std::vector<int>((this->matrixSize)+1, 0.0);
	this->colInd = new int[this->nzNum];
	this->val = new double[this->nzNum];
	this->rowPtrEGF = new int[(this->matrixSize)+1];
	this->colIndEGF = new int[this->nzNum];
	this->valEGF = new double[this->nzNum];*/
	this->crsO = new CRS(this->matrixSize, this->nzNum);
	this->crsEGF = new CRS(this->matrixSize, this->nzNum);
	for (int i = 0; i < this->matrixSize; i++){
		this->pos[i].z = (i/(this->unityCoarse.z*this->unityCoarse.y))*this->hCoarse;
		this->pos[i].y = ((i%(this->unityCoarse.z*this->unityCoarse.y))/this->unityCoarse.z)*this->hCoarse;
		this->pos[i].x = ((i%(this->unityCoarse.z*this->unityCoarse.y))%this->unityCoarse.z)*this->hCoarse;
		for (int j = 0; j < this->matrixSize; j++){
			if (i==j){
				this->crsO->createOrdered(i, j, 0);
				this->crsEGF->createOrdered(i, j, 1 + 6 * this->sigmaEGF);
				continue;
			}
			if ( (this->pos[i].x!=0 && i-j==1) || (this->pos[i].x!=domain.x && j-i==1) ){
				this->crsO->createOrdered(i, j, -this->sigmaO);
				this->crsEGF->createOrdered(i, j, -this->sigmaEGF);
				continue;
			}
			if ( (this->pos[i].y!=0 && i-j==this->unityCoarse.z) || (this->pos[i].y!=domain.y && j-i==this->unityCoarse.z) ){
				this->crsO->createOrdered(i,j, -this->sigmaO);
				this->crsEGF->createOrdered(i,j, -this->sigmaEGF);
				continue;
			}
			if ( (i-j==(this->unityCoarse.z*this->unityCoarse.y)) || (j-i==(this->unityCoarse.z*this->unityCoarse.y)) ){
				this->crsO->createOrdered(i,j, -this->sigmaO);
				this->crsEGF->createOrdered(i,j, -this->sigmaEGF);
			}
		}
	}
}

void Mesh3D::setBoundaryConditions(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy) {
	//Condição de Contorno oxigênio
	for(int i=0;i<this->matrixSize;i++){
		this->bO[i] = 0;
		this->bEGF[i] = 0;
		//Dirichlet
		if (this->pos[i].x == 0 || this->pos[i].x == domain.x || this->pos[i].y == 0 || this->pos[i].y==domain.y || this->pos[i].z == 0 || this->pos[i].z==domain.z){
			this->bO[i]= config->continuum.oxygen.oConsumptionBorder;
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++){
				if (this->crsO->colInd[j] == i){
					this->crsO->val[j]=1.0;
					this->crsEGF->val[j]=1.0;
				}
				else {
					this->crsO->val[j]=0.0;
					this->crsEGF->val[j]=0.0;
				}
			}
		}
	}
}

Mesh2D::Mesh2D(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy) {
	this->hCoarse = config->continuum.hCoarse;
	this->hRefined = config->continuum.hRefined;
	this->hRatio = this->hCoarse/this->hRefined;
	this->deltaT = config->continuum.deltaT;
	this->unityCoarse = Vector3i(domain.x/this->hCoarse + 1, domain.y/this->hCoarse + 1, 0);
	this->unityRefined = Vector3i(domain.x/this->hRefined + 1, domain.y/this->hRefined + 1, 0);
	this->matrixSize = this->unityCoarse.x * this->unityCoarse.y;
	this->nzNum = 4*3 + 2*(this->unityCoarse.x-2)*4 + 2*(this->unityCoarse.y-2)*4 + (this->unityCoarse.x-2)*(this->unityCoarse.y-2)*5; // GUSTAVO: ???

	for (size_t i = 0; i < std::max(this->unityRefined.x, this->unityRefined.y); i++) {
		this->refined.push_back(Vector3(i*this->hRefined, i*this->hRefined, 0.0));
	}

	this->pos = std::vector<Vector3>(this->matrixSize, Vector3());

	this->uO = std::vector<double>(this->matrixSize, 0.0);
	this->uEGF = std::vector<double>(this->matrixSize, 0.0);

	this->uChem.resize(chemotherapy->chemDrugs.size());
	for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++) {
		this->uChem[chemDrug].resize(this->matrixSize);
		for (int i = 0; i < this->matrixSize; i++)
			this->uChem[chemDrug][i] = 0.0;
	}

	// this->sigmaO = config->continuum.oxygen.oDiffusion/pow(this->hCoarse, 2); // Rocha, 2018
	this->sigmaO = this->deltaT*config->continuum.oxygen.oDiffusion/pow(this->hCoarse, 2); // Powathil, 2012
	this->sigmaEGF = this->deltaT*config->continuum.egf.egfDiffusion/pow(this->hCoarse, 2);

	this->sigmaChem = std::vector<double>(chemotherapy->chemDrugs.size(), 0.0);
	this->chemDecay = std::vector<double>(chemotherapy->chemDrugs.size(), 0.0);
	for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++) {
		if (chemotherapy->chemDrugs[chemDrug].compare("A") == 0) {
			this->sigmaChem[chemDrug] = this->deltaT*config->continuum.chemotherapy.CYC202.chemDiffusion/pow(this->hCoarse, 2);
			this->chemDecay[chemDrug] = config->continuum.chemotherapy.CYC202.chemDecay;
		} else if (chemotherapy->chemDrugs[chemDrug].compare("B") == 0) {
			this->sigmaChem[chemDrug] = this->deltaT*config->continuum.chemotherapy.cisplatin.chemDiffusion/pow(this->hCoarse, 2);
			this->chemDecay[chemDrug] = config->continuum.chemotherapy.cisplatin.chemDecay;
		} else if (chemotherapy->chemDrugs[chemDrug].compare("C") == 0) {
			this->sigmaChem[chemDrug] = this->deltaT*config->continuum.chemotherapy.taxotere.chemDiffusion/pow(this->hCoarse, 2);
			this->chemDecay[chemDrug] = config->continuum.chemotherapy.taxotere.chemDecay;
		} else if (chemotherapy->chemDrugs[chemDrug].compare("D") == 0) {
			this->sigmaChem[chemDrug] = this->deltaT*config->continuum.chemotherapy.taxol.chemDiffusion/pow(this->hCoarse, 2);
			this->chemDecay[chemDrug] = config->continuum.chemotherapy.taxol.chemDecay;
		}
	}

	this->bO = new double[this->matrixSize];
	this->bEGF = new double[this->matrixSize];

	this->bChem.resize(chemotherapy->chemDrugs.size());
	for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
		this->bChem[chemDrug] = new double[this->matrixSize];

	this->createCRS(domain, chemotherapy);
	this->setBoundaryConditions(domain, config, chemotherapy);
}

Mesh2D::~Mesh2D() {
	delete this->bO;
	delete this->bEGF;
	for (int i = 0; i < this->bChem.size(); i++)
		delete this->bChem[i];

	delete this->crsO;
	delete this->crsEGF;
	for (int i = 0; i < this->crsChem.size(); i++)
		delete this->crsChem[i];
}

void Mesh2D::createCRS(Vector3 domain, Chemotherapy* chemotherapy) {
	/*this->rowPtr = std::vector<int>(this->matrixSize + 1, 0.0);
	this->colInd = new int[this->nzNum];
	this->val = new double[this->nzNum];
	this->rowPtrEGF = new int[this->matrixSize + 1];
	this->colIndEGF = new int[this->nzNum];
	this->valEGF = new double[this->nzNum];*/
	this->crsO = new CRS(this->matrixSize, this->nzNum);
	this->crsEGF = new CRS(this->matrixSize, this->nzNum);

	this->crsChem.resize(chemotherapy->chemDrugs.size());
	for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
		this->crsChem[chemDrug] = new CRS(this->matrixSize, this->nzNum);

	for (int i = 0; i < this->matrixSize; i++) {
		this->pos[i].y = (i/this->unityCoarse.x)*this->hCoarse;
		this->pos[i].x = (i%this->unityCoarse.x)*this->hCoarse;
		for (int j = 0; j < this->matrixSize; j++) {
			if (i == j) {
				this->crsO->createOrdered(i, j, 0.0);
				this->crsEGF->createOrdered(i, j, 1.0 + 4.0*this->sigmaEGF);
				for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
					this->crsChem[chemDrug]->createOrdered(i, j, 0.0);
			} else if ((this->pos[i].x != 0 && i - j == 1) || (this->pos[i].x != domain.x && j - i == 1)) {
				this->crsO->createOrdered(i, j, -this->sigmaO);
				this->crsEGF->createOrdered(i, j, -this->sigmaEGF);
				for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
					this->crsChem[chemDrug]->createOrdered(i, j, -this->sigmaChem[chemDrug]);
			} else if ((this->pos[i].y != 0 && i - j == this->unityCoarse.x) || (this->pos[i].y != domain.y && j - i == this->unityCoarse.x)) {
				this->crsO->createOrdered(i, j, -this->sigmaO);
				this->crsEGF->createOrdered(i, j, -this->sigmaEGF);
				for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
					this->crsChem[chemDrug]->createOrdered(i, j, -this->sigmaChem[chemDrug]);
			}
		}
	}
}

void Mesh2D::setBoundaryConditions(Vector3 domain, ConfigHandler* config, Chemotherapy* chemotherapy) {
	//Condição de Contorno oxigênio
	for(int i = 0; i < this->matrixSize; i++) {
		this->bO[i] = 0.0;
		this->bEGF[i] = 0.0;
		for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
			this->bChem[chemDrug][i] = 0.0;

		//Dirichlet // Rocha, 2018
		if (this->pos[i].x == 0 || this->pos[i].x == domain.x || this->pos[i].y == 0 || this->pos[i].y == domain.y) {
			// this->bO[i] = config->continuum.oxygen.oConsumptionBorder;
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++) {
				if (this->crsO->colInd[j] == i) {
					// this->crsO->val[j] = 1.0;
					this->crsEGF->val[j] = 1.0;
				} else {
					// this->crsO->val[j] = 0.0;
					this->crsEGF->val[j] = 0.0;
				}
			}
		}

		//Neumann // Powathil, 2012
		if (this->pos[i].x == 0) {
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++) {
				if (i - this->crsO->colInd[j] == 1) {
					this->crsO->val[j] = 0.0;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = 0.0;
				} else if (this->crsO->colInd[j] - i == 1) {
					this->crsO->val[j] = -2.0*sigmaO;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = -2.0*sigmaChem[chemDrug];
				}
			}
		} else if (this->pos[i].x == domain.x) {
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++) {
				if (i - this->crsO->colInd[j] == 1) {
					this->crsO->val[j] = -2.0*sigmaO;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = -2.0*sigmaChem[chemDrug];
				} else if (this->crsO->colInd[j] - i == 1) {
					this->crsO->val[j] = 0.0;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = 0.0;
				}
			}
		}

		if (this->pos[i].y == 0) {
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++) {
				if (this->crsO->colInd[j] - i == this->unityCoarse.x) {
					this->crsO->val[j] = -2.0*sigmaO;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = -2.0*sigmaChem[chemDrug];
				}
			}
		} else if (this->pos[i].y == domain.y) {
			for(int j = this->crsO->rowPtr[i]; j < this->crsO->rowPtr[i+1]; j++) {
				if (i - this->crsO->colInd[j] == this->unityCoarse.x) {
					this->crsO->val[j] = -2.0*sigmaO;
					for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++)
						this->crsChem[chemDrug]->val[j] = -2.0*sigmaChem[chemDrug];
				}
			}
		}
	}
}
