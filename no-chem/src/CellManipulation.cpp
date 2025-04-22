#include "CellManipulation.hpp"

void CellManipulation::changeCellState(Cell* cell, int previousState, int state, int lifetime, double oConsumption, double egfSource, double chemUptake) {
	cell->previousState = previousState;
	cell->state = state;
	cell->lifetime = lifetime;
	cell->oConsumption = oConsumption;
	cell->egfSource = egfSource;
	cell->chemUptake = chemUptake;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Norma euclidiana
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline double CellManipulation::norma(Vector3 pos) {
	return sqrt(pow(pos.x, 2) + pow(pos.y, 2) + pow(pos.z, 2));
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Função phi
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation::func_var_phi(Vector3 normal, double actionRadius, int n) { 
	double r = norma(normal);
	if (r > 0.0 && r <= actionRadius) {
		double var = pow(1.0 - r/actionRadius, n + 1)/r;
		// return Vector3(var*normal.x, var*normal.y, var*normal.z);
		return normal*var;
	}

	return Vector3();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Função psi
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation::func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m) {
	double r = norma(normal);
	if (r > 0.0 && r < nucleusRadius) {
		double c = pow(1.0 - nucleusRadius/radius, m + 1) - M;
		double var = -(c*r/nucleusRadius + M)/r;
		// return Vector3(var*normal.x, var*normal.y, var*normal.z);
		return normal*var;
	}

	if (r >= nucleusRadius && r <= radius) {
		double var = -pow(1.0 - r/radius, m + 1)/r;
		// return Vector3(var*normal.x, var*normal.y, var*normal.z);
		return normal*var;
	}

	// if ( (r > radius) || (r == 0) )

	return Vector3();
}

void CellManipulation::updateInitialFrame(Frame* frame, ConfigHandler* config, Ran* ran) {
	if (config->output.prints)
		std::cout << "Updating frame..." << std::endl;

	// #pragma omp parallel for schedule(static) ordered
	for (int i = 0; i < frame->cells.size(); i++) {
		if (frame->cells[i].state > LAST_CELL_STATE)
			continue;

		frame->cells[i].nucleusRadius = config->agent.nucleusRadius;
		frame->cells[i].radius = config->agent.radius;
		frame->cells[i].actionRadius = config->agent.actionRadius;
		switch (frame->cells[i].state) {
			case QUIESCENT:
			frame->cells[i].oConsumption = config->agent.oConsumption;
			frame->cells[i].egfSource = config->agent.egfSource;
			break;

			case PROLIFERATIVE2:
			frame->cells[i].oConsumption = config->agent.oConsumption;
			frame->cells[i].egfSource = 0.0;
			break;

			case NORMOXIC:
			frame->cells[i].oConsumption = 0.0;
			frame->cells[i].egfSource = 0.0;
			break;
		}
	}
}

// Forças do movimento haptotático
double* CellManipulation::generateF_hap(Frame* frame, Mesh* mesh, Ran* ran) {
	int n = frame->domain.x/mesh->hCoarse + 1;
	int m = frame->domain.y/mesh->hCoarse + 1;
	
	double* F_hap = new double[n*m];
	for (int i = 0; i < n*m; i++)
		F_hap[i] = 0.5 + (1.0 - 0.5)*ran->doub();

	return F_hap;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//					Cálculo do gradiente do movimento haptotático
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CellManipulation::calculateGrad_F_hap(Frame* frame, Mesh* mesh, double* F_hap) {
	int n = frame->domain.x/mesh->hCoarse + 1;
	int m = frame->domain.y/mesh->hCoarse + 1;

	for (int i = 0; i < frame->cells.size(); i++) {
		if (frame->cells[i].state != MIGRATE)
			continue;

		// Calcula gradiente da força em cada ponto
		int xA = frame->cells[i].coordinates.x/mesh->hCoarse, xD = xA;
		int xB = xA + 1, xC = xB;
		int yA = frame->cells[i].coordinates.y/mesh->hCoarse, yB = yA;
		int yC = yA + 1, yD = yC;

		double F_hap_loc[4] = {F_hap[yA*n + xA], F_hap[yB*n + xB], F_hap[yC*n + xC], F_hap[yD*n + xD]};

		// Calcula concentração de oxigênio na célula
		double xCell = frame->cells[i].coordinates.x - xA*mesh->hCoarse;
		double yCell = frame->cells[i].coordinates.y - yA*mesh->hCoarse;

		if (xCell >= yCell) {
			frame->cells[i].grad_F_hap[0] = 1/pow(mesh->hCoarse, 2)*(-mesh->hCoarse)*F_hap_loc[0] + 1/pow(mesh->hCoarse, 2)*mesh->hCoarse*F_hap_loc[1];
			frame->cells[i].grad_F_hap[1] = 1/pow(mesh->hCoarse, 2)*(-mesh->hCoarse)*F_hap_loc[1] + 1/pow(mesh->hCoarse, 2)*mesh->hCoarse*F_hap_loc[2];
		} else {
			frame->cells[i].grad_F_hap[0] = 1/pow(mesh->hCoarse, 2)*mesh->hCoarse*F_hap_loc[2] + 1/pow(mesh->hCoarse, 2)*(-mesh->hCoarse)*F_hap_loc[3];
			frame->cells[i].grad_F_hap[1] = 1/pow(mesh->hCoarse, 2)*(-mesh->hCoarse)*F_hap_loc[0] + 1/pow(mesh->hCoarse, 2)*mesh->hCoarse*F_hap_loc[3];
		}
	}
}

// GUSTAVO: Conferir depois
void CellManipulation3D::calculateCellSigmas(Cell* cell, Mesh* mesh) {
	Vector3 pos((int)(cell->coordinates.x/mesh->hCoarse),
		(int)(cell->coordinates.y/mesh->hCoarse),
		(int)(cell->coordinates.z/mesh->hCoarse));
	Vector3 c   (((cell->coordinates.x - ((int)cell->coordinates.x))/mesh->hCoarse),
		((cell->coordinates.y - ((int)cell->coordinates.y))/mesh->hCoarse),
		((cell->coordinates.z - ((int)cell->coordinates.z))/mesh->hCoarse));

	double oSigma1 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + (int)pos.x];
	double oSigma2 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
	double oSigma3 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
	double oSigma4 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];
	double oSigma5 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + (int)pos.x];
	double oSigma6 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
	double oSigma7 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
	double oSigma8 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];

	double egfSigma1 = mesh->uEGF[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + (int)pos.x];
	double egfSigma2 = mesh->uEGF[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
	double egfSigma3 = mesh->uEGF[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
	double egfSigma4 = mesh->uEGF[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];
	double egfSigma5 = mesh->uEGF[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + (int)pos.x];
	double egfSigma6 = mesh->uEGF[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		(int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
	double egfSigma7 = mesh->uEGF[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
	double egfSigma8 = mesh->uEGF[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
		((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];
	cell->sigmaO = (1-c.x)*(1-c.y)*(1-c.z)*oSigma1 + c.x*(1-c.y)*(1-c.z)*oSigma2 + (1-c.x)*c.y*(1-c.z)*oSigma3 + c.x*c.y*(1-c.z)*oSigma4 + (1-c.x)*(1-c.y)*c.z*oSigma5 + c.x*(1-c.y)*c.z*oSigma6 + (1-c.x)*c.y*c.z*oSigma7 + c.x*c.y*c.z*oSigma8;
	cell->sigmaEGF = (1-c.x)*(1-c.y)*(1-c.z)*egfSigma1 + c.x*(1-c.y)*(1-c.z)*egfSigma2 + (1-c.x)*c.y*(1-c.z)*egfSigma3 + c.x*c.y*(1-c.z)*egfSigma4 + (1-c.x)*(1-c.y)*c.z*egfSigma5 + c.x*(1-c.y)*c.z*egfSigma6 + (1-c.x)*c.y*c.z*egfSigma7 + c.x*c.y*c.z*egfSigma8;
}

Cell CellManipulation3D::divide(Cell* cell, double rand1, double rand2) {
	double  theta = 2*M_PI*rand1,
	phi = M_PI*rand2,
	nucleusRadius = cell->nucleusRadius/pow(2,(1/3.0));
	Vector3 pos = cell->coordinates, dist(nucleusRadius*sin(phi)*cos(theta), nucleusRadius*sin(phi)*sin(theta), nucleusRadius*cos(phi));

	cell->coordinates.x = pos.x + dist.x;
	cell->coordinates.y = pos.y + dist.y;
	cell->coordinates.z = pos.z + dist.z;
	cell->radius = cell->radius/pow(2,(1/3.0));
	cell->nucleusRadius = nucleusRadius;
	cell->actionRadius = cell->actionRadius/pow(2,(1/3.0));
	cell->state = PROLIFERATIVE1;

	Cell child = *cell;
	child.coordinates.x = pos.x - dist.x;
	child.coordinates.y = pos.y - dist.y;
	child.coordinates.z = pos.z - dist.z;
	return child;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Calcula vetor normal
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation3D::normal(Vector3 coordinates, double domainRadius) {
	//Direção
	Vector3 u(coordinates.x - domainRadius, coordinates.y - domainRadius, coordinates.z - domainRadius);
	// double u1 = , u2 =  u3 =;
	double D = norma(u);
	//Vesor
	u /= D; // u1 = u1/D; u2 = u2/D; u3 = u3/D;
	//Ponto na casca esférica
	u *= domainRadius; // double X = u1*domainRadius, Y = u2*domainRadius, Z = u3*domainRadius;

	return coordinates - u; //Vector3(x_cel - u.x, y_cel - u.y, z_cel - u.z);
}

void CellManipulation3D::updateFrame(Frame* frame, ConfigHandler* config, Mesh* mesh, EGFRPathway* egfrPathway, CellCyclePathway* cellCyclePathway, Ran* ran, Chemotherapy* chemotherapy, PharmacologicModel* pharmacologic) {
	if (config->output.prints)
		std::cout << "Updating frame..." << std::endl;

	std::vector<int> cellsRemoved;
	std::vector<Cell> cellsInserted;
	std::vector<double> randomNumbers = std::vector<double>(frame->cells.size()+1, 0.0);
	// double randomNumbers[frame->cells.size()+1];
	for(int i = 0; i < frame->cells.size()+1; i++)
		randomNumbers[i] = ran->doub();
	// #pragma omp parallel for schedule(static) ordered
	for(int i = 0; i < frame->cells.size(); i++){
		calculateCellSigmas(&frame->cells[i], mesh);

		//Célula Normal
		if (frame->cells[i].state == NORMOXIC || frame->cells[i].state > LAST_CELL_STATE) continue;

		//EGFRPathway
		config->parameters.alphaP = 1;
		/*double ROC_PLC,ROC_ERK;
		if (frame->cells[i].state == QUIESCENT){
			egfrPathway->runge_kutta(frame->cells[i].sigmaEGF, &ROC_PLC, &ROC_ERK);
			if ((ROC_PLC < egfrPathway->ROC_PLC_th) && (ROC_ERK > egfrPathway->ROC_ERK_th)) config->parameters.alphaP = 1;
			else config->parameters.alphaP = 0;
		}*/

		if (frame->cells[i].state == QUIESCENT){
			//Transição Quiescente para Proliferativa
			double alpha_p = config->parameters.alphaP*(frame->cells[i].sigmaO-config->parameters.sigmaH)/(1 - config->parameters.sigmaH)*(1-frame->outCells/config->agent.maxOutCells);
			if (1-exp(-config->parameters.delta_tt*alpha_p) >= randomNumbers[i]){
				frame->cells[i].lifetime = frame->time;
				frame->cells[i].previousState = frame->cells[i].state;
				frame->cells[i].state = PROLIFERATIVE2;
				frame->cells[i].egfSource = 0.0;
			}
				//Transição Quiescente para Apoptótica
			else if (1-exp(-config->parameters.delta_tt*config->parameters.alphaA) >= randomNumbers[i]){
				frame->cells[i].lifetime = frame->time;
				frame->cells[i].previousState = frame->cells[i].state;
				frame->cells[i].state = APOPTOTIC;
				frame->cells[i].egfSource = 0.0;
			}
			continue;
		}

		//Transição para Hipóxica
		if (frame->cells[i].state != HYPOXIC && frame->cells[i].state != NECROTIC && frame->cells[i].sigmaO < config->parameters.sigmaH){
			frame->cells[i].lifetime = frame->time;
			frame->cells[i].previousState = frame->cells[i].state;
			frame->cells[i].state = HYPOXIC;
			frame->cells[i].egfSource = 0.0;
		}

		//Transição Hipóxica para Necrótica
		if (frame->cells[i].state == HYPOXIC){
			frame->cells[i].lifetime = frame->time;
			frame->cells[i].previousState = frame->cells[i].state;
			frame->cells[i].state = NECROTIC;
		//continue;
		}

		double tau = frame->time - frame->cells[i].lifetime;

		//Ocorre a Apoptose
		if (frame->cells[i].state == APOPTOTIC){
			if(tau >= config->parameters.tauA){
				// #pragma omp ordered
				cellsRemoved.push_back (i);
				//frame->cells.erase (frame->cells.begin() + i);
				frame->tumorCells -= 1;
				continue;
			}
			else{
				if (tau == 1)  frame->cells[i].calcification = frame->cells[i].radius/(config->parameters.tauA*0.5);
				if (tau > config->parameters.tauA*0.5){
					frame->cells[i].radius = frame->cells[i].radius - frame->cells[i].calcification;
					frame->cells[i].nucleusRadius = frame->cells[i].nucleusRadius - frame->cells[i].calcification/2;
				}
				continue;
			}
		}

		//Ocorre a Necrose
		if (frame->cells[i].state == NECROTIC){
			if(tau == 0) frame->cells[i].calcification = frame->cells[i].radius;
			if(tau <= config->parameters.tauNL){
				frame->cells[i].radius = frame->cells[i].calcification * pow(1+(config->parameters.fNS*tau)/config->parameters.tauNL,1/3.0);
				//frame->cells[i].actionRadius = 1.214*frame->cells[i].radius;
				continue;
			}
			if(tau > config->parameters.tauNL && tau <= config->parameters.tauC){
				if ((tau - config->parameters.tauNL)< 8) frame->cells[i].radius = frame->cells[i].radius - 0.1*frame->cells[i].radius;
				else frame->cells[i].radius = frame->cells[i].nucleusRadius;
				//frame->cells[i].actionRadius = 1.214*frame->cells[i].radius;
				frame->cells[i].calcification = tau/config->parameters.tauC;
				continue;
			}
			if (tau >= config->parameters.tauC){
				frame->cells[i].calcification = 1.0;
				continue;
			}
		}

		//Ocorre a divisão celular
		if (frame->cells[i].state == PROLIFERATIVE2 && tau == 9){
			//frame->cells.push_back(this->divide(&frame->cells[i], config->input.seed));
			// #pragma omp ordered
			cellsInserted.push_back(divide(&frame->cells[i], randomNumbers[i],randomNumbers[i+1]));
			frame->tumorCells++;
			continue;
		}

		//Ocorre a G1
		if ( frame->cells[i].state == PROLIFERATIVE1 && tau > (config->parameters.tauP-config->parameters.tauG1) ){
			if(tau < config->parameters.tauP){
				frame->cells[i].nucleusRadius = config->agent.nucleusRadius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
				frame->cells[i].radius        = config->agent.radius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
				frame->cells[i].actionRadius  = frame->cells[i].radius*1.214;
			}
			else{
				frame->cells[i].previousState = PROLIFERATIVE2;
				frame->cells[i].nucleusRadius = config->agent.nucleusRadius;
				frame->cells[i].radius        = config->agent.radius;
				frame->cells[i].actionRadius  = config->agent.actionRadius;
				frame->cells[i].state         = QUIESCENT;
				frame->cells[i].lifetime      = frame->time;
				frame->cells[i].egfSource     = config->agent.egfSource;
			}
		}
	}
	std::sort(cellsRemoved.begin(), cellsRemoved.end());
	//cellsRemoved Células
	for(int i = cellsRemoved.size()-1; i >=0; i--){
		frame->cells.erase (frame->cells.begin() + cellsRemoved[i]);
	}
	//Insere Células
	for(int i = 0; i < cellsInserted.size(); i++){
		frame->cells.push_back(cellsInserted[i]);
	}
}

void CellManipulation3D::force(Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel) {
	if (config->output.prints)
		std::cout << "Forces...";

	double iter = 0, domainRadius = frame->domain.y/2;
	int count=0;
	while (iter < 60){
		int n =1; int m = 1; double M = 1;
		for(int i = 0; i < frame->cells.size(); i++){
			for(int j = i+1; j < frame->cells.size(); j++){
				Vector3 r(
					frame->cells[j].coordinates.x - frame->cells[i].coordinates.x,
					frame->cells[j].coordinates.y - frame->cells[i].coordinates.y,
					frame->cells[j].coordinates.z - frame->cells[i].coordinates.z
					);
				double actionRadius = frame->cells[j].actionRadius + frame->cells[i].actionRadius;

				if (norma(r) > actionRadius) continue;

				Vector3 phi = func_var_phi(r, actionRadius, n),
				psi = func_var_psi( r,
					frame->cells[j].nucleusRadius + frame->cells[i].nucleusRadius,
					frame->cells[j].radius + frame->cells[i].radius,
					M, m);

				Vector3 f((phi*config->forces.c_cca) + (psi*config->forces.c_ccr));

				frame->cells[i].force += f;
				frame->cells[j].force -= f;
			}

			if ( norma(frame->cells[i].coordinates - domainRadius) < (domainRadius - frame->cells[i].actionRadius) ) continue;
			Vector3 _normal = normal(frame->cells[i].coordinates, domainRadius);
			Vector3 phi = func_var_phi(_normal, frame->cells[i].actionRadius, n),
			psi = func_var_psi(_normal, frame->cells[i].nucleusRadius, frame->cells[i].radius, M, m);

			Vector3 f((phi*config->forces.c_ct) + (psi*config->forces.c_rct));

			frame->cells[i].force += f;
		}

		//Calcula velocidade
		// #pragma omp parallel for
		for(int i = 0; i < frame->cells.size(); i++){
			frame->cells[i].speed.x = -0.5*(frame->cells[i].force.x);
			frame->cells[i].speed.y = -0.5*(frame->cells[i].force.y);
			frame->cells[i].speed.z = -0.5*(frame->cells[i].force.z);
		}

		//Norma da velocidade maior
		double v_max= 10e-10;
		// #pragma omp parallel for reduction(max : v_max)
		for(int i = 0; i < frame->cells.size(); i++){
			if (norma(frame->cells[i].speed) >= v_max) v_max = norma(frame->cells[i].speed);
		}
		double delta_tt = (1/v_max);

		// #pragma omp parallel for
		for(int i = 0; i < frame->cells.size(); i++){
			//Desloca Células
			frame->cells[i].coordinates.x += delta_tt * frame->cells[i].speed.x;
			frame->cells[i].coordinates.y += delta_tt * frame->cells[i].speed.y;
			frame->cells[i].coordinates.z += delta_tt * frame->cells[i].speed.z;

			//ZeactionRadius somatório de forças
			frame->cells[i].force = Vector3();

		}

		for(int i = 0; i < frame->cells.size(); i++){
			double dist = sqrt( pow(frame->cells[i].coordinates.x-domainRadius,2)+pow(frame->cells[i].coordinates.y-domainRadius,2)+pow(frame->cells[i].coordinates.z-domainRadius,2) );
			//Se passar do domínio exclui
			if (dist >= domainRadius){
				if (frame->cells[i].state != NORMOXIC) frame->tumorCells -= 1;
				frame->cells.erase (frame->cells.begin() + i);
				i--;
				frame->outCells += 1;
			}
		}

		iter+= delta_tt;
		count++;
	}

	if (config->output.prints)
		std::cout << " (" << count << " iterations)" << std::endl;
}

void CellManipulation2D::calculateCellSigmas(Cell* cell, Mesh* mesh) {
	Vector3 pos((int) (cell->coordinates.x/mesh->hCoarse), (int) (cell->coordinates.y/mesh->hCoarse), 0.0);
	Vector3 c((cell->coordinates.x - ((int) cell->coordinates.x))/mesh->hCoarse, (cell->coordinates.y - ((int) cell->coordinates.y))/mesh->hCoarse, 0.0);

	double max_uO = mesh->calculateMax_uO();

	double oSigma1 = mesh->uO[(int) pos.y*mesh->unityCoarse.y + (int) pos.x]/max_uO;
	double oSigma2 = mesh->uO[(int) pos.y*mesh->unityCoarse.y + ((int) pos.x + 1)]/max_uO;
	double oSigma3 = mesh->uO[((int) pos.y + 1)*mesh->unityCoarse.y + (int) pos.x]/max_uO;
	double oSigma4 = mesh->uO[((int) pos.y + 1)*mesh->unityCoarse.y + ((int) pos.x + 1)]/max_uO;

	double egfSigma1 = mesh->uEGF[(int) pos.y*mesh->unityCoarse.y + (int) pos.x];
	double egfSigma2 = mesh->uEGF[(int) pos.y*mesh->unityCoarse.y + ((int) pos.x + 1)];
	double egfSigma3 = mesh->uEGF[((int) pos.y + 1)*mesh->unityCoarse.y + (int) pos.x];
	double egfSigma4 = mesh->uEGF[((int) pos.y + 1)*mesh->unityCoarse.y + ((int) pos.x + 1)];

	cell->sigmaO = (1.0 - c.x)*(1.0 - c.y)*oSigma1 + c.x*(1.0 - c.y)*oSigma2 + (1.0 - c.x)*c.y*oSigma3 + c.x*c.y*oSigma4;
	cell->sigmaEGF = (1.0 - c.x)*(1.0 - c.y)*egfSigma1 + c.x*(1.0 - c.y)*egfSigma2 + (1.0 - c.x)*c.y*egfSigma3 + c.x*c.y*egfSigma4;

	for (int chemDrug = 0; chemDrug < mesh->uChem.size(); chemDrug++) {
		double max_uChem = mesh->calculateMax_uChem(chemDrug);

		double chemSigma1 = mesh->uChem[chemDrug][(int) pos.y*mesh->unityCoarse.y + (int) pos.x]/max_uChem;
		double chemSigma2 = mesh->uChem[chemDrug][(int) pos.y*mesh->unityCoarse.y + ((int) pos.x + 1)]/max_uChem;
		double chemSigma3 = mesh->uChem[chemDrug][((int) pos.y + 1)*mesh->unityCoarse.y + (int) pos.x]/max_uChem;
		double chemSigma4 = mesh->uChem[chemDrug][((int) pos.y + 1)*mesh->unityCoarse.y + ((int) pos.x + 1)]/max_uChem;

		cell->sigmaChem[chemDrug] = (1.0 - c.x)*(1.0 - c.y)*chemSigma1 + c.x*(1.0 - c.y)*chemSigma2 + (1.0 - c.x)*c.y*chemSigma3 + c.x*c.y*chemSigma4;
	}
}

Cell CellManipulation2D::divide(Cell* cell, double rand1, CellCyclePathway* cellCyclePathway) {
	double theta = 2.0*M_PI*rand1;
	double nucleusRadius = cell->nucleusRadius/sqrt(2);
	Vector3 pos = cell->coordinates;
	Vector3 dist(nucleusRadius*cos(theta), nucleusRadius*sin(theta), 0.0);

	cell->previousState = PROLIFERATIVE2;
	cell->state = PROLIFERATIVE1;
	cell->coordinates.x = pos.x + dist.x;
	cell->coordinates.y = pos.y + dist.y;
	cell->nucleusRadius = nucleusRadius;
	cell->radius /= sqrt(2);
	cell->actionRadius /= sqrt(2);
	cell->cellCycleConc[5] *= 0.5;

	Cell child = *cell;
	child.coordinates.x = pos.x - dist.x;
	child.coordinates.y = pos.y - dist.y;
	child.mu = cellCyclePathway->generateMu(rand1);

	return child;
}

Vector3 CellManipulation2D::normal(Vector3 coordinates, double domainRadius) {
	//Direção
	Vector3 u(coordinates.x - domainRadius, coordinates.y - domainRadius, 0.0);
	// double u1 = , u2 =  u3 =;
	double D = norma(u);
	//Vesor
	u /= D; // u1 = u1/D; u2 = u2/D; u3 = u3/D;
	//Ponto na casca esférica
	u *= domainRadius; // double X = u1*domainRadius, Y = u2*domainRadius, Z = u3*domainRadius;

	return coordinates - u; //Vector3(x_cel - u.x, y_cel - u.y, z_cel - u.z);
}

void CellManipulation2D::updateFrame(Frame* frame, ConfigHandler* config, Mesh* mesh, EGFRPathway* egfrPathway, CellCyclePathway* cellCyclePathway, Ran* ran, Chemotherapy* chemotherapy, PharmacologicModel* pharmacologic) {
	if (config->output.prints)
		std::cout << "Updating frame..." << std::endl;

	std::vector<int> cellsRemoved;
	std::vector<Cell> cellsInserted;
	std::vector<double> randomNumbers = std::vector<double>(frame->cells.size() + 1, 0.0);
	// double randomNumbers[frame->cells.size()+1];
	for (int i = 0; i < randomNumbers.size(); i++)
		randomNumbers[i] = ran->doub();

	// #pragma omp parallel for schedule(static) ordered
	for (int i = 0; i < frame->cells.size(); i++) {
		//Célula Normal
		if (frame->cells[i].state == NORMOXIC || frame->cells[i].state > LAST_CELL_STATE)
			continue;

		calculateCellSigmas(&frame->cells[i], mesh);

		//EGFRPathway
		double ROC_PLC, ROC_ERK;
		if (frame->cells[i].state == QUIESCENT || frame->cells[i].state == MIGRATE) {
			egfrPathway->runge_kutta(frame->cells[i].sigmaEGF, &ROC_PLC, &ROC_ERK);

			/*if (ROC_PLC < egfrPathway->ROC_PLC_th && ROC_ERK > egfrPathway->ROC_ERK_th)
				config->parameters.alphaP = 1;
			else {
				config->parameters.alphaP = 0;

				if (ROC_PLC > egfrPathway->ROC_PLC_th && ROC_ERK < egfrPathway->ROC_ERK_th)
					config->parameters.alphaM = 1;
				else
					config->parameters.alphaM = 0;
			}*/

			config->parameters.alphaP = (ROC_PLC < egfrPathway->ROC_PLC_th && ROC_ERK > egfrPathway->ROC_ERK_th) ? 1.0 : 0.0;
			config->parameters.alphaM = (ROC_PLC > egfrPathway->ROC_PLC_th && ROC_ERK < egfrPathway->ROC_ERK_th) ? 1.0 : 0.0;
		}

		//CellCyclePathway
		if (frame->cells[i].state == QUIESCENT || frame->cells[i].state == PROLIFERATIVE2 || frame->cells[i].state == HYPOXIC || frame->cells[i].state == PROLIFERATIVE1 || frame->cells[i].state == MIGRATE) {
			cellCyclePathway->runge_kutta(&frame->cells[i].cellCycleConc[0], frame->cells[i].HIF, frame->cells[i].mu);
		}

		if (frame->cells[i].state == QUIESCENT) {
			//Transição Quiescente para Proliferativa
			// double alpha_p = config->parameters.alphaP*(frame->cells[i].sigmaO - config->parameters.sigmaH)/(1.0 - config->parameters.sigmaH)*(1.0 - frame->outCells/config->agent.maxOutCells); // Rocha, 2018
			double alpha_p = config->parameters.alphaP*(frame->cells[i].sigmaO - config->parameters.sigmaH)/(1.0 - config->parameters.sigmaH); // Powathil, 2012
			if (1.0 - exp(-config->parameters.delta_tt*alpha_p) >= randomNumbers[i] && frame->cells[i].cellCycleConc[0] > cellCyclePathway->CycB_th) {
				changeCellState(&frame->cells[i], QUIESCENT, PROLIFERATIVE2, frame->time, config->agent.oConsumption, 0.0, config->agent.chemUptake);
				continue;
			}

			//Transição Quiescente para Migratória
			/*double alpha_m = config->parameters.alphaM*(frame->cells[i].sigmaO - config->parameters.sigmaH)/(1.0 - config->parameters.sigmaH);
			if (1.0 - exp(-config->parameters.delta_tt*alpha_m) >= randomNumbers[i]) {
				changeCellState(&frame->cells[i], QUIESCENT, MIGRATE, frame->time, config->agent.oConsumption, 0.0, config->agent.chemUptake);
				continue;
			}*/

			//Transição Quiescente para Apoptótica
			if (1.0 - exp(-config->parameters.delta_tt*config->parameters.alphaA) >= randomNumbers[i]) {
				changeCellState(&frame->cells[i], QUIESCENT, APOPTOTIC, frame->time, 0.0, 0.0, 0.0);
				frame->cells[i].calcification = frame->cells[i].radius/(config->parameters.tauA*0.5);
				continue;
			}
		}

		//Transição Migratória para Quiescente
		/*if (frame->cells[i].state == MIGRATE && frame->cells[i].lifetime < frame->time) {
			changeCellState(&frame->cells[i], MIGRATE, QUIESCENT, frame->time, config->agent.oConsumption, config->agent.egfSource, config->agent.chemUptake);
			continue;
		}*/

		//Transição para Hipóxica
		// frame->cells[i].state != HYPOXIC && frame->cells[i].state != NECROTIC && frame->cells[i].sigmaO < config->parameters.sigmaH
		if ((frame->cells[i].state == QUIESCENT || frame->cells[i].state == PROLIFERATIVE2 || frame->cells[i].state == PROLIFERATIVE1 || frame->cells[i].state == MIGRATE) && 
			frame->cells[i].sigmaO < config->parameters.sigmaH) {
			changeCellState(&frame->cells[i], frame->cells[i].state, HYPOXIC, frame->time, config->agent.oConsumption, 0.0, config->agent.chemUptake);
			frame->cells[i].HIF = 1;
			continue; // GUSTAVO: Teste
		}

		//Transição Hipóxica para Necrótica
		if (frame->cells[i].state == HYPOXIC) {
			changeCellState(&frame->cells[i], HYPOXIC, NECROTIC, frame->time, 0.0, 0.0, 0.0);
			frame->cells[i].calcification = frame->cells[i].radius;
			continue; // GUSTAVO: Teste
		}

		double tau = frame->time - frame->cells[i].lifetime;

		//Ocorre a Apoptose
		if (frame->cells[i].state == APOPTOTIC) {
			if (tau > config->parameters.tauA*0.5 && tau < config->parameters.tauA) {
				frame->cells[i].nucleusRadius -= frame->cells[i].calcification*0.5;
				frame->cells[i].radius -= frame->cells[i].calcification;
			} else if (tau >= config->parameters.tauA) {
				// #pragma omp ordered
				cellsRemoved.push_back(i);
				//}
				//frame->cells.erase (frame->cells.begin() + i);
				frame->tumorCells--;
			}

			continue;
		}

		//Ocorre a Necrose
		if (frame->cells[i].state == NECROTIC) {
			if (tau <= config->parameters.tauNL) {
				frame->cells[i].radius = frame->cells[i].calcification * pow(1.0 + config->parameters.fNS*tau/config->parameters.tauNL, 1.0/3.0);
				//frame->cells[i].actionRadius = 1.214*frame->cells[i].radius;
				// continue;
			} else if (tau <= config->parameters.tauC) {
				if (tau - config->parameters.tauNL < 8)
					frame->cells[i].radius *= 0.9; // frame->cells[i].radius - 0.1*frame->cells[i].radius;
				else
					frame->cells[i].radius = frame->cells[i].nucleusRadius;
				//frame->cells[i].actionRadius = 1.214*frame->cells[i].radius;
				frame->cells[i].calcification = tau/config->parameters.tauC;
				// continue;
			} else {
				frame->cells[i].calcification = 1.0;
				// continue;
			}

			continue;
		}

		if (frame->cells[i].state == PROLIFERATIVE2) {
			//Transição Proliferativa para Morta
			for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++) {
				double alpha_k = (frame->cells[i].sigmaChem[chemDrug] >= config->parameters.chemThreshold) ? (frame->cells[i].sigmaChem[chemDrug] - config->parameters.chemThreshold)/(1.0 - config->parameters.chemThreshold) : 0.0;
				if (tau < config->parameters.tauP - config->parameters.tauG1 && 
					(chemotherapy->chemDrugs[chemDrug].compare("C") == 0 || chemotherapy->chemDrugs[chemDrug].compare("D") == 0) &&
					1.0 - exp(-config->parameters.delta_tt*alpha_k) >= randomNumbers[i]) {
					changeCellState(&frame->cells[i], PROLIFERATIVE2, KILLED, frame->time, 0.0, 0.0, 0.0);
					continue;
				}
			}

			//Ocorre a divisão celular
			if (tau == config->parameters.tauP - config->parameters.tauG1) {
				//frame->cells.push_back(divide(&frame->cells[i], config->input.seed));
				// #pragma omp ordered
				cellsInserted.push_back(divide(&frame->cells[i], randomNumbers[i], cellCyclePathway));
				frame->tumorCells++;
				continue;
			}
		}

		//Ocorre a G1
		if (frame->cells[i].state == PROLIFERATIVE1 && tau > config->parameters.tauP - config->parameters.tauG1) {
			if (tau < config->parameters.tauP) {
				for (int chemDrug = 0; chemDrug < chemotherapy->chemDrugs.size(); chemDrug++) {
					double alpha_k = (frame->cells[i].sigmaChem[chemDrug] >= config->parameters.chemThreshold) ? (frame->cells[i].sigmaChem[chemDrug] - config->parameters.chemThreshold)/(1.0 - config->parameters.chemThreshold) : 0.0;
					if ((chemotherapy->chemDrugs[chemDrug].compare("A") == 0 || chemotherapy->chemDrugs[chemDrug].compare("B") == 0) &&
						1.0 - exp(-config->parameters.delta_tt*alpha_k) >= randomNumbers[i]) {
						changeCellState(&frame->cells[i], PROLIFERATIVE1, KILLED, frame->time, 0.0, 0.0, 0.0);
						continue;
					}
				}

				frame->cells[i].nucleusRadius = config->agent.nucleusRadius*pow(0.5*(1.0 + (config->parameters.tauG1 + tau - config->parameters.tauP)/config->parameters.tauG1), 1.0/3.0);
				frame->cells[i].radius = config->agent.radius*pow(0.5*(1.0 + (config->parameters.tauG1 + tau - config->parameters.tauP)/config->parameters.tauG1), 1.0/3.0);
				frame->cells[i].actionRadius = 1.214*frame->cells[i].radius;
				continue;
			}

			changeCellState(&frame->cells[i], PROLIFERATIVE1, QUIESCENT, frame->time, config->agent.oConsumption, config->agent.egfSource, config->agent.chemUptake);
			frame->cells[i].nucleusRadius = config->agent.nucleusRadius;
			frame->cells[i].radius = config->agent.radius;
			frame->cells[i].actionRadius = config->agent.actionRadius;
			continue;
		}

		//Ocorre a Morte por Quimioterapia
		if (frame->cells[i].state == KILLED) {
			// #pragma omp ordered
			cellsRemoved.push_back(i);
			frame->tumorCells--;
		}
	}

	std::sort(cellsRemoved.begin(), cellsRemoved.end());

	//cellsRemoved Células
	for (int i = cellsRemoved.size() - 1; i >= 0; i--)
		frame->cells.erase(frame->cells.begin() + cellsRemoved[i]);
	
	//Insere Células
	for (int i = 0; i < cellsInserted.size(); i++)
		frame->cells.push_back(cellsInserted[i]);
}

void CellManipulation2D::force(Frame* frame, ConfigHandler* config, std::vector<Vector3> bloodVessel) {
	if (config->output.prints)
		std::cout << "Forces...";

	double iter = 0.0, iterMax = 60.0;
	double domainRadius = frame->domain.y/2.0;

	int count = 0;
	while (iter < iterMax) {
		int n = 1;
		int m = 1;
		double M = 1.0;

		//// #pragma omp parallel for
		for (int i = 0; i < frame->cells.size(); i++) {
			for (int j = i + 1; j < frame->cells.size(); j++) {
				double F_cca_i[2] = {0.0, 0.0};
				double F_ccr_i[2] = {0.0, 0.0};
				Vector3 r(frame->cells[j].coordinates.x - frame->cells[i].coordinates.x,
					frame->cells[j].coordinates.y - frame->cells[i].coordinates.y,
					0.0);
				double actionRadius = frame->cells[j].actionRadius + frame->cells[i].actionRadius;
				double nucleusRadius = frame->cells[j].nucleusRadius + frame->cells[i].nucleusRadius;
				double radius = frame->cells[j].radius + frame->cells[i].radius;

				if (norma(r) > actionRadius)
					continue;

				// double phi_x,phi_y,phi_z,psi_x,psi_y,psi_z;
				Vector3 phi = func_var_phi(r, actionRadius, n);
				Vector3 psi = func_var_psi(r, nucleusRadius, radius, M, m);

				// double c_cca;
				// //if (frame->cells[i].state != NORMOXIC && frame->cells[j].state != NORMOXIC)
				// c_cca = 0.488836;
				// //else
				// //c_cca = 0.588836;
				//double c_ccr = 10;

				//Força adesão célula-célula (i)
				F_cca_i[0] = config->forces.c_cca*phi.x;
				F_cca_i[1] = config->forces.c_cca*phi.y;

				//Força repulsão célula-célula (i)
				F_ccr_i[0] = config->forces.c_ccr*psi.x;
				F_ccr_i[1] = config->forces.c_ccr*psi.y;

				//// #pragma omp critical
				//{
				//Somatório de forças
				frame->cells[i].force.x += (F_cca_i[0] + F_ccr_i[0]);
				frame->cells[i].force.y += (F_cca_i[1] + F_ccr_i[1]);

				frame->cells[j].force.x -= (F_cca_i[0] + F_ccr_i[0]);
				frame->cells[j].force.y -= (F_cca_i[1] + F_ccr_i[1]);
				//}
			}

			for (int j = 0; j < bloodVessel.size(); j++) {
				double F_cca_i[2] = {0.0, 0.0};
				double F_ccr_i[2] = {0.0, 0.0};
				Vector3 r(bloodVessel[j].x - frame->cells[i].coordinates.x,
					bloodVessel[j].y - frame->cells[i].coordinates.y,
					0.0);
				double actionRadius = config->bloodVessel.actionRadius + frame->cells[i].actionRadius;
				// double nucleusRadius = config->agent.nucleusRadius + frame->cells[i].nucleusRadius; // GUSTAVO: Verificar
				double radius = config->bloodVessel.radius + frame->cells[i].radius;

				if (norma(r) > actionRadius)
					continue;

				Vector3 phi = func_var_phi(r, actionRadius, n);
				Vector3 psi = func_var_psi(r, radius, radius, M, m);

				//Força adesão célula-célula (i)
				F_cca_i[0] = config->forces.c_cca*phi.x;
				F_cca_i[1] = config->forces.c_cca*phi.y;

				//Força repulsão célula-célula (i)
				F_ccr_i[0] = config->forces.c_ccr*psi.x;
				F_ccr_i[1] = config->forces.c_ccr*psi.y;

				//Somatório de forças
				frame->cells[i].force.x += (F_cca_i[0] + F_ccr_i[0]);
				frame->cells[i].force.y += (F_cca_i[1] + F_ccr_i[1]);
			}

			//Força de rigidez do tecido na célula (i)
			double F_ct[2] = {0.0, 0.0};
			double F_rct[2] = {0.0, 0.0};
			// double K = 0.1;
			// double c_ct = 10.0 * K;
			// double c_rct = 4.88836 * K;
			// double phi_x2,phi_y2,phi_z2,psi_x2,psi_y2,psi_z2;

			Vector3 Dif(frame->cells[i].coordinates.x - domainRadius, frame->cells[i].coordinates.y - domainRadius, 0.0);

			if (norma(Dif) < domainRadius - frame->cells[i].actionRadius)
				continue;

			Vector3 _normal = normal(frame->cells[i].coordinates, domainRadius);
			Vector3 phi2 = func_var_phi(_normal, frame->cells[i].actionRadius, n);
			Vector3 psi2 = func_var_psi(_normal, frame->cells[i].nucleusRadius, frame->cells[i].radius, M, m);

			//Adesão
			F_ct[0] = config->forces.c_ct*phi2.x;
			F_ct[1] = config->forces.c_ct*phi2.y;

			//Repulsão
			F_rct[0] = config->forces.c_rct*psi2.x;
			F_rct[1] = config->forces.c_rct*psi2.y;

			//Somatório de forças
			frame->cells[i].force.x += (F_ct[0] + F_rct[0]);
			frame->cells[i].force.y += (F_ct[1] + F_rct[1]);
		}

		//Calcula velocidade
		// #pragma omp parallel for
		for (int i = 0; i < frame->cells.size(); i++) {
			frame->cells[i].speed.x = -0.5*frame->cells[i].force.x;
			frame->cells[i].speed.y = -0.5*frame->cells[i].force.y;

			if (frame->cells[i].state == MIGRATE) {
				frame->cells[i].speed.x -= 0.5*(config->forces.c_hap*frame->cells[i].grad_F_hap[0]);
				frame->cells[i].speed.y -= 0.5*(config->forces.c_hap*frame->cells[i].grad_F_hap[1]);
			}
		}

		//Norma da velocidade maior
		double v_max = 1.0e-10;
		// #pragma omp parallel for reduction(max : v_max)
		for (int i = 0; i < frame->cells.size(); i++) {
			if (norma(frame->cells[i].speed) > v_max)
				v_max = norma(frame->cells[i].speed);
		}
		double delta_tt = 1.0/v_max;

		// #pragma omp parallel for
		for (int i = 0; i < frame->cells.size(); i++) {
			//Desloca Células
			frame->cells[i].coordinates.x += delta_tt*frame->cells[i].speed.x;
			frame->cells[i].coordinates.y += delta_tt*frame->cells[i].speed.y;

			//Zera somatório de forças
			frame->cells[i].force = Vector3();
		}

		for (int i = 0; i < frame->cells.size(); i++) {
			double dist = sqrt(pow(frame->cells[i].coordinates.x - domainRadius, 2.0) + pow(frame->cells[i].coordinates.y - domainRadius, 2.0));
			//Se passar do domínio exclui
			if (dist >= domainRadius) {
				if (frame->cells[i].state != NORMOXIC)
					frame->tumorCells--;
				frame->outCells++;
				frame->cells.erase(frame->cells.begin() + i);
				i--;
			}
		}

		iter += delta_tt;
		count++;
	}

	if (config->output.prints)
		std::cout << " (" << count << " iterations)" << std::endl;
}
