#include "FileFactory.hpp"

using namespace std;

/*void FileFactory::makeFrameFile(Frame* frame, ConfigHandler* config){
	char name[40];
	FILE *arq;

	if (config->input.fileFormat == BI_D) {
		sprintf(name,"%s/%s%d-%05d-2D.dat", config->output.paths.agent.c_str(), config->output.filenames.agent.c_str(), config->output.filenames.number, frame->time);
		//********** Saving the data **********
		arq = fopen(name, "w");
		fprintf(arq, "%lf %lf\n", frame->domain.x, frame->domain.y);
		fprintf(arq, "%ld %d\n", frame->cells.size(), frame->time);
		fprintf(arq, "%d %d\n", frame->outCells, frame->tumorCells);

		for(int i = 0; i < frame->cells.size(); i++){
			fprintf(arq,"%d\n",frame->cells[i].state);
			fprintf(arq,"%lf %lf\n",frame->cells[i].coordinates.x,frame->cells[i].coordinates.y);
			fprintf(arq,"%lf %lf %lf %d %d\n",frame->cells[i].nucleusRadius,frame->cells[i].radius,frame->cells[i].actionRadius,frame->cells[i].lifetime,frame->cells[i].previousState);
			if (frame->cells[i].calcification < 0 || frame->cells[i].calcification > 1 ){
				fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].egfSource,0.0,frame->cells[i].sigmaEGF,frame->cells[i].sigmaO);
			}
			else{
				fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].egfSource,frame->cells[i].calcification,frame->cells[i].sigmaEGF,frame->cells[i].sigmaO);
			}
			fprintf(arq,"%lf %lf\n", frame->cells[i].speed.x, frame->cells[i].speed.y);
		}
	} else {
		sprintf(name,"%s/%s%d-%05d-3D.dat", config->output.paths.agent.c_str(), config->output.filenames.agent.c_str(), config->output.filenames.number, frame->time);
		//********** Saving the data **********
		arq = fopen(name, "w");
		fprintf(arq, "%lf %lf %lf\n", frame->domain.x, frame->domain.y, frame->domain.z);
		fprintf(arq, "%ld %d\n", frame->cells.size(), frame->time);
		fprintf(arq, "%d %d\n", frame->outCells, frame->tumorCells);

		for(int i = 0; i < frame->cells.size(); i++){
			fprintf(arq,"%d\n",frame->cells[i].state);
			fprintf(arq,"%lf %lf %lf\n",frame->cells[i].coordinates.x,frame->cells[i].coordinates.y, frame->cells[i].coordinates.z);
			fprintf(arq,"%lf %lf %lf %d %d\n",frame->cells[i].nucleusRadius,frame->cells[i].radius,frame->cells[i].actionRadius,frame->cells[i].lifetime,frame->cells[i].previousState);
			if (frame->cells[i].calcification < 0 || frame->cells[i].calcification > 1 ){
				fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].egfSource,0.0,frame->cells[i].sigmaEGF,frame->cells[i].sigmaO);
			}
			else{
				fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].egfSource,frame->cells[i].calcification,frame->cells[i].sigmaEGF,frame->cells[i].sigmaO);
			}
			fprintf(arq,"%lf %lf %lf\n", frame->cells[i].speed.x, frame->cells[i].speed.y, frame->cells[i].speed.z);
		}
	}
	fclose(arq);
}*/

/*void FileFactory::makeFile(Frame* frame, ConfigHandler* config, Mesh* mesh, OutputMode mode){
	FILE *arq;
	char name[40];

	if (config->input.fileFormat == BI_D) {
		if (mode == EGF) {
			sprintf(name,"%s/%s%d-%05d2D.dat", config->output.paths.egf.c_str(), config->output.filenames.egf.c_str(), config->output.filenames.number, frame->time);
		} else {
			sprintf(name,"%s/%s%d-%05d2D.dat", config->output.paths.nut.c_str(), config->output.filenames.nut.c_str(), config->output.filenames.number, frame->time);
		}
		arq = fopen(name,"w");
		fprintf(arq,"# %lf %lf\n", frame->domain.x, frame->domain.y);
		fprintf(arq,"# %d %d\n", mesh->unityCoarse.x, mesh->unityCoarse.y);
		fprintf(arq,"# %d %d \n", 0, frame->time);
		if (mode == EGF) {
			for(int i = 0; i < mesh->matrixSize; i++){
				fprintf(arq,"%lf %lf %lf\n",    ( i % (int)mesh->unityCoarse.x ) * mesh->hCoarse,
					( i / (int)mesh->unityCoarse.x ) * mesh->hCoarse,
					mesh->uEGF[i]);
			}
		} else {
			for(int i = 0; i < mesh->matrixSize; i++){
				fprintf(arq,"%lf %lf %lf\n",    ( i % (int)mesh->unityCoarse.x ) * mesh->hCoarse,
					( i / (int)mesh->unityCoarse.x ) * mesh->hCoarse,
					mesh->uO[i]);
			}
		}
	} else {
		if (mode == EGF) {
			sprintf(name,"%s/%s%d-%05d3D.dat", config->output.paths.egf.c_str(), config->output.filenames.egf.c_str(), config->output.filenames.number, frame->time);
		} else {
			sprintf(name,"%s/%s%d-%05d3D.dat", config->output.paths.nut.c_str(), config->output.filenames.nut.c_str(), config->output.filenames.number, frame->time);
		}
		arq = fopen(name,"w");
		fprintf(arq,"# %lf %lf %lf \n", frame->domain.x, frame->domain.y, frame->domain.z);
		fprintf(arq,"# %d %d %d\n", mesh->unityCoarse.x, mesh->unityCoarse.y, mesh->unityCoarse.z);
		fprintf(arq,"# %d %d \n",0, frame->time);
		if (mode == EGF) {
			for(int i = 0; i < mesh->matrixSize; i++){
				fprintf(arq,"%lf %lf %lf %lf\n",    ( ( i % (int)( frame->domain.x * frame->domain.y ) ) % (int)frame->domain.x ) * mesh->hCoarse,
					( ( i % (int)( frame->domain.x * frame->domain.y ) ) / (int)frame->domain.x ) * mesh->hCoarse,
					( i / (int)( frame->domain.x * frame->domain.y ) ) * mesh->hCoarse,
					mesh->uEGF[i]);
			}
		} else {
			for(int i = 0; i < mesh->matrixSize; i++){
				fprintf(arq,"%lf %lf %lf %lf\n",    ( ( i % (int)( frame->domain.x * frame->domain.y ) ) % (int)frame->domain.x ) * mesh->hCoarse,
					( ( i % (int)( frame->domain.x * frame->domain.y ) ) / (int)frame->domain.x ) * mesh->hCoarse,
					( i / (int)( frame->domain.x * frame->domain.y ) ) * mesh->hCoarse,
					mesh->uO[i]);
			}
		}
	}

	fclose(arq);
}*/

void FileFactory::makeFrameFile(Frame* frame, ConfigHandler* config) {
	ostringstream formattedTime;
	ostringstream name;
	ofstream arq;

	formattedTime.fill('0');
	formattedTime.width(5);
	formattedTime << frame->time;
	if (config->input.fileFormat == BI_D) {
		name << config->output.paths.agent.c_str() << "/" << config->output.filenames.agent.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-2D.dat";

		//********** Saving the data **********
		arq.open(name.str());

		arq << frame->domain.x << " " << frame->domain.y << endl;
		arq << frame->cells.size() << " " << frame->time << endl;
		arq << frame->outCells << " " << frame->tumorCells << endl;

		for (int i = 0; i < frame->cells.size(); i++) {
			arq << frame->cells[i].state << endl;
			arq << frame->cells[i].coordinates.x << " " << frame->cells[i].coordinates.y << endl;
			arq << frame->cells[i].nucleusRadius << " " << frame->cells[i].radius << " " << frame->cells[i].actionRadius << " " << frame->cells[i].lifetime << " " << frame->cells[i].previousState << endl;

			if (frame->cells[i].calcification < 0 || frame->cells[i].calcification > 1)
				arq << frame->cells[i].oConsumption << " " << frame->cells[i].egfSource << " " << 0.0 << " " << frame->cells[i].sigmaEGF << " " << frame->cells[i].sigmaO << endl;
			else
				arq << frame->cells[i].oConsumption << " " << frame->cells[i].egfSource << " " << frame->cells[i].calcification << " " << frame->cells[i].sigmaEGF << " " << frame->cells[i].sigmaO << endl;

			arq << frame->cells[i].speed.x << " " << frame->cells[i].speed.y << endl;
		}
	} else {
		name << config->output.paths.agent.c_str() << "/" << config->output.filenames.agent.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-3D.dat";

		//********** Saving the data **********
		arq.open(name.str());

		arq << frame->domain.x << " " << frame->domain.y << " " << frame->domain.z << endl;
		arq << frame->cells.size() << " " << frame->time << endl;
		arq << frame->outCells << " " << frame->tumorCells << endl;

		for (int i = 0; i < frame->cells.size(); i++) {
			arq << frame->cells[i].state << endl;
			arq << frame->cells[i].coordinates.x << " " << frame->cells[i].coordinates.y << " " << frame->cells[i].coordinates.z << endl;
			arq << frame->cells[i].nucleusRadius << " " << frame->cells[i].radius << " " << frame->cells[i].actionRadius << " " << frame->cells[i].lifetime << " " << frame->cells[i].previousState << endl;

			if (frame->cells[i].calcification < 0 || frame->cells[i].calcification > 1)
				arq << frame->cells[i].oConsumption << " " << frame->cells[i].egfSource << " " << 0.0 << " " << frame->cells[i].sigmaEGF << " " << frame->cells[i].sigmaO << endl;
			else
				arq << frame->cells[i].oConsumption << " " << frame->cells[i].egfSource << " " << frame->cells[i].calcification << " " << frame->cells[i].sigmaEGF << " " << frame->cells[i].sigmaO << endl;

			arq << frame->cells[i].speed.x << " " << frame->cells[i].speed.y << " " << frame->cells[i].speed.z << endl;
		}
	}

	arq.close();
}

void FileFactory::makeFile(Frame* frame, ConfigHandler* config, Mesh* mesh, Chemotherapy* chemotherapy, OutputMode mode, int chemDrug) {
	ostringstream formattedTime;
	ostringstream name;
	ofstream arq;

	formattedTime.fill('0');
	formattedTime.width(5);
	formattedTime << frame->time;
	if (config->input.fileFormat == BI_D) {
		switch (mode) {
			case EGF:
			name << config->output.paths.egf.c_str() << "/" << config->output.filenames.egf.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-2D.dat";
			break;

			case NUT:
			name << config->output.paths.nut.c_str() << "/" << config->output.filenames.nut.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-2D.dat";
			break;

			case CHEM:
			name << config->output.paths.chem.c_str() << "/" << config->output.filenames.chem.c_str() << config->output.filenames.number << "-" << chemotherapy->chemDrugs[chemDrug] << "-" << formattedTime.str() << "-2D.dat";
			break;
		}

		arq.open(name.str());

		arq << "# " << frame->domain.x << " " << frame->domain.y << endl;
		arq << "# " << mesh->unityCoarse.x << " " << mesh->unityCoarse.y << endl;
		arq << "# " << 0 << " " << frame->time << endl;

		double max_uO, max_uChem;
		switch (mode) {
			case EGF:
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->uEGF[i] << endl;
			break;

			case NUT:
			max_uO = mesh->calculateMax_uO();
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->uO[i]/max_uO << endl;
			break;

			case CHEM:
			max_uChem = mesh->calculateMax_uChem(chemDrug);
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->uChem[chemDrug][i]/max_uChem << endl;
			break;
		}
	} else {
		switch (mode) {
			case EGF:
			name << config->output.paths.egf.c_str() << "/" << config->output.filenames.egf.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-3D.dat";
			break;

			case NUT:
			name << config->output.paths.nut.c_str() << "/" << config->output.filenames.nut.c_str() << config->output.filenames.number << "-" << formattedTime.str() << "-3D.dat";
			break;

			case CHEM:
			name << config->output.paths.chem.c_str() << "/" << config->output.filenames.chem.c_str() << config->output.filenames.number << "-" << chemotherapy->chemDrugs[chemDrug] << "-" << formattedTime.str() << "-3D.dat";
			break;
		}

		arq.open(name.str());

		arq << "# " << frame->domain.x << " " << frame->domain.y << " " << frame->domain.z << endl;
		arq << "# " << mesh->unityCoarse.x << " " << mesh->unityCoarse.y << " " << mesh->unityCoarse.z << endl;
		arq << "# " << 0 << " " << frame->time << endl;

		double max_uO, max_uChem;
		switch (mode) {
			case EGF:
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->pos[i].z << " " << mesh->uEGF[i] << endl;
			break;

			case NUT:
			max_uO = mesh->calculateMax_uO();
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->pos[i].z << " " << mesh->uO[i]/max_uO << endl;
			break;

			case CHEM:
			max_uChem = mesh->calculateMax_uChem(chemDrug);
			for (int i = 0; i < mesh->matrixSize; i++)
				arq << mesh->pos[i].x << " " << mesh->pos[i].y << " " << mesh->pos[i].z << " " << mesh->uChem[chemDrug][i]/max_uChem << endl;
			break;
		}
	}

	arq.close();
}

void FileFactory::makePKFile(PharmacologicModel* pharmacologic, ConfigHandler* config) {
	ostringstream name;
	ofstream arq;

	name << config->output.paths.agent.c_str() << "/" << "PK.dat";

	arq.open(name.str());

	for (int i = 0; i < pharmacologic->pkSize; i++)
		arq << pharmacologic->time[i] << " " << pharmacologic->pkCurve[i] << endl;

	arq.close();
}

void FileFactory::makePDFile(PharmacologicModel* pharmacologic, ConfigHandler* config) {
	ostringstream name;
	ofstream arq;

	name << config->output.paths.agent.c_str() << "/" << "PD.dat";
	
	arq.open(name.str());

	for (int i = 0; i < pharmacologic->pdSize; i++)
		arq << pharmacologic->concentration[i] << " " << pharmacologic->pdCurve[i] << endl;

	arq.close();
}
