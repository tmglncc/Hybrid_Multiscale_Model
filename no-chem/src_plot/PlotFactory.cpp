#include "PlotFactory.hpp"

using namespace std;

void PlotFactory::makePlot(char* inputFilename, char* outputFilename, const char* outputMode) {
	char comment;
	Vector3 frameDomain;
	Vector3i meshUnityCoarse;
	int frameCellsSize, frameTime;

	ifstream input(inputFilename);

	input >> comment >> frameDomain.x >> frameDomain.y;
	input >> comment >> meshUnityCoarse.x >> meshUnityCoarse.y;
	input >> comment >> frameCellsSize >> frameTime;

	input.close();

	char pltFilename[40];
	strcpy(pltFilename, outputFilename);
	strcat(pltFilename, ".plt");

	char epsFilename[40];
	strcpy(epsFilename, outputFilename);
	strcat(epsFilename, ".eps");

	ofstream output(pltFilename);

	output << "set terminal postscript eps enhanced color font 'Helvetica, 16'" << endl;
	output << "set colorsequence classic" << endl;
	output << "set encoding utf8" << endl;
	output << "set output '" << epsFilename << "'" << endl;

	// output << "set title '" << outputMode << " concentration - Agents: " << frameCellsSize << " Time step: " << frameTime << "'" << endl;
	output << "set label '{/Symbol m}m' at screen 0.775, screen 0.115" << endl;
	output << "set label '{/Symbol m}m' at screen 0.232, screen 0.900" << endl;
	// output << "set label 'nM/s' at screen 0.760,screen 0.893" << endl;
	output << "set cblabel '" << outputMode << " concentration' rotate by 270" << endl;
	output << "set autoscale" << endl;
	output << "set xtics 0, " << frameDomain.x/5 << ", " << frameDomain.x << endl;
	output << "set ytics 0, " << frameDomain.y/5 << ", " << frameDomain.y << endl;
	// output << "stats '" << inputFilename << "' using 3 nooutput" << endl;
	// output << "if(abs(STATS_max) == 0.0) set cbrange [0:1]" << endl;
	output << "set cbrange [0:1]" << endl;
	output << "set format cb '%3.1f'" << endl;
	output << "set border linewidth 2" << endl;

	output << "set dgrid " << meshUnityCoarse.y << ", " << meshUnityCoarse.x << endl;
	output << "set palette rgbformulae 33, 13, 10" << endl;
	output << "set view map" << endl;
	output << "set pm3d interpolate 0, 0" << endl;
	output << "set size ratio -1" << endl;
	output << "splot '" << inputFilename <<"' using 1:2:3 with pm3d notitle" << endl;

	output.close();
}

void PlotFactory::makeFramePlot(char* inputFilename, char* outputFilename, const char* bvFilename, double bvRadius) {
	Vector3 frameDomain;
	int frameCellsSize, frameTime;
	int frameOutCells, frameTumorCells;
	int numVessels;

	int cellState, cellPreviousState, cellLifetime;
	Vector3 cellCoordinates;
	double cellNucleusRadius, cellRadius, cellActionRadius;
	double cellOConsumption, cellEGFSource, cellCalcification, cellSigmaEGF, cellSigmaO;
	Vector3 cellSpeed;
	Vector3 bvCoordinates;

	char texFilename[40];
	strcpy(texFilename, outputFilename);
	strcat(texFilename, ".tex");

	ifstream input(inputFilename);
	ifstream bvInput(bvFilename);
	ofstream output(texFilename);

	input >> frameDomain.x >> frameDomain.y;
	input >> frameCellsSize >> frameTime;
	input >> frameOutCells >> frameTumorCells;
	bvInput >> numVessels;

	output << "\\documentclass[a4paper,10pt]{article}" << endl;

	output << "\\usepackage[usenames,dvipsnames]{xcolor}" << endl;
	output << "\\usepackage{pstricks,pst-plot}" << endl;
	output << "\\definecolor{MyGray}{HTML}{EEEEEE}" << endl;

	output << "\\begin{document}" << endl;
	output << "\\thispagestyle{empty}" << endl;

	output << "\\begin{figure}[!htbp]" << endl;
	output << "\\small " << endl;
	output << "\\centering" << endl;
	// output << "{Agents: " << frameCellsSize << " Tumor cells: " << frameTumorCells << " Outside cells: " << frameOutCells << " Time step: " << frameTime << "} \\par \\medskip \\medskip" << endl;
	output << "{$t = " << frameTime << "$ h} \\par \\medskip \\medskip" << endl;
	output << "\\psset{unit=0.015cm}" << endl;

	output << "\\begin{pspicture}(" << frameDomain.x << "," << frameDomain.y << ")" << endl;
	output << "\\psaxes[linecolor=black,tickcolor=black,Dx=" << frameDomain.x/5 << ", Dy=" << frameDomain.y/5 << ",labelFontSize=\\color{black}]{-}(0,0)(" << frameDomain.x << "," << frameDomain.y << ")[\\textcolor{black}{$\\mu m$},0][\\textcolor{black}{$\\mu m$},90]" << endl;
	output << "\\psarc[linewidth=1pt,linecolor=Gray,fillcolor=MyGray,fillstyle=solid](" << frameDomain.y/2. << "," << frameDomain.y/2. << "){" << frameDomain.y/2. << "}{0}{360}" << endl;

	double domainRadius = frameDomain.y/2.0, dist;
	for (int i = 0; i < numVessels; i++) {
		bvInput >> bvCoordinates.x >> bvCoordinates.y;

		dist = sqrt(pow(bvCoordinates.x - domainRadius, 2.0) + pow(bvCoordinates.y - domainRadius, 2.0));
		if (dist < domainRadius)
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=Mahogany,linewidth=0.1pt](" << bvCoordinates.x << "," << bvCoordinates.y << "){" << bvRadius << "}" << endl;
	}

	for (int i = 0; i < frameCellsSize; i++) {
		input >> cellState;
		input >> cellCoordinates.x >> cellCoordinates.y;
		input >> cellNucleusRadius >> cellRadius >> cellActionRadius >> cellLifetime >> cellPreviousState;
		input >> cellOConsumption >> cellEGFSource >> cellCalcification >> cellSigmaEGF >> cellSigmaO;
		input >> cellSpeed.x >> cellSpeed.y;

		switch (cellState) {
			case NECROTIC:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!10!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=Orange!"<< cellCalcification*100 <<"!yellow,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case QUIESCENT:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!90!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case PROLIFERATIVE2:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=green!75!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case HYPOXIC:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=gray!50!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case APOPTOTIC:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=red!90!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case PROLIFERATIVE1:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=green!25!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case NORMOXIC:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!90!MidnightBlue,linewidth=0.1pt,opacity = 0.2](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!90!MidnightBlue,linewidth=0.1pt,opacity = 0.2](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case MIGRATE:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=yellow!90!MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;

			case KILLED:
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=black!90!MidnightBlue,linewidth=0.1pt,opacity = 0.2](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellRadius << "}" << endl;
			output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=black!90!MidnightBlue,linewidth=0.1pt,opacity = 0.2](" << cellCoordinates.x << "," << cellCoordinates.y << "){" << cellNucleusRadius << "}" << endl;
			break;
		}
	}

	output << "\\end{pspicture}" << endl;
	output << "\\end{figure}" << endl;
	output << "\\end{document}" << endl;

	input.close();
	bvInput.close();
	output.close();
}

void PlotFactory::makePhenotypeData(char* outputFilename, int timeMax, int inputFilenamesNumber) {
	Vector3 frameDomain;
	int frameCellsSize, frameTime_;
	int frameOutCells, frameTumorCells;
	int necrotic, quiescent, proliferative2, hypoxic, apoptotic, proliferative1, normoxic, migrate, killed;

	int cellState, cellPreviousState, cellLifetime;
	Vector3 cellCoordinates;
	double cellNucleusRadius, cellRadius, cellActionRadius;
	double cellOConsumption, cellEGFSource, cellCalcification, cellSigmaEGF, cellSigmaO;
	Vector3 cellSpeed;

	ostringstream formattedTime;
	ostringstream inputFilename;

	ofstream output(outputFilename);

	for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
		formattedTime.fill('0');
		formattedTime.width(5);
		formattedTime << frameTime;
		inputFilename << "out" << inputFilenamesNumber << "-" << formattedTime.str() << "-2D.dat";

		cout << "Reading " << inputFilename.str() << "..." << endl;

		ifstream input(inputFilename.str());

		input >> frameDomain.x >> frameDomain.y;
		input >> frameCellsSize >> frameTime_;
		input >> frameOutCells >> frameTumorCells;

		necrotic = 0;
		quiescent = 0;
		proliferative2 = 0;
		hypoxic = 0;
		apoptotic = 0;
		proliferative1 = 0;
		normoxic = 0;
		migrate = 0;
		killed = 0;

		for (int i = 0; i < frameCellsSize; i++) {
			input >> cellState;
			input >> cellCoordinates.x >> cellCoordinates.y;
			input >> cellNucleusRadius >> cellRadius >> cellActionRadius >> cellLifetime >> cellPreviousState;
			input >> cellOConsumption >> cellEGFSource >> cellCalcification >> cellSigmaEGF >> cellSigmaO;
			input >> cellSpeed.x >> cellSpeed.y;

			switch (cellState) {
				case NECROTIC:
				necrotic++;
				break;

				case QUIESCENT:
				quiescent++;
				break;

				case PROLIFERATIVE2:
				proliferative2++;
				break;

				case HYPOXIC:
				hypoxic++;
				break;

				case APOPTOTIC:
				apoptotic++;
				break;

				case PROLIFERATIVE1:
				proliferative1++;
				break;

				case NORMOXIC:
				normoxic++;
				break;

				case MIGRATE:
				migrate++;
				break;

				case KILLED:
				killed++;
				break;
			}
		}

		input.close();

		formattedTime.str("");
		inputFilename.str("");

		output << frameTime << "\t";
		output << frameTumorCells << "\t";

		output << necrotic << "\t";
		output << quiescent << "\t";
		output << proliferative2 << "\t";
		output << hypoxic << "\t";
		output << apoptotic << "\t";
		output << proliferative1 << "\t";
		output << normoxic << "\t";
		output << migrate << "\t";
		output << killed << "\t";
		output << proliferative1 + proliferative2 << endl;
	}

	output.close();
}

void PlotFactory::makePhenotypePlot(char* inputFilename, char* outputFilename) {
	char pltFilename[40];
	strcpy(pltFilename, outputFilename);
	strcat(pltFilename, ".plt");

	char pdfFilename[40];
	strcpy(pdfFilename, outputFilename);
	strcat(pdfFilename, ".pdf");

	ofstream output(pltFilename);

	output << "set terminal pdfcairo enhanced color font 'Helvetica, 20'" << endl;
	output << "set colorsequence classic" << endl;
	output << "set encoding utf8" << endl;
	output << "set output '" << pdfFilename << "'" << endl;

	output << "set xlabel 'Time (h)'" << endl;
	output << "set ylabel 'Number of cells'" << endl;
	output << "set autoscale" << endl;
	output << "set xtic auto" << endl;
	output << "set ytic auto" << endl;
	// output << "set yrange [0:500]" << endl;
	output << "set key inside left top Left reverse font 'Helvetica, 16' width 1.5 height 0.5 box linewidth 2" << endl;
	output << "set border linewidth 2" << endl;
	output << "set size ratio -1" << endl;

	output << "set style line 1 linecolor rgb '#607D8B' dashtype 1 linewidth 2" << endl;
	output << "set style line 2 linecolor rgb '#F46E2B' dashtype 1 linewidth 2" << endl;
	output << "set style line 3 linecolor rgb '#06608F' dashtype 1 linewidth 2" << endl;
	output << "set style line 4 linecolor rgb '#00DB24' dashtype 1 linewidth 2" << endl;
	output << "set style line 5 linecolor rgb '#4E768F' dashtype 1 linewidth 2" << endl;
	output << "set style line 6 linecolor rgb '#E60B0F' dashtype 1 linewidth 2" << endl;
	output << "set style line 7 linecolor rgb '#00946D' dashtype 1 linewidth 2" << endl;
	output << "set style line 8 linecolor rgb '#DDEBF3' dashtype 2 linewidth 2" << endl;
	output << "set style line 9 linecolor rgb '#FABB15' dashtype 1 linewidth 2" << endl;
	output << "set style line 10 linecolor rgb '#2B292C' dashtype 1 linewidth 2" << endl;
	output << "set style line 11 linecolor rgb '#00B849' dashtype 1 linewidth 2" << endl;
	
	output << "plot '";
	// output << inputFilename << "' using 1:2 with lines linestyle 1 title 'Total', '";
	output << inputFilename << "' using 1:3 with lines linestyle 2 title 'N', '";
	output << inputFilename << "' using 1:4 with lines linestyle 3 title 'Q', '";
	// output << inputFilename << "' using 1:5 with lines linestyle 4 title 'P_{+}', '";
	// output << inputFilename << "' using 1:6 with lines linestyle 5 title 'H', '";
	// output << inputFilename << "' using 1:7 with lines linestyle 6 title 'A', '";
	// output << inputFilename << "' using 1:8 with lines linestyle 7 title 'P_{-}', '";
	// output << inputFilename << "' using 1:9 with lines linestyle 8 title 'Normal', '";
	output << inputFilename << "' using 1:10 with lines linestyle 9 title 'M', '";
	// output << inputFilename << "' using 1:11 with lines linestyle 10 title 'K', '";
	output << inputFilename << "' using 1:12 with lines linestyle 11 title 'P'";
	output << endl;

	output.close();
}

void PlotFactory::makeStochasticData(char* outputFilename, int seeds, int timeMax) {
	int frameTime_, frameTumorCells;
	int necrotic, quiescent, proliferative2, hypoxic, apoptotic, proliferative1, normoxic, migrate, killed, proliferative;

	vector<double> sumFrameTumorCells(timeMax + 1, 0.0);
	vector<double> sumNecrotic(timeMax + 1, 0.0);
	vector<double> sumQuiescent(timeMax + 1, 0.0);
	vector<double> sumProliferative2(timeMax + 1, 0.0);
	vector<double> sumHypoxic(timeMax + 1, 0.0);
	vector<double> sumApoptotic(timeMax + 1, 0.0);
	vector<double> sumProliferative1(timeMax + 1, 0.0);
	vector<double> sumNormoxic(timeMax + 1, 0.0);
	vector<double> sumMigrate(timeMax + 1, 0.0);
	vector<double> sumKilled(timeMax + 1, 0.0);
	vector<double> sumProliferative(timeMax + 1, 0.0);

	vector<double> meanFrameTumorCells(timeMax + 1, 0.0);
	vector<double> meanNecrotic(timeMax + 1, 0.0);
	vector<double> meanQuiescent(timeMax + 1, 0.0);
	vector<double> meanProliferative2(timeMax + 1, 0.0);
	vector<double> meanHypoxic(timeMax + 1, 0.0);
	vector<double> meanApoptotic(timeMax + 1, 0.0);
	vector<double> meanProliferative1(timeMax + 1, 0.0);
	vector<double> meanNormoxic(timeMax + 1, 0.0);
	vector<double> meanMigrate(timeMax + 1, 0.0);
	vector<double> meanKilled(timeMax + 1, 0.0);
	vector<double> meanProliferative(timeMax + 1, 0.0);

	vector<double> sdFrameTumorCells(timeMax + 1, 0.0);
	vector<double> sdNecrotic(timeMax + 1, 0.0);
	vector<double> sdQuiescent(timeMax + 1, 0.0);
	vector<double> sdProliferative2(timeMax + 1, 0.0);
	vector<double> sdHypoxic(timeMax + 1, 0.0);
	vector<double> sdApoptotic(timeMax + 1, 0.0);
	vector<double> sdProliferative1(timeMax + 1, 0.0);
	vector<double> sdNormoxic(timeMax + 1, 0.0);
	vector<double> sdMigrate(timeMax + 1, 0.0);
	vector<double> sdKilled(timeMax + 1, 0.0);
	vector<double> sdProliferative(timeMax + 1, 0.0);

	ostringstream inputFilename;

	for (int seed = 1; seed <= seeds; seed++) {
		inputFilename << "phenotype" << seed << ".dat";

		cout << "Reading " << inputFilename.str() << "..." << endl;

		ifstream input(inputFilename.str());

		for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
			input >> frameTime_;
			input >> frameTumorCells;

			input >> necrotic;
			input >> quiescent;
			input >> proliferative2;
			input >> hypoxic;
			input >> apoptotic;
			input >> proliferative1;
			input >> normoxic;
			input >> migrate;
			input >> killed;
			input >> proliferative;

			sumFrameTumorCells[frameTime] += frameTumorCells;
			sumNecrotic[frameTime] += necrotic;
			sumQuiescent[frameTime] += quiescent;
			sumProliferative2[frameTime] += proliferative2;
			sumHypoxic[frameTime] += hypoxic;
			sumApoptotic[frameTime] += apoptotic;
			sumProliferative1[frameTime] += proliferative1;
			sumNormoxic[frameTime] += normoxic;
			sumMigrate[frameTime] += migrate;
			sumKilled[frameTime] += killed;
			sumProliferative[frameTime] += proliferative;
		}

		input.close();

		inputFilename.str("");
	}

	for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
		meanFrameTumorCells[frameTime] = sumFrameTumorCells[frameTime]/seeds;
		meanNecrotic[frameTime] = sumNecrotic[frameTime]/seeds;
		meanQuiescent[frameTime] = sumQuiescent[frameTime]/seeds;
		meanProliferative2[frameTime] = sumProliferative2[frameTime]/seeds;
		meanHypoxic[frameTime] = sumHypoxic[frameTime]/seeds;
		meanApoptotic[frameTime] = sumApoptotic[frameTime]/seeds;
		meanProliferative1[frameTime] = sumProliferative1[frameTime]/seeds;
		meanNormoxic[frameTime] = sumNormoxic[frameTime]/seeds;
		meanMigrate[frameTime] = sumMigrate[frameTime]/seeds;
		meanKilled[frameTime] = sumKilled[frameTime]/seeds;
		meanProliferative[frameTime] = sumProliferative[frameTime]/seeds;
	}

	fill(sumFrameTumorCells.begin(), sumFrameTumorCells.end(), 0.0);
	fill(sumNecrotic.begin(), sumNecrotic.end(), 0.0);
	fill(sumQuiescent.begin(), sumQuiescent.end(), 0.0);
	fill(sumProliferative2.begin(), sumProliferative2.end(), 0.0);
	fill(sumHypoxic.begin(), sumHypoxic.end(), 0.0);
	fill(sumApoptotic.begin(), sumApoptotic.end(), 0.0);
	fill(sumProliferative1.begin(), sumProliferative1.end(), 0.0);
	fill(sumNormoxic.begin(), sumNormoxic.end(), 0.0);
	fill(sumMigrate.begin(), sumMigrate.end(), 0.0);
	fill(sumKilled.begin(), sumKilled.end(), 0.0);
	fill(sumProliferative.begin(), sumProliferative.end(), 0.0);

	for (int seed = 1; seed <= seeds; seed++) {
		inputFilename << "phenotype" << seed << ".dat";

		cout << "Reading " << inputFilename.str() << "..." << endl;

		ifstream input(inputFilename.str());

		for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
			input >> frameTime_;
			input >> frameTumorCells;

			input >> necrotic;
			input >> quiescent;
			input >> proliferative2;
			input >> hypoxic;
			input >> apoptotic;
			input >> proliferative1;
			input >> normoxic;
			input >> migrate;
			input >> killed;
			input >> proliferative;

			sumFrameTumorCells[frameTime] += pow(frameTumorCells - meanFrameTumorCells[frameTime], 2.0);
			sumNecrotic[frameTime] += pow(necrotic - meanNecrotic[frameTime], 2.0);
			sumQuiescent[frameTime] += pow(quiescent - meanQuiescent[frameTime], 2.0);
			sumProliferative2[frameTime] += pow(proliferative2 - meanProliferative2[frameTime], 2.0);
			sumHypoxic[frameTime] += pow(hypoxic - meanHypoxic[frameTime], 2.0);
			sumApoptotic[frameTime] += pow(apoptotic - meanApoptotic[frameTime], 2.0);
			sumProliferative1[frameTime] += pow(proliferative1 - meanProliferative1[frameTime], 2.0);
			sumNormoxic[frameTime] += pow(normoxic - meanNormoxic[frameTime], 2.0);
			sumMigrate[frameTime] += pow(migrate - meanMigrate[frameTime], 2.0);
			sumKilled[frameTime] += pow(killed - meanKilled[frameTime], 2.0);
			sumProliferative[frameTime] += pow(proliferative - meanProliferative[frameTime], 2.0);
		}

		input.close();

		inputFilename.str("");
	}

	for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
		sdFrameTumorCells[frameTime] = sqrt(sumFrameTumorCells[frameTime]/(seeds - 1));
		sdNecrotic[frameTime] = sqrt(sumNecrotic[frameTime]/(seeds - 1));
		sdQuiescent[frameTime] = sqrt(sumQuiescent[frameTime]/(seeds - 1));
		sdProliferative2[frameTime] = sqrt(sumProliferative2[frameTime]/(seeds - 1));
		sdHypoxic[frameTime] = sqrt(sumHypoxic[frameTime]/(seeds - 1));
		sdApoptotic[frameTime] = sqrt(sumApoptotic[frameTime]/(seeds - 1));
		sdProliferative1[frameTime] = sqrt(sumProliferative1[frameTime]/(seeds - 1));
		sdNormoxic[frameTime] = sqrt(sumNormoxic[frameTime]/(seeds - 1));
		sdMigrate[frameTime] = sqrt(sumMigrate[frameTime]/(seeds - 1));
		sdKilled[frameTime] = sqrt(sumKilled[frameTime]/(seeds - 1));
		sdProliferative[frameTime] = sqrt(sumProliferative[frameTime]/(seeds - 1));
	}

	ofstream output(outputFilename);

	for (int frameTime = 0; frameTime <= timeMax; frameTime++) {
		output << frameTime << "\t";
		output << meanFrameTumorCells[frameTime] << "\t" << sdFrameTumorCells[frameTime] << "\t";
		output << meanNecrotic[frameTime] << "\t" << sdNecrotic[frameTime] << "\t";
		output << meanQuiescent[frameTime] << "\t" << sdQuiescent[frameTime] << "\t";
		output << meanProliferative2[frameTime] << "\t" << sdProliferative2[frameTime] << "\t";
		output << meanHypoxic[frameTime] << "\t" << sdHypoxic[frameTime] << "\t";
		output << meanApoptotic[frameTime] << "\t" << sdApoptotic[frameTime] << "\t";
		output << meanProliferative1[frameTime] << "\t" << sdProliferative1[frameTime] << "\t";
		output << meanNormoxic[frameTime] << "\t" << sdNormoxic[frameTime] << "\t";
		output << meanMigrate[frameTime] << "\t" << sdMigrate[frameTime] << "\t";
		output << meanKilled[frameTime] << "\t" << sdKilled[frameTime] << "\t";
		output << meanProliferative[frameTime] << "\t" << sdProliferative[frameTime] << endl;
	}

	output.close();
}

void PlotFactory::makeStochasticPlot(char* inputFilename, char* outputFilename) {
	char pltFilename[40];
	strcpy(pltFilename, outputFilename);
	strcat(pltFilename, ".plt");

	char pdfFilename[40];
	strcpy(pdfFilename, outputFilename);
	strcat(pdfFilename, ".pdf");

	ofstream output(pltFilename);

	output << "set terminal pdfcairo enhanced color font 'Helvetica, 20'" << endl;
	output << "set colorsequence classic" << endl;
	output << "set encoding utf8" << endl;
	output << "set output '" << pdfFilename << "'" << endl;

	output << "set xlabel 'Time (h)'" << endl;
	output << "set ylabel 'Number of cells'" << endl;
	output << "set autoscale" << endl;
	output << "set xtic auto" << endl;
	output << "set ytic auto" << endl;
	output << "set yrange [0:500]" << endl;
	output << "set key inside left top Left reverse font 'Helvetica, 16' width 1.5 height 0.5 box linewidth 2" << endl;
	output << "set border linewidth 2" << endl;
	output << "set size ratio -1" << endl;

	output << "set style line 1 linecolor rgb '#607D8B' dashtype 1 linewidth 2" << endl;
	output << "set style line 2 linecolor rgb '#F46E2B' dashtype 1 linewidth 2" << endl;
	output << "set style line 3 linecolor rgb '#06608F' dashtype 1 linewidth 2" << endl;
	output << "set style line 4 linecolor rgb '#00DB24' dashtype 1 linewidth 2" << endl;
	output << "set style line 5 linecolor rgb '#4E768F' dashtype 1 linewidth 2" << endl;
	output << "set style line 6 linecolor rgb '#E60B0F' dashtype 1 linewidth 2" << endl;
	output << "set style line 7 linecolor rgb '#00946D' dashtype 1 linewidth 2" << endl;
	output << "set style line 8 linecolor rgb '#DDEBF3' dashtype 2 linewidth 2" << endl;
	output << "set style line 9 linecolor rgb '#FABB15' dashtype 1 linewidth 2" << endl;
	output << "set style line 10 linecolor rgb '#2B292C' dashtype 1 linewidth 2" << endl;
	output << "set style line 11 linecolor rgb '#00B849' dashtype 1 linewidth 2" << endl;
	output << "set style fill transparent solid 0.4" << endl;
	
	output << "plot '";
	// output << inputFilename << "' using 1:($2+$3):($2-$3) with filledcurve fill fc rgb '#607D8B' notitle, '' using 1:2 smooth mcspline linestyle 1 title 'Total', '";
	output << inputFilename << "' using 1:($4+$5):($4-$5) with filledcurve fill fc rgb '#F46E2B' notitle, '' using 1:4 smooth mcspline linestyle 2 title 'N', '";
	output << inputFilename << "' using 1:($6+$7):($6-$7) with filledcurve fill fc rgb '#06608F' notitle, '' using 1:6 smooth mcspline linestyle 3 title 'Q', '";
	// output << inputFilename << "' using 1:($8+$9):($8-$9) with filledcurve fill fc rgb '#00DB24' notitle, '' using 1:8 smooth mcspline linestyle 4 title 'P_{+}', '";
	// output << inputFilename << "' using 1:($10+$11):($10-$11) with filledcurve fill fc rgb '#4E768F' notitle, '' using 1:10 smooth mcspline linestyle 5 title 'H', '";
	// output << inputFilename << "' using 1:($12+$13):($12-$13) with filledcurve fill fc rgb '#E60B0F' notitle, '' using 1:12 smooth mcspline linestyle 6 title 'A', '";
	// output << inputFilename << "' using 1:($14+$15):($14-$15) with filledcurve fill fc rgb '#00946D' notitle, '' using 1:14 smooth mcspline linestyle 7 title 'P_{-}, '";
	// output << inputFilename << "' using 1:($16+$17):($16-$17) with filledcurve fill fc rgb '#DDEBF3' notitle, '' using 1:16 smooth mcspline linestyle 8 title 'Normal', '";
	// output << inputFilename << "' using 1:($18+$19):($18-$19) with filledcurve fill fc rgb '#FABB15' notitle, '' using 1:18 smooth mcspline linestyle 9 title 'M', '";
	// output << inputFilename << "' using 1:($20+$21):($20-$21) with filledcurve fill fc rgb '#2B292C' notitle, '' using 1:20 smooth mcspline linestyle 10 title 'K', '";
	output << inputFilename << "' using 1:($22+$23):($22-$23) with filledcurve fill fc rgb '#00B849' notitle, '' using 1:22 smooth mcspline linestyle 11 title 'P'";
	output << endl;

	output.close();
}
