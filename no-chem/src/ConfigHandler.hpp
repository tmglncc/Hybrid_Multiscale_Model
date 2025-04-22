#ifndef CONFIG_HANDLER
#define CONFIG_HANDLER

#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <vector>
#include "Util.hpp"

class ConfigHandler {
private:
	libconfig::Config cfg;
	bool created_;

	void getInputValues();
	void getOutputValues();
	void getContinuumValues();
	void getBloodVesselValues();
	void getAgentValues();
	void getForcesValues();
	void getParametersValues();
	void getPharmacologicValues();

public:
	ConfigHandler(std::string configFile = "config.cfg", int seed = 7);
	bool created();

	struct {
		std::string initialCondition;
		std::string bloodVessel;
		int seed;
		int timeMax;
		FileFormat fileFormat;
	} input;

	struct {
		struct {
			std::string agent;
			std::string egf;
			std::string nut;
			std::string chem;
		} paths;

		struct {
			std::string agent;
			std::string egf;
			std::string nut;
			std::string chem;
			int number;
		} filenames;

		bool egf;
		bool nut;
		bool chem;
		bool files;
		bool prints;
		bool justLastFile;
		bool justLastPrint;
		int filesFrequency;
	} output;

	struct {
		struct {
			double oDiffusion;
			double oConsumptionBg;
			double oConsumptionBorder;
		} oxygen;

		struct {
			double egfDiffusion;
			double egfSourceBg;
			double egfSourceBorder;
		} egf;

		struct {
			struct {
				double chemDiffusion;
				double chemDecay;
			} CYC202;

			struct {
				double chemDiffusion;
				double chemDecay;
			} cisplatin;

			struct {
				double chemDiffusion;
				double chemDecay;
			} taxotere;

			struct {
				double chemDiffusion;
				double chemDecay;
			} taxol;

			double chemUptakeBg;
			double chemUptakeBorder;
			std::vector<int> chemTimes;
			std::vector<int> chemDosages;
			std::vector<std::string> chemProtocol;
		} chemotherapy;
		
		double hCoarse;
		double hRefined;
		double deltaT;
	} continuum;

	struct {
		double radius;
		double actionRadius;
		double oProduction;
	} bloodVessel;

	struct {
		double nucleusRadius;
		double radius;
		double actionRadius;
		double oConsumption;
		double egfSource;
		double chemUptake;
		int maxOutCells;
	} agent;

	struct {
		double c_cca;
		double c_ccr;
		double K;
		double c_ct;
		double c_rct;
		double c_hap;
	} forces;

	struct {
		double tauA;
		double tauP;
		double tauG1;
		double tauC;
		double tauNL;
		double fNS;
		double delta_tt;
		double alphaP;
		double sigmaH;
		double alphaA;
		double alphaM;
		double chemThreshold;
	} parameters;

	struct {
		struct {
			double p1;
			double p2;
			double p3;
			double tau;
			double lambda_a;
			double lambda_d;
			double zeta;
		} pk;
		struct {
			double lowerThreshold;
			double upperThreshold;
			double maxConcentration;
			double q_u;
		} pd;
	} pharmacologic;
};

#endif /* end of include guard: CONFIG_HANDLER */
