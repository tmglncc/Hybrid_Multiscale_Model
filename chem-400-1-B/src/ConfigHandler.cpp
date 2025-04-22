#include "ConfigHandler.hpp"

using namespace libconfig;

void ConfigHandler::getInputValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["input"];

	try {
		// conf.lookupValue("seed", this->input.seed);
		conf.lookupValue("time-max", this->input.timeMax);
		conf.lookupValue("initial-condition", this->input.initialCondition);
		conf.lookupValue("blood-vessel", this->input.bloodVessel);
		
		std::string fileFormat;
		conf.lookupValue("file-format", fileFormat);
		if (fileFormat == "2d" || fileFormat == "2D") {
			this->input.fileFormat = BI_D;
		} else if (fileFormat == "3d" || fileFormat == "3D") {
			this->input.fileFormat = TRI_D;
		}
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Input Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getOutputValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["output"];

	try {
		conf["paths"].lookupValue("agent", this->output.paths.agent);
		conf["paths"].lookupValue("egf", this->output.paths.egf);
		conf["paths"].lookupValue("nut", this->output.paths.nut);
		conf["paths"].lookupValue("chem", this->output.paths.chem);
		
		conf["filenames"].lookupValue("number", this->output.filenames.number);
		conf["filenames"].lookupValue("agent", this->output.filenames.agent);
		conf["filenames"].lookupValue("egf", this->output.filenames.egf);
		conf["filenames"].lookupValue("nut", this->output.filenames.nut);
		conf["filenames"].lookupValue("chem", this->output.filenames.chem);
		
		conf.lookupValue("egf", this->output.egf);
		conf.lookupValue("nut", this->output.nut);
		conf.lookupValue("chem", this->output.chem);
		conf.lookupValue("files", this->output.files);
		conf.lookupValue("prints", this->output.prints);
		conf.lookupValue("just-last-file", this->output.justLastFile);
		conf.lookupValue("just-last-print", this->output.justLastPrint);
		conf.lookupValue("files-frequency", this->output.filesFrequency);
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Output Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getContinuumValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["continuum"];
	
	try {
		conf["oxygen"].lookupValue("o-diffusion", this->continuum.oxygen.oDiffusion);
		conf["oxygen"].lookupValue("o-consumption-bg", this->continuum.oxygen.oConsumptionBg);
		conf["oxygen"].lookupValue("o-consumption-border", this->continuum.oxygen.oConsumptionBorder);

		conf["egf"].lookupValue("egf-diffusion", this->continuum.egf.egfDiffusion);
		conf["egf"].lookupValue("egf-source-bg", this->continuum.egf.egfSourceBg);
		conf["egf"].lookupValue("egf-source-border", this->continuum.egf.egfSourceBorder);

		conf["chemotherapy"]["CYC202"].lookupValue("chem-diffusion", this->continuum.chemotherapy.CYC202.chemDiffusion);
		conf["chemotherapy"]["CYC202"].lookupValue("chem-decay", this->continuum.chemotherapy.CYC202.chemDecay);
		conf["chemotherapy"]["cisplatin"].lookupValue("chem-diffusion", this->continuum.chemotherapy.cisplatin.chemDiffusion);
		conf["chemotherapy"]["cisplatin"].lookupValue("chem-decay", this->continuum.chemotherapy.cisplatin.chemDecay);
		conf["chemotherapy"]["taxotere"].lookupValue("chem-diffusion", this->continuum.chemotherapy.taxotere.chemDiffusion);
		conf["chemotherapy"]["taxotere"].lookupValue("chem-decay", this->continuum.chemotherapy.taxotere.chemDecay);
		conf["chemotherapy"]["taxol"].lookupValue("chem-diffusion", this->continuum.chemotherapy.taxol.chemDiffusion);
		conf["chemotherapy"]["taxol"].lookupValue("chem-decay", this->continuum.chemotherapy.taxol.chemDecay);
		conf["chemotherapy"].lookupValue("chem-uptake-bg", this->continuum.chemotherapy.chemUptakeBg);
		conf["chemotherapy"].lookupValue("chem-uptake-border", this->continuum.chemotherapy.chemUptakeBorder);

		const Setting& confChemTimes = conf["chemotherapy"].lookup("chem-times");
		const Setting& confChemDosages = conf["chemotherapy"].lookup("chem-dosages");
		const Setting& confChemProtocol = conf["chemotherapy"].lookup("chem-protocol");
		for (int i = 0; i < confChemTimes.getLength(); i++) {
			this->continuum.chemotherapy.chemTimes.push_back(confChemTimes[i]);
			this->continuum.chemotherapy.chemDosages.push_back(confChemDosages[i]);
			this->continuum.chemotherapy.chemProtocol.push_back(confChemProtocol[i]);
		}

		conf.lookupValue("hCoarse", this->continuum.hCoarse);
		conf.lookupValue("hRefined", this->continuum.hRefined);
		conf.lookupValue("deltaT", this->continuum.deltaT);
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Continuum Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getBloodVesselValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["blood-vessel"];

	try {
		conf.lookupValue("radius", this->bloodVessel.radius);
		conf.lookupValue("action-radius", this->bloodVessel.actionRadius);
		conf.lookupValue("o-production", this->bloodVessel.oProduction);
		this->bloodVessel.actionRadius *= this->bloodVessel.radius;
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Blood Vessel Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getAgentValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["agent"];
	
	try {
		conf.lookupValue("nucleus-radius", this->agent.nucleusRadius);
		conf.lookupValue("radius", this->agent.radius);
		conf.lookupValue("action-radius", this->agent.actionRadius);
		conf.lookupValue("o-consumption", this->agent.oConsumption);
		conf.lookupValue("egf-source", this->agent.egfSource);
		conf.lookupValue("chem-uptake", this->agent.chemUptake);
		conf.lookupValue("max-out-cells", this->agent.maxOutCells);
		this->agent.actionRadius *= this->agent.radius;
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Agent Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getForcesValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["forces"];
	
	try {
		conf.lookupValue("c_cca", this->forces.c_cca);
		conf.lookupValue("c_ccr", this->forces.c_ccr);
		conf.lookupValue("K", this->forces.K);
		conf.lookupValue("c_ct", this->forces.c_ct);
		conf.lookupValue("c_rct", this->forces.c_rct);
		conf.lookupValue("c_hap", this->forces.c_hap);
		this->forces.c_ct *= this->forces.K;
		this->forces.c_rct *= this->forces.K;
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Forces Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getParametersValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["parameters"];
	
	try {
		double a1, a2;
		conf.lookupValue("tauA", this->parameters.tauA);
		conf.lookupValue("tauP", this->parameters.tauP);
		conf.lookupValue("tauG1", this->parameters.tauG1);
		conf.lookupValue("tauC", this->parameters.tauC);
		conf.lookupValue("tauNL", this->parameters.tauNL);
		conf.lookupValue("fNS", this->parameters.fNS);
		conf.lookupValue("delta_tt", this->parameters.delta_tt);
		conf.lookupValue("alphaP", a1);
		conf.lookupValue("sigmaH", this->parameters.sigmaH);
		conf.lookupValue("alphaA", a2);
		conf.lookupValue("chem-threshold", this->parameters.chemThreshold);
		this->parameters.alphaP = 1/a1;
		this->parameters.alphaA = 1/a2;
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Parameters Settings Not Found!" << std::endl;
	}
}

void ConfigHandler::getPharmacologicValues() {
	const Setting& conf = this->cfg.getRoot()["agents"]["pharmacologic"];
	
	try {
		conf["pharmacokinetic"].lookupValue("p1", this->pharmacologic.pk.p1);
		conf["pharmacokinetic"].lookupValue("p2", this->pharmacologic.pk.p2);
		conf["pharmacokinetic"].lookupValue("p3", this->pharmacologic.pk.p3);
		conf["pharmacokinetic"].lookupValue("tau", this->pharmacologic.pk.tau);
		conf["pharmacokinetic"].lookupValue("lambda_a", this->pharmacologic.pk.lambda_a);
		conf["pharmacokinetic"].lookupValue("lambda_d", this->pharmacologic.pk.lambda_d);
		conf["pharmacokinetic"].lookupValue("zeta", this->pharmacologic.pk.zeta);

		conf["pharmacodynamic"].lookupValue("lower-threshold", this->pharmacologic.pd.lowerThreshold);
		conf["pharmacodynamic"].lookupValue("upper-threshold", this->pharmacologic.pd.upperThreshold);
		conf["pharmacodynamic"].lookupValue("max-concentration", this->pharmacologic.pd.maxConcentration);
		conf["pharmacodynamic"].lookupValue("q_u", this->pharmacologic.pd.q_u);
	} catch(const SettingNotFoundException &nfex) {
		this->created_ = false;
		std::cerr << "Pharmacologic Settings Not Found!" << std::endl;
	}
}

ConfigHandler::ConfigHandler(std::string configFile, int seed) {
	this->created_ = false;
	
	try {
		this->cfg.readFile(configFile.c_str());
		this->created_ = true;
	} catch(const FileIOException &fioex) {
		std::cerr << "I/O error while reading file." << std::endl;
	} catch(const ParseException &pex) {
		std::cerr   << "Parse error at " << pex.getFile() << ":" << pex.getLine()
		<< " - " << pex.getError() << std::endl;
	}

	this->input.seed = seed;

	this->getInputValues();
	this->getOutputValues();
	this->getContinuumValues();
	this->getBloodVesselValues();
	this->getAgentValues();
	this->getForcesValues();
	this->getParametersValues();
	this->getPharmacologicValues();
}

bool ConfigHandler::created() {
	return this->created_;
}