#include "Chemotherapy.hpp"

Chemotherapy::Chemotherapy(std::vector<int> chemTimes, std::vector<int> chemDosages, std::vector<std::string> chemProtocol) {
	this->chemTimes = chemTimes;
	this->chemDosages = chemDosages;
	this->chemProtocol = chemProtocol;

	std::set<std::string> chemDrugsSet(this->chemProtocol.begin(), this->chemProtocol.end());
	this->chemDrugs = std::vector<std::string>(chemDrugsSet.begin(), chemDrugsSet.end());
	this->nextInjection = 0;
}

bool Chemotherapy::applyInjection(int frameTime) {
	if (this->nextInjection < this->chemProtocol.size() &&
		frameTime == this->chemTimes[this->nextInjection]) {
		this->nextInjection++;
		return true;
	}

	return false;
}
