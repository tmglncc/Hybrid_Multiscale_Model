#ifndef CHEMOTHERAPY
#define CHEMOTHERAPY

#include <iostream>
#include <string>
#include <vector>
#include <set>

class Chemotherapy {
private:
	std::vector<int> chemTimes;
	std::vector<int> chemDosages;
	std::vector<std::string> chemProtocol;

public:
	std::vector<std::string> chemDrugs;
	int nextInjection;

	Chemotherapy(std::vector<int> chemTimes, std::vector<int> chemDosages, std::vector<std::string> chemProtocol);
	bool applyInjection(int frameTime);
};

#endif