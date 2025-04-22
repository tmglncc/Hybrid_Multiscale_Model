#include "Frame.hpp"

void Frame::countCellStates() {
	this->necrotic = 0;
	this->quiescent = 0;
	this->proliferative2 = 0;
	this->hypoxic = 0;
	this->apoptotic = 0;
	this->proliferative1 = 0;
	this->migrate = 0;
	this->killed = 0;

	for (int i = 0; i < this->cells.size(); i++) {
		switch (this->cells[i].state) {
			case NECROTIC:
			this->necrotic++;
			break;

			case QUIESCENT:
			this->quiescent++;
			break;

			case PROLIFERATIVE2:
			this->proliferative2++;
			break;

			case HYPOXIC:
			this->hypoxic++;
			break;

			case APOPTOTIC:
			this->apoptotic++;
			break;

			case PROLIFERATIVE1:
			this->proliferative1++;
			break;

			case MIGRATE:
			this->migrate++;
			break;

			case KILLED:
			this->killed++;
			break;
		}
	}
}

Frame::Frame(int time, int tumorCells, int outCells, Vector3 domain, std::vector<Cell> cells): 
time(time), tumorCells(tumorCells), outCells(outCells), domain(domain), cells(cells) {}

std::string Frame::to_string() {
	this->countCellStates();

	std::string out = "============================================\n";
	out += "Time = " + std::to_string(this->time) + " h\n";
	out += "Tumor cells = " + std::to_string(this->tumorCells) + "\n";

	if (this->necrotic != 0)
		out += "- Necrotic cells = " + std::to_string(this->necrotic) + "\n";
	if (this->quiescent != 0)
		out += "- Quiescent (G0) cells = " + std::to_string(this->quiescent) + "\n";
	if (this->proliferative2 != 0)
		out += "- Proliferative (S-G2-M) cells = " + std::to_string(this->proliferative2) + "\n";
	if (this->hypoxic != 0)
		out += "- Hypoxic cells = " + std::to_string(this->hypoxic) + "\n";
	if (this->apoptotic != 0)
		out += "- Apoptotic cells = " + std::to_string(this->apoptotic) + "\n";
	if (this->proliferative1 != 0)
		out += "- Proliferative (G1) cells = " + std::to_string(this->proliferative1) + "\n";
	if (this->migrate != 0)
		out += "- Migrate cells = " + std::to_string(this->migrate) + "\n";
	if (this->killed != 0)
		out += "- Killed cells = " + std::to_string(this->killed) + "\n";

	out += "Normal cells = " + std::to_string(this->cells.size()-this->tumorCells) + "\n";
	out += "Total agents = " + std::to_string(this->cells.size()) + "\n";
	out += "Total cells left = " + std::to_string(this->outCells) + "\n";
	out += "============================================\n";

	return out;
}