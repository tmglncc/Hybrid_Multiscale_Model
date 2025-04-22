#ifndef FRAME
#define FRAME

#include <vector>
#include "Cell.hpp"

class Frame {
private:
	int necrotic, quiescent, proliferative2, hypoxic, apoptotic, proliferative1, migrate, killed;

	void countCellStates();

public:
	int time, tumorCells, outCells;
	Vector3 domain;
	std::vector<Cell> cells;

	Frame(int time = 0, int tumorCells = 0, int outCells = 0, Vector3 domain = Vector3(), std::vector<Cell> cells = std::vector<Cell>());
	std::string to_string();
};

#endif /* end of include guard: FRAME */