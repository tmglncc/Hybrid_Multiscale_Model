#ifndef CRS_H
#define CRS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class CRS {
public:
	int matrixSize, nzNum;
	int* rowPtr;
	int* colInd;
	double* val;

	CRS(int matrixSize, int nzNum);
	~CRS();

	void createOrdered(int row, int column, double value);
	void set(int row, int column, double newValue);
	double get(int row, int column);
	void add(int row, int column, double increment);
	void reset();
	void print();
	void printAsMatrix();
};

#endif /* end of include guard: CRS */
