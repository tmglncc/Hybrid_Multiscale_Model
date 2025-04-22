#include "CRS.hpp"

using namespace std;

/*void CRS::create_crs_ordened(int row, int column, int row_ptr[], int col_ind[], double val[], double value) {
	if (row_ptr[row+1] == 0)
		row_ptr[row+1] = row_ptr[row];

	row_ptr[row+1]++;
	col_ind[row_ptr[row+1]-1] = column;
	val[row_ptr[row+1]-1] = value;
}

void CRS::add_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value) {
	for(int i = row_ptr[row]; i < row_ptr[row+1]; i++) {
		if (column == col_ind[i])
			val[i] += new_value;
	}
}

void CRS::set_crs(int row, int column, int row_ptr[], int col_ind[], double val[], double new_value) {
	for(int i = row_ptr[row]; i < row_ptr[row+1]; i++) {
		if (column == col_ind[i])
			val[i] = new_value;
	}
}

void CRS::print_crs(int row_ptr[], int col_ind[], double val[], int num_nos, int nz_num) {
	for (int i = 0; i < num_nos + 1; i++)
		cout << row_ptr[i] << "\t";
	cout << endl;

	for (int i = 0; i < nz_num; i++)
		cout << col_ind[i] << "\t";
	cout << endl;

	for (int i = 0; i < nz_num; i++)
		cout << val[i] << "\t";
	cout << endl;
}

void CRS::zera_crs(int num_nos, int nz_num, int row_ptr[], int col_ind[], double val[]) {
	for (int i = 0; i < num_nos + 1; i++)
		row_ptr[i] = 0;

	for (int i = 0; i < nz_num; i++) {
		col_ind[i] = 0;
		val[i] = 0;
	}
}*/

CRS::CRS(int matrixSize, int nzNum) {
	this->matrixSize = matrixSize;
	this->nzNum = nzNum;

	this->rowPtr = new int[this->matrixSize + 1];
	this->colInd = new int[this->nzNum];
	this->val = new double[this->nzNum];

	this->reset();
}

CRS::~CRS() {
	delete this->rowPtr;
	delete this->colInd;
	delete this->val;
}

void CRS::createOrdered(int row, int column, double value) {
	if (this->rowPtr[row+1] == 0)
		this->rowPtr[row+1] = this->rowPtr[row];

	this->rowPtr[row+1]++;
	this->colInd[this->rowPtr[row+1]-1] = column;
	this->val[this->rowPtr[row+1]-1] = value;
}

void CRS::set(int row, int column, double newValue) {
	for (int i = this->rowPtr[row]; i < this->rowPtr[row+1]; i++) {
		if (column == this->colInd[i])
			this->val[i] = newValue;
	}
}

double CRS::get(int row, int column) {
	for (int i = this->rowPtr[row]; i < this->rowPtr[row+1]; i++) {
		if (column == this->colInd[i])
			return this->val[i];
	}

	return 0.0;
}

void CRS::add(int row, int column, double increment) {
	for (int i = this->rowPtr[row]; i < this->rowPtr[row+1]; i++) {
		if (column == this->colInd[i])
			this->val[i] += increment;
	}
}

void CRS::reset() {
	for (int i = 0; i < this->matrixSize + 1; i++)
		this->rowPtr[i] = 0;

	for (int i = 0; i < this->nzNum; i++) {
		this->colInd[i] = 0;
		this->val[i] = 0.0;
	}
}

void CRS::print() {
	for (int i = 0; i < this->matrixSize + 1; i++)
		cout << this->rowPtr[i] << "\t";
	cout << endl;

	for (int i = 0; i < this->nzNum; i++)
		cout << this->colInd[i] << "\t";
	cout << endl;

	for (int i = 0; i < this->nzNum; i++)
		cout << this->val[i] << "\t";
	cout << endl;
}

void CRS::printAsMatrix() {
	for (int i = 0; i < this->matrixSize; i++) {
		for (int j = 0; j < this->matrixSize; j++)
			cout << this->get(i, j) << "\t";
		cout << endl;
	}
}
