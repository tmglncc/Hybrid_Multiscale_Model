#ifndef MGMRES_H
#define MGMRES_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>

class MGMRES {
public:
	static void atx_cr(int n, int nz_num, double a[], int ia[], int ja[], double x[], double w[]);
	static void atx_st(int n, int nz_num, double a[], int ia[], int ja[], double x[], double w[]);
	static void ax_cr(int n, int nz_num, double a[], int ia[], int ja[], double x[], double w[]);
	static void ax_st(int n, int nz_num, double a[], int ia[], int ja[], double x[], double w[]);
	static void diagonal_pointer_cr(int n, int nz_num, int ia[], int ja[], int ua[]);
	static void ilu_cr(int n, int nz_num, int ia[], int ja[], double a[], int ua[], double l[]);
	static void lus_cr(int n, int nz_num, int ia[], int ja[], double l[], int ua[], double r[], double z[]);
	static void mgmres_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel);
	static void mult_givens(double c, double s, int k, double g[]);
	static void pmgmres_ilu_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel);
	static double r8vec_dot(int n, double a1[], double a2[]);
	static double* r8vec_uniform_01(int n, int *seed);
	static void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[]);
	static void timestamp();
};

#endif /* end of include guard: MGMRES */
