
#include <mex.h>
#include <mat.h>
#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <vector> /* For STL */
#include <complex>
#include "runge_kutta.h"

#define BUFSIZE 256

void testSolver(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void readY(std::vector< std::complex<double> >&  y)
{
	MATFile *pmat;
	mxArray * pa1[1];
	const char *file = "/Users/tanay/Dropbox/dev/cftd_mex_test/y1_fp.mat";

	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
	}

	pa1[0] = matGetVariable(pmat, "y");
	if (pa1[0] == NULL) {
		printf("Error reading existing matrix noise_vec\n");
	}

	const size_t N = mxGetDimensions(pa1[0])[0];

	double * re = (double *)mxGetPr(pa1[0]);
	double * im = (double *)mxGetPi(pa1[0]);

	for (size_t i = 0; i < N; i++) {
		y.push_back(std::complex<double>(re[i], im[i]));
	}

	mxDestroyArray(pa1[0]);

	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n", file);
	}
}

int main() {
	MATFile *pmat;
	mxArray * pa1[7];
	mxArray * pa2[2];
	std::vector<int> myInts;
	myInts.push_back(1);
	myInts.push_back(2);
	printf("Accessing a STL vector: %d\n", myInts[1]);

	const char *file = "/Users/tanay/Dropbox/dev/cftd_mex_test/sd_data.mat";

	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
		return(EXIT_FAILURE);
	}

	//save('sd_data.out', 'S_coredata', 'pump_pwr', 't_initial', 't_final', 'noise_vec', 'basis_type');

	pa1[0] = matGetVariable(pmat, "S_coredata");
	if (pa1[0] == NULL) {
		printf("Error reading existing matrix S_coredata\n");
		return(EXIT_FAILURE);
	}

	pa1[1] = matGetVariable(pmat, "pump_pwr");
	if (pa1[1] == NULL) {
		printf("Error reading existing matrix pump_pwr\n");
		return(EXIT_FAILURE);
	}

	pa1[2] = matGetVariable(pmat, "t_initial");
	if (pa1[2] == NULL) {
		printf("Error reading existing matrix t_initial\n");
		return(EXIT_FAILURE);
	}

	pa1[3] = matGetVariable(pmat, "t_final");
	if (pa1[3] == NULL) {
		printf("Error reading existing matrix t_final\n");
		return(EXIT_FAILURE);
	}

	pa1[4] = matGetVariable(pmat, "noise_vec");
	if (pa1[4] == NULL) {
		printf("Error reading existing matrix noise_vec\n");
		return(EXIT_FAILURE);
	}

	pa1[5] = matGetVariable(pmat, "basis_type");
	if (pa1[5] == NULL) {
		printf("Error reading existing matrix basis_type\n");
		return(EXIT_FAILURE);
	}

	pa1[6] = matGetVariable(pmat, "ode_options");
	if (pa1[6] == NULL) {
		printf("Error reading existing matrix basis_type\n");
		return(EXIT_FAILURE);
	}

	const mxArray *const_pa[7] = { pa1[0], pa1[1], pa1[2], pa1[3], pa1[4], pa1[5], pa1[6] };
	//testSolver(2, pa2, 6, const_pa);
	mexFunction(2, pa2, 7, const_pa);

	/* clean up before exit */
	for (int i = 0; i < 7; i++)
		mxDestroyArray(pa1[i]);


	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n", file);
		return(EXIT_FAILURE);
	}

	printf("Done\n");
	return(EXIT_SUCCESS);
}

void testSolver(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{

	// Initialize logger
	FILELog::ReportingLevel() = logDEBUG;

	Parameters params;
	State state;

	FILE_LOG(logINFO) << "beginning of mexfunction";

	params.g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
	params.g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
	params.k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
	char * temp_basis_type = mxArrayToString(prhs[5]);
	std::string basis_type = std::string(temp_basis_type);
	mxFree(temp_basis_type);
	params.N = mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
	double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
	size_t dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
	double * ptB = mxGetPr(mxGetField(prhs[0], 0, "A"));
	size_t dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
	double * CFvalsReal = (double *)mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
	double * CFvalsImag = (double *)mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
	params.CFvals.resize(params.N);

	for (int i = 0; i < params.N; i++) {
		params.CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
	}

	params.A.resize(dimA, dimA);

	for (int i = 0; i < dimA; i++) {
		for (int j = 0; j < dimA; j++) {
			params.A(i, j) = ptA[i + j * dimA];
		}
	}

	params.B.resize(dimB, dimB);

	for (int i = 0; i < dimB; i++) {
		for (int j = 0; j < dimB; j++) {
			params.B(i, j) = ptB[i + j * dimB];
		}
	}

	auto pump_pt = mxGetPr(prhs[1]);
	params.pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

	for (int i = 0; i < params.pump_pwr.size(); i++)
		params.pump_pwr[i] = pump_pt[i];

	double re = *mxGetPr(mxGetField(prhs[0], 0, "n"));
	double im = *mxGetPi(mxGetField(prhs[0], 0, "n"));
	params.n = std::complex<double>(re, im);
	auto noise_pt = mxGetPr(prhs[4]);

	size_t outvec_size = params.N * (params.N + 2);
	params.noise_vec.resize(outvec_size);

	for (int i = 0; i < params.noise_vec.size(); i++)
		params.noise_vec[i] = noise_pt[i];

	// test for one iteration
	std::vector< std::complex<double> >  y, dy;
	readY(y);
	dy.resize(y.size());
	TDSSolvers_FP fp(params);
	fp(y, dy, 0);

	printf("Done\n");
}
