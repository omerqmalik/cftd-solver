#include <mex.h>
#include <matrix.h>
#include <complex>
#include <iostream>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/numeric/odeint.hpp>
#include "output.h"
#include "log.h"
#include "runge_kutta.h"

// Entry point for MATLAB
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	// Redirect cout to MEX printf
	mstream mout;
	std::streambuf *outbuf = std::cout.rdbuf(&mout);

	// Initialize logger
	FILELog::ReportingLevel() = logDEBUG;

    //Set up holders for input data and transitionary states
	Parameters params;
	State state;

	FILE_LOG(logINFO) << "beginning of mexfunction";
    
    //Get g_per from S_coredata
    //g_per is a value that is used to calculate the next state
	params.g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    
    
    
    //Get g_par from S_coredata
    //g_par is a value that is used to calculate the next state
	params.g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    
    
    
    //Get k_a from S_coredata
    //k_a is a value that is used to calculate the next state
	params.k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    
    
    
    //Get basis_type from inputs
    //basis_type is a string that represents the type of Solver to be used
    //It is either FP, UCF, or RING
    
    //Get temporary char * pointer using mxArrayToString
	char * temp_basis_type = mxArrayToString(prhs[5]);
    
	std::string basis_type = std::string(temp_basis_type);
    
    //free the temporary pointer
	mxFree(temp_basis_type);
    
    
    
    
    //get N from S_coredata
    //N is the size of CFvals
	params.N = mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    
    
    
    //get A from S_coredata
    //A is a matrix used to calculate the next state
    //A should be N^2 by N^2
    
    //get a pointer to the beginning of A
	double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    
    //get the length of A
	size_t dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    
    //set the size of A, guaranteed to be square
    params.A.resize(dimA, dimA);
    
    //iterate from pointer to fill array
	for (int i = 0; i < dimA; i++) {
		for (int j = 0; j < dimA; j++) {
			params.A(i, j) = ptA[i + j * dimA];
		}
	}
    
    
    
	//get the matrix B from S_coredata
    //B is a matrix used to calculate the next state
    //B is only used in the UCF Solver
    //B should be N by N
    
    //get a pointer to the beginning of B
	double * ptB = mxGetPr(mxGetField(prhs[0], 0, "B"));
    
    //get the length of B
	size_t dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
    
    //set the size of B, guaranteed to be square
    params.B.resize(dimB, dimB);
    
    //iterate from pointer to fill array
	for (int i = 0; i < dimB; i++) {
		for (int j = 0; j < dimB; j++) {
			params.B(i, j) = ptB[i + j * dimB];
		}
	}
    
    
    
    //get CFvals from S_coredata
    //CFvals is a complex double 
    
    //get both a real and imaginary pointer to CFvals
	double * CFvalsReal = (double *)mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
	double * CFvalsImag = (double *)mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
    
    //Resize CFvals
	params.CFvals.resize(params.N);
    
    //fill CFvals with complex numbers from the real and imaginary pointers
	for (int i = 0; i < params.N; i++) {
		params.CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
	}
    
    
    
    
    //get the pump power from input
    
    //get a pointer to pump_pwr
	auto pump_pt = mxGetPr(prhs[1]);
    
    //resize pump_pwr
	params.pump_pwr.resize(mxGetDimensions(prhs[1])[0]);
    
    //iterate from pointer to fill pump_pwr
	for (int i = 0; i < params.pump_pwr.size(); i++)
		params.pump_pwr[i] = pump_pt[i];
    
    
    
    
    //get n (little n) from S_coredata
    
    //get the real and imaginary part of n
	double re = *mxGetPr(mxGetField(prhs[0], 0, "n"));
	double im = *mxGetPi(mxGetField(prhs[0], 0, "n"));
    
    //combine the real and imaginary parts to form n
	params.n = std::complex<double>(re, im);
    
    
    
    
    //get noise_vec from input
    
    //get a pointer to noise_vec
	auto noise_pt = mxGetPr(prhs[4]);
    
	size_t outvec_size = params.N * (params.N + 2);
    
	params.noise_vec.resize(outvec_size);

	for (int i = 0; i < params.noise_vec.size(); i++)
		params.noise_vec[i] = noise_pt[i];

	auto t_initial_pt = mxGetPr(prhs[2]);
	auto t_final_pt = mxGetPr(prhs[3]);

	// Read options structure
	mxArray* mxA = mxGetField(prhs[6], 0, "odeint_const");
	bool ode_const = *mxGetLogicals(mxA);
	double dt = *mxGetPr(mxGetField(prhs[6], 0, "dt"));
	double abs = *mxGetPr(mxGetField(prhs[6], 0, "abs_error"));
	double rel = *mxGetPr(mxGetField(prhs[6], 0, "rel_error"));

	char * plog_level = mxArrayToString(mxGetField(prhs[6], 0, "log_level"));
	FILELog::ReportingLevel() = getLogLevel(plog_level);
	mxFree(plog_level);

	FILE_LOG(logINFO) << "Before Integrate";

	Observer observer(state);
	std::unique_ptr<TDSSolver> solver = solverFactory(basis_type, params);
	OdeIntegrator integrator(ode_const, abs, rel);
	integrator.integrate(*solver, params.noise_vec, *t_initial_pt, *t_final_pt, dt, boost::ref(observer));

	FILE_LOG(logINFO) << "After Integrate";
	FILE_LOG(logDEBUG4) << "before assigning output";
	FILE_LOG(logDEBUG4) << state.num;

	int index = 0;
	double * outT = (double *)mxMalloc(state.num * sizeof(double));

	FILE_LOG(logDEBUG4) << outT;
	FILE_LOG(logDEBUG4) << state.tout.size();
	FILE_LOG(logDEBUG4) << "after malloc";

	double * outYr = (double *)mxCalloc(state.num * outvec_size, sizeof(double));
	double * outYi = (double *)mxCalloc(state.num * outvec_size, sizeof(double));

	for (index = 0; index < state.num; index++) {
		FILE_LOG(logDEBUG4) << "in loop " << index << std::endl;
		outT[index] = state.tout[index];

		for (int j = 0; j < outvec_size; j++) {
			if (j > outvec_size - 2 && index > state.num - 2) {
				FILE_LOG(logDEBUG) << "Values Re: " <<
					state.yout[index][j].real() << " Im: " <<
					state.yout[index][j].imag();
			}

			outYr[index * outvec_size + j] = state.yout[index][j].real();
			outYi[index * outvec_size + j] = state.yout[index][j].imag();
		}
	}

	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetPr(plhs[0], outT);
	mxSetM(plhs[0], state.num);
	mxSetN(plhs[0], 1);

	FILE_LOG(logDEBUG4) << "in between assigning output";

	plhs[1] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
	mxSetPr(plhs[1], outYr);
	mxSetPi(plhs[1], outYi);

	// Needs a transpose in Matlab!
	mxSetN(plhs[1], state.num);
	mxSetM(plhs[1], outvec_size);

	FILE_LOG(logINFO) << "End of mexfunction";

	// Replace redirection
	std::cout.rdbuf(outbuf);
}
