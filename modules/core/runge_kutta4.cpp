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

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;
using namespace boost::assign;

typedef vector<std::complex<double>> state_type;
enum solver_type {FP, RING, UCF};
const std::complex<double> halfI = std::complex<double>(0, 0.5);

// Parameters passed 
struct Parameters {
    
public:
	size_t N = 0;
	double g_per = 0.0;
	double g_par = 0.0;
	double k_a = 0.0;

	std::complex<double> n;
	vector <double> pump_pwr;
	matrix <double> A;
	matrix <double> B;
	state_type noise_vec;
	state_type CFvals;

private:

};

// Output holder
struct State {
	int num = 0;
	std::vector<double> tout;
	std::vector<state_type> yout;
};

// Collects intermediate state
struct Observer {

private:
	State& state;

public:
    Observer(State& state_) : state(state_)
    {
        FILE_LOG(logINFO) << "Initializing Observer";
    }

    void operator()(const state_type &x, const double t)
    {
        FILE_LOG(logDEBUG4) << "begin of operator";
        state.tout.push_back(t);
		state.yout.push_back(x);

        //for(int j = 0; j < N * (N + 2); j++) {
        //    if (j > N*N) {
        //     FILE_LOG(logDEBUG) << "Observer Vals Re: " << 
        //             x[j].real() << " Im: " << x[j].imag();
        //    }
        //    yout_real.push_back(x[j].real());
        //    yout_imag.push_back(x[j].imag());
        //}       
        
        state.num++;
    }
};

// ODE system solver
struct TDSSolvers_RING {
private:
	state_type a_vec;
	state_type b_vec1;
	state_type D_vec;
	state_type aD_vec;
	state_type Q;
	state_type NL_term;

	Parameters& params;
	
	void resize(size_t N) 
	{
		a_vec.resize(N);
		b_vec1.resize(N);
		D_vec.resize(N * N);
		aD_vec.resize(N);
		Q.resize(N * N);
		NL_term.resize(N * N);
	}

public:
    
    TDSSolvers_RING(Parameters& params_) : params(params_)
    {
        FILE_LOG(logDEBUG) << "Init RING";
		resize(params.N);
    }

    void operator()(const state_type& y, state_type& dy, const double t)
    {
		size_t N = params.N;
        //a_vec  = y(1:N);
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N), y.size());
        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                size_t a = j / N;
                size_t b = j % N;
                NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(j * N + i) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
        state_type & dy_a = dy;
        state_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vec1 / (params.n * params.n));
        state_type dy2 = -params.g_per * b_vec1 - halfI * 2 * params.g_per * a_vecTD_mat;
        state_type dy3 = -params.g_par * (D_vec - params.pump_pwr) + halfI * params.g_par * NL_term;

        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

// ODE System Solver
struct TDSSolvers_UCF {

private:
	state_type a_vec;
	state_type b_vec1;
	state_type D_vec;
	state_type aD_vec;
	state_type Q;
	state_type NL_term;

	Parameters& params;

	void resize(size_t N)
	{
		a_vec.resize(N);
		b_vec1.resize(N);
		D_vec.resize(N * N);
		aD_vec.resize(N);
		Q.resize(N * N);
		NL_term.resize(N * N);
	}

public:

	TDSSolvers_UCF(Parameters& params_) : params(params_)
	{
        FILE_LOG(logDEBUG) << "Init UCF";
		resize(params.N);
    }

    void operator()(const state_type& y, state_type& dy, const double t)
    {
		size_t N = params.N;
        FILE_LOG(logDEBUG4) << "begin UCF";
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N), y.size());
        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                size_t a = j / N;
				size_t b = j % N;
                NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecT_D_mat(N);
        state_type b_vecT_B(N);
		a_vecT_D_mat.clear();
		b_vecT_B.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecT_D_mat[i] += D_vec(j * N + i) * a_vec(i);
                b_vecT_B[i] += params.B(i, j) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
        state_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vecT_B);
        state_type dy2 = -params.g_per * b_vec1 - halfI * 2 * params.g_per * a_vecT_D_mat;
        state_type dy3 = -params.g_par * (D_vec - params.pump_pwr) + 2 * halfI * params.g_par * NL_term;

        state_type & dy_a = dy;
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

// ODE System Solver
struct TDSSolvers_FP {

private:
	state_type a_vec;
	state_type b_vec1;
	state_type D_vec;
	state_type aD_vec;
	state_type Q;
	state_type NL_term;

	Parameters& params;

	void resize(size_t N)
	{
		a_vec.resize(N);
		b_vec1.resize(N);
		D_vec.resize(N * N);
		aD_vec.resize(N);
		Q.resize(N * N);
		NL_term.resize(N * N);
	}

public:

	TDSSolvers_FP(Parameters& params_) : params(params_)
	{
        FILE_LOG(logDEBUG) << "Init FP";
		resize(params.N);
    }
    
    void operator()(const state_type& y, state_type& dy, const double t)
    {
		size_t N = params.N;
        FILE_LOG(logDEBUG4) << "begin FP";
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, 2 * N , y.size());
        matrix<std::complex<double>> Q = imag(outer_prod(conj(b_vec1), a_vec));

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
				size_t a = j / N;
				size_t b = j % N ;
                NL_term[i] += params.A(i, j) * Q(a, b);
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(j * N + i) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
        state_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vec1 / (params.n * params.n));
        state_type dy2 = -params.g_per * b_vec1 - halfI * 2 * params.g_per * a_vecTD_mat;
        state_type dy3 = -params.g_par * (D_vec - params.pump_pwr) - params.g_par * NL_term;

		//std::ostringstream oss;
		//for (int i = 0; i < N*(N + 2); i++) {
		//	if (i < N)
		//		oss << dy1(i);
		//	else if (i < 2 * N)
		//		oss << dy2(i - N);
		//	else
		//		oss << dy3(i - (2*N));
		//}
		//FILE_LOG(logDEBUG) << oss.str();

        state_type & dy_a = dy;
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);

		//std::ostringstream osr;
		//for (int i = 0; i < dy.size(); i++)
		//	osr << dy(i);
		//FILE_LOG(logDEBUG) << osr.str();
    }
};

// Entry point for MATLAB
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Redirect cout to MEX printf
    mstream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);

    // Initialize logger
    FILELog::ReportingLevel() = logINFO;

	Parameters params;
	State state;
    
    FILE_LOG(logINFO) << "beginning of mexfunction";
	params.g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.1";
	params.g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.2";
	params.k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.31";
    char * temp_basis_type = mxArrayToString(prhs[5]);
    FILE_LOG(logDEBUG4) << "checkpoint 0.32";
    std::string basis_type = std::string(temp_basis_type);
    FILE_LOG(logDEBUG4) << "checkpoint 0.33";
    mxFree(temp_basis_type);
    FILE_LOG(logDEBUG4) << "checkpoint 0.4";
	params.N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.5";
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.6";
    size_t dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.7";
    double * ptB = mxGetPr(mxGetField(prhs[0], 0, "A"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.8";
    size_t dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.9";
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.10";
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
    FILE_LOG(logDEBUG4) << "checkpoint 1";
	params.CFvals.resize(params.N);
    FILE_LOG(logDEBUG4) << "checkpoint 1.01";

    for(int i = 0; i < params.N; i++) {
		params.CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
    }

    FILE_LOG(logDEBUG4) << "checkpoint 1.1";
	params.A.resize(dimA, dimA);
    FILE_LOG(logDEBUG4) << "checkpoint 1.2";

    for(int i = 0; i < dimA; i++) {
        for(int j = 0; j < dimA; j++) {
			params.A(i, j) = ptA[i + j * dimA];
        }
    }

    FILE_LOG(logDEBUG4) << "checkpoint 1.3";
	params.B.resize(dimB, dimB);
    FILE_LOG(logDEBUG4) << "checkpoint 1.4";

    for(int i = 0; i < dimB; i++) {
        for(int j = 0; j < dimB; j++) {
			params.B(i, j) = ptB[i + j * dimB];
        }
    }

    FILE_LOG(logDEBUG4) << "checkpoint 2";
    auto pump_pt = mxGetPr(prhs[1]);
	params.pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

    for(int i = 0; i < params.pump_pwr.size(); i++)
		params.pump_pwr[i] = pump_pt[i];

    double re = *mxGetPr(mxGetField(prhs[0], 0, "n"));
	double im = *mxGetPi(mxGetField(prhs[0], 0, "n"));
	params.n = std::complex<double>(re, im);
    FILE_LOG(logDEBUG4) << "checkpoint 3";
    auto noise_pt = mxGetPr(prhs[4]);
    //noise_vec.resize(mxGetDimensions(prhs[4])[0]);

	size_t outvec_size = params.N * (params.N + 2);
	params.noise_vec.resize(outvec_size);

    for(int i = 0; i < params.noise_vec.size(); i++)
		params.noise_vec[i] = noise_pt[i];

    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);
    FILE_LOG(logDEBUG4) << "Resizing all structures";
    
    FILE_LOG(logINFO) << "Before Integrate";

	typedef runge_kutta4<vector<double>> rk4;
	//typedef runge_kutta_dopri5<vector<double>> stepper_type;
	typedef runge_kutta_dopri5 <state_type> stepper_type;
    
	Observer observer(state);
    if(basis_type.compare("RING") == 0) {
		TDSSolvers_RING ring(params);
        FILE_LOG(logDEBUG4) << "Before calling RING integrator";
        integrate_const(stepper_type(), ring, params.noise_vec,
                        *t_initial_pt, *t_final_pt, 0.061, observer);
    }

    else if(basis_type.compare("UCF") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling UCF integrator";
		TDSSolvers_UCF ucf(params);
        integrate_const(stepper_type(),
                        ucf, params.noise_vec,
                        *t_initial_pt, *t_final_pt, 0.061, observer);
    }

    else if(basis_type.compare("FP") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling FP integrator";
		TDSSolvers_FP fp(params);
        integrate_const(stepper_type(),
						fp, params.noise_vec,
                        *t_initial_pt, *t_final_pt, 0.061, observer);
    }

    else {
        exit(EXIT_FAILURE);
    }

    FILE_LOG(logINFO) << "After Integrate";
    FILE_LOG(logDEBUG) << "before assigning output";
    FILE_LOG(logDEBUG) << state.num;
    int index = 0;
    double * outT = (double *) mxMalloc(state.num * sizeof(double));
    FILE_LOG(logDEBUG4) << outT;
    FILE_LOG(logDEBUG4) << state.tout.size();
    FILE_LOG(logDEBUG4) << "after malloc";

    FILE_LOG(logDEBUG4) << "after first loop";
    double * outYr = (double *) mxCalloc(state.num * outvec_size, sizeof(double));
    double * outYi = (double *) mxCalloc(state.num * outvec_size, sizeof(double));
    
    for(index = 0; index < state.num; index++) {
        FILE_LOG(logDEBUG4) << "in loop " << index << std::endl;
        outT[index] = state.tout[index]; 
    
        for(int j = 0; j < outvec_size; j++) {
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
    
    FILE_LOG(logINFO) << "End of Mex Function";
    
    // Replace redirection
    std::cout.rdbuf(outbuf);
}
