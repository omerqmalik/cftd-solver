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
namespace bnu = boost::numeric::ublas;

// std vector is used for making it work with steppers
typedef std::vector<std::complex<double>> state_type;

// ublas vector is used for efficient matrix algebra
typedef bnu::vector<std::complex<double>> ustate_type;

const std::complex<double> halfI = std::complex<double>(0, 0.5);
const std::complex<double> oneI = std::complex<double>(0, 1.0);

enum solver_type { FP, RING, UCF };

// Parameters passed 
struct Parameters {

public:
	size_t N = 0;
	double g_per = 0.0;
	double g_par = 0.0;
	double k_a = 0.0;

	std::complex<double> n;
	bnu::vector <double> pump_pwr;
	bnu::matrix <double> A;
	bnu::matrix <double> B;
	state_type noise_vec;
	ustate_type CFvals;

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
        FILE_LOG(logDEBUG) << "Initializing Observer";
    }

    void operator()(const state_type &x, const double t)
    {
        FILE_LOG(logDEBUG4) << "begin of observer call";
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
	ustate_type a_vec;
	ustate_type b_vec1;
	ustate_type D_vec;
	ustate_type aD_vec;
	ustate_type Q;
	ustate_type NL_term;

	Parameters& params;
	
public:
    
    TDSSolvers_RING(Parameters& params_) : params(params_),
		a_vec(params_.N), b_vec1(params_.N), D_vec(params_.N*params_.N),
		aD_vec(params_.N), Q(params_.N*params_.N), NL_term(params_.N*params_.N)
    {
        FILE_LOG(logDEBUG) << "Init RING";
    }

    void operator()(const state_type& y, state_type& dy, const double t)
    {		
		int N = params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());

		bnu::matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                size_t a = j / N;
                size_t b = j % N;
                NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        ustate_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(i * N + j) * a_vec(j);
            }
        }

        bnu::scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
		ustate_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vec1 / (params.n * params.n));
		ustate_type dy2 = -params.g_per * b_vec1 - oneI * params.g_per * a_vecTD_mat;
		ustate_type dy3 = -params.g_par * (D_vec - params.pump_pwr) + halfI * params.g_par * NL_term;

		state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

// ODE System Solver
struct TDSSolvers_UCF {

private:
	ustate_type a_vec;
	ustate_type b_vec1;
	ustate_type D_vec;
	ustate_type aD_vec;
	ustate_type Q;
	ustate_type NL_term;

	Parameters& params;

public:

	TDSSolvers_UCF(Parameters& params_) : params(params_),
		a_vec(params_.N), b_vec1(params_.N), D_vec(params_.N*params_.N),
		aD_vec(params_.N), Q(params_.N*params_.N), NL_term(params_.N*params_.N)
	{
        FILE_LOG(logDEBUG) << "Init UCF";
    }

    void operator()(const state_type& y, state_type& dy, const double t)
    {
		int N = params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());
		bnu::matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                size_t a = j / N;
				size_t b = j % N;
                NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

		ustate_type a_vecT_D_mat(N);
		ustate_type b_vecT_B(N);
		a_vecT_D_mat.clear();
		b_vecT_B.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecT_D_mat[i] += D_vec(i * N + j) * a_vec(j);
                b_vecT_B[i] += params.B(j, i) * b_vec1(j);
            }
        }

        bnu::scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
		ustate_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vecT_B);
		ustate_type dy2 = -params.g_per * b_vec1 - oneI * params.g_per * a_vecT_D_mat;
		ustate_type dy3 = -params.g_par * (D_vec - params.pump_pwr) + oneI * params.g_par * NL_term;

		state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

// ODE System Solver
struct TDSSolvers_FP {

private:
	ustate_type a_vec;
	ustate_type b_vec1;
	ustate_type D_vec;
	ustate_type aD_vec;
	ustate_type Q;
	ustate_type NL_term;

	Parameters& params;

public:

	TDSSolvers_FP(Parameters& params_) : params(params_), 
		a_vec(params_.N), b_vec1(params_.N), D_vec(params_.N*params_.N), 
		aD_vec(params_.N), Q(params_.N*params_.N), NL_term(params_.N*params_.N)
	{
        FILE_LOG(logDEBUG) << "Init FP";
    }
    
    void operator()(const state_type& y, state_type& dy, const double t)
    {
		int N = params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());

		bnu::matrix<std::complex<double>> Q = imag(outer_prod(conj(b_vec1), a_vec));

		NL_term.clear();
        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
				size_t a = j / N;
				size_t b = j % N ;
                NL_term[i] += params.A(i, j) * Q(a, b);
            }
        }

        ustate_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(i * N + j) * a_vec(j);
            }
        }

		bnu::scalar_vector<std::complex<double>> k_avec(a_vec.size(), params.k_a);
        ustate_type dy1 = element_prod((k_avec - (element_prod(params.CFvals, params.CFvals) / params.k_a)), (halfI * a_vec)) + (halfI * params.k_a * b_vec1 / (params.n * params.n));
        ustate_type dy2 = -params.g_per * b_vec1 - oneI * params.g_per * a_vecTD_mat;
        ustate_type dy3 = -params.g_par * (D_vec - params.pump_pwr) - params.g_par * NL_term;

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

        //ustate_type & dy_a = dy;
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy.begin());
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
	params.N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    size_t dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    double * ptB = mxGetPr(mxGetField(prhs[0], 0, "A"));
    size_t dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
	params.CFvals.resize(params.N);

    for(int i = 0; i < params.N; i++) {
		params.CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
    }

	params.A.resize(dimA, dimA);

    for(int i = 0; i < dimA; i++) {
        for(int j = 0; j < dimA; j++) {
			params.A(i, j) = ptA[i + j * dimA];
        }
    }

	params.B.resize(dimB, dimB);

    for(int i = 0; i < dimB; i++) {
        for(int j = 0; j < dimB; j++) {
			params.B(i, j) = ptB[i + j * dimB];
        }
    }

    auto pump_pt = mxGetPr(prhs[1]);
	params.pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

    for(int i = 0; i < params.pump_pwr.size(); i++)
		params.pump_pwr[i] = pump_pt[i];

    double re = *mxGetPr(mxGetField(prhs[0], 0, "n"));
	double im = *mxGetPi(mxGetField(prhs[0], 0, "n"));
	params.n = std::complex<double>(re, im);
    auto noise_pt = mxGetPr(prhs[4]);

	size_t outvec_size = params.N * (params.N + 2);
	params.noise_vec.resize(outvec_size);

    for(int i = 0; i < params.noise_vec.size(); i++)
		params.noise_vec[i] = noise_pt[i];

    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);
    
    FILE_LOG(logINFO) << "Before Integrate";

	//typedef runge_kutta_dopri5 <state_type> stepper_type;
	//typedef dense_output_runge_kutta<controlled_runge_kutta< runge_kutta_dopri5<state_type> > > stepper_type;

	auto rkd = runge_kutta_dopri5<state_type>{};
	auto stepper = make_dense_output(1.0e-6, 1.0e-3, rkd);

	double dt = 0.1;

	//std::vector<double> times(91);
	//for (size_t i = 0; i<times.size(); ++i)
	//	times[i] = 1.0 + dt*i;

	Observer observer(state);
    if(basis_type.compare("RING") == 0) {
		TDSSolvers_RING ring(params);
        FILE_LOG(logDEBUG4) << "Before calling RING integrator";
		integrate_const(stepper,
			ring, params.noise_vec,
			*t_initial_pt, *t_final_pt, dt, observer);
        //integrate_times(stepper, ring, params.noise_vec,
		//	boost::begin(times), boost::end(times), dt, observer);
    }

    else if(basis_type.compare("UCF") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling UCF integrator";
		TDSSolvers_UCF ucf(params);
        integrate_const(stepper,
                        ucf, params.noise_vec,
                        *t_initial_pt, *t_final_pt, dt, observer);
    }

    else if(basis_type.compare("FP") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling FP integrator";
		TDSSolvers_FP fp(params);
        integrate_const(stepper,
						fp, params.noise_vec,
                        *t_initial_pt, *t_final_pt, dt, observer);
    }

    else {
        exit(EXIT_FAILURE);
    }

    FILE_LOG(logINFO) << "After Integrate";
    FILE_LOG(logDEBUG4) << "before assigning output";
    FILE_LOG(logDEBUG4) << state.num;

    int index = 0;
    double * outT = (double *) mxMalloc(state.num * sizeof(double));

    FILE_LOG(logDEBUG4) << outT;
    FILE_LOG(logDEBUG4) << state.tout.size();
    FILE_LOG(logDEBUG4) << "after malloc";

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
    
    FILE_LOG(logINFO) << "End of mexfunction";
    
    // Replace redirection
    std::cout.rdbuf(outbuf);
}
