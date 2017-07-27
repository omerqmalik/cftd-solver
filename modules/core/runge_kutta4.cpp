#include <mex.h>
#include <complex>
#include <iostream>
#include <matrix.h>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/assign/list_inserter.hpp>

#include "output.h"
#include "log.h"

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;
using namespace boost::assign;

typedef vector<std::complex<double>> state_type;

int N = 0;
int num = 0;

double g_per = 0.0;
double g_par = 0.0;
double k_a = 0.0;

const std::complex<double> halfI = std::complex<double>(0, 0.5);
std::complex<double> n;

vector <double> pump_pwr;
vector <double> noise_vec;

state_type CFvals;

state_type a_vec;
state_type b_vec1;
state_type b_vec2;
state_type Q;
state_type NL_term;
state_type aD_vec;
state_type D_vec;
state_type holder;
state_type bB_vec;
state_type obs_temp;

matrix<std::complex<double>> D_mat;

matrix <double> A;
matrix <double> B;

enum solver_type {FP, RING, UCF};

struct observer {

private:
    std::vector<double>& tout;
    std::vector<double>& yout_real;
    std::vector<double>& yout_imag;
    std::vector<state_type> & yout;

public:
    observer(std::vector<double>& tout_,
                       std::vector<double>& real_,
                       std::vector<double>& imag_,
                       std::vector<state_type> & yout_) :
        tout(tout_), yout_real(real_), yout_imag(imag_), yout(yout_)
    {
        FILE_LOG(logINFO) << "Initializing observer";
    }

    void operator()(const state_type &x, const double t)
    {
        FILE_LOG(logDEBUG4) << "begin of operator";
        tout.push_back(t);

        
        for(int j = 0; j < N * (N + 2); j++) {
            if (j > N*N) {
             FILE_LOG(logDEBUG4) << "Observer Vals Re: " << 
                     x[j].real() << " Im: " << x[j].imag();
            }
            yout_real.push_back(x[j].real());
            yout_imag.push_back(x[j].imag());
        }
        
        yout.push_back(x);
        
        num++;
    }

};

struct TDSSolvers_RING {
    
    TDSSolvers_RING()
    {
        FILE_LOG(logDEBUG) << "Init FP";
    }

    void operator()(const state_type& y, const state_type& dy, const double t)
    {
        //a_vec  = y(1:N);
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N), y.size());
        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N;
                NL_term[i] += A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(i * (N - 1) + j) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), k_a);
        state_type & dy_a = *(state_type*)&dy;
        state_type dy1 = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vec1 / n * n);
        state_type dy2 = -(g_per * b_vec1 - halfI * 2 * g_per * a_vecTD_mat);
        state_type dy3 = - (g_par * (D_vec - pump_pwr) + halfI * g_par * NL_term);
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

struct TDSSolvers_UCF {
    
    TDSSolvers_UCF() 
    {
        FILE_LOG(logDEBUG) << "Init FP";
    }

    void operator()(const state_type& y, const state_type& dy, const double t)
    {
        FILE_LOG(logDEBUG4) << "begin UCF";
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N), y.size());
        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N;
                NL_term[i] += A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecT_D_mat(N);
        state_type b_vecT_B(N);

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecT_D_mat[i] += D_vec(i * (N - 1) + j) * a_vec(i);
                b_vecT_B[i] += B(i, j) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), k_a);
        state_type dy1 = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vecT_B);
        state_type dy2 = -(g_per * b_vec1 - halfI * 2 * g_per * a_vecT_D_mat);
        state_type dy3 = -(g_par * (D_vec - pump_pwr) + 2 * halfI * g_par * NL_term);
        state_type & dy_a = *(state_type*)&dy;
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

struct TDSSolvers_FP {

    TDSSolvers_FP() 
    {
        FILE_LOG(logDEBUG) << "Init FP";
    }
    
    void operator()(const state_type& y, const state_type& dy, const double t)
    {
        FILE_LOG(logDEBUG4) << "begin FP";
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, 2 * N , y.size());
        matrix<std::complex<double>> Q = imag(outer_prod(conj(b_vec1), a_vec));

        for(int i = 0; i < N * N; i++) {
            for(int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N ;
                NL_term[i] += A(i, j) * Q(a, b);
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                a_vecTD_mat[i] += D_vec(i * (N - 1) + j) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), k_a);
        state_type dy1 = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vec1 / n * n);
        state_type dy2 = -(g_per * b_vec1 - halfI * 2 * g_per * a_vecTD_mat);
        state_type dy3 = -(g_par * (D_vec - pump_pwr) - g_par * NL_term);
        state_type & dy_a = *(state_type*)&dy;
        state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
        oi = std::copy_n(dy2.begin(), dy2.size(), oi);
        std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Redirect cout to MEX printf
    mstream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    // Initialize logger
    FILELog::ReportingLevel() = logDEBUG;
    
    FILE_LOG(logINFO) << "beginning of mexfunction";
    g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.1";
    g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.2";
    k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.31";
    char * temp_basis_type = mxArrayToString(prhs[5]);
    FILE_LOG(logDEBUG4) << "checkpoint 0.32";
    std::string basis_type = std::string(temp_basis_type);
    FILE_LOG(logDEBUG4) << "checkpoint 0.33";
    mxFree(temp_basis_type);
    FILE_LOG(logDEBUG4) << "checkpoint 0.4";
    N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.5";
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.6";
    int dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.7";
    double * ptB = mxGetPr(mxGetField(prhs[0], 0, "A"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.8";
    int dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
    FILE_LOG(logDEBUG4) << "checkpoint 0.9";
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    FILE_LOG(logDEBUG4) << "checkpoint 0.10";
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
    FILE_LOG(logDEBUG4) << "checkpoint 1";
    CFvals.resize(N);
    FILE_LOG(logDEBUG4) << "checkpoint 1.01";

    for(int i = 0; i < N; i++) {
        CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
    }

    FILE_LOG(logDEBUG4) << "checkpoint 1.1";
    A.resize(dimA, dimA);
    FILE_LOG(logDEBUG4) << "checkpoint 1.2";

    for(int i = 0; i < dimA; i++) {
        for(int j = 0; j < dimA; j++) {
            A(i, j) = ptA[i + j * dimA];
        }
    }

    FILE_LOG(logDEBUG4) << "checkpoint 1.3";
    B.resize(dimB, dimB);
    FILE_LOG(logDEBUG4) << "checkpoint 1.4";

    for(int i = 0; i < dimB; i++) {
        for(int j = 0; j < dimB; j++) {
            B(i, j) = ptB[i + j * dimB];
        }
    }

    FILE_LOG(logDEBUG4) << "checkpoint 2";
    auto pump_pt = mxGetPr(prhs[1]);
    pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

    for(int i = 0; i < pump_pwr.size(); i++)
        pump_pwr[i] = pump_pt[i];

    n = *mxGetPr(mxGetField(prhs[0], 0, "n"));
    FILE_LOG(logDEBUG4) << "checkpoint 3";
    auto noise_pt = mxGetPr(prhs[4]);
    //noise_vec.resize(mxGetDimensions(prhs[4])[0]);
    noise_vec.resize(N * (N + 2));

    for(int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = noise_pt[i];

    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);
    FILE_LOG(logDEBUG4) << "Resizing all structures";
    /*
    const int dt_steps = 5089;
    double dt = (*t_final_pt - *t_initial_pt) / dt_steps;
    outT.resize(dt_steps + 1, 1);
    outY.resize(dt_steps + 1, N * (N + 2));
    */
    FILE_LOG(logDEBUG4) << "checkpoint 4";
    holder.resize(N * (N + 2));
    a_vec.resize(N);
    b_vec1.resize(N);
    b_vec2.resize(N);
    D_mat.resize(N, N);
    D_vec.resize(N * N);
    aD_vec.resize(N);
    Q.resize(N * N);
    NL_term.resize(N * N);
    obs_temp.resize(N * (N + 2));
    num = 0;
    
    FILE_LOG(logDEBUG4) << "checkpoint 5";
    FILE_LOG(logDEBUG4) << "checkpoint 6";

    FILE_LOG(logDEBUG4) << "checkpoint 7";
    FILE_LOG(logDEBUG4) << "checkpoint 8";
    
    std::vector<double> tout;
    std::vector<double> yout_real;
    std::vector<double> yout_imag;
    std::vector<state_type> yout;
    
    FILE_LOG(logINFO) << "Before Integrate";

    TDSSolvers_FP fp;
	typedef runge_kutta4<vector<double>> rk4; //unused
	typedef runge_kutta_dopri5<vector<double>> stepper_type;
	//typedef runge_kutta_dopri5<state_type> stepper_type; //does not work?
    
    if(basis_type.compare("RING") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling RING integrator";
        integrate_const(make_dense_output(1.0e-6, 1.0e-6, stepper_type()),
                        TDSSolvers_RING(), noise_vec,
                        *t_initial_pt, *t_final_pt, 0.1,
                        observer(tout, yout_real, yout_imag, yout));
    }

    else if(basis_type.compare("UCF") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling UCF integrator";
        integrate_const(make_dense_output(1.0e-6, 1.0e-3, stepper_type()),
                        TDSSolvers_UCF(), noise_vec,
                        *t_initial_pt, *t_final_pt, 0.1,
                        observer(tout, yout_real, yout_imag, yout));
    }

    else if(basis_type.compare("FP") == 0) {
        FILE_LOG(logDEBUG4) << "Before calling FP integrator";
        integrate_const(make_dense_output(1.0e-6, 1.0e-3, stepper_type()),
                        fp, noise_vec,
                        *t_initial_pt, *t_final_pt, 0.1,
                        observer(tout, yout_real, yout_imag, yout));
    }

    else {
        exit(EXIT_FAILURE);
    }

    FILE_LOG(logINFO) << "After Integrate";
    FILE_LOG(logDEBUG) << "before assigning output";
    FILE_LOG(logDEBUG) << num;
    int index = 0;
    double * outT = (double *) mxMalloc(num * sizeof(double));
    FILE_LOG(logDEBUG4) << outT;
    FILE_LOG(logDEBUG4) << tout.size();
    FILE_LOG(logDEBUG4) << "after malloc";

    FILE_LOG(logDEBUG4) << "after first loop";
    double * outYr = (double *) mxCalloc(num * N * (N + 2), sizeof(double));
    double * outYi = (double *) mxCalloc(num * N * (N + 2), sizeof(double));
    
    for(index = 0; index < num; index++) {
        FILE_LOG(logDEBUG4) << "in loop " << index << std::endl;
        outT[index] = tout[index]; 
    
        for(int j = 0; j < N * (N + 2); j++) {
            if (j > N*N && index > num-2) {
                FILE_LOG(logDEBUG) << "Values Re: " << 
                        yout[index][j].real() << " Im: " << 
                        yout[index][j].imag();
            }
            
            outYr[index *N * (N + 2) + j] = yout[index][j].real();
            outYi[index *N * (N + 2) + j] = yout[index][j].imag();
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL); 
    mxSetPr(plhs[0], outT);
    mxSetM(plhs[0], num);
    mxSetN(plhs[0], 1);
    
    FILE_LOG(logDEBUG4) << "in between assigning output";
    
    plhs[1] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); 
    mxSetPr(plhs[1], outYr);
    mxSetPi(plhs[1], outYi);
    mxSetM(plhs[1], num);
    mxSetN(plhs[1], N * (N + 2));
    
    FILE_LOG(logINFO) << "End of Mex Function";
    
    // Replace redirection
    std::cout.rdbuf(outbuf);
}
