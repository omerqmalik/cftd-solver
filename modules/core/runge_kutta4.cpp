#include <mex.h>
#include <matrix.h>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <iostream>
#include <boost/assign/list_inserter.hpp>

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

typedef runge_kutta4<vector<double>> rk4;
typedef runge_kutta_dopri5<vector<double>> stepper_type;

enum solver_type {FP, RING, UCF};

//struct observer {
    
    std::vector<double> tout;
    std::vector<double> yout_real;
    std::vector<double> yout_imag;

    void observing (const state_type &x, const double t)
    {
        ////std::cout << "begin of operator" << std::endl;
        tout.push_back(t);

        for(int j = 0; j < N * (N + 2); j++) {
            yout_real.push_back(x[j].real());
            yout_imag.push_back(x[j].imag());
        }
        
        num++;
    }

//};

struct TDSSolvers_RING {
    //bool first_time = true;
    
    void operator () (const state_type& y, const state_type& dy, const double t)
    {
        ////std::cout << "begin RING" << std::endl;
        //if (first_time) {
            ////std::cout << "Size of y:" << y.size() << " dy :" << dy.size() << std::endl;
            //first_time = false;
        //}
         
        //a_vec  = y(1:N);
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N + 1), y.size());

        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

        for (int i = 0; i < N * N; i++) {
            for (int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N - 1;
                NL_term[i] += A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int a = j / N;
                int b = j % N - 1;
                a_vecTD_mat[i] += D_vec(i *(N - 1) + j) * a_vec(i);
            }
        }
       
		scalar_vector<std::complex<double>> k_avec(a_vec.size(), k_a);
		state_type & dy_a = *(state_type*)&dy;
		state_type dy1 = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vec1 / n*n);
		state_type dy2 = -(g_per * b_vec1 - halfI * 2 * g_per * a_vecTD_mat);
		state_type dy3 = - (g_par * (D_vec - pump_pwr) + halfI * g_par * NL_term);
		
		state_type::iterator oi = std::copy_n(dy1.begin(), dy1.size(), dy_a.begin());
		oi = std::copy_n(dy2.begin(), dy2.size(), oi);
		std::copy_n(dy3.begin(), dy3.size(), oi);
    }
};

struct TDSSolvers_UCF {
    void operator () (const state_type& y, const state_type& dy, const double t)
    {
        ////std::cout << "begin UCF" << std::endl;
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2 * N);
        D_vec = subrange(y, (2 * N + 1), y.size());

        matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

        for (int i = 0; i < N * N; i++) {
            for (int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N - 1;
                NL_term[i] += A(i, j) * (Q(a, b) - Q(b, a));
            }
        }

        state_type a_vecT_D_mat(N);
        state_type b_vecT_B(N);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int a = j / N;
                int b = j % N - 1;
                a_vecT_D_mat[i] += D_vec(i *(N - 1) + j) * a_vec(i);
                b_vecT_B[i] += B(i,j) * a_vec(i);
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
    void operator () (const state_type& y, const state_type& dy, const double t)
    {
        ////std::cout << "begin FP" << std::endl;
        a_vec = subrange(y, 0, N);
        b_vec1 = subrange(y, N, 2*N);
        D_vec = subrange(y, (2*N+1), y.size());

        matrix<std::complex<double>> Q = imag(outer_prod(conj(b_vec1), a_vec));

        for (int i = 0; i < N * N; i++) {
            for (int j = 0; j < N * N; j++) {
                int a = j / N;
                int b = j % N - 1;
                NL_term[i] += A(i, j) * Q(a, b);
            }
        }

        state_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int a = j / N;
                int b = j % N - 1;
                a_vecTD_mat[i] += D_vec(i *(N-1) + j) * a_vec(i);
            }
        }

        scalar_vector<std::complex<double>> k_avec(a_vec.size(), k_a);
        state_type dy1 = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vec1 / n*n);
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
    //std::cout << "beginning of mexfunction" << std::endl;
    
    g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    
    //std::cout << "checkpoint 0.1" << std::endl;
    
    g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    
    //std::cout << "checkpoint 0.2" << std::endl;
    
    k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    
    ////std::cout << "checkpoint 0.31" << std::endl;
    
    char * temp_basis_type = mxArrayToString(prhs[5]);
    
    ////std::cout << "checkpoint 0.32" << std::endl;
    
    std::string basis_type = std::string(temp_basis_type);
    
    ////std::cout << "checkpoint 0.33" << std::endl;
    
    mxFree(temp_basis_type);
    
    ////std::cout << "checkpoint 0.4" << std::endl;
    
    N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    
    ////std::cout << "checkpoint 0.5" << std::endl;
    
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    
    ////std::cout << "checkpoint 0.6" << std::endl;
    
    int dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    
    ////std::cout << "checkpoint 0.7" << std::endl;
    
    double * ptB = mxGetPr(mxGetField(prhs[0], 0, "A"));
    
    ////std::cout << "checkpoint 0.8" << std::endl;
    
    int dimB = mxGetDimensions(mxGetField(prhs[0], 0, "B"))[0];
    
    ////std::cout << "checkpoint 0.9" << std::endl;
    
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    
    ////std::cout << "checkpoint 0.10" << std::endl;
    
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
    
    ////std::cout << "checkpoint 1" << std::endl;
    
    CFvals.resize(N);
    
    ////std::cout << "checkpoint 1.01" << std::endl;
    
    for(int i = 0; i < N; i++) {
        CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
    }
    
    ////std::cout << "checkpoint 1.1" << std::endl;
    
    A.resize(dimA, dimA);

    ////std::cout << "checkpoint 1.2" << std::endl;
    
    for(int i = 0; i < dimA; i++) {
        for(int j = 0; j < dimA; j++) {
            A(i, j) = ptA[i + j * dimA];
        }
    }
    
    ////std::cout << "checkpoint 1.3" << std::endl;
    
    B.resize(dimB, dimB);

    ////std::cout << "checkpoint 1.4" << std::endl;
    
    for(int i = 0; i < dimB; i++) {
        for(int j = 0; j < dimB; j++) {
            B(i, j) = ptB[i + j * dimB];
        }
    }
    
    //std::cout << "checkpoint 2" << std::endl;
    
    auto pump_pt = mxGetPr(prhs[1]);
    pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

    for(int i = 0; i < pump_pwr.size(); i++)
        pump_pwr[i] = pump_pt[i];

    n = *mxGetPr(mxGetField(prhs[0], 0, "n"));
    
    //std::cout << "checkpoint 3" << std::endl;
    
    auto noise_pt = mxGetPr(prhs[4]);
    //noise_vec.resize(mxGetDimensions(prhs[4])[0]);
    noise_vec.resize(N*(N+2));
    
    for(int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = noise_pt[i];

    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);

    //////std::cout << "Resizing all structures" << std::endl;
    
    /*
    const int dt_steps = 5089;
    double dt = (*t_final_pt - *t_initial_pt) / dt_steps;
    outT.resize(dt_steps + 1, 1);
    outY.resize(dt_steps + 1, N * (N + 2));
    */
    
    ////std::cout << "checkpoint 4" << std::endl;
    
    holder.resize(N * (N + 2));
    a_vec.resize(N);
    b_vec1.resize(N);
    b_vec2.resize(N);
    D_mat.resize(N, N);
    D_vec.resize(N * N);
    aD_vec.resize(N);
    Q.resize(N * N);
    NL_term.resize(N * N);
    obs_temp.resize(N*(N+2));
    
    num = 0;
    
    ////std::cout << "checkpoint 5" << std::endl;
    
    //std::cout << "checkpoint 6" << std::endl;
    
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
    
    //std::cout << "checkpoint 7" << std::endl;
    
    //observer observing;
    
    //std::cout << "checkpoint 8" << std::endl;

    if(basis_type.compare("RING") == 0) {
        ////std::cout << "Before calling RING integrator" << std::endl;
        integrate_const(make_controlled(1.0e-6, 1.0e-6, stepper_type()), TDSSolvers_RING(), noise_vec, *t_initial_pt, *t_final_pt, 0.1, observing);     
    }

    else if(basis_type.compare("UCF") == 0) {
        ////std::cout << "Before calling UCF integrator" << std::endl;
        integrate_const(make_controlled(1.0e-6, 1.0e-6, stepper_type()), TDSSolvers_UCF(), noise_vec, *t_initial_pt, *t_final_pt, 0.1, observing);      
    }

    else if(basis_type.compare("FP") == 0) {
        ////std::cout << "Before calling FP integrator" << std::endl;
        integrate_const(make_controlled(1.0e-6, 1.0e-6, stepper_type()), TDSSolvers_FP(), noise_vec, *t_initial_pt, *t_final_pt, 0.1, observing);       
    }

    else {
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "before assigning output" << std::endl;
    
    int index = 0;
    
    //std::cout << num << std::endl;
    
    
    double * outT = (double *) mxMalloc(num * sizeof(double));
    
    //std::cout << outT << std::endl;
    
    //std::cout << tout.size() << std::endl;
    
    //std::cout << "after malloc" << std::endl;
    
    for ( index = 0; index < num; index++ ) {
        //std::cout << "in loop " << index << std::endl;
        outT[index] = tout[index];
    }
    
    //std::cout << "after first loop" << std::endl;
    
    double * outYr = (double *) mxMalloc(num * N * (N+2) * sizeof(double));
    double * outYi = (double *) mxMalloc(num * N * (N+2) * sizeof(double));
    for ( index = 0; index < num * N * (N+2); index++ ) {
        outYr[index] = yout_real[index];
        outYi[index] = yout_imag[index];
    }
            
    mxSetPr(plhs[0], outT);
    mxSetM(plhs[0], num);
    mxSetN(plhs[0], 1);
    
    //std::cout << "in between assigning output" << std::endl;
    
    mxSetPr(plhs[1], outYr);
    mxSetPi(plhs[1], outYi);
    mxSetM(plhs[1], num);
    mxSetN(plhs[1], N*(N+2));
    
    //std::cout << "after assigning output" << std::endl;
}


