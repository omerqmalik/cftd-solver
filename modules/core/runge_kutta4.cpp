#include <mex.h>
#include <matrix.h>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <iostream>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

typedef vector<std::complex<double>> state_type;

int mindex = 0;
int N = 0;

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

matrix<std::complex<double>> D_mat;

matrix <double> A;
matrix <double> B;
matrix <double> outT;

matrix<std::complex<double>> outY;

typedef runge_kutta4<vector<double>> rk4;

enum solver_type {FP, RING, UCF};


void observer(const state_type &x, const double t)
{
    //std::cout << "beginning of observer, mindex = " << mindex << std::endl;
    outT(mindex, 0) = t;

    for(int j = 0; j < N * (N + 2); j++) {
        //std::cout << "j = " << j << " x[j] = " << x[j] << std::endl;
        outY(mindex, j) = x[j];
    }

    /*
    matrix_row<matrix<std::complex<double> > > mr(outY, mindex);
    std::copy(x.begin(), x.end(), mr.begin());
    */
    mindex++;
    //std::cout << "after observer" << std::endl;
}

struct TDSSolvers_RING {
    void operator () (const state_type& y, const state_type& dy, const double t) const
    {
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
        *(state_type*)&dy = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + (halfI * k_a * b_vec1 / n*n) 
            - (g_per * b_vec1 - halfI * 2 * g_per * a_vecTD_mat) 
            - (g_par * (D_vec - pump_pwr) + halfI * g_par * NL_term);
    }
};

struct TDSSolvers_UCF {
    void operator () (const state_type& y, const state_type& dy, const double t) const
    {
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
        *(state_type*)&dy = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) +
            (halfI * k_a * b_vecT_B) -
            (g_per * b_vec1 - halfI * 2 * g_per * a_vecT_D_mat) -
            (g_par * (D_vec - pump_pwr) + 2 * halfI * g_par * NL_term);
    }
};

struct TDSSolvers_FP {
    void operator () (const state_type& y, const state_type& dy, const double t) const
    {
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
        *(state_type*)&dy = element_prod((k_avec - (element_prod(CFvals, CFvals) / k_a)), (halfI * a_vec)) + 
            (halfI * k_a * b_vec1 / n*n) - 
            (g_per * b_vec1 - halfI * 2 * g_per * a_vecTD_mat) - 
            (g_par * (D_vec - pump_pwr) - g_par * NL_term);
    }
};

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    //std::cout << "beginning of mexfunction" << std::endl;
    g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    std::string basis_type = mxArrayToString(prhs[5]);
    N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    int dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));

    for(int i = 0; i < N; i++) {
        CFvals[i] = std::complex<double>(CFvalsReal[i], CFvalsImag[i]);
    }

    A.resize(dimA, dimA);

    for(int i = 0; i < dimA; i++) {
        for(int j = 0; j < dimA; j++) {
            A(i, j) = ptA[i + j * dimA];
        }
    }

    auto pump_pt = mxGetPr(prhs[1]);
    pump_pwr.resize(mxGetDimensions(prhs[1])[0]);

    for(int i = 0; i < pump_pwr.size(); i++)
        pump_pwr[i] = pump_pt[i];

    n = *mxGetPr(mxGetField(prhs[0], 0, "n"));
    auto noise_pt = mxGetPr(prhs[4]);
    noise_vec.resize(mxGetDimensions(prhs[4])[0]);

    for(int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = noise_pt[i];

    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);
    mindex = 0;

    std::cout << "Resizing all structures" << std::endl;

    const int dt_steps = 5089;
    double dt = (*t_final_pt - *t_initial_pt) / dt_steps;
    outT.resize(dt_steps + 1, 1);
    outY.resize(dt_steps + 1, N * (N + 2));
    holder.resize(N * (N + 2));
    a_vec.resize(N);
    b_vec1.resize(N);
    b_vec2.resize(N);
    D_mat.resize(N, N);
    D_vec.resize(N * N);
    aD_vec.resize(N);
    Q.resize(N * N);
    NL_term.resize(N * N);

    if(basis_type.compare("RING") == 0) {
        std::cout << "Before calling RING integrator" << std::endl;
        integrate_const(rk4(), TDSSolvers_RING(), noise_vec, *t_initial_pt, *t_final_pt, dt, observer);
    }

    else if(basis_type.compare("UCF") == 0) {
        std::cout << "Before calling UCF integrator" << std::endl;
        integrate_const(rk4(), TDSSolvers_UCF(), noise_vec, *t_initial_pt, *t_final_pt, dt, observer);
    }

    else if(basis_type.compare("FP") == 0) {
        std::cout << "Before calling FP integrator" << std::endl;
        integrate_const(rk4(), TDSSolvers_FP(), noise_vec, *t_initial_pt, *t_final_pt, dt, observer);
    }

    else {
        exit(EXIT_FAILURE);
    }

    //outT.resize (mindex, 1);
    //outY.resize (mindex, N*(N+2));
    std::cout << "after integrate" << std::endl;
    //std::cout << "before assigning output" << std::endl;
    plhs[0] = mxCreateDoubleMatrix((mwSize) outT.size1(), (mwSize) outT.size2(), mxREAL);
    //std::cout << "before setting field 1" << std::endl;
    double* outputMatrix1 = (double *)mxGetData(plhs[0]);

    //std::cout << "before output assignment loop" << std::endl;

    for(int col = 0; col < outT.size2(); col++) {
        for(int row = 0; row < outT.size1(); row++) {
            int i = row + col * outT.size1();
            //std::cout << "i = " << i << std::endl;
            //std::cout << "outT value = " << outT(col,row) << std::endl;
            outputMatrix1[i] = outT(row, col);
        }
    }

    //std::cout << "in the middle of assigning output" << std::endl;
    plhs[1] = mxCreateDoubleMatrix((mwSize) outY.size1(), (mwSize) outY.size2(), mxCOMPLEX);
    double * opReal = (double *) mxGetPr(plhs[1]);
    double * opImag = (double *) mxGetPi(plhs[1]);

    for(int col = 0; col < outY.size2(); col++) {
        for(int row = 0; row < outY.size1(); row++) {
            int i = row + col * outY.size1();
            //std::cout << "i = " << i << std::endl;
            //std::cout << "outY value = " << outY(col,row) << std::endl;
            opReal[i] = outY(row, col).real();
            opImag[i] = outY(row, col).imag();
        }
    }

    std::cout << "after assigning output" << std::endl;
}


