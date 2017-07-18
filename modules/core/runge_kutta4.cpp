#include <mex.h>
#include <matrix.h>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <iostream>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

typedef std::vector<double> state_type;

int N;
double g_per;
double g_par;
double k_a;
vector <std::complex<double>> CFvals;
vector <double> pump_pwr;
matrix <double> A;
matrix <double> B;
std::complex <double> n;

/*
void TDSSolvers_RING(const state_type &y, state_type &dy, const double t) {
    
    a_vec  = y(1:N);
    b_vec  = y(N+1:2*N);
    D_vec  = y((2*N+1):end);
    D_mat  = reshape(D_vec,[N,N]);
    
    Q = conj(b_vec)*a_vec.';
    NL_term = A*reshape(Q-Q',[N^2,1]);

    dy = [(k_a-CFvals.^2/k_a)*1i*0.5.*a_vec + 1i*k_a*0.5*b_vec/n^2; 
          -g_per*b_vec - 1i*g_per*(a_vec.'*D_mat).'; 
          -g_par*(D_vec - pump_pwr) + 1i*g_par*0.5*NL_term];
}
*/

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    
    
    
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    int dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    
    A.resize(dimA,dimA);
    for(int i=0; i<dimA; i++) {
      for(int j=0; j<dimA; j++) {
         A(i,j) = ptA[i+j*dimA];
      }
    }
    
    auto pump_pt = mxGetPr(prhs[1]);
    
    pump_pwr.resize(mxGetDimensions(prhs[1])[0]);
    for (int i = 0; i < pump_pwr.size(); i++)
        pump_pwr[i] = pump_pt[i];
    
    
    n = *mxGetPr(mxGetField(prhs[0], 0, "n"));
    
    auto bleh = mxGetPr(prhs[3]);
    vector <double> noise_vec(mxGetDimensions(prhs[3])[0]);
    for (int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = bleh[i];
    
    double *T;
    double *Y;
    
    
    
    
    
}


