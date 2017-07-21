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

int mindex;
vector<std::complex<double>> noise_vec;
int N;
double g_per;
double g_par;
double k_a;
vector <std::complex<double>> CFvals;
vector <double> outT;
matrix<std::complex<double>> outY;
vector <double> pump_pwr;
matrix <double> A;
matrix <double> B;
std::complex <double> n;

void TDSSolvers_RING(const state_type &y, const state_type &dy, const double t) {
    
    state_type a_vec;
    state_type b_vec1;
    state_type b_vec2;
    matrix<std::complex<double>> D_mat;
    state_type Q;
    state_type NL_term;
    state_type aD_vec;
    state_type D_vec;
    state_type holder;
    
    holder.resize(N*(N+2));
    a_vec.resize(N);
    b_vec1.resize(N);
    b_vec2.resize(N);
    D_mat.resize(N,N);
    D_vec.resize(N*N);
    aD_vec.resize(N);
    Q.resize(N*N);
    NL_term.resize(N*N);
    
    
    for (int i = 0; i < N; i++)
       a_vec[i] = y[i];
    
    
    for (int i = 0; i < N; i++) {
       b_vec1[i] = y[i+N];
       b_vec2[i] = std::conj(y[i+N]);
    }
    
    
    for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
           const std::complex<double> & x = y[((i+2)*N) + j];
           D_mat.insert_element(i,j,x);
           D_vec[(i*N) + j] = x;
       }
    }
    
    
    for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
           const std::complex<double> & x = b_vec2[i] * a_vec[j];
           const std::complex<double> & y = b_vec2[j] * a_vec[i];
           Q[(i*N) + j] = (x-y);
       }
    }
    
    for (int i = 0; i < N*N; i++) {
       for (int j = 0; j < N*N; j++) {
           NL_term[i] += A.at_element(i,j) * Q[j];
       }
    }
    
    const std::complex<double> halfI = std::complex<double>{0,0.5};
    
    for (int i = 0; i < N; i++) {
        CFvals[i] = a_vec[i] * (k_a - ((CFvals[i] * CFvals[i]) / k_a)) * (halfI);
    }
    
    for (int i = 0; i < N; i++) {
        holder[i] = CFvals[i] + (halfI * k_a * b_vec1[i] / (n * n));
    }
    
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aD_vec[i] += a_vec[j] * D_mat.at_element(j,i);
        }
    }
    
    for (int i = N; i < N+N; i++) {
        holder[i] = (-1 * g_per * b_vec1[i-N]) - (halfI * 2 * g_per * aD_vec[i-N]);
    }
    
    for (int i = N+N; i < N * (N+2); i++) {
        holder[i] = (-1 * g_par * (D_vec[i-(N+N)] - pump_pwr[i-(N+N)])) 
                + (halfI * g_par * NL_term[i-(N+N)]);
    }
    
    * (state_type *) &dy = holder;
}

void observer( const state_type &x , const double t )
{
    outT[mindex] = t;
    for (int j = 0; j < N; j++) {
        outY.insert_element(mindex,j,x[j]);
    }
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*
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
    
    auto noise_pt = mxGetPr(prhs[4]);
    noise_vec.resize(mxGetDimensions(prhs[4])[0]);
    for (int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = noise_pt[i];
    
    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[4]);
    
    mindex = 0;
    /*
    integrate(TDSSolvers_RING , noise_vec , *t_initial_pt , *t_final_pt , 0.050 , observer );
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) mindex, (mwSize) 1, mxREAL);
    mxSetData(plhs[0], &outT);
    
    plhs[1] = mxCreateDoubleMatrix((mwSize) mindex, (mwSize) (N*(N+2)), mxCOMPLEX);
    mxSetData(plhs[1], &outY);
    */
}


