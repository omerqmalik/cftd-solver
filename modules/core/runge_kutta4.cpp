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

std::complex<double> n;

vector <double> pump_pwr(100);
vector <double> noise_vec(100);

state_type CFvals(1000);

matrix <double> A;
matrix <double> B;
matrix <double> outT;

matrix<std::complex<double>> outY;


void TDSSolvers_RING(const state_type &y, const state_type &dy, const double t) {
    //std::cout << "start of solver" << std::endl;
    
    state_type a_vec;
    state_type b_vec1;
    state_type b_vec2;
    state_type Q;
    state_type NL_term;
    state_type aD_vec;
    state_type D_vec;
    state_type holder;
    
    matrix<std::complex<double>> D_mat;
    
    //std::cout << "after setup" << std::endl;
    
    holder.resize(N*(N+2));
    a_vec.resize(N);
    b_vec1.resize(N);
    b_vec2.resize(N);
    D_mat.resize(N,N);
    D_vec.resize(N*N);
    aD_vec.resize(N);
    Q.resize(N*N);
    NL_term.resize(N*N);
    
    //std::cout << "after resizing" << std::endl;
    
    for (int i = 0; i < N; i++)
       a_vec[i] = y[i];
    
    //std::cout << "after loop 1" << std::endl;
    
    for (int i = 0; i < N; i++) {
       b_vec1[i] = y[i+N];
       b_vec2[i] = std::conj(y[i+N]);
    }
    
    //std::cout << "after loop 2" << std::endl;
    
    for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
           const std::complex<double> & x = y[((i+2)*N) + j];
           D_mat.insert_element(i,j,x);
           D_vec[(i*N) + j] = x;
       }
    }
    
    //std::cout << "after loop 3" << std::endl;
    
    for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
           const std::complex<double> & x = b_vec2[i] * a_vec[j];
           const std::complex<double> & y = b_vec2[j] * a_vec[i];
           Q[(i*N) + j] = (x-y);
       }
    }
    
    //std::cout << "after loop 4" << std::endl;
    
    for (int i = 0; i < N*N; i++) {
       for (int j = 0; j < N*N; j++) {
           NL_term[i] += A.at_element(i,j) * Q[j];
       }
    }
    
    //std::cout << "after loop 5" << std::endl;
    
    const std::complex<double> halfI = std::complex<double>(0,0.5);
    
    for (int i = 0; i < N; i++) {
        CFvals(i) = a_vec[i] * (k_a - ((CFvals(i) * CFvals(i)) / k_a)) * (halfI);
    }
    
    //std::cout << "after loop 6" << std::endl;
    
    for (int i = 0; i < N; i++) {
        holder[i] = CFvals(i) + (halfI * k_a * b_vec1[i] / (n * n));
    }
    
    //std::cout << "after loop 7" << std::endl;
    
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aD_vec[i] += a_vec[j] * D_mat.at_element(j,i);
        }
    }
    
    //std::cout << "after loop 8" << std::endl;
    
    for (int i = N; i < N+N; i++) {
        holder[i] = (-1 * g_per * b_vec1[i-N]) - (halfI * 2 * g_per * aD_vec[i-N]);
    }
    
    //std::cout << "after loop 9" << std::endl;
    
    for (int i = N+N; i < N * (N+2); i++) {
        holder[i] = (-1 * g_par * (D_vec[i-(N+N)] - pump_pwr[i-(N+N)])) 
                + (halfI * g_par * NL_term[i-(N+N)]);
    }
    
    //std::cout << "after loop 10" << std::endl;
    
    * (state_type *) &dy = holder;
    
    //std::cout << "after final" << std::endl;
}

void observer( const state_type &x , const double t )
{
    //std::cout << "beginning of observer, mindex = " << mindex << std::endl;
    
    //std::cout << "t = " << t << std::endl;
    
    outT.resize(std::max(outT.size1(), (unsigned long) (mindex+100)*2), 1);
    
    outT(mindex, 0) = t;
    
    //std::cout << "before loop, N = " << N << std::endl;
    
    outY.resize(std::max(outY.size1(), (unsigned long) (mindex+100)*2), N*(N+2));
    
    for (int j = 0; j < N*(N+2); j++) {
        
        //std::cout << "j = " << j << " x[j] = " << x[j] << std::endl;
        
        outY(mindex,j) = x[j];
    }
    
    mindex++;
    
    //std::cout << "after observer" << std::endl;
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    
    
    //std::cout << "beginning of mexfunction" << std::endl;
    
    g_per = *mxGetPr(mxGetField(prhs[0], 0, "g_per"));
    g_par = *mxGetPr(mxGetField(prhs[0], 0, "g_par"));
    k_a = *mxGetPr(mxGetField(prhs[0], 0, "k_a"));
    N =  mxGetDimensions(mxGetField(prhs[0], 0, "CFvals"))[0];
    double * ptA = mxGetPr(mxGetField(prhs[0], 0, "A"));
    int dimA = mxGetDimensions(mxGetField(prhs[0], 0, "A"))[0];
    
    double * CFvalsReal = (double *) mxGetPr(mxGetField(prhs[0], 0, "CFvals"));
    
    double * CFvalsImag = (double *) mxGetPi(mxGetField(prhs[0], 0, "CFvals"));
    
    for (int i = 0; i < N; i++) {
        CFvals[i] = std::complex<double>(CFvalsReal[i],CFvalsImag[i]);
    }
    
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
    
    
    
    auto noise_pt = mxGetPr(prhs[4]);
    noise_vec.resize(mxGetDimensions(prhs[4])[0]);
    
    for (int i = 0; i < noise_vec.size(); i++)
        noise_vec[i] = noise_pt[i];
    
    auto t_initial_pt = mxGetPr(prhs[2]);
    auto t_final_pt = mxGetPr(prhs[3]);
    
    mindex = 0;
    outT.resize(10000,1);
    outY.resize(10000,N*(N+2));
    
    //std::cout << "before integrate" << std::endl;
    
    integrate(TDSSolvers_RING , noise_vec , *t_initial_pt , *t_final_pt , 0.050 , observer);
    
    //std::cout << "after integrate" << std::endl;
    
    //std::cout << "before assigning output" << std::endl;
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) outT.size1(), (mwSize) outT.size2(), mxREAL);
    
    //std::cout << "before setting field 1" << std::endl;
    
    double* outputMatrix1 = (double *)mxGetData(plhs[0]);
    
    //std::cout << "before output assignment loop" << std::endl;
    
    for (int col=0; col < outT.size2(); col++) {
        for (int row=0; row < outT.size1(); row++) {
            
            int i = row + col*outT.size1();
            
            //std::cout << "i = " << i << std::endl;
            
            //std::cout << "outT value = " << outT(col,row) << std::endl;
            
            outputMatrix1[i] = outT(row,col);
        }
    }
    
    
    //std::cout << "in the middle of assigning output" << std::endl;
    
    plhs[1] = mxCreateDoubleMatrix((mwSize) outY.size1(), (mwSize) outY.size2(), mxCOMPLEX);
    
    double * opReal = (double *) mxGetPr(plhs[1]);
    
    double * opImag = (double *) mxGetPi(plhs[1]);
    
    
    for (int col=0; col < outY.size2(); col++) {
        for (int row=0; row < outY.size1(); row++) {
            
            int i = row + col*outY.size1();
            
            //std::cout << "i = " << i << std::endl;
            
            //std::cout << "outY value = " << outY(col,row) << std::endl;
            
            opReal[i] = outY(row,col).real();
            
            opImag[i] = outY(row,col).imag();
        }
    }
    
    
    std::cout << "after assigning output" << std::endl;
    
}


