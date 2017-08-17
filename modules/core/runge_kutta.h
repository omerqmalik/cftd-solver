#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

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
        
        // Uncomment for debugging
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

struct TDSSolver {
public:
	virtual ~TDSSolver() {}

	virtual void operator()(const state_type& y, state_type& dy, const double t) = 0;
};

// ODE system solver
struct TDSSolvers_RING : public TDSSolver {

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

	virtual ~TDSSolvers_RING() {}

	void operator()(const state_type& y, state_type& dy, const double t)
	{
		int N = (int) params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());

		bnu::matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
		for (int i = 0; i < N * N; i++) {
			for (int j = 0; j < N * N; j++) {
				size_t a = j / N;
				size_t b = j % N;
				NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
			}
		}

		ustate_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
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
struct TDSSolvers_UCF : public TDSSolver {

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

	virtual ~TDSSolvers_UCF() {}

	void operator()(const state_type& y, state_type& dy, const double t)
	{
		int N = (int) params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());
		bnu::matrix<std::complex<double>> Q = outer_prod(conj(b_vec1), a_vec);

		NL_term.clear();
		for (int i = 0; i < N * N; i++) {
			for (int j = 0; j < N * N; j++) {
				size_t a = j / N;
				size_t b = j % N;
				NL_term[i] += params.A(i, j) * (Q(a, b) - Q(b, a));
			}
		}

		ustate_type a_vecT_D_mat(N);
		ustate_type b_vecT_B(N);
		a_vecT_D_mat.clear();
		b_vecT_B.clear();

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
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
struct TDSSolvers_FP : public TDSSolver {

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

	virtual ~TDSSolvers_FP() {}

	void operator()(const state_type& y, state_type& dy, const double t)
	{
		int N = (int) params.N;
		state_type::const_iterator slice1 = y.begin() + N;
		state_type::const_iterator slice2 = slice1 + N;

		std::copy(y.begin(), slice1, a_vec.begin());
		std::copy(slice1, slice2, b_vec1.begin());
		std::copy(slice2, y.end(), D_vec.begin());

		bnu::matrix<std::complex<double>> Q = imag(outer_prod(conj(b_vec1), a_vec));

		NL_term.clear();
		for (int i = 0; i < N * N; i++) {
			for (int j = 0; j < N * N; j++) {
				size_t a = j / N;
				size_t b = j % N;
				NL_term[i] += params.A(i, j) * Q(a, b);
			}
		}

		ustate_type a_vecTD_mat(N); // = (a_vec.' * D_mat).'
		a_vecTD_mat.clear();

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
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

inline std::unique_ptr<TDSSolver> solverFactory(std::string& type, Parameters& p)
{
	if (type.compare("RING") == 0) {
		return std::unique_ptr<TDSSolvers_RING>(new TDSSolvers_RING(p));
	}
	else if (type.compare("UCF") == 0) {
		return std::unique_ptr<TDSSolvers_UCF>(new TDSSolvers_UCF(p));
	}
	else if (type.compare("FP") == 0) {
		return std::unique_ptr<TDSSolvers_FP>(new TDSSolvers_FP(p));
	}
	else {
		throw std::runtime_error("Unknown solver: " + type);
	}
}

struct OdeIntegrator {

public:
	OdeIntegrator(bool ode_const) : use_const(ode_const), abs_error(1.03 - 6), rel_error(1.0e-3)
	{
	}

	OdeIntegrator(bool ode_const, double aerr, double rerr) : use_const(ode_const), abs_error(aerr), rel_error(rerr)
	{
	}

	void integrate(TDSSolver& solver, state_type initial_state, double start, double end, double dt, Observer& observer)
	{
		auto rkd = runge_kutta_dopri5<state_type>{};
		auto stepper = make_dense_output(abs_error, rel_error, rkd);

		if (use_const) {
			integrate_const(boost::ref(stepper), boost::ref(solver), initial_state, start, end, dt, boost::ref(observer));
		}
		else {
			integrate_adaptive(boost::ref(stepper), boost::ref(solver), initial_state, start, end, dt, boost::ref(observer));
		}
	}

private:
	bool use_const;
	double abs_error;
	double rel_error;
};

#endif /* RUNGE_KUTTA_H */