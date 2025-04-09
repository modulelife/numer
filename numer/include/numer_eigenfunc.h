//	numer_eigenfunc.h	-*- C_++20 -*-
#pragma once
//@brief: quantum mechanics eigenfunction header
//
//@functions:
//	
//	
//@description:
//	
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <cmath>
#include <numer_complex.h>
#include <numer_common.h>




namespace numer {


	namespace detail_ {

		inline double factorial__(unsigned N_) {

			static const double kFactorials[21] = {
				1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.,
				39916800., 479001600., 6227020800., 87178291200., 1307674368000.,
				20922789888000., 355687428096000., 6402373705728000.,
				121645100408832000., 2432902008176640000.
			};

			if (N_ <= 20) return kFactorials[N_];
			double x = static_cast<double>(N_);
			return exp(0.5 * log(2.0 * Pi * x) + x * log(x) - x + 1.0 / (12.0 * x));
		}

	}//namespace detail end


	class HermiPolyno {
	private:
		const unsigned n_;
	public:
		HermiPolyno(unsigned Order_) : n_(Order_) {}

		double operator()(double X_) {
			if (n_ == 0) return 1.0;
			if (n_ == 1) return 2.0 * X_;

			double hn_2 = 1.0, hn_1 = 2.0 * X_, hn;
			for (unsigned n = 2; n <= n_; ++n) {
				hn = 2.0 * (X_ * hn_1 - (n - 1) * hn_2);
				hn_2 = hn_1;
				hn_1 = hn;
			}
			return hn;
		}
	};


	class HermiFunc {
	private:
		HermiPolyno hn_;
		double beta_;
		double coef_;
	public:
		HermiFunc(unsigned N_, double Beta_)
			: hn_(N_), beta_(Beta_)
		{
			coef_ = pow(beta_ * beta_ / Pi, 0.25) / sqrt(pow(2.0, N_) * detail_::factorial__(N_));
		}

		double operator()(double X_) {
			return coef_ * exp(-beta_ * beta_ * X_ * X_ / 2.0) * hn_(beta_ * X_);
		}
	};


	class CoherentState1D {
	private:
		Complex alpha_;
		double beta_;
		Complex coef_;
		double p_over_hbar_;
		double center_;
	public:
		CoherentState1D(Complex Alpha_, double Beta_)
			:alpha_(Alpha_), beta_(Beta_) 
		{
			coef_ = pow(beta_ * beta_ / Pi, 0.25) * exp((Alpha_.conj() * Alpha_.conj() - Alpha_ * Alpha_) / 4.0);
			p_over_hbar_ = sqrt(2.0) * beta_ * alpha_.im();
			center_ = sqrt(2.0) * alpha_.re() / beta_;
		}

		Complex operator()(double X_) {
			double x_shifted = X_ - center_;
			return coef_ * Complex::expi(p_over_hbar_ * X_) * exp(-beta_ * beta_ * x_shifted * x_shifted / 2.0);
		}
	};


}//namespace numer end