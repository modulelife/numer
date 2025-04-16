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

		

	}//namespace detail end


	class HermiPolyno {
	private:
		const unsigned n_;
	public:
		HermiPolyno(unsigned Order_) : n_(Order_) {}

		double operator()(double X_) const {
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
			coef_ = pow(beta_ * beta_ / Pi, 0.25) / sqrt(pow(2.0, N_) * factorial(N_));
		}

		double operator()(double X_) const {
			return coef_ * exp(-beta_ * beta_ * X_ * X_ / 2.0) * hn_(beta_ * X_);
		}
	};


	class HermiGaussMode {
	private:
		HermiFunc hfx_;
		HermiFunc hfy_;

	public:
		HermiGaussMode(unsigned Nx_, double Beta_x_, unsigned Ny_, double Beta_y_)
			: hfx_(Nx_, Beta_x_), hfy_(Ny_, Beta_y_)
		{}

		double operator()(double X_, double Y_) const {
			return hfx_(X_) * hfy_(Y_);
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

		Complex operator()(double X_) const {
			double x_shifted = X_ - center_;
			return coef_ * Complex::expi(p_over_hbar_ * X_) * exp(-beta_ * beta_ * x_shifted * x_shifted / 2.0);
		}
	};

	class CoherentState2D {
	private:
		CoherentState1D phiax_;
		CoherentState1D phiay_;

	public:
		CoherentState2D(Complex Alpha_x_, double Beta_x_, Complex Alpha_y_, double Beta_y_)
			:phiax_(Alpha_x_, Beta_x_), phiay_(Alpha_y_, Beta_y_)
		{}

		Complex operator()(double X_, double Y_) const {
			return phiax_(X_) * phiay_(Y_);
		}
	};


}//namespace numer end