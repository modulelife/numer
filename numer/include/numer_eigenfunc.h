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
#include <exception>
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
			for (unsigned i = 2; i <= n_; ++i) {
				double n = static_cast<double>(i);
				hn = 2.0 * (X_ * hn_1 - (n - 1.0) * hn_2);
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

	class CoherentState3D {
	private:
		CoherentState1D phiax_;
		CoherentState1D phiay_;
		CoherentState1D phiaz_;

	public:
		CoherentState3D(Complex Alpha_x_, double Beta_x_, Complex Alpha_y_, double Beta_y_, Complex Alpha_z_, double Beta_z_)
			:phiax_(Alpha_x_, Beta_x_), phiay_(Alpha_y_, Beta_y_), phiaz_(Alpha_z_, Beta_z_)
		{
		}

		Complex operator()(double X_, double Y_, double Z_) const {
			return phiax_(X_) * phiay_(Y_) * phiaz_(Z_);
		}
	};


	class AssoLaguerrePolyno {
	private:
		const double a_;
		const unsigned n_;

	public:
		AssoLaguerrePolyno(unsigned Alpha_, unsigned Order_) : a_(static_cast<double>(Alpha_)), n_(Order_) {}

		double operator()(double X_) const {
			if (n_ == 0) return 1.0;
			if (n_ == 1) return a_ + 1.0 - X_;

			double lak_2 = 1.0, lak_1 = a_ + 1.0 - X_, lak;
			for (unsigned i = 2; i <= n_; ++i) {
				double k = static_cast<double>(i);
				lak = ((2.0 * k + a_ - 1.0 - X_) * lak_1 - (k - 1.0 + a_) * lak_2) / k;
				lak_2 = lak_1;
				lak_1 = lak;
			}
			return lak;
		}
	};


	class HydrogenRadical {
	private:
		AssoLaguerrePolyno lan_;
		unsigned l_;
		double scale_;
		double coef_;

	public:
		HydrogenRadical(unsigned N_, unsigned L_, double Radius_a0_) 
			: lan_(2 * L_ + 1, N_ - L_ - 1), l_(L_) {
			if (N_ == 0) throw std::invalid_argument("energy level must be at least 1.");
			if (N_ <= L_) throw std::invalid_argument("energy level must be greater than angular quantum number.");

			double n = static_cast<double>(N_);
			scale_ = 1.0 / (n * Radius_a0_);
			coef_ = sqrt(8.0 * scale_ * scale_ * scale_ / (2.0 * n) * (factorial(N_ - L_ - 1) / factorial(N_ + L_)));
		}

		double operator()(double R_) const {
			return coef_ * exp(-R_ * scale_) * pow(2.0 * R_ * scale_, l_) * lan_(2.0 * R_ * scale_);
		}
	};


	class LegendrePolyno {
	private:
		const unsigned l_;

	public:
		LegendrePolyno(unsigned L_) : l_(L_) {}

		double operator()(double X_) const {
			if (l_ == 0) return 1.0;
			if (l_ == 1) return X_;

			double pl_2 = 1.0, pl_1 = X_, pl;
			for (unsigned i = 2; i <= l_; ++i) {
				double l = static_cast<double>(i);
				pl = ((2.0 * l - 1.0) * X_ * pl_1 - (l - 1.0) * pl_2) / l;
				pl_2 = pl_1;
				pl_1 = pl;
			}
			return pl;
		}
	};

	class AssoLegendrePolyno {
	private:
		LegendrePolyno pl_;
		const unsigned l_;
		const unsigned m_;


	public:
		AssoLegendrePolyno(unsigned L_, unsigned M_)
			: pl_(L_), l_(L_), m_(M_) {
			if (L_ < M_) throw std::invalid_argument("L must be greater than or equal to M for associated Legendre function.");
		}

		double operator()(double X_) const {
			constexpr double epsilon = std::numeric_limits<double>::epsilon() * 10;

			if (std::abs(X_) >= 1.0 - epsilon) {
				return (m_ == 0) ? pl_(X_) : 0.0;  // m>0 Ê±·µ»Ø0
			}

			double plm_2 = pl_(X_);
			if (m_ == 0) return plm_2;

			const double l = static_cast<double>(l_);
			LegendrePolyno pl_1(l_ - 1);
			double denom = sqrt(1.0 - X_ * X_);
			double plm_1 = l * (X_ * plm_2 - pl_1(X_)) / denom;
			if (m_ == 1) return plm_1;

			double plm;
			for (unsigned i = 2; i <= m_; ++i) {
				double m = static_cast<double>(i);
				plm = -(l + m - 1.0) * (l - m + 2.0) * plm_2 - 2.0 * (m - 1.0) * X_ * plm_1 / denom;
				plm_2 = plm_1;
				plm_1 = plm;
			}
			return plm;
		}
	};

	class SphericalHarmonic {
	private:
		AssoLegendrePolyno plm_;
		const double m_;
		double coef_;

	public:
		SphericalHarmonic(unsigned L_, int M_)
			:plm_(L_, abs(M_)), m_(static_cast<double>(M_)) {
			double l = static_cast<double>(L_);
			unsigned m = abs(M_);
			coef_ = sqrt((2.0 * l + 1.0) / (4.0 * Pi) * (factorial(L_ - m) / factorial(L_ + m)));
			if (M_ >= 0) {
				coef_ *= parity_u(m);
			}
		}

		Complex operator()(double Theta_, double Phi_) const {
			return coef_ * plm_(cos(Theta_)) * Complex::expi(m_ * Phi_);
		}
	};


	class HydrogenState {
	private:
		HydrogenRadical Rnl_;
		SphericalHarmonic Ylm_;

	public:
		HydrogenState(unsigned N_, unsigned L_, int M_, double Radius_a0_)
			: Rnl_(N_, L_, Radius_a0_), Ylm_(L_, M_) {
		}

		Complex operator()(double R_, double Theta_, double Phi_) const {
			return Rnl_(R_) * Ylm_(Theta_, Phi_);
		}
	};






}//namespace numer end