//	numer_complex.h	-*- C_++20 -*-
#pragma once
//@brief: Complex number class header
//
//@classes:
//	numer::Complex
//								------------ " Complex number class, double accuracy"
//	
//@description:
//	Complex number & basic arithmetic & primary functions
//	any initial phase beyond [0, 2Pi) will not be reserved
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <cmath>
#include <numer_constant.h>




namespace numer {

	
	//Complex number type
	//any initial phase beyond [0, 2Pi) will not be reserved
	class Complex {
	private:
		double re_{ 0.0 };
		double im_{ 0.0 };

	public:
		constexpr Complex() {}
		constexpr Complex(double Real_) : re_(Real_) {}
		constexpr Complex(double Real_, double Imag_) : re_(Real_), im_(Imag_) {}
		

		//creat imaginary unit
		static constexpr Complex i() {
			return Complex(0.0, 1.0);
		}
		//creat Complex number exp[i*phase]
		static Complex expi(double Phase_rad_) {
			return Complex(cos(Phase_rad_), sin(Phase_rad_));
		}
		//creat Complex number 1+0i
		static constexpr Complex identity() {
			return Complex(1.0, 0.0);
		}
		//creat Complex number 0+0i
		static constexpr Complex zero() {
			return Complex();
		}

		void setRe(double Real_) { re_ = Real_; }
		void setIm(double Imag_) { im_ = Imag_; }
		constexpr double re() const { return re_; }
		constexpr double im() const { return im_; }
		//a.k.a. modulus lenght
		double amplitude() const { return sqrt(re_ * re_ + im_ * im_); }
		//a.k.a. squared modulus lenght
		constexpr double sqrdAmp() const { return re_ * re_ + im_ * im_; }
		//a.k.a. argument, in range [0, 2Pi)
		double phase() const {
			double phase = atan2(im_, re_);
			return phase >= 0.0 ? phase : 2.0 * Pi + phase;
		}
		//Complex conjugate
		constexpr Complex conj() const {
			return Complex(re_, -im_);
		}

		constexpr Complex operator-() const {
			return Complex(-re_, -im_);
		}

		constexpr Complex& operator+=(const Complex& Right_) {
			re_ += Right_.re();
			im_ += Right_.im();
			return *this;
		}

		constexpr Complex& operator-=(const Complex& Right_) {
			re_ -= Right_.re();
			im_ -= Right_.im();
			return *this;
		}

		constexpr Complex& operator*=(const Complex& Right_) {
			double re = re_ * Right_.re() - im_ * Right_.im();
			double im = re_ * Right_.im() + im_ * Right_.re();
			re_ = re;
			im_ = im;
			return *this;
		}

		constexpr Complex& operator*=(const double& Right_) {
			re_ *= Right_;
			im_ *= Right_;
			return *this;
		}

		constexpr Complex& operator/=(const Complex& Right_) {
			double re = re_ * Right_.re() + im_ * Right_.im();
			double im = -re_ * Right_.im() + im_ * Right_.re();
			double denom = Right_.re() * Right_.re() + Right_.im() * Right_.im();
			re_ = re / denom;
			im_ = im / denom;
			return *this;
		}

		constexpr Complex& operator/=(const double& Right_) {
			re_ /= Right_;
			im_ /= Right_;
			return *this;
		}
	};


	inline constexpr bool operator==(const Complex& Left_, const Complex& Right_);
	inline constexpr bool operator!=(const Complex& Left_, const Complex& Right_);
	inline constexpr Complex operator+(const Complex& Left_, const Complex& Right_);
	inline constexpr Complex operator-(const Complex& Left_, const Complex& Right_);
	inline bool approx_eq(const Complex& Left_, const Complex& Right_, double Tolerance_);
	inline bool approx_neq(const Complex& Left_, const Complex& Right_, double Tolerance_);
	inline constexpr Complex operator*(const Complex& Left_, const Complex& Right_);
	inline constexpr Complex operator*(const double& Left_, const Complex& Right_);
	inline constexpr Complex operator*(const Complex& Left_, const double& Right_);
	inline constexpr Complex operator/(const Complex& Left_, const Complex& Right_);
	inline constexpr Complex operator/(const double& Left_, const Complex& Right_);
	inline constexpr Complex operator/(const Complex& Left_, const double& Right_);


}//namespace numer end




inline constexpr
bool numer::operator==(const numer::Complex& Left_, const numer::Complex& Right_) {
	return Left_.re() == Right_.re() && Left_.im() == Right_.im();
}

inline constexpr
bool numer::operator!=(const numer::Complex& Left_, const numer::Complex& Right_) {
	return Left_.re() != Right_.re() || Left_.im() != Right_.im();
}

inline constexpr
numer::Complex numer::operator+(const numer::Complex& Left_, const numer::Complex& Right_) {
	return numer::Complex(Left_.re() + Right_.re(), Left_.im() + Right_.im());
}

inline constexpr
numer::Complex numer::operator-(const numer::Complex& Left_, const numer::Complex& Right_) {
	return numer::Complex(Left_.re() - Right_.re(), Left_.im() - Right_.im());
}

inline
bool numer::approx_eq(const numer::Complex& Left_, const numer::Complex& Right_, double Tolerance_) {
	return (Left_ - Right_).amplitude() < Tolerance_;
}

inline
bool numer::approx_neq(const numer::Complex& Left_, const numer::Complex& Right_, double Tolerance_) {
	return (Left_ - Right_).amplitude() >= Tolerance_;
}

inline constexpr
numer::Complex numer::operator*(const numer::Complex& Left_, const numer::Complex& Right_) {
	double re = Left_.re() * Right_.re() - Left_.im() * Right_.im();
	double im = Left_.re() * Right_.im() + Left_.im() * Right_.re();
	return numer::Complex(re, im);
}

inline constexpr
numer::Complex numer::operator*(const double& Left_, const numer::Complex& Right_) {
	return numer::Complex(Right_.re() * Left_, Right_.im() * Left_);
}

inline constexpr
numer::Complex numer::operator*(const numer::Complex& Left_, const double& Right_) {
	return numer::Complex(Left_.re() * Right_, Left_.im() * Right_);
}

inline constexpr
numer::Complex numer::operator/(const numer::Complex& Left_, const numer::Complex& Right_) {
	double re = Left_.re() * Right_.re() + Left_.im() * Right_.im();
	double im = -Left_.re() * Right_.im() + Left_.im() * Right_.re();
	double denom = Right_.re() * Right_.re() + Right_.im() * Right_.im();
	return numer::Complex(re / denom, im / denom);
}

inline constexpr
numer::Complex numer::operator/(const double& Left_, const numer::Complex& Right_) {
	double re = Left_ * Right_.re();
	double im = -Left_ * Right_.im();
	double denom = Right_.re() * Right_.re() + Right_.im() * Right_.im();
	return numer::Complex(re / denom, im / denom);
}

inline constexpr
numer::Complex numer::operator/(const numer::Complex& Left_, const double& Right_) {
	return numer::Complex(Left_.re() / Right_, Left_.im() / Right_);
}


inline constexpr double Re(const numer::Complex& X_) {
	return X_.re();
}

inline constexpr double Im(const numer::Complex& X_) {
	return X_.im();
}

inline double Arg(const numer::Complex& X_) {
	return X_.phase();
}

inline constexpr numer::Complex conj(const numer::Complex& X_) {
	return X_.conj();
}

inline double abs(const numer::Complex& X_) {
	return X_.amplitude();
}

inline double norm(const numer::Complex& X_) {
	return X_.sqrdAmp();
}

inline numer::Complex exp(const numer::Complex& X_) {
	return numer::Complex::expi(X_.im()) * exp(X_.re());
}

inline numer::Complex pow(const numer::Complex& X_, double Y_) {
	return numer::Complex::expi(X_.phase() * Y_) * pow(X_.amplitude(), Y_);
}

inline numer::Complex sin(const numer::Complex& X_) {
	constexpr numer::Complex i = numer::Complex::i();
	return (exp(X_ * i) - exp(-X_ * i)) / (2 * i);
}

inline numer::Complex cos(const numer::Complex& X_) {
	constexpr numer::Complex i = numer::Complex::i();
	return (exp(X_ * i) + exp(-X_ * i)) / 2;
}

inline numer::Complex tan(const numer::Complex& X_) {
	return sin(X_) / cos(X_);
}

inline numer::Complex sinh(const numer::Complex& X_) {
	return (exp(X_) - exp(-X_)) / 2;
}

inline numer::Complex cosh(const numer::Complex& X_) {
	return (exp(X_) + exp(-X_)) / 2;
}

inline numer::Complex tanh(const numer::Complex& X_) {
	return sinh(X_) / cosh(X_);
}