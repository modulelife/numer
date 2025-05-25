//	numer_grid.h	-*- C_++20 -*-
#pragma once
//@brief: sampler, coordinate transform & indexer header
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

#include <cstddef>
#include <cmath>
#include <utility>




namespace numer {


	struct RangeSpec {
		double start;
		double end;
		size_t length;
	};

	class RangeSampler {
	private:
		double start_;
		double diff_;
		double idm_;

		void setTo_(const RangeSpec& Specifier_) {
			start_ = Specifier_.start;
			diff_ = Specifier_.end - Specifier_.start;
			idm_ = static_cast<double>(Specifier_.length);
		}

	public:
		RangeSampler() : start_(0.0), diff_(0.0), idm_(1.0) {}

		explicit RangeSampler(const RangeSpec& Specifier_) {
			setTo_(Specifier_);
		}

		void setRange(const RangeSpec& Specifier_) {
			setTo_(Specifier_);
		}

		double operator()(size_t Idx_) const {
			if (idm_ == 0.0) return start_;
			double ratio = static_cast<double>(Idx_) / idm_;
			return start_ + diff_ * ratio;
		}

		double step() const {
			return diff_ / idm_;
		}

		bool verifyAndIndex(double Coord_, size_t& Idx_var_) const {
			if(diff_ == 0.0) return false;
			double off = Coord_ - start_;
			if (off < 0.0 || off >= diff_) return false;
			double idx = idm_ * off / diff_;
			Idx_var_ = static_cast<size_t>(idx);
			return true;
		}
	};


	struct Cartes2_To_Polar {
		static double r(double x, double y){
			return sqrt(x * x + y * y);
		}

		static double phi(double x, double y) {
			double result = atan2(y, x);
			return result >= 0.0 ? result : 2.0 * Pi + result;
		}
	private:
		Cartes2_To_Polar() = delete;
	};

	struct Polar_To_Cartes2 {
		static double x(double r, double phi) {
			return r * cos(phi);
		}

		static double y(double r, double phi) {
			return r * sin(phi);
		}
	private:
		Polar_To_Cartes2() = delete;
	};

	struct Cartes3_To_Spheric {
		static double r(double x, double y, double z) {
			return sqrt(x * x + y * y + z * z);
		}

		static double theta(double x, double y, double z) {
			double rxy = sqrt(x * x + y * y);
			return atan2(rxy, z);
		}

		static double phi(double x, double y, double z) {
			double result = atan2(y, x);
			return result >= 0.0 ? result : 2.0 * Pi + result;
		}
	private:
		Cartes3_To_Spheric() = delete;
	};

	struct Spheric_To_Cartes3 {
		static double x(double r, double theta, double phi) {
			return r * sin(theta) *cos(phi);
		}

		static double y(double r, double theta, double phi) {
			return r * sin(theta) * sin(phi);
		}

		static double z(double r, double theta, double phi) {
			return r * cos(theta);
		}
	private:
		Spheric_To_Cartes3() = delete;
	};


	template<class Func>
	class Cartes2Adp {
	private:
		Func func_;

	public:
		Cartes2Adp(const Func& BiFunc_polar_) : func_(BiFunc_polar_) {}

		auto operator()(double Xcoord_, double Ycoord_) -> decltype(func_(std::declval<double>(), std::declval<double>())) const {
			return func_(Cartes2_To_Polar::r(Xcoord_, Ycoord_), Cartes2_To_Polar::phi(Xcoord_, Ycoord_));
		}
	};

	template<class Func>
	class PolarAdp {
	private:
		Func func_;

	public:
		PolarAdp(const Func& BiFunc_cartes2_) : func_(BiFunc_cartes2_) {}

		auto operator()(double Rcoord_, double Phi_coord_) -> decltype(func_(std::declval<double>(), std::declval<double>())) const {
			return func_(Polar_To_Cartes2::x(Rcoord_, Phi_coord_), Polar_To_Cartes2::y(Rcoord_, Phi_coord_));
		}
	};

	template<class Func>
	class Cartes3Adp {
	private:
		Func func_;

	public:
		Cartes3Adp(const Func& TrFunc_spheric_) : func_(TrFunc_spheric_) {}

		auto operator()(double Xcoord_, double Ycoord_, double Zcoord_) -> decltype(func_(std::declval<double>(), std::declval<double>(), std::declval<double>())) const {
			return func_(Cartes3_To_Spheric::r(Xcoord_, Ycoord_, Zcoord_), Cartes3_To_Spheric::theta(Xcoord_, Ycoord_, Zcoord_), Cartes3_To_Spheric::phi(Xcoord_, Ycoord_, Zcoord_));
		}
	};

	template<class Func>
	class SphericAdp {
	private:
		Func func_;

	public:
		SphericAdp(const Func& TrFunc_cartes3_) : func_(TrFunc_cartes3_) {}

		auto operator()(double Rcoord_, double Theta_coord_, double Phi_coord_) -> decltype(func_(std::declval<double>(), std::declval<double>(), std::declval<double>())) const {
			return func_(Spheric_To_Cartes3::x(Rcoord_, Theta_coord_, Phi_coord_), Spheric_To_Cartes3::y(Rcoord_, Theta_coord_, Phi_coord_), Spheric_To_Cartes3::z(Rcoord_, Theta_coord_, Phi_coord_));
		}
	};




	template<typename Func>
	class UnaryFuncSampler {
	private:
		Func func_;
		RangeSampler to_domain_;

	public:
		UnaryFuncSampler(const Func& UFunc_, const RangeSpec& Specifier_)
			:func_(UFunc_), to_domain_(Specifier_)
		{}

		auto operator()(size_t Idx_) -> decltype(func_(std::declval<double>())) const {
			return func_(to_domain_(Idx_));
		}

		double stepOfDim(unsigned Which_dim_) const { 
			if (Which_dim_ == 0) return to_domain_.step();
			return 0.0;
		}

		double diffElem() const {
			return to_domain_.step();
		}
	};

	template<typename Func>
	class BinaryFuncSampler {
	private:
		Func func_;
		RangeSampler to_domain1_;
		RangeSampler to_domain2_;

	public:
		BinaryFuncSampler(const Func& BiFunc_, const RangeSpec& Specifier1_, const RangeSpec& Specifier2_)
			:func_(BiFunc_), to_domain1_(Specifier1_), to_domain2_(Specifier2_)
		{}

		auto operator()(size_t Idx1_, size_t Idx2_) -> decltype(func_(std::declval<double>(), std::declval<double>())) const {
			return func_(to_domain1_(Idx1_), to_domain2_(Idx2_));
		}

		double stepOfDim(unsigned Which_dim_) const {
			if (Which_dim_ == 0) return to_domain1_.step();
			if (Which_dim_ == 1) return to_domain2_.step();
			return 0.0;
		}

		double diffElem() const {
			return to_domain1_.step() * to_domain2_.step();
		}
	};

	template<typename Func>
	class TernaryFuncSampler {
	private:
		Func func_;
		RangeSampler to_domain1_;
		RangeSampler to_domain2_;
		RangeSampler to_domain3_;

	public:
		TernaryFuncSampler(const Func& TrFunc_, const RangeSpec& Specifier1_, const RangeSpec& Specifier2_, const RangeSpec& Specifier3_)
			:func_(TrFunc_), to_domain1_(Specifier1_), to_domain2_(Specifier2_), to_domain3_(Specifier3_)
		{}

		auto operator()(size_t Idx1_, size_t Idx2_, size_t Idx3_) -> decltype(func_(std::declval<double>(), std::declval<double>(), std::declval<double>())) const {
			return func_(to_domain1_(Idx1_), to_domain2_(Idx2_), to_domain3_(Idx3_));
		}

		double stepOfDim(unsigned Which_dim_) const {
			if (Which_dim_ == 0) return to_domain1_.step();
			if (Which_dim_ == 1) return to_domain2_.step();
			if (Which_dim_ == 2) return to_domain3_.step();
			return 0.0;
		}

		double diffElem() const {
			return to_domain1_.step() * to_domain2_.step() * to_domain3_.step();
		}
	};



}//namespace numer end