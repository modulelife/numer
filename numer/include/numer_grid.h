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
		size_t max_index;
	};

	class RangeSampler {
	private:
		double start_;
		double diff_;
		double idm_;

	public:
		explicit RangeSampler(const RangeSpec& Specifier_) {

			if (Specifier_.max_index == 0) {
				start_ = (Specifier_.start + Specifier_.end) / 2.0;
				diff_ = 0.0;
				idm_ = 1.0;
			}
			else {
				start_ = Specifier_.start;
				diff_ = Specifier_.end - Specifier_.start;
				idm_ = static_cast<double>(Specifier_.max_index);
			}
		}

		double operator()(size_t Idx_) const {
			double ratio = static_cast<double>(Idx_) / idm_;
			return start_ + diff_ * ratio;
		}

		double step() const {
			return diff_ / idm_;
		}

		bool verifyAndIndex(double Coord_, size_t& Idx_Var_) const {
			if(diff_ == 0.0) return false;
			double off = Coord_ - start_;
			if (off < 0.0 || off > diff_) return false;
			double idx = idm_ * off / diff_;
			Idx_Var_ = static_cast<size_t>(idx + 0.5);
			return true;
		}
	};


	template<typename Func>
	class UnaryFuncSampler {
	private:
		Func func_;
		RangeSampler to_domain_;

	public:
		UnaryFuncSampler(Func UFunc_, RangeSpec Specifier_)
			:func_(UFunc_), to_domain_(Specifier_)
		{}

		auto operator()(size_t Idx_) -> decltype(func_(std::declval<double>())) const {
			return func_(to_domain_(Idx_));
		}

		double stepOfDim(unsigned Which_Dim_) const { 
			if (Which_Dim_ == 0) return to_domain_.step();
			return 0.0;
		}

		double lineElem() const {
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
		BinaryFuncSampler(Func BiFunc_, RangeSpec Specifier1_, RangeSpec Specifier2_)
			:func_(BiFunc_), to_domain1_(Specifier1_), to_domain2_(Specifier2_)
		{}

		auto operator()(size_t Idx1_, size_t Idx2_) -> decltype(func_(std::declval<double>(), std::declval<double>())) const {
			return func_(to_domain1_(Idx1_), to_domain2_(Idx2_));
		}

		double stepOfDim(unsigned Which_Dim_) const {
			if (Which_Dim_ == 0) return to_domain1_.step();
			if (Which_Dim_ == 1) return to_domain2_.step();
			return 0.0;
		}

		double panlElem() const {
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
		TernaryFuncSampler(Func TrFunc_, RangeSpec Specifier1_, RangeSpec Specifier2_, RangeSpec Specifier3_)
			:func_(TrFunc_), to_domain1_(Specifier1_), to_domain2_(Specifier2_), to_domain3_(Specifier3_)
		{}

		auto operator()(size_t Idx1_, size_t Idx2_, size_t Idx3_) -> decltype(func_(std::declval<double>(), std::declval<double>(), std::declval<double>())) const {
			return func_(to_domain1_(Idx1_), to_domain2_(Idx2_), to_domain3_(Idx3_));
		}

		double stepOfDim(unsigned Which_Dim_) const {
			if (Which_Dim_ == 0) return to_domain1_.step();
			if (Which_Dim_ == 1) return to_domain2_.step();
			if (Which_Dim_ == 2) return to_domain3_.step();
			return 0.0;
		}

		double volElem() const {
			return to_domain1_.step() * to_domain2_.step() * to_domain3_.step();
		}
	};



}//namespace numer end