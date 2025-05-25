//	numer_common.h	-*- C_++20 -*-
#pragma once
//@brief: common utilities header
//
//@functions:
//	bool is_pow2(Uint Integer_)
//	Uint to_pow2_down(Uint Integer_)
//	Uint to_pow2_up(Uint Integer_)
//	...in coming
//	
//@description:
//	function like their own names
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <limits>




namespace numer {


	template<typename Uint>
	inline constexpr bool is_pow2(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		return Integer_ != 0 && (Integer_ & (Integer_ - 1)) == 0;
	}

	template<typename Uint>
	inline constexpr Uint to_pow2_down(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		if (Integer_ == 0) return 0;
		for (unsigned shift = 1; shift < sizeof(Uint) * 8; shift <<= 1)
			Integer_ |= Integer_ >> shift;
		return (Integer_ + 1) >> 1;
	}

	template<typename Uint>
	inline constexpr Uint to_pow2_up(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		if (Integer_ == 0) return 1;
		--Integer_;
		for (unsigned shift = 1; shift < sizeof(Uint) * 8; shift <<= 1) {
			Integer_ |= Integer_ >> shift;
		}
		return ++Integer_;
	}

	inline double factorial(unsigned N_) {

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

	template<typename Uint>
	inline constexpr double parity_u(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		if (Integer_ % 2U == 0) return 1.0;
		return -1.0;
	}

}//namespace numer end
