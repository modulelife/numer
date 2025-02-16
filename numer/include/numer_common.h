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
	inline bool is_pow2(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		return Integer_ != 0 && (Integer_ & (Integer_ - 1)) == 0;
	}

	template<typename Uint>
	inline Uint to_pow2_down(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		if (Integer_ == 0) return 0;
		for (unsigned shift = 1; shift < sizeof(Uint) * 8; shift <<= 1)
			Integer_ |= Integer_ >> shift;
		return (Integer_ + 1) >> 1;
	}

	template<typename Uint>
	inline Uint to_pow2_up(Uint Integer_) {
		static_assert(!std::numeric_limits<Uint>::is_signed, "must pass an unsigned type");
		if (Integer_ == 0) return 1;
		--Integer_;
		for (unsigned shift = 1; shift < sizeof(Uint) * 8; shift <<= 1) {
			Integer_ |= Integer_ >> shift;
		}
		return ++Integer_;
	}


}//namespace numer end
