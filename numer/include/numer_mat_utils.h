//	numer_mat_utils.h		-*- C++20 -*-
#pragma once
//@brief: 2D array container template class mat<Ty, Alloc> header
//
//@classes:
//	
// 
//	
//
//@description:
//	
//
//@usage:
//
//----------------------------content begins----------------------------




namespace numer {
	
	template<typename Ty>
	inline constexpr Ty chooseSecond(const Ty& First, const Ty& Second) {
		return Second;
	}

	template<typename Ty, typename Tx>
	inline constexpr auto multiply(const Ty& First, const Tx& Second) {
		return First * Second;
	}

}//namespace numer end