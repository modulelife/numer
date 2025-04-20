//	numer_indexer.h		-*- C++20 -*-
#pragma once
//@brief: 
//
//@classes:
//	
// 
//
//@description:
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <cstddef>




namespace numer {


	template<class VecContainer>
	class VecIndexer {
	private:
		VecContainer& base_;

	public:
		VecIndexer(VecContainer& Vec_base_) : base_(Vec_base_) {}

		decltype(auto) operator()(size_t Pos_) {
			return base_[Pos_];
		}

		decltype(auto) operator()(size_t Pos_) const {
			return base_[Pos_];
		}
	};

	template<class MatContainer>
	class MatIndexer {
	private:
		MatContainer& base_;

	public:
		MatIndexer(MatContainer& Mat_base_) : base_(Mat_base_) {}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_) {
			return base_[Pos1_][Pos2_];
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_) const {
			return base_[Pos1_][Pos2_];
		}
	};

	template<class CubeContainer>
	class CubeIndexer {
	private:
		CubeContainer& base_;

	public:
		CubeIndexer(CubeContainer& Mat_base_) : base_(Mat_base_) {}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_, size_t Pos3_) {
			return base_[Pos1_][Pos2_][Pos3_];
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_, size_t Pos3_) const {
			return base_[Pos1_][Pos2_][Pos3_];
		}
	};


	namespace id_order {

		class original {
		private:
			size_t idm_;

		public:
			original(size_t Max_id_) : idm_(Max_id_) {}
			size_t operator()(size_t Idx_) const {
				return Idx_;
			}
		};

		class reverse {
		private:
			size_t idm_;

		public:
			reverse(size_t Max_id_) : idm_(Max_id_) {}
			size_t operator()(size_t Idx_) const {
				return idm_ - Idx_;
			}
		};

		class cyclic {
		private:
			size_t idm_;

		public:
			cyclic(size_t Max_id_) : idm_(Max_id_) {}
			size_t operator()(size_t Idx_) const {
				return Idx_ % idm_;
			}
		};
	}


	template<class AbstractMat>
	class ImageOrientation {
	private:
		AbstractMat& base_;
		id_order::original xorder_;
		id_order::reverse yorder_;

	public:
		ImageOrientation(AbstractMat& Mat_base_, size_t Rows_, size_t Cols_)
			: base_(Mat_base_), xorder_(Rows_), yorder_(Cols_) {
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_) {
			return base_(yorder_(Pos2_), xorder_(Pos1_));
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_) const {
			return base_(yorder_(Pos2_), xorder_(Pos1_));
		}
	};




}//namespace numer end



