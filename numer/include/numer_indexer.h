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
		CubeIndexer(CubeContainer& Cube_base_) : base_(Cube_base_) {}

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
			size_t off_;

		public:
			cyclic(size_t Max_id_) : idm_(Max_id_), off_(0) {}
			cyclic(size_t Max_id_, size_t Offset_) : idm_(Max_id_), off_(Offset_) {}
			size_t operator()(size_t Idx_) const {
				return (Idx_ + off_) % idm_;
			}
		};
	}


	template<class Id_Order, class VecContainer>
	class VecReIndexer {
	private:
		VecContainer& base_;
		Id_Order order_;

	public:
		VecReIndexer(VecContainer& Vec_base_, Id_Order&& Order_)
			: base_(Vec_base_), order_(Order_){
		}

		decltype(auto) operator()(size_t Pos_) {
			return base_[order_(Pos_)];
		}

		decltype(auto) operator()(size_t Pos_) const {
			return base_[order_(Pos_)];
		}
	};


	template<class AbstractMat>
	class ImageOrientation {
	private:
		AbstractMat& base_;
		id_order::reverse xorder_;
		id_order::original yorder_;

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


	template<class AbstractCube>
	class DimRearranger {
	private:
		AbstractCube& base_;
		unsigned permutation_[3];

	public:
		DimRearranger(AbstractCube& Mat_base_, unsigned Dim1_id_, unsigned Dim2_id_, unsigned Dim3_id_)
			: base_(Mat_base_) {
			permutation_[0] = Dim1_id_;
			permutation_[1] = Dim2_id_;
			permutation_[2] = Dim3_id_;
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_, size_t Pos3_) {
			size_t indice[3] = { Pos1_, Pos2_, Pos3_ };
			return base_(indice[permutation_[0]], indice[permutation_[1]], indice[permutation_[2]]);
		}

		decltype(auto) operator()(size_t Pos1_, size_t Pos2_, size_t Pos3_) const {
			size_t indice[3] = { Pos1_, Pos2_, Pos3_ };
			return base_(indice[permutation_[0]], indice[permutation_[1]], indice[permutation_[2]]);
		}
	};




}//namespace numer end



