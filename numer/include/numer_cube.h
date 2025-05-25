//	numer_cube.h		-*- C++20 -*-
#pragma once
//@brief: 3D array container template class Cube<Ty, Alloc> header
//
//@classes:
//	numer::Cube<Ty, Alloc>
//								------------ " container template class "
// 
//	
// 
//
//@description:
//
//@usage:
//
//----------------------------content begins----------------------------

#include <cstddef>
#include <numer_mat.h>
#include <vector>
#include <functional>




namespace numer {


	template<class MyCube_>
	class Cube_Depth_iterator_ {
	public:
		friend MyCube_;

		using iterator_concept						= std::contiguous_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyCube_::value_type;
		using difference_type						= typename MyCube_::difference_type;
		using pointer								= typename MyCube_::pointer;
		using reference								= typename MyCube_::reference;

		using Base_iter_							= typename MyCube_::Mat_Vec_iter_;

		using self									= Cube_Depth_iterator_;

	private:
		Base_iter_			base_iter_;
		size_t				row_;
		size_t				col_;


		Cube_Depth_iterator_(Base_iter_ Iter_, size_t Row_, size_t Col_)
			: base_iter_(Iter_), row_(Row_), col_(Col_) {
		}

	public:
		Cube_Depth_iterator_() : base_iter_(nullptr), row_(0), col_(0) {}

		reference operator*() const noexcept {
			return *base_iter_[row_][col_];
		}

		pointer operator->() const noexcept {
			return &(*base_iter_[row_][col_]);
		}

		self& operator++() noexcept {
			++base_iter_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++base_iter_;
			return tmp;
		}

		self& operator--() noexcept {
			--base_iter_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--base_iter_;
			return tmp;
		}

		self& operator+=(const difference_type Off_) noexcept {
			base_iter_ += Off_;
			return *this;
		}

		self operator+(const difference_type Off_) const noexcept {
			return self(base_iter_ + Off_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(const difference_type Off_) noexcept {
			base_iter_ -= Off_;
			return *this;
		}

		self operator-(const difference_type Off_) const noexcept {
			return self(base_iter_ - Off_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(base_iter_ - Right_.base_iter_);
		}

		reference operator[](const difference_type Off_) const noexcept {
			return base_iter_[Off_][row_][col_];
		}

		bool operator==(const self& Right_) const noexcept {
			return base_iter_ == Right_.base_iter_;
		}

		bool operator!=(const self& Right_) const noexcept {
			return base_iter_ != Right_.base_iter_;
		}

		bool operator<(const self& Right_) const noexcept {
			return base_iter_ < Right_.base_iter_;
		}

		bool operator<=(const self& Right_) const noexcept {
			return base_iter_ <= Right_.base_iter_;
		}

		bool operator>(const self& Right_) const noexcept {
			return base_iter_ > Right_.base_iter_;
		}

		bool operator>=(const self& Right_) const noexcept {
			return base_iter_ >= Right_.base_iter_;
		}

	};



	template<class MyCube_>
	class Cube_Depth_const_iterator_ {
	public:
		friend MyCube_;

		using iterator_concept						= std::contiguous_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyCube_::value_type;
		using difference_type						= typename MyCube_::difference_type;
		using pointer								= typename MyCube_::const_pointer;
		using reference								= typename MyCube_::const_reference;

		using Base_iter_							= typename MyCube_::Mat_Vec_citer_;

		using self									= Cube_Depth_const_iterator_;

	private:
		Base_iter_			base_iter_;
		size_t				row_;
		size_t				col_;

		Cube_Depth_const_iterator_(Base_iter_ Iter_, size_t Row_, size_t Col_)
			: base_iter_(Iter_), row_(Row_), col_(Col_) {
		}

	public:
		Cube_Depth_const_iterator_() : base_iter_(nullptr), row_(0), col_(0) {}

		reference operator*() const noexcept {
			return *base_iter_[row_][col_];
		}

		pointer operator->() const noexcept {
			return &(*base_iter_[row_][col_]);
		}

		self& operator++() noexcept {
			++base_iter_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++base_iter_;
			return tmp;
		}

		self& operator--() noexcept {
			--base_iter_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--base_iter_;
			return tmp;
		}

		self& operator+=(const difference_type Off_) noexcept {
			base_iter_ += Off_;
			return *this;
		}

		self operator+(const difference_type Off_) const noexcept {
			return self(base_iter_ + Off_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(const difference_type Off_) noexcept {
			base_iter_ -= Off_;
			return *this;
		}

		self operator-(const difference_type Off_) const noexcept {
			return self(base_iter_ - Off_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(base_iter_ - Right_.base_iter_);
		}

		reference operator[](const difference_type Off_) const noexcept {
			return base_iter_[Off_][row_][col_];
		}

		bool operator==(const self& Right_) const noexcept {
			return base_iter_ == Right_.base_iter_;
		}
		
		bool operator!=(const self& Right_) const noexcept {
			return base_iter_ != Right_.base_iter_;
		}

		bool operator<(const self & Right_) const noexcept {
			return base_iter_ < Right_.base_iter_;
		}

		bool operator<=(const self& Right_) const noexcept {
			return base_iter_ <= Right_.base_iter_;
		}

		bool operator>(const self& Right_) const noexcept {
			return base_iter_ > Right_.base_iter_;
		}

		bool operator>=(const self& Right_) const noexcept {
			return base_iter_ >= Right_.base_iter_;
		}

	};


	template<class MyCube_>
	class Cube_Layer_const_proxy_ {
	public:
		friend MyCube_;
		using Mat									= typename MyCube_::Layer_Mat_T_;

	private:
		const Mat&			layer_;

		Cube_Layer_const_proxy_(const Mat& Layer_) : layer_(Layer_) {}


	public:
		decltype(auto) operator[](size_t Row_id_) const noexcept
		{
			return layer_[Row_id_];
		}

		//careful with this
		const Mat& get() const noexcept {
			return layer_;
		}

	};


	template<class MyCube_>
	class Cube_Layer_proxy_ {
	public:
		friend MyCube_;
		using Mat									= typename MyCube_::Layer_Mat_T_;

	private:
		Mat&				layer_;

		Cube_Layer_proxy_(Mat& Layer_) : layer_(Layer_) {}


	public:
		decltype(auto) operator[](size_t Row_id_) const noexcept
		{
			return layer_[Row_id_];
		}

		Mat& get() noexcept {
			return layer_;
		}

		const Mat& get() const noexcept {
			return layer_;
		}

	};


	template<typename Ty, class Alloc = std::allocator<Ty>>
	class cube final {
	public:
		using value_type							= Ty;
		using allocator_type						= Alloc;
		using pointer								= Ty*;
		using const_pointer							= const Ty*;
		using reference								= Ty&;
		using const_reference						= const Ty&;
		using size_type								= size_t;
		using difference_type						= ptrdiff_t;

		using Layer_Mat_T_							= mat<Ty, Alloc>;
		using Mat_Vec_								= std::vector<Layer_Mat_T_>;
		using Mat_Vec_iter_							= typename Mat_Vec_::iterator;
		using Mat_Vec_citer_						= typename Mat_Vec_::const_iterator;

		using layer_proxy							= Cube_Layer_proxy_<cube>;
		using layer_const_proxy						= Cube_Layer_const_proxy_<cube>;
		using depth_iterator						= Cube_Depth_iterator_<cube>;
		using depth_const_iterator					= Cube_Depth_const_iterator_<cube>;

	private:
		size_t			depth_;
		size_t			height_;
		size_t			width_;
		Mat_Vec_		data_;

	public:

		cube() : depth_(0), height_(0), width_(0), data_() {}
		cube(size_t Depth_, size_t Height_, size_t Width_)
			: depth_(Depth_), height_(Height_), width_(Width_), data_(Depth_)
		{
			for (size_t i = 0; i < Depth_; ++i) {
				Layer_Mat_T_ layer(Height_, Width_);
				data_[i] = std::move(layer);
			}
		}
		cube(size_t Depth_, size_t Height_, size_t Width_, const Ty& Value_)
			: depth_(Depth_), height_(Height_), width_(Width_), data_(Depth_)
		{
			for (size_t i = 0; i < Depth_; ++i) {
				Layer_Mat_T_ layer(Height_, Width_, Value_);
				data_[i] = std::move(layer);
			}
		}

		cube(const cube&) = default;
		cube& operator=(const cube&) = default;
		cube(cube&&) = default;
		cube& operator=(cube&&) = default;

		void clear() {
			data_.clear();
		}

		size_type depth() const { return depth_; }
		size_type height() const { return height_; }
		size_type width() const { return width_; }

		bool empty() const { return data_.empty(); }

		layer_proxy operator[](size_t Layer_id_)
		{
			return layer_proxy(data_[Layer_id_]);
		}

		layer_const_proxy operator[](size_t Layer_id_) const
		{
			return layer_const_proxy(data_[Layer_id_]);
		}

		depth_iterator begin_at(size_t Row_, size_t Col_) {
			return depth_iterator(data_.begin(), Row_, Col_);
		}

		depth_const_iterator cbegin_at(size_t Row_, size_t Col_) const {
			return depth_const_iterator(data_.cbegin(), Row_, Col_);
		}

		depth_iterator end_at(size_t Row_, size_t Col_) {
			return depth_iterator(data_.end(), Row_, Col_);
		}

		depth_const_iterator cend_at(size_t Row_, size_t Col_) const {
			return depth_const_iterator(data_.cend(), Row_, Col_);
		}


		template<class EntrywiseGenerator>
		static cube<Ty, Alloc> creat_par(const size_t Depth_, const size_t Height_, const size_t Width_, EntrywiseGenerator&& Generate_);

	};

	template<typename Ty, class Alloc>
	template<class EntrywiseGenerator>
	cube<Ty, Alloc> cube<Ty, Alloc>::creat_par(const size_t Depth_, const size_t Height_, const size_t Width_, EntrywiseGenerator&& Generate_)
	{
		std::vector<mat<Ty, Alloc>> data(Depth_);
		for (unsigned i = 0; i < Depth_; ++i) {
			mat<Ty, Alloc> layer = mat<Ty, Alloc>::creat_par(Height_, Width_, std::bind_front(Generate_, i));
			data[i] = std::move(layer);
		}
		cube<Ty, Alloc> result;
		result.depth_ = Depth_;
		result.height_ = Height_;
		result.width_ = Width_;
		result.data_ = std::move(data);
		return result;
	}


}//namespace numer end