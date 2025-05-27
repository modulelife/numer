//	numer_mat.h		-*- C++20 -*-
#pragma once
//@brief: 2D array container template class mat<Ty, Alloc> header
//
//@classes:
//	numer::mat<Ty, Alloc>
//								------------ " container template class "
// 
//	numer::Mat_iterator_<mat>
//	numer::Mat_const_iterator_<mat>
//	numer::Mat_Row_iterator_<mat>
//	numer::Mat_Row_const_iterator_<mat>
//	numer::Mat_Row_cyc_iterator_<mat>
//	numer::Mat_Row_const_cyc_iterator_<mat>
//	numer::Mat_Col_iterator_<mat>
//	numer::Mat_Col_const_iterator_<mat>
//	numer::Mat_Col_cyc_iterator_<mat>
//	numer::Mat_Col_const_cyc_iterator_<mat>
//								------------ " subsidiary iterators "
// 
//	numer::Mat_row_view_<mat>
//	numer::Mat_col_view_<mat>
//								------------ " auxiliary access controller "
//
//@description:
//	2D array in compact memory layout  ( resizing is not supported !!! )
//	designated for organizing large amount of data that is logically 2 dimensional
//	but requires 1 dimensional continuous memory layout, e.g. images
//
//@usage:
//
//----------------------------content begins----------------------------

#include <cstddef>
#include <stdexcept>
#include <memory>
#include <utility>
#include <algorithm>
#include <execution>
#include <vector>

#include <numer_mat_utils.h>




namespace numer {


	template<class MyMat_>
	class Mat_iterator_ {
	public:
		friend MyMat_;

		using iterator_concept						= std::contiguous_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_::value_type;
		using difference_type						= typename MyMat_::difference_type;
		using pointer								= typename MyMat_::pointer;
		using reference								= typename MyMat_::reference;

		using self									= Mat_iterator_;

	private:
		pointer ptr_;


		Mat_iterator_(pointer Ptr_) : ptr_(Ptr_) {}

	public:
		Mat_iterator_() : ptr_(nullptr) {}

		reference operator*() const noexcept {
			return *ptr_;
		}

		pointer operator->() const noexcept {
			return ptr_;
		}

		self& operator++() noexcept {
			++ptr_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++ptr_;
			return tmp;
		}

		self& operator--() noexcept {
			--ptr_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--ptr_;
			return tmp;
		}

		self& operator+=(const difference_type Off_) noexcept {
			ptr_ += Off_;
			return *this;
		}

		self operator+(const difference_type Off_) const noexcept {
			return self(ptr_ + Off_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(const difference_type Off_) noexcept {
			ptr_ -= Off_;
			return *this;
		}

		self operator-(const difference_type Off_) const noexcept {
			return self(ptr_ - Off_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(ptr_ - Right_.ptr_);
		}

		reference operator[](const difference_type Off_) const noexcept {
			return *(ptr_ + Off_);
		}

		bool operator==(const self& Right_) const noexcept {
			return ptr_ == Right_.ptr_;
		}

		bool operator!=(const self& Right_) const noexcept {
			return ptr_ != Right_.ptr_;
		}

		bool operator<(const self& Right_) const noexcept {
			return ptr_ < Right_.ptr_;
		}

		bool operator<=(const self& Right_) const noexcept {
			return ptr_ <= Right_.ptr_;
		}

		bool operator>(const self& Right_) const noexcept {
			return ptr_ > Right_.ptr_;
		}

		bool operator>=(const self& Right_) const noexcept {
			return ptr_ >= Right_.ptr_;
		}

	};



	template<class MyMat_>
	class Mat_const_iterator_ {
	public:
		friend MyMat_;

		using iterator_concept						= std::contiguous_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_::value_type;
		using difference_type						= typename MyMat_::difference_type;
		using pointer								= typename MyMat_::const_pointer;
		using reference								= typename MyMat_::const_reference;

		using self									= Mat_const_iterator_;

	private:
		pointer ptr_;

		Mat_const_iterator_(pointer Ptr_) : ptr_(Ptr_) {}

	public:
		Mat_const_iterator_() : ptr_(nullptr) {}

		reference operator*() const noexcept {
			return *ptr_;
		}

		pointer operator->() const noexcept {
			return ptr_;
		}

		self& operator++() noexcept {
			++ptr_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++ptr_;
			return tmp;
		}

		self& operator--() noexcept {
			--ptr_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--ptr_;
			return tmp;
		}

		self& operator+=(const difference_type Off_) noexcept {
			ptr_ += Off_;
			return *this;
		}

		self operator+(const difference_type Off_) const noexcept {
			return self(ptr_ + Off_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(const difference_type Off_) noexcept {
			ptr_ -= Off_;
			return *this;
		}

		self operator-(const difference_type Off_) const noexcept {
			return self(ptr_ - Off_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(ptr_ - Right_.ptr_);
		}

		reference operator[](const difference_type Off_) const noexcept {
			return *(ptr_ + Off_);
		}

		bool operator==(const self& Right_) const noexcept {
			return ptr_ == Right_.ptr_;
		}
		
		bool operator!=(const self& Right_) const noexcept {
			return ptr_ != Right_.ptr_;
		}

		bool operator<(const self & Right_) const noexcept {
			return ptr_ < Right_.ptr_;
		}

		bool operator<=(const self& Right_) const noexcept {
			return ptr_ <= Right_.ptr_;
		}

		bool operator>(const self& Right_) const noexcept {
			return ptr_ > Right_.ptr_;
		}

		bool operator>=(const self& Right_) const noexcept {
			return ptr_ >= Right_.ptr_;
		}

	};



	template<class MyMat_rv_>
	class Mat_Row_iterator_ {
	public:
		friend MyMat_rv_;

		using iterator_concept						= std::random_access_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_rv_::value_type;
		using difference_type						= typename MyMat_rv_::difference_type;
		using pointer								= typename MyMat_rv_::pointer;
		using reference								= typename MyMat_rv_::reference;

		using self									= Mat_Row_iterator_;

	private:
		pointer row_begin_;
		difference_type curr_id_;
		difference_type max_id_;

		Mat_Row_iterator_(pointer Row_begin_, size_t Pos_, size_t Row_len_)
			: row_begin_(Row_begin_), curr_id_(Pos_), max_id_(Row_len_)
		{
		}

	public:
		Mat_Row_iterator_() : row_begin_(nullptr), curr_id_(0), max_id_(0) {}

		reference operator*() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return *(row_begin_ + curr_id_);
		}

		pointer operator->() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return row_begin_ + curr_id_;
		}

		self& operator++() {
			++curr_id_;
			return *this;
		}

		self operator++(int) {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() {
			--curr_id_;
			return *this;
		}

		self operator--(int) {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		self& operator+=(difference_type Off_) const {
			curr_id_ += Off_;
			return *this;
		}

		self operator+(difference_type Off_) const {
			return self(row_begin_, curr_id_ + Off_, max_id_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(difference_type Off_) const {
			curr_id_ += Off_;
			return *this;
		}

		self operator-(difference_type Off_) const {
			return self(row_begin_, curr_id_ - Off_, max_id_);
		}

		difference_type operator-(const self& Right_) const {
			return static_cast<difference_type>(curr_id_ - Right_.curr_id_);
		}

		reference operator[](const difference_type Off_) const {
			return *(*this + Off_);
		}

		bool operator==(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ == Right_.curr_id_;
		}

		bool operator!=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ != Right_.curr_id_;
		}

		bool operator<(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ < Right_.curr_id_;
		}

		bool operator<=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ <= Right_.curr_id_;
		}

		bool operator>(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ > Right_.curr_id_;
		}

		bool operator>=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ >= Right_.curr_id_;
		}

	};



	
	template<class MyMat_rv_>
	class Mat_Row_const_iterator_ {
	public:
		friend MyMat_rv_;

		using iterator_concept						= std::random_access_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_rv_::value_type;
		using difference_type						= typename MyMat_rv_::difference_type;
		using pointer								= typename MyMat_rv_::const_pointer;
		using reference								= typename MyMat_rv_::const_reference;

		using self									= Mat_Row_const_iterator_;

	private:
		pointer row_begin_;
		difference_type curr_id_;
		difference_type max_id_;

		Mat_Row_const_iterator_(pointer Row_begin_, size_t Pos_, size_t Row_len_)
			: row_begin_(Row_begin_), curr_id_(Pos_), max_id_(Row_len_)
		{
		}

	public:
		Mat_Row_const_iterator_() : row_begin_(nullptr), curr_id_(0), max_id_(0) {}

		reference operator*() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return *(row_begin_ + curr_id_);
		}

		pointer operator->() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return row_begin_ + curr_id_;
		}

		self& operator++() {
			++curr_id_;
			return *this;
		}

		self operator++(int) {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() {
			--curr_id_;
			return *this;
		}

		self operator--(int) {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		self& operator+=(difference_type Off_) const {
			curr_id_ += Off_;
			return *this;
		}

		self operator+(difference_type Off_) const {
			return self(row_begin_, curr_id_ + Off_, max_id_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(difference_type Off_) const {
			curr_id_ += Off_;
			return *this;
		}

		self operator-(difference_type Off_) const {
			return self(row_begin_, curr_id_ - Off_, max_id_);
		}

		difference_type operator-(const self& Right_) const {
			return static_cast<difference_type>(curr_id_ - Right_.curr_id_);
		}

		reference operator[](const difference_type Off_) const noexcept {
			return *(*this + Off_);
		}

		bool operator==(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ == Right_.curr_id_;
		}

		bool operator!=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ != Right_.curr_id_;
		}

		bool operator<(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ < Right_.curr_id_;
		}

		bool operator<=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ <= Right_.curr_id_;
		}

		bool operator>(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ > Right_.curr_id_;
		}

		bool operator>=(const self& Right_) const noexcept(false) {
			if (row_begin_ != Right_.row_begin_) throw std::logic_error("comparision of two iterators not in the same row");
			return curr_id_ >= Right_.curr_id_;
		}

	};



	template<class MyMat_rv_>
	class Mat_Row_cyc_iterator_ {
	public:
		friend MyMat_rv_;

		using iterator_concept						= std::bidirectional_iterator_tag;

		using iterator_category						= std::bidirectional_iterator_tag;
		using value_type							= typename MyMat_rv_::value_type;
		using difference_type						= typename MyMat_rv_::difference_type;
		using pointer								= typename MyMat_rv_::pointer;
		using reference								= typename MyMat_rv_::reference;

		using self									= Mat_Row_cyc_iterator_;

	private:
		pointer row_begin_;
		difference_type curr_id_;
		difference_type max_id_;

		Mat_Row_cyc_iterator_(pointer Row_begin_, size_t Pos_, size_t Row_len_)
			: row_begin_(Row_begin_), curr_id_(Pos_), max_id_(Row_len_)
		{
		}

	public:
		Mat_Row_cyc_iterator_()
			: row_begin_(nullptr), curr_id_(0), max_id_(0)
		{
		}

		reference operator*() const noexcept {
			return *(row_begin_ + (curr_id_ % max_id_ + max_id_) % max_id_);
		}

		pointer operator->() const noexcept {
			return row_begin_ + (curr_id_ % max_id_ + max_id_) % max_id_;
		}

		self& operator++() {
			++curr_id_;
			return *this;
		}

		self operator++(int) {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() {
			--curr_id_;
			return *this;
		}

		self operator--(int) {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		bool operator!=(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator==(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator<(const self&) const = delete;
		bool operator<=(const self&) const = delete;
		bool operator>(const self&) const = delete;
		bool operator>=(const self&) const = delete;
	};



	template<class MyMat_rv_>
	class Mat_Row_const_cyc_iterator_ {
	public:
		friend MyMat_rv_;

		using iterator_concept						= std::bidirectional_iterator_tag;

		using iterator_category						= std::bidirectional_iterator_tag;
		using value_type							= typename MyMat_rv_::value_type;
		using difference_type						= typename MyMat_rv_::difference_type;
		using pointer								= typename MyMat_rv_::const_pointer;
		using reference								= typename MyMat_rv_::const_reference;

		using self									= Mat_Row_const_cyc_iterator_;

	private:
		pointer row_begin_;
		difference_type curr_id_;
		difference_type max_id_;

		Mat_Row_const_cyc_iterator_(pointer Row_begin_, size_t Pos_, size_t Row_len_)
			: row_begin_(Row_begin_), curr_id_(Pos_), max_id_(Row_len_)
		{
		}

	public:
		Mat_Row_const_cyc_iterator_()
			: row_begin_(nullptr), curr_id_(0), max_id_(0)
		{
		}

		reference operator*() const noexcept {
			return *(row_begin_ + (curr_id_ % max_id_ + max_id_) % max_id_);
		}

		pointer operator->() const noexcept {
			return row_begin_ + (curr_id_ % max_id_ + max_id_) % max_id_;
		}

		self& operator++() {
			++curr_id_;
			return *this;
		}

		self operator++(int) {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() {
			--curr_id_;
			return *this;
		}

		self operator--(int) {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		bool operator!=(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator==(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator<(const self&) const = delete;
		bool operator<=(const self&) const = delete;
		bool operator>(const self&) const = delete;
		bool operator>=(const self&) const = delete;
	};



	template<class MyMat_>
	class Mat_row_view_ {
	public:
		friend MyMat_;

		using value_type							= typename MyMat_::value_type;
		using allocator_type						= typename MyMat_::allocator_type;
		using pointer								= typename MyMat_::pointer;
		using const_pointer							= typename MyMat_::const_pointer;
		using reference								= typename MyMat_::reference;
		using const_reference						= typename MyMat_::const_reference;
		using size_type								= typename MyMat_::size_type;
		using difference_type						= typename MyMat_::difference_type;

		using iterator								= Mat_Row_iterator_<Mat_row_view_>;
		using const_iterator						= Mat_Row_const_iterator_<Mat_row_view_>;
		using reverse_iterator						= std::reverse_iterator<iterator>;
		using const_reverse_iterator				= std::reverse_iterator<const_iterator>;

		using cyc_iterator							= Mat_Row_cyc_iterator_<Mat_row_view_>;
		using const_cyc_iterator					= Mat_Row_const_cyc_iterator_<Mat_row_view_>;
		using reverse_cyc_iterator					= std::reverse_iterator<cyc_iterator>;
		using const_reverse_cyc_iterator			= std::reverse_iterator<const_cyc_iterator>;

	private:
		const size_type row_id_;
		value_type* const& vdata_;
		const size_type& max_cid_;

		explicit Mat_row_view_(const size_type Row_id_, value_type* const& Data_, const size_type& Col_count_)
			: row_id_(Row_id_), vdata_(Data_), max_cid_(Col_count_)
		{
		}

	public:
		Mat_row_view_() = delete;
		Mat_row_view_(const Mat_row_view_&) = delete;
		Mat_row_view_& operator=(const Mat_row_view_&) = delete;

		bool operator!=(const Mat_row_view_&) const = delete;
		bool operator==(const Mat_row_view_&) const = delete;
		bool operator<(const Mat_row_view_&) const = delete;
		bool operator<=(const Mat_row_view_&) const = delete;
		bool operator>(const Mat_row_view_&) const = delete;
		bool operator>=(const Mat_row_view_&) const = delete;

		reference operator[](size_t Pos_) noexcept(false) {
			if (Pos_ >= max_cid_) throw std::out_of_range("index out of range");
			return *(vdata_ + row_id_ * max_cid_ + Pos_);
		}

		const_reference operator[](size_t Pos_) const noexcept(false) {
			if (Pos_ >= max_cid_) throw std::out_of_range("index out of range");
			return *(vdata_ + row_id_ * max_cid_ + Pos_);
		}

		size_type size() const {
			return max_cid_;
		}

		iterator begin() {
			return iterator(vdata_ + row_id_ * max_cid_, 0, max_cid_);
		}

		iterator end() {
			return iterator(vdata_ + row_id_ * max_cid_, max_cid_, max_cid_);
		}

		const_iterator cbegin() const {
			return const_iterator(vdata_ + row_id_ * max_cid_, 0, max_cid_);
		}

		const_iterator cend() const {
			return const_iterator(vdata_ + row_id_ * max_cid_, max_cid_, max_cid_);
		}

		reverse_iterator rbegin() {
			return reverse_iterator(iterator(vdata_ + row_id_ * max_cid_, max_cid_, max_cid_));
		}

		reverse_iterator rend() {
			return reverse_iterator(iterator(vdata_ + row_id_ * max_cid_, 0, max_cid_));
		}

		const_reverse_iterator crbegin() const {
			return const_reverse_iterator(const_iterator(vdata_ + row_id_ * max_cid_, max_cid_, max_cid_));
		}

		const_reverse_iterator crend() const {
			return const_reverse_iterator(const_iterator(vdata_ + row_id_ * max_cid_, 0, max_cid_));
		}

		cyc_iterator cycle_from(ptrdiff_t Pos_) {
			return cyc_iterator(vdata_ + row_id_ * max_cid_, Pos_, max_cid_);
		}

		const_cyc_iterator ccycle_from(ptrdiff_t Pos_) const {
			return const_cyc_iterator(vdata_ + row_id_ * max_cid_, Pos_, max_cid_);
		}

		reverse_cyc_iterator rcycle_from(ptrdiff_t Pos_) {
			return reverse_cyc_iterator(
				cyc_iterator(vdata_ + row_id_ * max_cid_, Pos_ + 1, max_cid_));
		}

		const_reverse_cyc_iterator crcycle_from(ptrdiff_t Pos_) const {
			return const_reverse_cyc_iterator(
				const_cyc_iterator(vdata_ + row_id_ * max_cid_, Pos_ + 1, max_cid_));
		}

	};



	template<class MyMat_cv_>
	class Mat_Col_iterator_ {
	public:
		friend MyMat_cv_;

		using iterator_concept						= std::random_access_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_cv_::value_type;
		using difference_type						= typename MyMat_cv_::difference_type;
		using pointer								= typename MyMat_cv_::pointer;
		using reference								= typename MyMat_cv_::reference;

		using self									= Mat_Col_iterator_;

	private:
		pointer col_begin_;
		difference_type curr_id_;
		difference_type max_id_;
		difference_type stride_;

		Mat_Col_iterator_(pointer Col_begin_, size_t Pos_, size_t Col_len_, size_t Row_len_)
			: col_begin_(Col_begin_), curr_id_(Pos_), max_id_(Col_len_), stride_(Row_len_)
		{
		}

	public:
		Mat_Col_iterator_()
			: col_begin_(nullptr), curr_id_(0), max_id_(0), stride_(0)
		{
		}

		reference operator*() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return *(col_begin_ + curr_id_ * stride_);
		}

		pointer operator->() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return col_begin_ + curr_id_ * stride_;
		}

		self& operator++() noexcept {
			++curr_id_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() noexcept {
			--curr_id_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		self& operator+=(difference_type Off_) const noexcept {
			curr_id_ += Off_;
			return *this;
		}

		self operator+(difference_type Off_) const noexcept{
			return self(col_begin_, curr_id_ + Off_, max_id_, stride_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(difference_type Off_) const noexcept {
			curr_id_ -= Off_;
			return *this;
		}

		self operator-(difference_type Off_) const noexcept {
			return self(col_begin_, curr_id_ - Off_, max_id_, stride_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(curr_id_ - Right_.curr_id_);
		}

		reference operator[](const difference_type Off_) const {
			return *(*this + Off_);
		}

		bool operator==(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ == Right_.curr_id_;
		}

		bool operator!=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ != Right_.curr_id_;
		}

		bool operator<(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ < Right_.curr_id_;
		}

		bool operator<=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ <= Right_.curr_id_;
		}

		bool operator>(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ > Right_.curr_id_;
		}

		bool operator>=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ >= Right_.curr_id_;
		}

	};



	template<class MyMat_cv_>
	class Mat_Col_const_iterator_ {
	public:
		friend MyMat_cv_;

		using iterator_concept						= std::random_access_iterator_tag;

		using iterator_category						= std::random_access_iterator_tag;
		using value_type							= typename MyMat_cv_::value_type;
		using difference_type						= typename MyMat_cv_::difference_type;
		using pointer								= typename MyMat_cv_::const_pointer;
		using reference								= typename MyMat_cv_::const_reference;

		using self									= Mat_Col_const_iterator_;

	private:
		pointer col_begin_;
		difference_type curr_id_;
		difference_type max_id_;
		difference_type stride_;

		Mat_Col_const_iterator_(pointer Col_begin_, size_t Pos_, size_t Col_len_, size_t Row_len_)
			: col_begin_(Col_begin_), curr_id_(Pos_), max_id_(Col_len_), stride_(Row_len_)
		{
		}

	public:
		Mat_Col_const_iterator_()
			: col_begin_(nullptr), curr_id_(0), max_id_(0), stride_(0)
		{
		}

		reference operator*() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return *(col_begin_ + curr_id_ * stride_);
		}

		pointer operator->() const noexcept(false) {
			if (curr_id_ < 0 || curr_id_ >= max_id_) throw std::out_of_range("iterator is out of range");
			return col_begin_ + curr_id_ * stride_;
		}

		self& operator++() noexcept {
			++curr_id_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() noexcept {
			--curr_id_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		self& operator+=(difference_type Off_) const noexcept {
			curr_id_ += Off_;
			return *this;
		}

		self operator+(difference_type Off_) const noexcept {
			return self(col_begin_, curr_id_ + Off_, max_id_, stride_);
		}

		friend self operator+(const difference_type Off_, self Next_) noexcept {
			Next_ += Off_;
			return Next_;
		}

		self& operator-=(difference_type Off_) const noexcept {
			curr_id_ -= Off_;
			return *this;
		}

		self operator-(difference_type Off_) const noexcept {
			return self(col_begin_, curr_id_ - Off_, max_id_, stride_);
		}

		difference_type operator-(const self& Right_) const noexcept {
			return static_cast<difference_type>(curr_id_ - Right_.curr_id_);
		}

		reference operator[](const difference_type Off_) const {
			return *(*this + Off_);
		}

		bool operator==(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ == Right_.curr_id_;
		}

		bool operator!=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ != Right_.curr_id_;
		}

		bool operator<(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ < Right_.curr_id_;
		}

		bool operator<=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ <= Right_.curr_id_;
		}

		bool operator>(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ > Right_.curr_id_;
		}

		bool operator>=(const self& Right_) const noexcept(false) {
			if (col_begin_ != Right_.col_begin_) throw std::logic_error("comparision of two iterators not in the same column");
			return curr_id_ >= Right_.curr_id_;
		}
	};



	template<class MyMat_cv_>
	class Mat_Col_cyc_iterator_ {
	public:
		friend MyMat_cv_;

		using iterator_concept						= std::bidirectional_iterator_tag;

		using iterator_category						= std::bidirectional_iterator_tag;
		using value_type							= typename MyMat_cv_::value_type;
		using difference_type						= typename MyMat_cv_::difference_type;
		using pointer								= typename MyMat_cv_::pointer;
		using reference								= typename MyMat_cv_::reference;

		using self									= Mat_Col_cyc_iterator_;

	private:
		pointer col_begin_;
		difference_type curr_id_;
		difference_type max_id_;
		difference_type stride_;

		Mat_Col_cyc_iterator_(pointer Col_begin_, size_t Pos_, size_t Col_len_, size_t Row_len_)
			: col_begin_(Col_begin_), curr_id_(Pos_), max_id_(Col_len_), stride_(Row_len_)
		{
		}

	public:

		Mat_Col_cyc_iterator_()
			: col_begin_(nullptr), curr_id_(0), max_id_(0), stride_(0)
		{
		}

		reference operator*() const noexcept {
			return *(col_begin_ + ((curr_id_ % max_id_ + max_id_) % max_id_) * stride_);
		}

		pointer operator->() const noexcept {
			return col_begin_ + ((curr_id_ % max_id_ + max_id_) % max_id_) * stride_;
		}

		self& operator++() noexcept {
			++curr_id_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() noexcept {
			--curr_id_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		bool operator!=(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator==(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator<(const self&) const = delete;
		bool operator<=(const self&) const = delete;
		bool operator>(const self&) const = delete;
		bool operator>=(const self&) const = delete;
	};



	template<class MyMat_cv_>
	class Mat_Col_const_cyc_iterator_ {
	public:
		friend MyMat_cv_;

		using iterator_concept						= std::bidirectional_iterator_tag;

		using iterator_category						= std::bidirectional_iterator_tag;
		using value_type							= typename MyMat_cv_::value_type;
		using difference_type						= typename MyMat_cv_::difference_type;
		using pointer								= typename MyMat_cv_::const_pointer;
		using reference								= typename MyMat_cv_::const_reference;

		using self = Mat_Col_const_cyc_iterator_;

	private:
		pointer col_begin_;
		difference_type curr_id_;
		difference_type max_id_;
		difference_type stride_;

		Mat_Col_const_cyc_iterator_(pointer Col_begin_, size_t Pos_, size_t Col_len_, size_t Row_len_)
			: col_begin_(Col_begin_), curr_id_(Pos_), max_id_(Col_len_), stride_(Row_len_)
		{
		}

	public:

		Mat_Col_const_cyc_iterator_()
			: col_begin_(nullptr), curr_id_(0), max_id_(0), stride_(0)
		{
		}

		reference operator*() const noexcept {
			return *(col_begin_ + ((curr_id_ % max_id_ + max_id_) % max_id_) * stride_);
		}

		pointer operator->() const noexcept {
			return col_begin_ + ((curr_id_ % max_id_ + max_id_) % max_id_) * stride_;
		}

		self& operator++() noexcept {
			++curr_id_;
			return *this;
		}

		self operator++(int) noexcept {
			self tmp = *this;
			++curr_id_;
			return tmp;
		}

		self& operator--() noexcept {
			--curr_id_;
			return *this;
		}

		self operator--(int) noexcept {
			self tmp = *this;
			--curr_id_;
			return tmp;
		}

		bool operator!=(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator==(const self& Right_) const { return curr_id_ != Right_.curr_id_; }
		bool operator<(const self&) const = delete;
		bool operator<=(const self&) const = delete;
		bool operator>(const self&) const = delete;
		bool operator>=(const self&) const = delete;
	};



	template<class MyMat_>
	class Mat_col_view_ {
	public:
		friend MyMat_;

		using value_type							= typename MyMat_::value_type;
		using allocator_type						= typename MyMat_::allocator_type;
		using pointer								= typename MyMat_::pointer;
		using const_pointer							= typename MyMat_::const_pointer;
		using reference								= typename MyMat_::reference;
		using const_reference						= typename MyMat_::const_reference;
		using size_type								= typename MyMat_::size_type;
		using difference_type						= typename MyMat_::difference_type;

		using iterator								= Mat_Col_iterator_<Mat_col_view_>;
		using const_iterator						= Mat_Col_const_iterator_<Mat_col_view_>;
		using reverse_iterator						= std::reverse_iterator<iterator>;
		using const_reverse_iterator				= std::reverse_iterator<const_iterator>;

		using cyc_iterator							= Mat_Col_cyc_iterator_<Mat_col_view_>;
		using const_cyc_iterator					= Mat_Col_const_cyc_iterator_<Mat_col_view_>;
		using reverse_cyc_iterator					= std::reverse_iterator<cyc_iterator>;
		using const_reverse_cyc_iterator			= std::reverse_iterator<const_cyc_iterator>;

	private:
		const size_type col_id_;
		value_type* const& vdata_;
		const size_type& max_cid_;
		const size_type& max_rid_;

		explicit Mat_col_view_(const size_type Col_id_, value_type* const& Data_, const size_type& Col_count_, const size_type& Row_count_)
			: col_id_(Col_id_), vdata_(Data_), max_cid_(Col_count_), max_rid_(Row_count_)
		{
		}

	public:
		Mat_col_view_() = delete;
		Mat_col_view_(const Mat_col_view_&) = delete;
		Mat_col_view_& operator=(const Mat_col_view_&) = delete;

		bool operator!=(const Mat_col_view_&) const = delete;
		bool operator==(const Mat_col_view_&) const = delete;
		bool operator<(const Mat_col_view_&) const = delete;
		bool operator<=(const Mat_col_view_&) const = delete;
		bool operator>(const Mat_col_view_&) const = delete;
		bool operator>=(const Mat_col_view_&) const = delete;

		reference operator[](size_type Pos_) noexcept(false) {
			if (Pos_ >= max_rid_) throw std::out_of_range("index out of range");
			return *(vdata_ + col_id_ + max_cid_ * Pos_);
		}

		const_reference operator[](size_type Pos_) const noexcept(false) {
			if (Pos_ >= max_rid_) throw std::out_of_range("index out of range");
			return *(vdata_ + col_id_ + max_cid_ * Pos_);
		}

		size_type size() const {
			return max_rid_;
		}

		iterator begin() {
			return iterator(vdata_ + col_id_, 0, max_rid_, max_cid_);
		}

		iterator end() {
			return iterator(vdata_ + col_id_, max_rid_, max_rid_, max_cid_);
		}

		const_iterator cbegin() const {
			return const_iterator(vdata_ + col_id_, 0, max_rid_, max_cid_);
		}

		const_iterator cend() const {
			return const_iterator(vdata_ + col_id_, max_rid_, max_rid_, max_cid_);
		}

		reverse_iterator rbegin() {
			return reverse_iterator(
				iterator(vdata_ + col_id_, max_rid_, max_rid_, max_cid_));
		}

		reverse_iterator rend() {
			return reverse_iterator(
				iterator(vdata_ + col_id_, 0, max_rid_, max_cid_));
		}

		const_reverse_iterator crbegin() const {
			return const_reverse_iterator(
				const_iterator(vdata_ + col_id_, max_rid_, max_rid_, max_cid_));
		}

		const_reverse_iterator crend() const {
			return const_reverse_iterator(
				const_iterator(vdata_ + col_id_, 0, max_rid_, max_cid_));
		}

		cyc_iterator cycle_from(ptrdiff_t Pos_) {
			return cyc_iterator(vdata_ + col_id_, Pos_, max_rid_, max_cid_);
		}

		const_cyc_iterator ccycle_from(ptrdiff_t Pos_) const {
			return const_cyc_iterator(vdata_ + col_id_, Pos_, max_rid_, max_cid_);
		}

		reverse_cyc_iterator rcycle_from(ptrdiff_t Pos_) {
			return reverse_cyc_iterator(
				cyc_iterator(vdata_ + col_id_, Pos_ + 1, max_rid_, max_cid_));
		}

		const_reverse_cyc_iterator crcycle_from(ptrdiff_t Pos_) const {
			return const_reverse_cyc_iterator(
				const_cyc_iterator(vdata_ + col_id_, Pos_ + 1, max_rid_, max_cid_));
		}

	};




	/**
	*2D array in compact memory layout  ( resizing is not supported !!! )
	*designated for organizing large amount of data that is logically 2 dimensional
	*but requires 1 dimensional continuous memory layout, e.g. images
	*/
	template<typename Ty, class Alloc = std::allocator<Ty> >
	class mat final {
	private:
		Ty*				data_;
		size_t			nrows_;
		size_t			ncols_;
		Alloc			allocator_;


		Ty* allocate_(size_t n_) {
			if (n_ == 0) return nullptr;
			Ty* p = allocator_.allocate(n_);
			return p;
		}

		Ty* allocate_and_construct_(size_t n_) {
			if (n_ == 0) return nullptr;

			Ty* p = allocator_.allocate(n_);
			try {
				std::uninitialized_default_construct_n(p, n_);
			}
			catch (...) {
				allocator_.deallocate(p, n_);
				throw;
			}
			return p;
		}

		void destroy_and_deallocate_() {
			if (data_ != nullptr) {
				std::destroy_n(data_, count_());
				allocator_.deallocate(data_, count_());
				data_ = nullptr;
			}
		}

		Ty* allocate_and_fill_(size_t n_, const Ty& Value_) {
			if (n_ == 0) return nullptr;
			Ty* p = allocator_.allocate(n_);
			std::uninitialized_fill_n(p, n_, Value_);
			return p;
		}

		Ty* mem_begin_() const { return data_; }
		Ty* mem_end_() const { return data_ + nrows_ * ncols_; }
		size_t count_() const { return nrows_ * ncols_; }

	public:
		template<typename Tx, class Allocx> friend class mat;

		using value_type							= Ty;
		using allocator_type						= Alloc;
		using pointer								= Ty*;
		using const_pointer							= const Ty*;
		using reference								= Ty&;
		using const_reference						= const Ty&;
		using size_type								= size_t;
		using difference_type						= ptrdiff_t;

		using iterator								= Mat_iterator_<mat>;
		using const_iterator						= Mat_const_iterator_<mat>;
		using reverse_iterator						= std::reverse_iterator<iterator>;
		using const_reverse_iterator				= std::reverse_iterator<const_iterator>;

		using row_view								= Mat_row_view_<mat>;
		using col_view								= Mat_col_view_<mat>;



		//default ctor
		mat() : nrows_(0), ncols_(0), data_(nullptr), allocator_() {}

		//ctor specifying initial size
		mat(const size_t Rows_, const size_t Cols_)
			: nrows_(Rows_), ncols_(Cols_), data_(nullptr), allocator_()
		{
			data_ = allocate_and_construct_(count_());
		}

		//ctor specifying initial size & element value
		mat(const size_t Rows_, const size_t Cols_, const Ty& Value_)
			:nrows_(Rows_), ncols_(Cols_), data_(nullptr), allocator_()
		{
			data_ = allocate_and_fill_(count_(), Value_);
		}

		//copy ctor
		mat(const mat<Ty, Alloc>& Right_) noexcept
			:nrows_(Right_.nrows()), ncols_(Right_.ncols()), data_(nullptr), allocator_(Right_.allocator_)
		{
			data_ = allocate_(count_());
			std::uninitialized_copy_n(Right_.mem_begin_(), count_(), mem_begin_());
		}

		//dtor
		~mat() {
			destroy_and_deallocate_();
			nrows_ = 0;
			ncols_ = 0;
		}

		//assignment
		mat<Ty, Alloc>& operator=(const mat<Ty, Alloc>& Right_) {
			if (this == &Right_) return *this;
			Ty* new_data = allocate_(Right_.count_());
			try {
				std::uninitialized_copy_n(Right_.mem_begin_(), Right_.count_(),new_data);
			}
			catch (...) {
				allocator_.deallocate(new_data, Right_.count_());
				throw;
			}

			destroy_and_deallocate_();
			nrows_ = Right_.nrows_;
			ncols_ = Right_.ncols_;
			data_ = new_data;
			return *this;
		}

		//move ctor
		mat(mat<Ty, Alloc>&& Right_) noexcept
			:nrows_(Right_.nrows_), ncols_(Right_.ncols_), data_(Right_.data_), allocator_(std::move(Right_.allocator_))
		{
			Right_.data_ = nullptr;
			Right_.nrows_ = 0;
			Right_.ncols_ = 0;
		}

		//move assignment
		mat<Ty, Alloc>& operator=(mat<Ty, Alloc>&& Right_) noexcept {
			if (this == &Right_) return *this; 
			destroy_and_deallocate_();
			nrows_ = Right_.nrows_;
			ncols_ = Right_.ncols_;
			data_ = Right_.data_;
			allocator_ = std::move(Right_.allocator_);

			Right_.data_ = nullptr;
			Right_.nrows_ = 0;
			Right_.ncols_ = 0;
			return *this;
		}

		//clear content
		void clear() {
			destroy_and_deallocate_();
			nrows_ = 0;
			ncols_ = 0;
		}

		//static swap
		static void swap(mat<Ty, Alloc>& Mat_X_, mat<Ty, Alloc>& Mat_Y_) {
			using std::swap;
			swap(Mat_X_.nrows_, Mat_Y_.nrows_);
			swap(Mat_X_.ncols_, Mat_Y_.ncols_);
			swap(Mat_X_.data_, Mat_Y_.data_);
			swap(Mat_X_.allocator_, Mat_Y_.allocator_);
		}

		//swap
		void swap(mat<Ty, Alloc>& Right_) {
			using std::swap;
			swap(nrows_, Right_.nrows_);
			swap(ncols_, Right_.ncols_);
			swap(data_, Right_.data_);
			swap(allocator_, Right_.allocator_);
		}

		//number of rows
		size_type nrows() const { return nrows_; }

		//number of columns
		size_type ncols() const { return ncols_; }

		//availability test
		explicit operator bool() const { return data_ != nullptr; }

		//empty test
		bool empty() const { return data_ == nullptr; }

		//total number of entries
		size_type size() const {
			return nrows_ * ncols_;
		}

		//not equal
		bool operator!=(const mat<Ty, Alloc>& Right_) const;

		//equal
		bool operator==(const mat<Ty, Alloc>& Right_) const;

		bool operator<(const mat<Ty, Alloc>&) const = delete;
		bool operator<=(const mat<Ty, Alloc>&) const = delete;
		bool operator>(const mat<Ty, Alloc>&) const = delete;
		bool operator>=(const mat<Ty, Alloc>&) const = delete;

		//get all data in 1D array : A[rows*cols]
		iterator begin()
		{
			return iterator(data_);
		}

		//const get all data in 1D array : A[rows*cols]
		const_iterator cbegin() const
		{
			return const_iterator(data_);
		}

		//get end of all data in 1D array : A[rows*cols]
		iterator end()
		{
			return iterator(data_ + count_());
		}

		//const get end of all data in 1D array : A[rows*cols]
		const_iterator cend() const
		{
			return const_iterator(data_ + count_());
		}

		//get all data in 1D array : A[rows*cols]
		reverse_iterator rbegin()
		{
			return reverse_iterator(iterator(data_ + count_()));
		}

		//const get all data in 1D array : A[rows*cols]
		const_reverse_iterator crbegin() const
		{
			return const_reverse_iterator(const_iterator(data_ + count_()));
		}

		//get end of all data in 1D array : A[rows*cols]
		reverse_iterator rend()
		{
			return reverse_iterator(iterator(data_));
		}

		//const get end of all data in 1D array : A[rows*cols]
		const_reverse_iterator crend() const
		{
			return const_reverse_iterator(const_iterator(data_));
		}

		//static 2d-array style access row
		pointer operator[](size_t Row_id_) noexcept
		{
			return data_ + ncols_ * Row_id_;
		}

		//classic 2d-array const access row
		const_pointer operator[](size_t Row_id_) const noexcept
		{
			return data_ + ncols_ * Row_id_;
		}

		//access row
		row_view row(size_t Row_id_) noexcept(false)
		{
			if (Row_id_ >= nrows()) throw std::out_of_range("index out of range");
			return row_view(Row_id_, data_, ncols_);
		}

		//const access row
		const row_view row(size_t Row_id_) const noexcept(false)
		{
			if (Row_id_ >= nrows()) throw std::out_of_range("index out of range");
			return row_view(Row_id_, data_, ncols_);
		}

		//access column
		col_view col(size_t Col_id_) noexcept(false)
		{
			if (Col_id_ >= ncols()) throw std::out_of_range("index out of range");
			return col_view(Col_id_, data_, ncols_, nrows_);
		}

		//const access column
		const col_view col(size_t Col_id_) const noexcept(false)
		{
			if (Col_id_ >= ncols()) throw std::out_of_range("index out of range");
			return col_view(Col_id_, data_, ncols_, nrows_);
		}

		//factory method accepting a Generator : Ty (size_t i, size_t j), i/j stands for row/col index
		template<class EntrywiseGenerator>
		static mat<Ty, Alloc> creat(const size_t Rows_, const size_t Cols_, EntrywiseGenerator&& Generate_);

		//factory method accepting mat of another type and a Converter : Ty (const Tx&)
		template<class EntrywiseConverter, typename Tx, class Allocx>
		static mat<Ty, Alloc> creat(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_);

		//parapllel execution ver. of factory method accepting a Generator : Ty (size_t i, size_t j), i/j stands for row/col index
		template<class EntrywiseGenerator>
		static mat<Ty, Alloc> creat_par(const size_t Rows_, const size_t Cols_, EntrywiseGenerator&& Generate_);

		//parapllel execution ver. of factory method accepting mat of another type and a  : Ty Convert_(const Tx&)
		template<class EntrywiseConverter, typename Tx, class Allocx>
		static mat<Ty, Alloc> creat_par(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_);

		//refill the mat with a new value
		void set_all_to(const Ty& Value_);

		//refill the mat with a Generate
		template<class EntrywiseGenerator>
		void refill(EntrywiseGenerator&& Generate_);

		//refill the mat with another mat of the same size
		template<class EntrywiseConverter, typename Tx, class Allocx>
		void refill(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_);

		//refill the mat with a Generate
		template<class EntrywiseGenerator>
		void refill_par(EntrywiseGenerator&& Generate_);

		//refill the mat with another mat of the same size
		template<class EntrywiseConverter, typename Tx, class Allocx>
		void refill_par(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_);

		//select entries in range [R1, R2) X [C1, C2), returning a new mat
		mat<Ty, Alloc> select_range(size_t R1_, size_t R2_, size_t C1_, size_t C2_) const noexcept(false);

		//modify all elements with a Modifer : void (Ty&)
		template<class Modifier>
		void modify(Modifier&& Modify_);

		//parallel execution ver. modify all elements with a Modifer : void (Ty&)
		template<class Modifier>
		void modify_par(Modifier&& Modify_);

		//put a Patch_ on caller mat at (Row_Pos, Col_Pos), aligning upper-left corner of the Patch_
		//if range is invalid, no changes will be made, thus "noexcept"
		//Mixer should be any callable object with signature " Ty (const Ty&, const Ty&) "
		template<class Mixer>
		mat& overlay(const mat<Ty, Alloc>& Patch_, size_t Row_pos_, size_t Col_pos_, Mixer&& Mix_) noexcept;

		//overlay method for all compatible types, Mixer : any callable " Ty (const Ty&, const Tx&) "
		template<class Mixer, typename Tx, class Allocx>
		mat& overlay(const mat<Tx, Allocx>& Patch_, size_t Row_pos_, size_t Col_pos_, Mixer&& Mix_) noexcept;

		//overlay method dedicated for patch of the same shape
		template<class Mixer>
		void overlay(const mat<Ty, Alloc>& Mask_, Mixer&& Mix_) noexcept(false);

		//overlay method dedicated for patch of the same shape, for all compatible types
		template<class Mixer, typename Tx, class Allocx>
		void overlay(const mat<Tx, Allocx>& Mask_, Mixer&& Mix_) noexcept(false);

		//parallel execution ver. of overlay method dedicated for patch of the same shape
		template<class Mixer>
		void overlay_par(const mat<Ty, Alloc>& Mask_, Mixer&& Mix_) noexcept(false);

		//parallel execution ver. of overlay method dedicated for patch of the same shape, for all compatible types
		template<class Mixer, typename Tx, class Allocx>
		void overlay_par(const mat<Tx, Allocx>& Mask_, Mixer&& Mix_) noexcept(false);

	};


	//mat member functions def

	template<typename Ty, class Alloc>
	template<class EntrywiseGenerator>
	mat<Ty, Alloc> mat<Ty, Alloc>::creat(const size_t Rows_, const size_t Cols_, EntrywiseGenerator&& Generate_)
	{
		mat<Ty, Alloc> result(Rows_, Cols_);
		Ty* iter = result.mem_begin_();
		for (size_t i = 0; i < Rows_; ++i)
		{
			for (size_t j = 0; j < Cols_; ++j) {
				*iter++ = Generate_(i, j);
			}
		}
		return result;
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseConverter, typename Tx, class Allocx>
	mat<Ty, Alloc> mat<Ty, Alloc>::creat(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_)
	{
		mat<Ty, Alloc> result(Src_.nrows(), Src_.ncols());
		Ty* iter = result.mem_begin_();
		for (const Tx* oiter = Src_.mem_begin_(), * oend = Src_.mem_end_();
			oiter != oend;
			++iter, ++oiter)
		{
			*iter = Convert_(*oiter);
		}
		return result;
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseGenerator>
	mat<Ty, Alloc> mat<Ty, Alloc>::creat_par(const size_t Rows_, const size_t Cols_, EntrywiseGenerator&& Generate_)
	{
		using args__ = struct { Ty* row_begin__; size_t i__; };

		mat<Ty, Alloc> result(Rows_, Cols_);
		
		Ty* row_begin = result.mem_begin_();
		std::vector<args__> args_list(Rows_);
		for (size_t i = 0; i < Rows_; ++i) {
			args_list[i] = { row_begin, i};
			row_begin += Cols_;
		}

		std::for_each_n(std::execution::par, args_list.begin(), Rows_,
			[&](args__ args) {
				for (size_t j = 0; j < Cols_; ++j)
				{
					*args.row_begin__++ = Generate_(args.i__, j);
				}
			});
		return result;
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseConverter, typename Tx, class Allocx>
	mat<Ty, Alloc> mat<Ty, Alloc>::creat_par(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_)
	{
		mat<Ty, Alloc> result(Src_.nrows(), Src_.ncols());
		std::transform(std::execution::par_unseq, Src_.mem_begin_(), Src_.mem_end_(), result.mem_begin_(), Convert_);
		return result;
	}

	template<typename Ty, class Alloc>
	bool mat<Ty, Alloc>::operator!=(const mat<Ty, Alloc>& Right_) const
	{
		if (this->empty() || Right_.empty()) return true;
		if (nrows() != Right_.nrows() || ncols() != Right_.ncols()) return true;

		for (const Ty* iter = mem_begin_(), *oiter = Right_.mem_begin_(), *oend = Right_.mem_end_();
			oiter != oend;
			++iter, ++oiter)
		{
			if (*iter != *oiter) return true;
		}
		return false;
	}

	template<typename Ty, class Alloc>
	bool mat<Ty, Alloc>::operator==(const mat<Ty, Alloc>& Right_) const
	{
		if (this->empty() || Right_.empty()) return false;
		if (nrows() != Right_.nrows() || ncols() != Right_.ncols()) return false;

		for (const Ty* iter = mem_begin_(), *oiter = Right_.mem_begin_(), *oend = Right_.mem_end_();
			oiter != oend;
			++iter, ++oiter)
		{
			if (*iter != *oiter) return false;
		}
		return true;
	}

	template<typename Ty, class Alloc>
	void mat<Ty, Alloc>::set_all_to(const Ty& Value_)
	{
		Ty* iter = mem_begin_();
		for (size_t i = 0; i < count_(); ++i)
		{
			*iter++ = Value_;
		}
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseGenerator>
	void mat<Ty, Alloc>::refill(EntrywiseGenerator&& Generate_)
	{
		Ty* iter = mem_begin_();
		const size_t rows = nrows(), cols = ncols();
		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < cols; ++j) {
				*iter++ = Generate_(i, j);
			}
		}
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseConverter, typename Tx, class Allocx>
	void mat<Ty, Alloc>::refill(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_)
	{
		if (nrows() != Src_.nrows() || ncols() != Src_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		Ty* iter = mem_begin_();
		for (const Tx* oiter = Src_.mem_begin_(), *oend = Src_.mem_end_();
			oiter != oend;
			++iter, ++oiter)
		{
			*iter = Convert_(*oiter);
		}
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseGenerator>
	void mat<Ty, Alloc>::refill_par(EntrywiseGenerator&& Generate_)
	{
		using args__ = struct { Ty* row_begin__; size_t i__; };

		Ty* row_begin = mem_begin_();
		const size_t rows = nrows(), cols = ncols();

		std::vector<args__> args_list(rows);
		for (size_t i = 0; i < rows; ++i) {
			args_list[i] = { row_begin, i };
			row_begin += cols;
		}

		std::for_each_n(std::execution::par, args_list.begin(), rows,
			[&](args__ args) {
				for (size_t j = 0; j < cols; ++j)
				{
					*args.row_begin__++ = Generate_(args.i__, j);
				}
			});
	}

	template<typename Ty, class Alloc>
	template<class EntrywiseConverter, typename Tx, class Allocx>
	void mat<Ty, Alloc>::refill_par(const mat<Tx, Allocx>& Src_, EntrywiseConverter&& Convert_)
	{
		if (nrows() != Src_.nrows() || ncols() != Src_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		std::transform(std::execution::par_unseq, Src_.mem_begin_(), Src_.mem_end_(), mem_begin_(), Convert_);
	}

	template<typename Ty, class Alloc>
	mat<Ty, Alloc> mat<Ty, Alloc>::select_range(size_t R1_, size_t R2_, size_t C1_, size_t C2_) const noexcept(false)
	{
		if (R1_ >= R2_ || C1_ >= C2_ || R2_ > nrows() || C2_ > ncols()) {
			throw std::range_error("invalid range");
		}

		size_t rows = R2_ - R1_;
		size_t cols = C2_ - C1_;
		mat<Ty> result(rows, cols);

		for (size_t i = 0; i < rows; ++i)
		{
			std::uninitialized_copy_n(
				mem_begin_() + (i + R1_) * ncols() + C1_,
				(size_t)cols,
				result.mem_begin_() + i * cols);
		}
		return result;
	}

	template<typename Ty, class Alloc>
	template<class Modifier>
	void mat<Ty, Alloc>::modify(Modifier&& Modify_)
	{
		std::for_each_n(mem_begin_(), count_(), Modify_);
	}

	template<typename Ty, class Alloc>
	template<class Modifier>
	void mat<Ty, Alloc>::modify_par(Modifier&& Modify_)
	{
		std::for_each_n(std::execution::par_unseq, mem_begin_(), count_(), Modify_);
	}

	template<typename Ty, class Alloc>
	template<class Mixer>
	mat<Ty, Alloc>& mat<Ty, Alloc>::overlay(const mat<Ty, Alloc>& Patch_, size_t Row_pos_, size_t Col_pos_, Mixer&& Mix_) noexcept
	{
		size_t rows = Patch_.nrows();
		size_t cols = Patch_.ncols();

		if (Row_pos_ + rows > nrows() || Col_pos_ + cols > ncols()) return *this;

		const Ty* patch_iter = Patch_.mem_begin_();
		for (size_t i = 0; i < rows; ++i)
		{
			Ty* iter = mem_begin_() + (i + Row_pos_) * ncols() + Col_pos_;
			for (size_t j = 0; j < cols; ++j)
			{
				*iter = Mix_(*iter, *patch_iter++);
				++iter;
			}
		}
		return *this;
	}

	template<typename Ty, class Alloc>
	template<class Mixer, typename Tx, class Allocx>
	mat<Ty, Alloc>& mat<Ty, Alloc>::overlay(const mat<Tx, Allocx>& Patch_, size_t Row_pos_, size_t Col_pos_, Mixer&& Mix_) noexcept
	{
		size_t rows = Patch_.nrows();
		size_t cols = Patch_.ncols();

		if (Row_pos_ + rows > nrows() || Col_pos_ + cols > ncols()) return *this;

		const Tx* patch_iter = Patch_.mem_begin_();
		for (size_t i = 0; i < rows; ++i)
		{
			Ty* iter = mem_begin_() + (i + Row_pos_) * ncols() + Col_pos_;
			for (size_t j = 0; j < cols; ++j)
			{
				*iter = Mix_(*iter, *patch_iter++);
				++iter;
			}
		}
		return *this;
	}

	template<typename Ty, class Alloc>
	template<class Mixer>
	void mat<Ty, Alloc>::overlay(const mat<Ty, Alloc>& Mask_, Mixer&& Mix_) noexcept(false)
	{
		if (nrows() != Mask_.nrows() || ncols() != Mask_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		Ty* iter = mem_begin_();
		for (const Ty* miter = Mask_.mem_begin_(), * mend = Mask_.mem_end_();
			miter != mend;
			++iter, ++miter)
		{
			*iter = Mix_(*iter, *miter);
		}
	}

	template<typename Ty, class Alloc>
	template<class Mixer, typename Tx, class Allocx>
	void mat<Ty, Alloc>::overlay(const mat<Tx, Allocx>& Mask_, Mixer&& Mix_) noexcept(false)
	{
		if (nrows() != Mask_.nrows() || ncols() != Mask_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		Ty* iter = mem_begin_();
		for (const Tx* miter = Mask_.mem_begin_(), *mend = Mask_.mem_end_();
			miter != mend;
			++iter, ++miter)
		{
			*iter = Mix_(*iter, *miter);
		}
	}

	template<typename Ty, class Alloc>
	template<class Mixer>
	void mat<Ty, Alloc>::overlay_par(const mat<Ty, Alloc>& Mask_, Mixer&& Mix_) noexcept(false)
	{
		if (nrows() != Mask_.nrows() || ncols() != Mask_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		std::transform(std::execution::par_unseq, mem_begin_(), mem_end_(), Mask_.mem_begin_(), mem_begin_(), Mix_);
	}

	template<typename Ty, class Alloc>
	template<class Mixer, typename Tx, class Allocx>
	void mat<Ty, Alloc>::overlay_par(const mat<Tx, Allocx>& Mask_, Mixer&& Mix_) noexcept(false)
	{
		if (nrows() != Mask_.nrows() || ncols() != Mask_.ncols()) throw std::range_error("the shapes of the two mats do not match");

		std::transform(std::execution::par_unseq, mem_begin_(), mem_end_(), Mask_.mem_begin_(), mem_begin_(), Mix_);
	}


}// namespace numer end