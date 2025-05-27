//	numer_file.h		-*- C++20 -*-
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
#include <numer_mat.h>
#include <numer_cube.h>
#include <type_traits>
#include <vector>
#include <string>
#include <fstream>




namespace numer {


	using Seq_File_Header = struct {
		size_t length;
		short elem_size;
		short info_size;
		unsigned reserved__;
	};

	using Mat_File_Header = struct {
		size_t rows;
		size_t cols;
		short elem_size;
		short info_size;
		unsigned reserved__;
	};

	using Cube_File_Header = struct {
		size_t depth;
		size_t height;
		size_t width;
		short elem_size;
		short info_size;
		unsigned reserved__;
	};

	template<typename Ty, class Alloc__>
	bool writeMat(std::string Path_no_suffix_, const mat<Ty, Alloc__>& Src_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		if (Src_.empty()) return false;

		const Mat_File_Header head{ Src_.nrows(), Src_.ncols(), sizeof(Ty), 0, 0};
		const size_t bytes = head.rows * head.cols * head.elem_size;
		Path_no_suffix_ += ".nmmat";

		std::ofstream nmmat_out(Path_no_suffix_.c_str(), std::ios::binary | std::ios::trunc);
		if (!nmmat_out.is_open()) return false;

		nmmat_out.write((char*)&head, sizeof(Mat_File_Header));
		nmmat_out.write((char*)(Src_[0]), bytes);
		nmmat_out.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool writeMat(std::string Path_no_suffix_, const mat<Ty, Alloc__>& Src_, const Info_Header Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");
		if (Src_.empty()) return false;

		const Mat_File_Header head{ Src_.nrows(), Src_.ncols(), sizeof(Ty), sizeof(Info_Header), 0};
		const size_t bytes = head.rows * head.cols * head.elem_size;
		Path_no_suffix_ += ".nmmat";

		std::ofstream nmmat_out(Path_no_suffix_.c_str(), std::ios::binary | std::ios::trunc);
		if (!nmmat_out.is_open()) return false;

		nmmat_out.write((char*)&head, sizeof(Mat_File_Header));
		nmmat_out.write((char*)&Info_, sizeof(Info_Header));
		nmmat_out.write((char*)(Src_[0]), bytes);
		nmmat_out.close();
		return true;
	}

	template<typename Ty, class Alloc__>
	bool readMat(std::string Path_no_suffix_, mat<Ty, Alloc__>& Dst_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");

		Path_no_suffix_ += ".nmmat";
		std::ifstream nmmat_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmmat_in.is_open()) return false;

		Mat_File_Header head{};
		nmmat_in.read((char*)&head, sizeof(Mat_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) nmmat_in.seekg(head.info_size, std::ios::cur);

		mat<Ty, Alloc__> buffer(head.rows, head.cols);
		const size_t bytes = head.rows * head.cols * head.elem_size;
		nmmat_in.read((char*)buffer[0], bytes);
		Dst_ = std::move(buffer);
		nmmat_in.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool readMat_weak(std::string Path_no_suffix_, mat<Ty, Alloc__>& Dst_, Info_Header& Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");

		Path_no_suffix_ += ".nmmat";
		std::ifstream nmmat_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmmat_in.is_open()) return false;

		Mat_File_Header head{};
		nmmat_in.read((char*)&head, sizeof(Mat_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) {
			if (head.info_size == sizeof(Info_Header)) nmmat_in.read((char*)&Info_, sizeof(Info_Header));
			else nmmat_in.seekg(head.info_size, std::ios::cur);
		}

		mat<Ty, Alloc__> buffer(head.rows, head.cols);
		const size_t bytes = head.rows * head.cols * head.elem_size;
		nmmat_in.read((char*)buffer[0], bytes);
		Dst_ = std::move(buffer);
		nmmat_in.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool readMat_strong(std::string Path_no_suffix_, mat<Ty, Alloc__>& Dst_, Info_Header& Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");

		Path_no_suffix_ += ".nmmat";
		std::ifstream nmmat_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmmat_in.is_open()) return false;

		Mat_File_Header head{};
		nmmat_in.read((char*)&head, sizeof(Mat_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) {
			if (head.info_size == sizeof(Info_Header)) nmmat_in.read((char*)&Info_, sizeof(Info_Header));
			else return false;
		}
		else {
			return false;
		}

		mat<Ty, Alloc__> buffer(head.rows, head.cols);
		const size_t bytes = head.rows * head.cols * head.elem_size;
		nmmat_in.read((char*)buffer[0], bytes);
		Dst_ = std::move(buffer);
		nmmat_in.close();
		return true;
	}

	template<typename Ty, class Alloc__>
	bool writeCube(std::string Path_no_suffix_, const cube<Ty, Alloc__>& Src_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		if (Src_.empty()) return false;

		const Cube_File_Header head{ Src_.depth(), Src_.height(), Src_.width(), sizeof(Ty), 0, 0 };
		const size_t bytes = head.height * head.width * head.elem_size;
		Path_no_suffix_ += ".nmcube";

		std::ofstream nmcube_out(Path_no_suffix_.c_str(), std::ios::binary | std::ios::trunc);
		if (!nmcube_out.is_open()) return false;

		nmcube_out.write((char*)&head, sizeof(Cube_File_Header));
		for (size_t i = 0; i < head.depth; ++i) {
			nmcube_out.write((char*)(Src_[i][0]), bytes);
		}
		nmcube_out.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool writeCube(std::string Path_no_suffix_, const cube<Ty, Alloc__>& Src_, const Info_Header Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");
		if (Src_.empty()) return false;

		const Cube_File_Header head{ Src_.depth(), Src_.height(), Src_.width(), sizeof(Ty), sizeof(Info_Header), 0 };
		const size_t bytes = head.height * head.width * head.elem_size;
		Path_no_suffix_ += ".nmcube";

		std::ofstream nmcube_out(Path_no_suffix_.c_str(), std::ios::binary | std::ios::trunc);
		if (!nmcube_out.is_open()) return false;

		nmcube_out.write((char*)&head, sizeof(Cube_File_Header));
		nmcube_out.write((char*)&Info_, sizeof(Info_Header));
		for (size_t i = 0; i < head.depth; ++i) {
			nmcube_out.write((char*)(Src_[i][0]), bytes);
		}
		nmcube_out.close();
		return true;
	}

	template<typename Ty, class Alloc__>
	bool readCube(std::string Path_no_suffix_, cube<Ty, Alloc__>& Dst_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");

		Path_no_suffix_ += ".nmcube";
		std::ifstream nmcube_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmcube_in.is_open()) return false;

		Cube_File_Header head{};
		nmcube_in.read((char*)&head, sizeof(Cube_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) nmcube_in.seekg(head.info_size, std::ios::cur);

		cube<Ty, Alloc__> buffer(head.depth, head.height, head.width);
		const size_t bytes = head.height * head.width * head.elem_size;
		for (size_t i = 0; i < head.depth; ++i) {
			nmcube_in.read((char*)buffer[i][0], bytes);
		}
		Dst_ = std::move(buffer);
		nmcube_in.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool readCube_weak(std::string Path_no_suffix_, cube<Ty, Alloc__>& Dst_, Info_Header& Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");

		Path_no_suffix_ += ".nmcube";
		std::ifstream nmcube_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmcube_in.is_open()) return false;

		Cube_File_Header head{};
		nmcube_in.read((char*)&head, sizeof(Cube_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) {
			if (head.info_size == sizeof(Info_Header)) nmcube_in.read((char*)&Info_, sizeof(Info_Header));
			else nmcube_in.seekg(head.info_size, std::ios::cur);
		}

		cube<Ty, Alloc__> buffer(head.depth, head.height, head.width);
		const size_t bytes = head.height * head.width * head.elem_size;
		for (size_t i = 0; i < head.depth; ++i) {
			nmcube_in.read((char*)buffer[i][0], bytes);
		}
		Dst_ = std::move(buffer);
		nmcube_in.close();
		return true;
	}

	template<typename Ty, class Alloc__, typename Info_Header>
	bool readCube_strong(std::string Path_no_suffix_, cube<Ty, Alloc__>& Dst_, Info_Header& Info_)
	{
		static_assert(std::is_trivially_copyable<Ty>::value, "Element must be plain data type");
		static_assert(std::is_trivially_copyable<Info_Header>::value, "Information header must be plain data type");

		Path_no_suffix_ += ".nmcube";
		std::ifstream nmcube_in(Path_no_suffix_.c_str(), std::ios::binary | std::ios::in);
		if (!nmcube_in.is_open()) return false;

		Cube_File_Header head{};
		nmcube_in.read((char*)&head, sizeof(Cube_File_Header));
		if (head.elem_size != sizeof(Ty)) return false;
		if (head.info_size != 0) {
			if (head.info_size == sizeof(Info_Header)) nmcube_in.read((char*)&Info_, sizeof(Info_Header));
			else return false;
		}
		else {
			return false;
		}

		cube<Ty, Alloc__> buffer(head.depth, head.height, head.width);
		const size_t bytes = head.height * head.width * head.elem_size;
		for (size_t i = 0; i < head.depth; ++i) {
			nmcube_in.read((char*)buffer[i][0], bytes);
		}
		Dst_ = std::move(buffer);
		nmcube_in.close();
		return true;
	}

}//namespace numer end