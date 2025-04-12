//	numer_visualize.h	-*- C++20 -*-
#pragma once
//
//
//----------------------------content begins----------------------------

#include <cstdint>
#include <numer_constant.h>
#include <concepts>
#include <functional>

namespace numer { class Complex; }




namespace numer {


	struct RGB { uint8_t R, G, B; };
	struct RGBA { uint8_t R, G, B, A; };

	constexpr RGBA attachAlpha(uint8_t A_, RGB Clr_) {
		return RGBA{ Clr_.R, Clr_.G, Clr_.B, A_ };
	}

	constexpr RGB removeAlpha(RGBA Clr_) {
		return RGB{ Clr_.R, Clr_.G, Clr_.B};
	}

	constexpr RGBA mixAlpha(RGBA Lower_, RGBA Upper_) {

		const float upperA = static_cast<float>(Upper_.A) / 255.0f;
		const float lowerA = static_cast<float>(Lower_.A) / 255.0f;

		const float resultA = upperA + lowerA * (1.0f - upperA);

		const auto computeChannel = [&](uint8_t upperC, uint8_t lowerC) -> uint8_t {
			const float upperC_norm = static_cast<float>(upperC) / 255.0f;
			const float lowerC_norm = static_cast<float>(lowerC) / 255.0f;

			const float blended =
				(upperC_norm * upperA) +
				(lowerC_norm * lowerA) * (1.0f - upperA);

			const float final = (resultA > 0.0f) ? (blended / resultA) : 0.0f;

			return static_cast<uint8_t>(final * 255.0f + 0.5f);
		};

		const uint8_t R = computeChannel(Upper_.R, Lower_.R);
		const uint8_t G = computeChannel(Upper_.G, Lower_.G);
		const uint8_t B = computeChannel(Upper_.B, Lower_.B);

		const uint8_t A = static_cast<uint8_t>(resultA * 255.0f + 0.5f);

		return RGBA{ R, G, B, A };
	}


	//X should stay within [0, 1]
	using NormalizedColorMap = std::function<RGB(double X_)>;
	//brt should largely stay within [-1, 1]
	using BrightnessModifier = std::function<void(RGB& Clr_, double Brt_)>;


	template <class ClrImpl, typename AnyType>
	concept AnyColorizer = requires(ClrImpl Clr_impl_, AnyType x_) {
		{ Clr_impl_(x_) } -> std::same_as<RGB>;
	};

	template <class ClrImpl>
	concept RealNumColorizer = requires(ClrImpl Clr_impl_, double x_) {
		{ Clr_impl_(x_) } -> std::same_as<RGB>;
	};

	template <class ClrImpl>
	concept ComplxNumColorizer = requires(ClrImpl Clr_impl_, Complex x_) {
		{ Clr_impl_(x_) } -> std::same_as<RGB>;
	};


	struct Color {
		static NormalizedColorMap Cividis();
		static NormalizedColorMap Coolwarm();
		static NormalizedColorMap Glacier();
		static NormalizedColorMap GrayScale();
		static NormalizedColorMap Hue();
		static NormalizedColorMap Inferno();
		static NormalizedColorMap Plasma();
		static NormalizedColorMap Rainbow();
		static NormalizedColorMap Thermo();
		static NormalizedColorMap Vaporwave();
		static NormalizedColorMap Viridis();
		static NormalizedColorMap Zone_even();
		static NormalizedColorMap Zone_odd();
		static NormalizedColorMap Zone_half();
	private:
		Color() = delete;
	};

	class ReverseColor {
	private:
		NormalizedColorMap cl_map_;
	public:
		ReverseColor(NormalizedColorMap Color_map_) : cl_map_(Color_map_) {}

		RGB operator()(double x_) const {
			return cl_map_(1.0 - x_);
		}
	};


	class GrayScale {
	private:
		double min_thrs_;
		double max_thrs_;
	public:
		GrayScale() : min_thrs_(-10.0), max_thrs_(10.0) {}
		GrayScale(double Minimum_, double Maximum_)
			: min_thrs_(Minimum_), max_thrs_(Maximum_)
		{}

		double getMinThreshold() const { return min_thrs_; }
		double getMaxThreshold() const { return max_thrs_; }
		void setMinThreshold(double Minimum_) { min_thrs_ = Minimum_; }
		void setMaxThreshold(double Maximum_) { max_thrs_ = Maximum_; }
		uint8_t operator()(double x_) const;
		double operator()(uint8_t ch_) const;
	};

	class LinearHeatMap {
	private:
		double min_thrs_;
		double max_thrs_;
		NormalizedColorMap cl_map_;
	public:
		LinearHeatMap();
		LinearHeatMap(double Minimum_, double Maximum_, NormalizedColorMap Color_map_)
			: min_thrs_(Minimum_), max_thrs_(Maximum_), cl_map_(Color_map_)
		{}

		double getMinThreshold() const { return min_thrs_; }
		double getMaxThreshold() const { return max_thrs_; }
		void setMinThreshold(double Minimum_) { min_thrs_ = Minimum_; }
		void setMaxThreshold(double Maximum_) { max_thrs_ = Maximum_; }
		RGB operator()(double x_) const;
	};

	class LogthmHeatMap {
	private:
		double min_thrs_;
		double max_thrs_;
		NormalizedColorMap cl_map_;
	public:
		LogthmHeatMap();
		LogthmHeatMap(double Minimum_, double Maximum_, NormalizedColorMap Color_map_)
			: min_thrs_(Minimum_), max_thrs_(Maximum_), cl_map_(Color_map_)
		{}

		double getMinThreshold() const { return min_thrs_; }
		double getMaxThreshold() const { return max_thrs_; }
		void setMinThreshold(double Minimum_) { min_thrs_ = Minimum_; }
		void setMaxThreshold(double Maximum_) { max_thrs_ = Maximum_; }
		RGB operator()(double x_) const;
	};

	class CompressedHeatMap {
	private:
		double min_thrs_;
		double max_thrs_;
		NormalizedColorMap cl_map_;
	public:
		CompressedHeatMap();
		CompressedHeatMap(double Minimum_, double Maximum_, NormalizedColorMap Color_map_)
			: min_thrs_(Minimum_), max_thrs_(Maximum_), cl_map_(Color_map_)
		{
		}

		double getMinThreshold() const { return min_thrs_; }
		double getMaxThreshold() const { return max_thrs_; }
		void setMinThreshold(double Minimum_) { min_thrs_ = Minimum_; }
		void setMaxThreshold(double Maximum_) { max_thrs_ = Maximum_; }
		RGB operator()(double x_) const;
	};

	class ComplxRainbowClr{
	private:
		double max_amp_;
	public:
		ComplxRainbowClr(): max_amp_(10.0) 
		{}
		ComplxRainbowClr(double Max_amplitude_)
			: max_amp_(Max_amplitude_)
		{}

		double getMinThreshold() const { return 0.0; }
		double getMaxThreshold() const { return max_amp_; }
		void setMinThreshold(double Minimum_) {}
		void setMaxThreshold(double Max_amplitude_) { max_amp_ = Max_amplitude_; }
		RGB operator()(Complex C_) const;
	};

	class ComplxPhaseClr {
	private:
		NormalizedColorMap cl_map_;
	public:
		ComplxPhaseClr();
		ComplxPhaseClr(NormalizedColorMap Color_map_) : cl_map_(Color_map_)
		{}

		double getMinThreshold() const { return 0.0; }
		double getMaxThreshold() const { return 2.0 * Pi; }
		void setMinThreshold(double Minimum_) {}
		void setMaxThreshold(double Maximum_) {}
		RGB operator()(Complex C_) const;
	};



}//namespace numer end