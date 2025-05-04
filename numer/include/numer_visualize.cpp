//	numer_visualize.cpp		-*- C++20 -*-
//
//
//----------------------------content begins----------------------------

#include <numer_visualize.h>

#include <cmath>
#include <numer_complex.h>







namespace numer {


	namespace {


		RGB clrCividis__(double X_);
		RGB clrCoolwarm__(double X_);
		RGB clrGlacier__(double X_);
		RGB clrGrayScale__(double X_);
		RGB clrInferno__(double X_);
		RGB clrPlasma__(double X_);
		RGB clrRainbow__(double X_);
		RGB clrThermo__(double X_);
		RGB clrVaporwave__(double X_);
		RGB clrViridis__(double X_);

		RGB HSL_to_RGB__(double h_, double l_, double s_);
		RGB clrHue__(double X_);
		RGB clrWhiteBlack_cyc__(double X_);
		RGB clrBlackWhite_cyc__(double X_);
		RGB clrBlackWhite__(double X_);

		void linearBrightness__(RGB& clr_, double bri_);
		void tanhBrightness__(RGB& clr_, double bri_);


		//======================================================================
		//							Color Functions
		//======================================================================

		inline double nGamma__(double gamma_, double X_) {
			return pow(X_, gamma_);
		}

		inline double rGamma__(double gamma_, double X_) {
			return 1.0 - pow(1.0 - X_, gamma_);
		}

		inline double curvS__(double gamma_, double X_)
		{
			if (X_ < 0.5)
				return 0.5 * pow(2.0 * X_, gamma_);
			else
				return 1.0 - 0.5 * pow(2.0 * (1.0 - X_), gamma_);
		}

		inline double bump__(double X_)
		{
			return (1.0 - cos(2.0 * Pi * X_)) / 2.0;
		}

		inline RGB linearMix__(double ratio_1, const RGB& first, const RGB& second)
		{
			return RGB{
				static_cast<uint8_t>(ratio_1 * first.R + (1.0 - ratio_1) * second.R),
				static_cast<uint8_t>(ratio_1 * first.G + (1.0 - ratio_1) * second.G),
				static_cast<uint8_t>(ratio_1 * first.B + (1.0 - ratio_1) * second.B)
			};
		}

		//---------------------------- basic colormaps ----------------------------
		/* for convenience, inside any NormalizedColorMap

			constexpr double t1 = 0.0, t2 = 0.2, t3 = 0.4, t4 = 0.6, t5 = 0.8, t6 = 1.0;
			constexpr RGB
				Clr1{ 0x, 0x, 0x },
				Clr2{ 0x, 0x, 0x },
				Clr3{ 0x, 0x, 0x },
				Clr4{ 0x, 0x, 0x },
				Clr5{ 0x, 0x, 0x },
				Clr6{ 0x, 0x, 0x };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__((t4 - X_) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__((t5 - X_) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__((t6 - X_) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
		*/

		RGB clrCividis__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.5, t3 = 1.0;
			constexpr RGB
				Clr1{ 0x00, 0x22, 0x50 },
				Clr2{ 0x6b, 0x6a, 0x70 },
				Clr3{ 0xff, 0xe8, 0x3a };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else return Clr3;
		}

		RGB clrCoolwarm__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.5, t3 = 1.0;
			constexpr RGB
				Clr1{ 0x4a, 0x91, 0xcc },
				Clr2{ 0xcf, 0xcf, 0xcf },
				Clr3{ 0xd3, 0x62, 0x4b };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else return Clr3;
		}

		RGB clrGlacier__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.4, t3 = 1.0;
			constexpr RGB
				Clr1{ 0x10, 0x20, 0x3b },
				Clr2{ 0x27, 0x51, 0x79 },
				Clr3{ 0x9c, 0xb5, 0xd2 };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else return Clr3;
		}

		RGB clrGrayScale__(double X_) {
			if (X_ < 0.0) return { 0, 0, 0 };
			if (X_ > 1.0) return { 255, 255, 255 };
			return RGB{
				static_cast<uint8_t>(255.0 * X_),
				static_cast<uint8_t>(255.0 * X_),
				static_cast<uint8_t>(255.0 * X_)
			};
		}

		RGB clrInferno__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.33, t3 = 0.55, t4 = 0.8, t5 = 1.0;
			constexpr RGB
				Clr1{ 0x04, 0x00, 0x00 },
				Clr2{ 0x78, 0x1c, 0x71 },
				Clr3{ 0xd0, 0x4f, 0x3c },
				Clr4{ 0xfb, 0xb6, 0x1a },
				Clr5{ 0xff, 0xf5, 0xa4 };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__((t4 - X_) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__((t5 - X_) / (t5 - t4), Clr4, Clr5);
			else return Clr5;
		}

		RGB clrPlasma__(double X_) {


			constexpr double t1 = 0.0, t2 = 0.45, t3 = 0.6, t4 = 0.85, t5 = 1.0;
			constexpr RGB
				Clr1{ 0x0e, 0x07, 0x89 },
				Clr2{ 0xb3, 0x30, 0x8d },
				Clr3{ 0xee, 0x78, 0x50 },
				Clr4{ 0xf9, 0xd6, 0x24 },
				Clr5{ 0xf3, 0xea, 0x84 };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__((t4 - X_) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__((t5 - X_) / (t5 - t4), Clr4, Clr5);
			else return Clr5;
		}

		RGB clrRainbow__(double X_) {
			if (X_ < 0.0) return { 68, 51, 154 };
			if (X_ > 1.0) return { 231, 51, 54 };
			RGB result{ 68, 51, 154 };
			result.R += static_cast<uint8_t>(163 * curvS__(3.0, X_));
			result.G += static_cast<uint8_t>(144 * bump__(X_));
			result.B += static_cast<uint8_t>(-100 * curvS__(3.0, X_));
			return result;
		}

		RGB clrThermo__(double X) {

			constexpr double t0 = 0.0, t1 = 0.4, t2 = 0.55, t3 = 0.7, t4 = 0.8, t5 = 0.9, t6 = 1.0;
			constexpr RGB
				Clr0{ 0x00, 0x00, 0x00 },
				Clr1{ 0x71, 0x21, 0x0c },
				Clr2{ 0xb3, 0x3b, 0x1b },
				Clr3{ 0xe6, 0x98, 0x42 },
				Clr4{ 0xe6, 0xd4, 0x65 },
				Clr5{ 0xea, 0xe6, 0xbc },
				Clr6{ 0xf1, 0xfd, 0xff };


			if (X < t0) return Clr0;
			else if (X < t1) return linearMix__((t1 - X) / (t1 - t0), Clr0, Clr1);
			else if (X < t2) return linearMix__((t2 - X) / (t2 - t1), Clr1, Clr2);
			else if (X < t3) return linearMix__((t3 - X) / (t3 - t2), Clr2, Clr3);
			else if (X < t4) return linearMix__((t4 - X) / (t4 - t3), Clr3, Clr4);
			else if (X < t5) return linearMix__((t5 - X) / (t5 - t4), Clr4, Clr5);
			else if (X < t6) return linearMix__((t6 - X) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
		}

		RGB clrVaporwave__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.2, t3 = 0.4, t4 = 0.6, t5 = 0.8, t6 = 1.0;
			constexpr RGB
				Clr1{ 0x00, 0x00, 0x00 },
				Clr2{ 0x12, 0x27, 0x3f },
				Clr3{ 0x44, 0x42, 0x85 },
				Clr4{ 0x71, 0x46, 0xc1 },
				Clr5{ 0xd9, 0x5f, 0x9e },
				Clr6{ 0xf7, 0xdc, 0xa2 };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__((t4 - X_) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__((t5 - X_) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__((t6 - X_) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
		}

		RGB clrViridis__(double X_) {

			constexpr double t1 = 0.0, t2 = 0.25, t3 = 0.5, t4 = 0.75, t5 = 1.0;
			constexpr RGB
				Clr1{ 0x42, 0x01, 0x53 },
				Clr2{ 0x3e, 0x53, 0x8c },
				Clr3{ 0x1e, 0xa3, 0x84 },
				Clr4{ 0x91, 0xd7, 0x43 },
				Clr5{ 0xfe, 0xe7, 0x1f };


			if (X_ < t1) return Clr1;
			else if (X_ < t2) return linearMix__((t2 - X_) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__((t3 - X_) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__((t4 - X_) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__((t5 - X_) / (t5 - t4), Clr4, Clr5);
			else return Clr5;
		}


		//---------------------------- phase color ----------------------------
		inline double floatingMod__(double num, double modu) {
			if (num >= 0.0) {
				if (num < modu)return num;
				if (num > modu)return num - modu * std::floor(num / modu);
				return 0.0;
			}
			else {
				if (-num < modu)return modu + num;
				if (-num > modu)return num + modu * std::floor(-num / modu);
				return 0.0;
			}
		}

		RGB HSL_to_RGB__(double h_, double l_, double s_) {

			h_ /= 360.0;
			l_ /= 100.0;
			s_ /= 100.0;

			double c = (1 - std::abs(2.0 * l_ - 1.0)) * s_;
			double x = c * (1 - std::abs(floatingMod__(h_ * 6.0, 2.0) - 1.0));
			double m = l_ - c / 2;

			double r, g, b;
			if (0 <= h_ && h_ < 1 / 6.0) {
				r = c; g = x; b = 0;
			}
			else if (1 / 6.0 <= h_ && h_ < 1 / 3.0) {
				r = x; g = c; b = 0;
			}
			else if (1 / 3.0 <= h_ && h_ < 1 / 2.0) {
				r = 0; g = c; b = x;
			}
			else if (1 / 2.0 <= h_ && h_ < 2 / 3.0) {
				r = 0; g = x; b = c;
			}
			else if (2 / 3.0 <= h_ && h_ < 5 / 6.0) {
				r = x; g = 0; b = c;
			}
			else {
				r = c; g = 0; b = x;
			}

			uint8_t r_final = static_cast<uint8_t>(std::round((r + m) * 255));
			uint8_t g_final = static_cast<uint8_t>(std::round((g + m) * 255));
			uint8_t b_final = static_cast<uint8_t>(std::round((b + m) * 255));

			return RGB{ r_final, g_final, b_final };
		}



		RGB clrHue__(double X_) {
			if (X_ < 0.0) return HSL_to_RGB__(0.0, 50.0, 100.0);
			if (X_ > 1.0) return HSL_to_RGB__(270.0, 50.0, 100.0);
			return HSL_to_RGB__(X_ * 270.0, 50.0, 100.0);
		}

		RGB clrWhiteBlack_cyc__(double X_) {
			if (X_ < 0.0 || X_ > 1.0) return { 255, 255, 255 };
			RGB result{};
			result.R = static_cast<uint8_t>(255 * (1.0 + cos(2.0 * Pi * X_)) / 2.0);
			result.G = static_cast<uint8_t>(255 * (1.0 + cos(2.0 * Pi * X_)) / 2.0);
			result.B = static_cast<uint8_t>(255 * (1.0 + cos(2.0 * Pi * X_)) / 2.0);
			return result;
		}

		RGB clrBlackWhite_cyc__(double X_) {
			if (X_ < 0.0 || X_ > 1.0) return { 0, 0, 0 };
			RGB result{};
			result.R = static_cast<uint8_t>(255 * (1.0 - cos(2.0 * Pi * X_)) / 2.0);
			result.G = static_cast<uint8_t>(255 * (1.0 - cos(2.0 * Pi * X_)) / 2.0);
			result.B = static_cast<uint8_t>(255 * (1.0 - cos(2.0 * Pi * X_)) / 2.0);
			return result;
		}

		RGB clrBlackWhite__(double X_) {
			if (X_ < 0.0) return { 0, 0, 0 };
			if (X_ > 1.0) return { 255, 255, 255 };
			RGB result{};
			result.R = static_cast<uint8_t>(255 * (1.0 - cos(Pi * X_)) / 2.0);
			result.G = static_cast<uint8_t>(255 * (1.0 - cos(Pi * X_)) / 2.0);
			result.B = static_cast<uint8_t>(255 * (1.0 - cos(Pi * X_)) / 2.0);
			return result;
		}

		//---------------------------- brightness modifier ----------------------------

		inline uint8_t clampChannl__(int32_t channl) {
			if (channl < 0) return 0U;
			if (channl > 255) return 255U;
			return static_cast<uint8_t>(channl);
		}


		void linearBrightness__(RGB& clr_, double bri_) {
			clr_.R = clampChannl__(static_cast<int32_t>(clr_.R) + static_cast<int32_t>(255.0 * bri_));
			clr_.G = clampChannl__(static_cast<int32_t>(clr_.G) + static_cast<int32_t>(255.0 * bri_));
			clr_.B = clampChannl__(static_cast<int32_t>(clr_.B) + static_cast<int32_t>(255.0 * bri_));
		}

		void tanhBrightness__(RGB& clr_, double bri_) {
			clr_.R = clampChannl__(static_cast<int32_t>(clr_.R) + static_cast<int32_t>(255.0 * tanh(bri_)));
			clr_.G = clampChannl__(static_cast<int32_t>(clr_.G) + static_cast<int32_t>(255.0 * tanh(bri_)));
			clr_.B = clampChannl__(static_cast<int32_t>(clr_.B) + static_cast<int32_t>(255.0 * tanh(bri_)));
		}


	}//anonymous namespace end





	NormalizedColorMap Color::Cividis() {
		return clrCividis__;
	}
	NormalizedColorMap Color::Coolwarm() {
		return clrCoolwarm__;
	}
	NormalizedColorMap Color::Glacier() {
		return clrGlacier__;
	}
	NormalizedColorMap Color::GrayScale() {
		return clrGrayScale__;
	}
	NormalizedColorMap Color::Hue() {
		return clrHue__;
	}
	NormalizedColorMap Color::Inferno() {
		return clrInferno__;
	}
	NormalizedColorMap Color::Plasma() {
		return clrPlasma__;
	}
	NormalizedColorMap Color::Rainbow() {
		return clrRainbow__;
	}
	NormalizedColorMap Color::Thermo() {
		return clrThermo__;
	}
	NormalizedColorMap Color::Vaporwave() {
		return clrVaporwave__;
	}
	NormalizedColorMap Color::Viridis() {
		return clrViridis__;
	}
	NormalizedColorMap Color::Zone_even() {
		return clrWhiteBlack_cyc__;
	}
	NormalizedColorMap Color::Zone_odd() {
		return clrBlackWhite_cyc__;
	}
	NormalizedColorMap Color::Zone_half() {
		return clrBlackWhite__;
	}

	uint8_t GrayScale::operator()(double x_) const
	{
		if (x_ < min_thrs_) return 0U;
		if (x_ < max_thrs_)
		{
			double y = (x_ - min_thrs_) / (max_thrs_ - min_thrs_);
			return static_cast<uint8_t>(255.0 * y);
		}
		return 255U;
	}

	double GrayScale::operator()(uint8_t ch_) const {
		double volume = static_cast<double>(ch_) / 255.0;
		return min_thrs_ + volume * (max_thrs_ - min_thrs_);
	}


	LinearHeatMap::LinearHeatMap()
		:min_thrs_(-10.0), max_thrs_(10.0), cl_map_(clrViridis__)
	{
	}

	RGB LinearHeatMap::operator()(double x_) const {
		if (x_ < min_thrs_) return cl_map_(0.0);
		if (x_ < max_thrs_){
			double X = (x_ - min_thrs_) / (max_thrs_ - min_thrs_);
			return cl_map_(X);
		}
		return cl_map_(1.0);
	}


	LogthmHeatMap::LogthmHeatMap()
		:min_thrs_(-10.0), max_thrs_(10.0), cl_map_(clrViridis__)
	{
	}

	RGB LogthmHeatMap::operator()(double x_) const {
		double y = log(x_);
		if (y < min_thrs_) return cl_map_(0.0);
		if (y < max_thrs_) {
			double X = (y - min_thrs_) / (max_thrs_ - min_thrs_);
			return cl_map_(X);
		}
		return cl_map_(1.0);
	}


	CompressedHeatMap::CompressedHeatMap()
		:min_thrs_(-10.0), max_thrs_(10.0), cl_map_(clrViridis__)
	{
	}

	RGB CompressedHeatMap::operator()(double x_) const {
		double diff = max_thrs_ - min_thrs_;
		double y = 3.0 * (x_ - 0.5 * diff - min_thrs_) / diff;
		double X = (1.0 + tanh(y)) / 2.0;
		return cl_map_(X);
	}


	RGB ComplxRainbowClr::operator()(Complex C_) const {
		double H = 360.0 * C_.phase() / (2.0 * Pi);
		double L = 100.0 * 0.5 * tanh(C_.amplitude() / max_amp_);
		return HSL_to_RGB__(H, L, 100.0);
	}


	ComplxPhaseClr::ComplxPhaseClr() : cl_map_(clrGlacier__) {}

	RGB ComplxPhaseClr::operator()(Complex C_) const {
		return cl_map_(C_.phase() / (2.0 * Pi));
	}

	


}//namespace numer end

