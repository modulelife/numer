//	numer_visualize.cpp		-*- C++20 -*-
//
//
//----------------------------content begins----------------------------

#include <numer_visualize.h>

#include <cmath>
#include <numer_complex.h>







namespace numer {


	namespace {


		RGB clrAqueousBlue__(double X_);
		RGB clrCividis__(double X_);
		RGB clrCoolwarm__(double X_);
		RGB clrCoolTech__(double X_);
		RGB clrEmerald__(double X_);
		RGB clrGlacier__(double X_);
		RGB clrGoldenBlue__(double X_);
		RGB clrGrayScale__(double X_);
		RGB clrInferno__(double X_);
		RGB clrMist__(double X_);
		RGB clrMulberryTea__(double X_);
		RGB clrPlasma__(double X_);
		RGB clrRainbow__(double X_);
		RGB clrSandstone__(double X_);
		RGB clrThermo__(double X_);
		RGB clrTwilight__(double X_);
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

		inline double cubic__(double inflec, double a, double X_)
		{
			return a * X_ * X_ * X_ - 3.0 * a * inflec * X_ * X_ + (3.0 * a * inflec + 1.0 - a) * X_;
		}

		inline double quart__(double a, double b, double c, double X_)
		{
			return a * X_ * X_ * X_ * X_ + b * X_ * X_ * X_ + c * X_ * X_ + (1.0 - a - b - c) * X_;
		}

		inline double bump__(double k, double h, double l, double X_)
		{
			return X_ * (1.0 - X_) * (k * X_ * X_ + h * X_ + l);
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

		RGB clrAqueousBlue__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.18, t2 = 0.35, t3 = 0.52, t4 = 0.69, t5 = 0.86, t6 = 1.0;

			constexpr numer::RGB
				Clr0{ 0xfd, 0xfe, 0xff },
				Clr1{ 0xee, 0xf4, 0xfa },
				Clr2{ 0xd6, 0xe6, 0xef },
				Clr3{ 0xb3, 0xcc, 0xdc },
				Clr4{ 0x86, 0xa9, 0xbf },
				Clr5{ 0x57, 0x82, 0x99 },
				Clr6{ 0x32, 0x55, 0x64 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__(1.0 - (X_ - t5) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
		}

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

		RGB clrCoolTech__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.2, t2 = 0.35, t3 = 0.5,
				t4 = 0.65, t5 = 0.8, t6 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x0b, 0x0b, 0x2b },
				Clr1{ 0x2d, 0x0b, 0x64 },
				Clr2{ 0x4a, 0x16, 0x96 },
				Clr3{ 0x6a, 0x34, 0xc8 },
				Clr4{ 0x29, 0x9c, 0xda },
				Clr5{ 0x45, 0xe1, 0xd9 },
				Clr6{ 0xa5, 0xff, 0xf0 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__(1.0 - (X_ - t5) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
		}

		RGB clrEmerald__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.25, t2 = 0.5, t3 = 0.75, t4 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x3c, 0x7d, 0x8c },
				Clr1{ 0x4a, 0x9d, 0x76 },
				Clr2{ 0x60, 0xb4, 0x7b },
				Clr3{ 0x9a, 0xbf, 0x7b },
				Clr4{ 0xd1, 0xd9, 0x80 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else return Clr4;
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

		RGB clrGoldenBlue__(double X_) {
			constexpr double max_thrs = 1.0, min_thrs = 0.00;

			RGB result = { 0, 0, 0 };

			if (X_ < min_thrs) return result;
			if (X_ < max_thrs)
			{
				float y = curvS__(0.7, (X_ - min_thrs) / (max_thrs - min_thrs));
				result.R += 255 * cubic__(0.51, -2.1, y);
				result.G += 255 * cubic__(0.49, 2.2, y);
				result.B += 255 * cubic__(0.50, 7.0, y);
				return result;
			}
			return { 255, 255, 255 };
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

		RGB clrMist__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.2, t2 = 0.4, t3 = 0.6, t4 = 0.8, t5 = 1.0;

			constexpr numer::RGB
				Clr0{ 0xfb, 0xfc, 0xf9 },
				Clr1{ 0xe6, 0xf2, 0xeb },
				Clr2{ 0xc4, 0xe3, 0xd6 },
				Clr3{ 0x94, 0xc7, 0xb0 },
				Clr4{ 0x5a, 0xa5, 0x8e },
				Clr5{ 0x34, 0x65, 0x54 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else return Clr5;
		}

		RGB clrMulberryTea__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.2, t2 = 0.45, t3 = 0.7, t4 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x6a, 0x2a, 0x3d },
				Clr1{ 0x9e, 0x4b, 0x5e },
				Clr2{ 0xc7, 0x7e, 0x7f },
				Clr3{ 0xdc, 0xa8, 0x82 },
				Clr4{ 0xe9, 0xc1, 0x7a };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else return Clr4;
		}

		RGB clrPlasma__(double X_) {

			constexpr double t0 = 0.0, t1 = 0.11, t2 = 0.22, t3 = 0.33,
				t4 = 0.44, t5 = 0.55, t6 = 0.66, t7 = 0.77, t8 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x06, 0x00, 0x1d },
				Clr1{ 0x23, 0x02, 0x5d },
				Clr2{ 0x3c, 0x04, 0xa6 },
				Clr3{ 0x80, 0x00, 0x86 },
				Clr4{ 0xb7, 0x15, 0x40 },
				Clr5{ 0xe3, 0x3b, 0x12 },
				Clr6{ 0xff, 0x78, 0x00 },
				Clr7{ 0xff, 0xc4, 0x00 },
				Clr8{ 0xff, 0xfe, 0xd4 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__(1.0 - (X_ - t5) / (t6 - t5), Clr5, Clr6);
			else if (X_ < t7) return linearMix__(1.0 - (X_ - t6) / (t7 - t6), Clr6, Clr7);
			else if (X_ < t8) return linearMix__(1.0 - (X_ - t7) / (t8 - t7), Clr7, Clr8);
			else return Clr8;
		}

		RGB clrRainbow__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.14, t2 = 0.27, t3 = 0.39, t4 = 0.52,
				t5 = 0.64, t6 = 0.73, t7 = 0.82, t8 = 0.91, t9 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x0a, 0x00, 0x27 },
				Clr1{ 0x26, 0x00, 0x7a },
				Clr2{ 0x33, 0x0b, 0xc8 },
				Clr3{ 0x00, 0x7f, 0xf4 },
				Clr4{ 0x00, 0xd9, 0xf0 },
				Clr5{ 0x3a, 0xfc, 0xb7 },
				Clr6{ 0xa8, 0xff, 0x60 },
				Clr7{ 0xf1, 0xef, 0x1c },
				Clr8{ 0xff, 0x87, 0x00 },
				Clr9{ 0xff, 0x00, 0x53 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__(1.0 - (X_ - t5) / (t6 - t5), Clr5, Clr6);
			else if (X_ < t7) return linearMix__(1.0 - (X_ - t6) / (t7 - t6), Clr6, Clr7);
			else if (X_ < t8) return linearMix__(1.0 - (X_ - t7) / (t8 - t7), Clr7, Clr8);
			else if (X_ < t9) return linearMix__(1.0 - (X_ - t8) / (t9 - t8), Clr8, Clr9);
			else return Clr9;
		}

		RGB clrSandstone__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.12, t2 = 0.3, t3 = 0.5, t4 = 0.7, t5 = 0.85, t6 = 1.0;

			constexpr numer::RGB
				Clr0{ 0xff, 0xfd, 0xfa },
				Clr1{ 0xfd, 0xf6, 0xe9 },
				Clr2{ 0xfa, 0xec, 0xd7 },
				Clr3{ 0xf2, 0xdc, 0xbf },
				Clr4{ 0xdb, 0xc0, 0xa2 },
				Clr5{ 0xb1, 0x98, 0x86 },
				Clr6{ 0x82, 0x6c, 0x61 };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else if (X_ < t5) return linearMix__(1.0 - (X_ - t4) / (t5 - t4), Clr4, Clr5);
			else if (X_ < t6) return linearMix__(1.0 - (X_ - t5) / (t6 - t5), Clr5, Clr6);
			else return Clr6;
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

		RGB clrTwilight__(double X_) {
			constexpr double t0 = 0.0, t1 = 0.25, t2 = 0.5, t3 = 0.75, t4 = 1.0;

			constexpr numer::RGB
				Clr0{ 0x2A, 0x1B, 0x4D },
				Clr1{ 0x5D, 0x3F, 0x9D },
				Clr2{ 0x9A, 0x5A, 0xA8 },
				Clr3{ 0xDF, 0x9E, 0x5B },
				Clr4{ 0xFE, 0xD8, 0x9E };

			if (X_ < t0) return Clr0;
			else if (X_ < t1) return linearMix__(1.0 - (X_ - t0) / (t1 - t0), Clr0, Clr1);
			else if (X_ < t2) return linearMix__(1.0 - (X_ - t1) / (t2 - t1), Clr1, Clr2);
			else if (X_ < t3) return linearMix__(1.0 - (X_ - t2) / (t3 - t2), Clr2, Clr3);
			else if (X_ < t4) return linearMix__(1.0 - (X_ - t3) / (t4 - t3), Clr3, Clr4);
			else return Clr4;
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





	NormalizedColorMap Color::AqueousBlue() {
		return clrAqueousBlue__;
	}
	NormalizedColorMap Color::Cividis() {
		return clrCividis__;
	}
	NormalizedColorMap Color::Coolwarm() {
		return clrCoolwarm__;
	}
	NormalizedColorMap Color::CoolTech() {
		return clrCoolTech__;
	}
	NormalizedColorMap Color::Emerald() {
		return clrEmerald__;
	}
	NormalizedColorMap Color::Glacier() {
		return clrGlacier__;
	}
	NormalizedColorMap Color::GoldenBlue() {
		return clrGoldenBlue__;
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
	NormalizedColorMap Color::Mist() {
		return clrMist__;
	}
	NormalizedColorMap Color::MulberryTea() {
		return clrMulberryTea__;
	}
	NormalizedColorMap Color::Plasma() {
		return clrPlasma__;
	}
	NormalizedColorMap Color::Rainbow() {
		return clrRainbow__;
	}
	NormalizedColorMap Color::Sandstone() {
		return clrSandstone__;
	}
	NormalizedColorMap Color::Thermo() {
		return clrThermo__;
	}
	NormalizedColorMap Color::Twilight() {
		return clrTwilight__;
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

