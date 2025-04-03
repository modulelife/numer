#include "mandelbrot_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_complex.h>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <stbi.h>//requires stbi library: https://github.com/modulelife/stbi



using namespace numer;


class ComplxPlainSampler {
private:
	double re_min_;
	double im_min_;
	double re_range_;
	double im_range_;
	double imax_;
	double jmax_;
public:
	ComplxPlainSampler(double Re_min, double Im_min, double Re_max, double Im_max, size_t imax, size_t jmax)
		: re_min_(Re_min), im_min_(Im_min), re_range_(Re_max - re_min_), im_range_(Im_max - im_min_),
		imax_((double)imax), jmax_((double)jmax)
	{
	}

	Complex operator()(size_t i, size_t j) const {
		double re = re_min_ + re_range_ * ((double)j / jmax_);
		double im = im_min_ + im_range_ * ((imax_ - (double)i) / imax_);
		return Complex(re, im);
	}
};


constexpr double cx = 0.0, cy = 0.0, d = 2.0, r = 1.2;
constexpr double lm = cx - d * r, rm = cx + d * r, dm = cy - d, um = cy + d;
constexpr size_t height = 14460;
constexpr size_t width = height * r;
constexpr unsigned max_iteration = 128;
constexpr double norm_sup = 2.4;




void MandelbrotTest::run()
{

	std::cout << "\n>..MandelbrotTest: begin" << std::endl;
	BENCHMARK_BEGIN(mandelbrot_test);

	auto c_plain = mat<Complex>::creat_par(height, width, ComplxPlainSampler(lm, dm, rm, um, height, width));

	auto iteration = mat<unsigned>::creat_par(c_plain,
		[](Complex c) {
			Complex z = 0.0;
			unsigned count = 0;
			while (z.amplitude() <= norm_sup && count < max_iteration) {
				z = z * z + c;
				++count;
			}
			return count;
		});

	auto fractal_img = mat<RGB>::creat_par(iteration, LinearHeatMap(0, max_iteration, Color::GrayScale()));

	stbi::ImageWriter<stbi::format::PNG> writer;

	writer.writeInto("./image/mandelbrot/mandelbrot", fractal_img[0], fractal_img.ncols(), fractal_img.nrows());


	std::cout << "\n>..MandelbrotTest: end" << std::endl;
	BENCHMARK_END(mandelbrot_test);
}
