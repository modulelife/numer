#include "complx_mat_test.h"

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
	{}

	Complex operator()(size_t i, size_t j) const {
		double re = re_min_ + re_range_ * ((double)j / jmax_);
		double im = im_min_ + im_range_ * ((imax_ - (double)i) / imax_);
		return Complex(re, im);
	}
};







constexpr size_t height = 2048, width = 2048;
constexpr double L = 1.0 * Pi;

void ComplxMatTest::run()
{

	std::cout << "\n>..ComplxMatTest: begin" << std::endl;
	BENCHMARK_BEGIN(complxmat_test);


	BENCHMARK_BEGIN(cplx_field);
	auto c_field = mat<Complex>::creat(height, width, ComplxPlainSampler(-L, -L, L, L, height, width));
	auto pha_shift = Complex::expi(-Pi / 3.0);
	c_field.modify_par([pha_shift](Complex& c) {
		constexpr Complex i{ 0.0, 1.0 };
		c = pha_shift * (5.0 * (c * c * c * c - i)) / (c * c * c * c * c - 1.0);
		});
	BENCHMARK_END(cplx_field);

	BENCHMARK_BEGIN(amp_field);
	auto amp_field = mat<double>::creat_par(c_field, [](const Complex& c) {return c.amplitude(); });
	BENCHMARK_END(amp_field);

	BENCHMARK_BEGIN(amp_img);
	auto amp_img = mat<RGB>::creat_par(
		amp_field, CompressedHeatMap(0.0, 10.0, Color::Vaporwave()));
	BENCHMARK_END(amp_img);

	BENCHMARK_BEGIN(complxrainbow_img);
	auto complxrainbow_img = mat<RGB>::creat_par(
		c_field, ComplxRainbowClr(5.0));
	BENCHMARK_END(complxrainbow_img);

	BENCHMARK_BEGIN(complxphase_img);
	auto complxphase_img = mat<RGB>::creat_par(
		c_field, ComplxPhaseClr(Color::GrayScale()));
	BENCHMARK_END(complxphase_img);

	stbi::ImageWriter<stbi::format::PNG> writer;

	BENCHMARK_BEGIN(write_img_file);
	writer.writeInto("./image/complx_mat/amplitude", amp_img[0], amp_img.ncols(), amp_img.nrows());
	BENCHMARK_END(write_img_file);

	writer.writeInto("./image/complx_mat/complex", complxrainbow_img[0], complxrainbow_img.ncols(), complxrainbow_img.nrows());
	writer.writeInto("./image/complx_mat/phase", complxphase_img[0], complxphase_img.ncols(), complxphase_img.nrows());

	std::cout << "\n>..ComplxMatTest: end" << std::endl;
	BENCHMARK_END(complxmat_test);
}
