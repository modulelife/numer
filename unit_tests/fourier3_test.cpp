#include "fourier3_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_complex.h>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_fourier.h>
#include <stbi.h>//requires stbi library: https://github.com/modulelife/stbi

using namespace numer;




static inline double GetRReal(const numer::RGB& Pixel) {
	return (double)Pixel.R / 255.0;
}

static inline double GetGReal(const numer::RGB& Pixel) {
	return (double)Pixel.G / 255.0;
}

static inline double GetBReal(const numer::RGB& Pixel) {
	return (double)Pixel.B / 255.0;
}

static inline numer::RGB SetRChannl(numer::RGB Pixel, uint8_t ch) {
	Pixel.R = ch;
	return Pixel;
}

static inline numer::RGB SetGChannl(numer::RGB Pixel, uint8_t ch) {
	Pixel.G = ch;
	return Pixel;
}

static inline numer::RGB SetBChannl(numer::RGB Pixel, uint8_t ch) {
	Pixel.B = ch;
	return Pixel;
}


static inline double GetAmplitude(const numer::Complex& Complx) {
	return Complx.amplitude();
}

static inline double GetSquaredAmp(const numer::Complex& Complx) {
	return Complx.sqrdAmp();
}

static inline numer::Complex ComplxRealMutip(const numer::Complex& Complx, const double& Real) {
	return Complx * Real;
}

static inline numer::Complex ComplxAdd(const numer::Complex& C1, const numer::Complex& C2) {
	return C1 + C2;
}

static inline numer::Complex ComplxMutip(const numer::Complex& C1, const numer::Complex& C2) {
	return C1 * C2;
}

static inline numer::Complex GetConjugate(const numer::Complex& Complx) {
	return Complx.conj();
}

static inline double GetPhase(const numer::Complex& Complx) {
	return Complx.phase();
}

static mat<Complex> centralize(const mat<Complex>& Spectr) {

	size_t X = Spectr.ncols();
	size_t Y = Spectr.nrows();
	mat<Complex> xcentered(Y, X);
	for (size_t i = 0; i < Y; ++i) {
		auto src = Spectr.row(i).ccycle_from(X / 2);
		auto des = xcentered.row(i).begin();
		for (size_t j = 0; j < X; ++j) {
			*des++ = *src++;
		}
	}
	mat<Complex> centered(Y, X);
	for (size_t j = 0; j < X; ++j) {
		auto src = xcentered.col(j).ccycle_from(Y / 2);
		auto des = centered.col(j).begin();
		for (size_t i = 0; i < Y; ++i) {
			*des++ = *src++;
		}
	}

	return centered;
}


void Fourier3Test::run()
{
	stbi::ImageWriter<stbi::format::PNG> writer;
	stbi::ImageLoader loader;

	std::cout << "\n>..Fourier3Test: begin" << std::endl;
	BENCHMARK_BEGIN(fourier3_test);

	loader.load("./image/fourier/testimg/neon.jpg", stbi::LOAD_RGB);
	mat<RGB> img(loader.height(), loader.width());
	loader.putInto(img.begin());

	writer.writeInto("./image/fourier/original", img.cbegin(), img.ncols(), img.nrows());

	auto imgr_num = mat<double>::creat_par(img, GetRReal);
	mat<Complex> imgr_cplx(img.nrows(), img.ncols(), Complex::identity());
	imgr_cplx.overlay_par(imgr_num, ComplxRealMutip);

	auto imgg_num = mat<double>::creat_par(img, GetGReal);
	mat<Complex> imgg_cplx(img.nrows(), img.ncols(), Complex::identity());
	imgg_cplx.overlay_par(imgg_num, ComplxRealMutip);

	auto imgb_num = mat<double>::creat_par(img, GetBReal);
	mat<Complex> imgb_cplx(img.nrows(), img.ncols(), Complex::identity());
	imgb_cplx.overlay_par(imgb_num, ComplxRealMutip);

	BENCHMARK_BEGIN(ft);
	fft2d_ortho_par(imgr_cplx);
	fft2d_ortho_par(imgg_cplx);
	fft2d_ortho_par(imgb_cplx);
	BENCHMARK_END(ft);
	
	auto center_rspectr = centralize(imgr_cplx);
	auto center_gspectr = centralize(imgg_cplx);
	auto center_bspectr = centralize(imgb_cplx);

	BENCHMARK_BEGIN(ift);
	ifft2d_ortho_par(imgr_cplx);
	ifft2d_ortho_par(imgg_cplx);
	ifft2d_ortho_par(imgb_cplx);
	BENCHMARK_END(ift);


	auto amp_rspectr = mat<double>::creat_par(center_rspectr, GetAmplitude);
	auto amp_gspectr = mat<double>::creat_par(center_gspectr, GetAmplitude);
	auto amp_bspectr = mat<double>::creat_par(center_bspectr, GetAmplitude);

	auto pha_rspectr = mat<double>::creat_par(center_rspectr, GetPhase);
	auto pha_gspectr = mat<double>::creat_par(center_gspectr, GetPhase);
	auto pha_bspectr = mat<double>::creat_par(center_bspectr, GetPhase);

	auto imgr_re_num = mat<double>::creat_par(imgr_cplx, GetAmplitude);
	auto imgg_re_num = mat<double>::creat_par(imgg_cplx, GetAmplitude);
	auto imgb_re_num = mat<double>::creat_par(imgb_cplx, GetAmplitude);


	auto ramp = mat<uint8_t>::creat_par(amp_rspectr, GrayScale(0.0, 1.0));
	auto gamp = mat<uint8_t>::creat_par(amp_gspectr, GrayScale(0.0, 1.0));
	auto bamp = mat<uint8_t>::creat_par(amp_bspectr, GrayScale(0.0, 1.0));

	auto rpha = mat<uint8_t>::creat_par(pha_rspectr, GrayScale(0.0, 2.0 * Pi));
	auto gpha = mat<uint8_t>::creat_par(pha_gspectr, GrayScale(0.0, 2.0 * Pi));
	auto bpha = mat<uint8_t>::creat_par(pha_bspectr, GrayScale(0.0, 2.0 * Pi));

	auto rch_re = mat<uint8_t>::creat_par(imgr_re_num, GrayScale(0.0, 1.0));
	auto gch_re = mat<uint8_t>::creat_par(imgg_re_num, GrayScale(0.0, 1.0));
	auto bch_re = mat<uint8_t>::creat_par(imgb_re_num, GrayScale(0.0, 1.0));



	mat<RGB> amp_img(img.nrows(), img.ncols());
	amp_img.overlay(ramp, SetRChannl);
	amp_img.overlay(gamp, SetGChannl);
	amp_img.overlay(bamp, SetBChannl);
	
	mat<RGB> pha_img(img.nrows(), img.ncols());
	pha_img.overlay(rpha, SetRChannl);
	pha_img.overlay(gpha, SetGChannl);
	pha_img.overlay(bpha, SetBChannl);
	
	mat<RGB> img_re(img.nrows(), img.ncols());
	img_re.overlay(rch_re, SetRChannl);
	img_re.overlay(gch_re, SetGChannl);
	img_re.overlay(bch_re, SetBChannl);


	writer.writeInto("./image/fourier/amp_spectr3", amp_img.cbegin(), amp_img.ncols(), amp_img.nrows());
	writer.writeInto("./image/fourier/pha_spectr3", pha_img.cbegin(), pha_img.ncols(), amp_img.nrows());
	writer.writeInto("./image/fourier/reconstructed3", img_re.cbegin(), img_re.ncols(), img_re.nrows());




	std::cout << "\n>..Fourier3Test: end" << std::endl;
	BENCHMARK_END(fourier3_test);
}