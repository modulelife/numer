//	numer_fourier.h	-*- C_++20 -*-
#pragma once
//@brief: fourier transform algorithm header
//
//@functions:
//  void fft_ortho(...)
//  void ifft_ortho(...)
//								------------ " symmetrical normalization FT/IFT, physists like this"
// 
//  void fft(...)
//  void ifft(...)
//								------------ " well-known normalization FT/IFT, one track minded"
//  
//	
//@description:
//	Single sequense fft algorithm implementation.
//	In-place ver. can be used for any Array(requires operator[]) of numer::Complex.
//  Non-mutating ver. however, taking in an iterator, can fit for any sequence of integers, floating
//  points or numer::Complex, it then ouputs a sequence of numer::Complex as a result.
//
//@usage:
//
//----------------------------content begins----------------------------

#include <numer_constant.h>
#include <numer_complex.h>
#include <numer_common.h>
#include <numer_mat.h>
#include <bit>
#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>
#include <execution>


namespace numer {


    //you shouldn't use these
	namespace detail_{


		template <typename Uint>
        inline Uint bit_reverse__(Uint x_, int low_bits_) {
            Uint reversed = 0;
            for (int i = 0; i < low_bits_; ++i) {
                reversed = (reversed << 1) | (x_ & 1);
                x_ >>= 1;
            }
            return reversed;
        }


        //hey! you shouldn't use this
		template<bool Inverse_, class Array>
		inline void fft_radix2__(Array&& seq_cplx_, const size_t len_)
        {
            constexpr double sign = Inverse_ ? 1.0 : -1.0;

            thread_local std::vector<size_t> pos(len_);
            thread_local bool init = false;

            if (pos.size() != len_) init = false;
            if (pos.size() < len_) pos.resize(len_);
            if (!init) {
                const int log2n = std::countr_zero(len_);
                for (size_t i = 0; i < len_; ++i) {
                    pos[i] = bit_reverse__(i, log2n);
                }
            }

            // in-place rearrangement, good for cache hit
            for (size_t i = 0; i < len_; ++i) {
                size_t rev = pos[i];
                if (i < rev) std::swap(seq_cplx_[i], seq_cplx_[rev]);
            }

            // staged butterfly operation
            for (size_t m = 2; m <= len_; m *= 2) {
                const size_t m2 = m / 2;
                const Complex wm = Complex::expi(-2.0 * Pi * sign / m);

                for (size_t k = 0; k < len_; k += m) {
                    Complex w = 1.0;
#pragma omp simd
                    for (size_t j = 0; j < m2; ++j) {

                        const Complex t = w * seq_cplx_[k + j + m2];
                        const Complex u = seq_cplx_[k + j];

                        seq_cplx_[k + j] = u + t;
                        seq_cplx_[k + j + m2] = u - t;

                        w *= wm;
                    }
                }
            }
        }


        //hey! you shouldn't use this
        template<bool Inverse_, class Array>
        inline void fft_bluestein__(Array&& seq_cplx_, const size_t len_)
        {
            constexpr double sign = Inverse_ ? -1.0 : 1.0;
            const size_t len_ex = to_pow2_up(2 * len_ - 1);
            const double N = static_cast<double>(len_);
            const double M = static_cast<double>(len_ex);

            thread_local std::vector<Complex> a(len_ex);
            thread_local std::vector<Complex> b(len_ex);
            thread_local std::vector<Complex> w(len_ex);
            thread_local bool init = false;

            if(w.size() != len_ex) init = false;
            if (w.size() < len_ex) {
                a.resize(len_ex);
                b.resize(len_ex);
                w.resize(len_ex);
            }
            if (!init) {
                for (size_t i = 0; i < len_; ++i) {
                    const double n = static_cast<double>(i);
                    w[i] = Complex::expi(sign * Pi * n * n / N);
                }
                for (size_t i = len_; i < len_ex; ++i) {
                    const double n = static_cast<double>(i);
                    w[i] = Complex::expi(sign * Pi * (M - n) * (M - n) / N);
                }
                init = true;
            }

#pragma omp simd
            for (size_t i = 0; i < len_; ++i) {
                a[i] = seq_cplx_[i] * w[i].conj();
                b[i] = w[i];
            }
#pragma omp simd
            for (size_t i = len_; i < len_ex; ++i) {
                a[i] = 0.0;
                b[i] = w[i];
            }

            fft_radix2__<false>(a, len_ex);
            fft_radix2__<false>(b, len_ex);
#pragma omp simd
            for (size_t i = 0; i < len_ex; ++i) {
                a[i] *= b[i];
            }
            fft_radix2__<true>(a, len_ex);
            double normalizer = static_cast<double>(len_ex);
#pragma omp simd
            for (size_t i = 0; i < len_; ++i) {
                a[i] /= normalizer;
            }

#pragma omp simd
            for (size_t i = 0; i < len_; ++i) {
                double k = static_cast<double>(i);
                seq_cplx_[i] = a[i] * w[i].conj();
            }
        }


	}//namespace detail_ end


    //in-place FFT
    //normalize factor N^0.5 for FT and IFT
    template<class Array>
    inline void fft_ortho(Array&& Seq_complx_, const size_t Len_)
    {
        if (is_pow2(Len_)) detail_::fft_radix2__<false>(Seq_complx_, Len_);
        else detail_::fft_bluestein__<false>(Seq_complx_, Len_);
        double normalizer = sqrt(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            Seq_complx_[i] /= normalizer;
        }
    }

    //in-place FFT
    //normalize factor N^0.5 for FT and IFT
    template<class Array>
    inline void ifft_ortho(Array&& Seq_complx_, const size_t Len_)
    {
        if (is_pow2(Len_)) detail_::fft_radix2__<true>(Seq_complx_, Len_);
        else detail_::fft_bluestein__<true>(Seq_complx_, Len_);
        double normalizer = sqrt(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            Seq_complx_[i] /= normalizer;
        }
    }


    //in-place FFT
    //normalize factor 1 for FT and N for IFT
    template<class Array>
    inline void fft(Array&& Seq_complx_, const size_t Len_)
    {
        if (is_pow2(Len_)) detail_::fft_radix2__<false>(Seq_complx_, Len_);
        else detail_::fft_bluestein__<false>(Seq_complx_, Len_);
    }

    //in-place FFT
    //normalize factor 1 for FT and N for IFT
    template<class Array>
    inline void ifft(Array&& Seq_complx_, const size_t Len_)
    {
        if (is_pow2(Len_)) detail_::fft_radix2__<true>(Seq_complx_, Len_);
        else detail_::fft_bluestein__<true>(Seq_complx_, Len_);
        double normalizer = static_cast<double>(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            Seq_complx_[i] /= normalizer;
        }
    }


    //non-mutating FFT
    //normalize factor N^0.5 for FT and IFT
    template<class InIt, class OutIt>
    inline void fft_ortho(InIt First_, const size_t Len_, OutIt Dest_complx_)
    {
        std::vector<Complex> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = *First_++;
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<false>(seq_buf, Len_);
        else detail_::fft_bluestein__<false>(seq_buf, Len_);
        double normalizer = sqrt(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] /= normalizer;
        }

        for (size_t i = 0; i < Len_; ++i) {
            *Dest_complx_++ = seq_buf[i];
        }
    }

    //non-mutating FFT
    //normalize factor N^0.5 for FT and IFT
    template<class InIt, class OutIt>
    inline void ifft_ortho(InIt First_, const size_t Len_, OutIt Dest_complx_)
    {
        std::vector<Complex> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = *First_++;
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<true>(seq_buf, Len_);
        else detail_::fft_bluestein__<true>(seq_buf, Len_);
        double normalizer = sqrt(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] /= normalizer;
        }

        for (size_t i = 0; i < Len_; ++i) {
            *Dest_complx_++ = seq_buf[i];
        }
    }


    //non-mutating FFT
    //normalize factor 1 for FT and N for IFT
    template<class InIt, class OutIt>
    inline void fft(InIt First_, const size_t Len_, OutIt Dest_complx_)
    {
        std::vector<Complex> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = *First_++;
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<false>(seq_buf, Len_);
        else detail_::fft_bluestein__<false>(seq_buf, Len_);

        for (size_t i = 0; i < Len_; ++i) {
            *Dest_complx_++ = seq_buf[i];
        }
    }

    //non-mutating FFT
    //normalize factor 1 for FT and N for IFT
    template<class InIt, class OutIt>
    inline void ifft(InIt First_, const size_t Len_, OutIt Dest_complx_)
    {
        std::vector<Complex> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = *First_++;
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<true>(seq_buf, Len_);
        else detail_::fft_bluestein__<true>(seq_buf, Len_);
        double normalizer = static_cast<double>(Len_);
#pragma omp simd
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] /= normalizer;
        }

        for (size_t i = 0; i < Len_; ++i) {
            *Dest_complx_++ = seq_buf[i];
        }
    }

	
	//in-place FFT for 2d complex field in mat<Complex>
    //normalize factor N^0.5 for FT and IFT
    template<class Alloc__>
    inline void fft2d_ortho_par(mat<Complex, Alloc__>& Field_complx_)
    {
        size_t X = Field_complx_.ncols();
        size_t Y = Field_complx_.nrows();
        std::vector<size_t> argsx(Y);
        std::vector<size_t> argsy(X);

        for (size_t i = 0; i < Y; ++i) {
            argsx[i] = i;
        }
        for (size_t j = 0; j < X; ++j) {
            argsy[j] = j;
        }

        std::for_each_n(std::execution::par, argsx.begin(), Y,
            [&](size_t id) { fft_ortho(Field_complx_[id], X); });
        std::for_each_n(std::execution::par, argsy.begin(), X,
            [&](size_t id) { fft_ortho(Field_complx_.col(id).begin(), Y); });
    }

    //in-place FFT for 2d complex field in mat<Complex>
    //normalize factor N^0.5 for FT and IFT
    template<class Alloc__>
    inline void ifft2d_ortho_par(mat<Complex, Alloc__>& Field_complx_)
    {
        size_t X = Field_complx_.ncols();
        size_t Y = Field_complx_.nrows();
        std::vector<size_t> argsx(Y);
        std::vector<size_t> argsy(X);

        for (size_t i = 0; i < Y; ++i) {
            argsx[i] = i;
        }
        for (size_t j = 0; j < X; ++j) {
            argsy[j] = j;
        }

        std::for_each_n(std::execution::par, argsx.begin(), Y,
            [&](size_t id) { ifft_ortho(Field_complx_[id], X); });
        std::for_each_n(std::execution::par, argsy.begin(), X,
            [&](size_t id) { ifft_ortho(Field_complx_.col(id).begin(), Y); });
    }

    //in-place FFT for 2d complex field in mat<Complex>
    //normalize factor 1 for FT and N for IFT
    template<class Alloc__>
    inline void fft2d_par(mat<Complex, Alloc__>& Field_complx_)
    {
        size_t X = Field_complx_.ncols();
        size_t Y = Field_complx_.nrows();
        std::vector<size_t> argsx(Y);
        std::vector<size_t> argsy(X);

        for (size_t i = 0; i < Y; ++i) {
            argsx[i] = i;
        }
        for (size_t j = 0; j < X; ++j) {
            argsy[j] = j;
        }

        std::for_each_n(std::execution::par, argsx.begin(), Y,
            [&](size_t id) { fft(Field_complx_[id], X); });
        std::for_each_n(std::execution::par, argsy.begin(), X,
            [&](size_t id) { fft(Field_complx_.col(id).begin(), Y); });
    }

    //in-place FFT for 2d complex field in mat<Complex>
    //normalize factor 1 for FT and N for IFT
    template<class Alloc__>
    inline void ifft2d_par(mat<Complex, Alloc__>& Field_complx_)
    {
        size_t X = Field_complx_.ncols();
        size_t Y = Field_complx_.nrows();
        std::vector<size_t> argsx(Y);
        std::vector<size_t> argsy(X);

        for (size_t i = 0; i < Y; ++i) {
            argsx[i] = i;
        }
        for (size_t j = 0; j < X; ++j) {
            argsy[j] = j;
        }

        std::for_each_n(std::execution::par, argsx.begin(), Y,
            [&](size_t id) { ifft(Field_complx_[id], X); });
        std::for_each_n(std::execution::par, argsy.begin(), X,
            [&](size_t id) { ifft(Field_complx_.col(id).begin(), Y); });
    }


}//namespace numer end


