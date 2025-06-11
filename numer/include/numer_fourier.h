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
#include <numer_cube.h>
#include <bit>
#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>
#include <execution>


namespace numer {


    class FFTFreq {
    private:
        long long seq_len_;
        double range_len_;

    public:
        FFTFreq(unsigned Seq_length_, double Range_lenght_)
            : seq_len_(Seq_length_), range_len_(Range_lenght_) {
        }

        double operator()(unsigned Idx_) const {
            long long id = static_cast<long long>(Idx_);
            long long half_seql = seq_len_ / 2;
            if (id <= half_seql) return 2.0 * Pi * static_cast<double>(id) / range_len_;
            else if (id < seq_len_) return 2.0 * Pi * static_cast<double>(id - seq_len_) / range_len_;
            else return 0.0;
        }
    };


    //you shouldn't use these
	namespace detail_{


		template<typename Uint>
        inline Uint bit_reverse__(Uint x_, int low_bits_) {
            Uint reversed = 0;
            for (int i = 0; i < low_bits_; ++i) {
                reversed = (reversed << 1) | (x_ & 1);
                x_ >>= 1;
            }
            return reversed;
        }

        template<typename T>
        concept has_cplx_multipy_assignment = requires(T obj, Complex c) {
            obj.operator*=(c);
        };


        //hey! you shouldn't use this
		template<bool Inverse_, class Array>
		inline void fft_radix2__(Array&& seq_cplx_, const size_t len_)
        {
            using Complex_Absorbing_T = std::decay_t<decltype(seq_cplx_[0])>;
            using Product_T = std::decay_t<decltype(std::declval<Complex_Absorbing_T>() * std::declval<Complex>())>;
            static_assert(std::is_same<Complex_Absorbing_T, Product_T>::value, "Element type can not carry out complex scalar multipication");


            constexpr double sign = Inverse_ ? -1.0 : 1.0;

            thread_local size_t prev_len = 0;
            thread_local std::vector<size_t> pos;
            thread_local bool init = false;

            if (prev_len != len_) init = false;
            if (pos.size() < len_) pos.resize(len_);
            if (!init) {
                const int log2n = std::countr_zero(len_);
                for (size_t i = 0; i < len_; ++i) {
                    pos[i] = bit_reverse__(i, log2n);
                }
                prev_len = len_;
                init = true;
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

                    for (size_t j = 0; j < m2; ++j) {

                        const Complex_Absorbing_T t = w * seq_cplx_[k + j + m2];
                        const Complex_Absorbing_T u = seq_cplx_[k + j];

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
            using Complex_Absorbing_T = std::decay_t<decltype(seq_cplx_[0])>;
            using Product_T = std::decay_t<decltype(std::declval<Complex_Absorbing_T>()* std::declval<Complex>())>;
            static_assert(std::is_same<Complex_Absorbing_T, Product_T>::value, "Element type can not carry out complex scalar multipication");
            static_assert(has_cplx_multipy_assignment<Complex_Absorbing_T>, "Element type can not carry out complex scalar assignment multipication");


            constexpr double sign = Inverse_ ? -1.0 : 1.0;

            thread_local size_t prev_len = 0;
            thread_local std::vector<Complex_Absorbing_T> a;
            thread_local std::vector<Complex> b, w;
            thread_local bool init = false;

            const size_t len_ex = to_pow2_up(2 * len_ - 1);
            const double N = static_cast<double>(len_);
            const double M = static_cast<double>(len_ex);

            if(prev_len != len_) init = false;
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
                prev_len = len_;
                init = true;
            }

            for (size_t i = 0; i < len_; ++i) {
                a[i] = seq_cplx_[i] * w[i].conj();
                b[i] = w[i];
            }

            for (size_t i = len_; i < len_ex; ++i) {
                a[i] = Complex_Absorbing_T{};
                b[i] = w[i];
            }

            fft_radix2__<false>(a, len_ex);
            fft_radix2__<false>(b, len_ex);

            for (size_t i = 0; i < len_ex; ++i) {
                a[i] *= b[i];
            }
            fft_radix2__<true>(a, len_ex);
            double normalizer = static_cast<double>(len_ex);

            for (size_t i = 0; i < len_; ++i) {
                a[i] /= normalizer;
            }

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

        for (size_t i = 0; i < Len_; ++i) {
            Seq_complx_[i] /= normalizer;
        }
    }


    //non-mutating FFT
    //normalize factor N^0.5 for FT and IFT
    template<class InIt, class OutIt>
    inline void fft_ortho(InIt First_, const size_t Len_, OutIt Dest_complx_)
    {
        using Complex_Absorbing_T = std::decay_t<decltype(*Dest_complx_)>;
        std::vector<Complex_Absorbing_T> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = (*First_++) * Complex::identity();
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<false>(seq_buf, Len_);
        else detail_::fft_bluestein__<false>(seq_buf, Len_);
        double normalizer = sqrt(Len_);

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
        using Complex_Absorbing_T = std::decay_t<decltype(*Dest_complx_)>;
        std::vector<Complex_Absorbing_T> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = (*First_++) * Complex::identity();
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<true>(seq_buf, Len_);
        else detail_::fft_bluestein__<true>(seq_buf, Len_);
        double normalizer = sqrt(Len_);

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
        using Complex_Absorbing_T = std::decay_t<decltype(*Dest_complx_)>;
        std::vector<Complex_Absorbing_T> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = (*First_++) * Complex::identity();
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
        using Complex_Absorbing_T = std::decay_t<decltype(*Dest_complx_)>;
        std::vector<Complex_Absorbing_T> seq_buf(Len_);
        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] = (*First_++) * Complex::identity();
        }

        if (is_pow2(Len_)) detail_::fft_radix2__<true>(seq_buf, Len_);
        else detail_::fft_bluestein__<true>(seq_buf, Len_);
        double normalizer = static_cast<double>(Len_);

        for (size_t i = 0; i < Len_; ++i) {
            seq_buf[i] /= normalizer;
        }

        for (size_t i = 0; i < Len_; ++i) {
            *Dest_complx_++ = seq_buf[i];
        }
    }

	
	//in-place FFT for 2d complex field in mat<Complex>
    //normalize factor N^0.5 for FT and IFT
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void fft2d_ortho_par(mat<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.ncols();
        const size_t Y = Field_complx_.nrows();
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
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void ifft2d_ortho_par(mat<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.ncols();
        const size_t Y = Field_complx_.nrows();
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
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void fft2d_par(mat<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.ncols();
        const size_t Y = Field_complx_.nrows();
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
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void ifft2d_par(mat<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.ncols();
        const size_t Y = Field_complx_.nrows();
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

    template<typename Complex_Absorbing_T, class Alloc__>
    inline void centralize(mat<Complex_Absorbing_T, Alloc__>& Spectr_) {

        const size_t X = Spectr_.ncols();
        const size_t Y = Spectr_.nrows();
        mat<Complex_Absorbing_T, Alloc__> xcentered(Y, X);
        for (size_t i = 0; i < Y; ++i) {
            auto src = Spectr_.row(i).ccycle_from(X / 2);
            auto des = xcentered.row(i).begin();
            for (size_t j = 0; j < X; ++j) {
                *des++ = *src++;
            }
        }

        for (size_t j = 0; j < X; ++j) {
            auto src = xcentered.col(j).ccycle_from(Y / 2);
            auto des = Spectr_.col(j).begin();
            for (size_t i = 0; i < Y; ++i) {
                *des++ = *src++;
            }
        }
    }

    //in-place FFT for 3d complex field in cube<Complex>
    //normalize factor N^0.5 for FT and IFT
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void fft3d_ortho_par(cube<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.depth();
        const size_t Y = Field_complx_.height();
        const size_t Z = Field_complx_.width();

        for (size_t i = 0; i < X; ++i) {
            fft2d_ortho_par(Field_complx_[i].get());
        }

        using args__ = struct { size_t i__; size_t j__; };
        mat<args__> arg_list = mat<args__>::creat(Y, Z, [](size_t i, size_t j) {return args__{ i, j }; });

        const auto depth_fft = [&](const args__& arg) -> void {
            auto dep_ary = Field_complx_.begin_at(arg.i__, arg.j__);
            fft_ortho(dep_ary, X);
            };

        std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), depth_fft);
    }

    //in-place FFT for 3d complex field in cube<Complex>
    //normalize factor N^0.5 for FT and IFT
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void ifft3d_ortho_par(cube<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.depth();
        const size_t Y = Field_complx_.height();
        const size_t Z = Field_complx_.width();

        for (size_t i = 0; i < X; ++i) {
            ifft2d_ortho_par(Field_complx_[i].get());
        }

        using args__ = struct { size_t i__; size_t j__; };
        mat<args__> arg_list = mat<args__>::creat(Y, Z, [](size_t i, size_t j) {return args__{ i, j }; });

        const auto depth_fft = [&](const args__& arg) -> void {
            auto dep_ary = Field_complx_.begin_at(arg.i__, arg.j__);
            ifft_ortho(dep_ary, X);
            };

        std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), depth_fft);
    }

    //in-place FFT for 3d complex field in cube<Complex>
    //normalize factor 1 for FT and N for IFT
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void fft3d_par(cube<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.depth();
        const size_t Y = Field_complx_.height();
        const size_t Z = Field_complx_.width();

        for (size_t i = 0; i < X; ++i) {
            fft2d_par(Field_complx_[i].get());
        }

        using args__ = struct { size_t i__; size_t j__; };
        mat<args__> arg_list = mat<args__>::creat(Y, Z, [](size_t i, size_t j) {return args__{ i, j }; });

        const auto depth_fft = [&](const args__& arg) -> void {
            auto dep_ary = Field_complx_.begin_at(arg.i__, arg.j__);
            fft(dep_ary, X);
            };

        std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), depth_fft);
    }

    //in-place FFT for 3d complex field in cube<Complex>
    //normalize factor 1 for FT and N for IFT
    template<typename Complex_Absorbing_T, class Alloc__>
    inline void ifft3d_par(cube<Complex_Absorbing_T, Alloc__>& Field_complx_)
    {
        const size_t X = Field_complx_.depth();
        const size_t Y = Field_complx_.height();
        const size_t Z = Field_complx_.width();

        for (size_t i = 0; i < X; ++i) {
            ifft2d_par(Field_complx_[i].get());
        }

        using args__ = struct { size_t i__; size_t j__; };
        mat<args__> arg_list = mat<args__>::creat(Y, Z, [](size_t i, size_t j) {return args__{ i, j }; });

        const auto depth_fft = [&](const args__& arg) -> void {
            auto dep_ary = Field_complx_.begin_at(arg.i__, arg.j__);
            ifft(dep_ary, X);
            };

        std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), depth_fft);
    }


    template<typename Complex_Absorbing_T, class Alloc__>
    inline void centralize(cube<Complex_Absorbing_T, Alloc__>& Spectr_) {

        const size_t X = Spectr_.depth();

        for (size_t i = 0; i < X; ++i) {
            centralize(Spectr_[i].get());
        }
        const size_t off = X / 2;
        std::vector<mat<Complex_Absorbing_T, Alloc__>> temp(X);
        for (size_t i = 0; i < X; ++i) {
            std::swap(Spectr_[i].get(), temp[i]);
        }
        for (size_t i = 0; i < X; ++i) {
            std::swap(Spectr_[i].get(), temp[(i + off) % X]);
        }
    }

}//namespace numer end


