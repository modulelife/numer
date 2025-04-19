//	numer_qmkit.h	-*- C_++20 -*-
#pragma once
//@brief: quantum mechanic tool kit header
//
//@functions:
//	
//	
//@description:
//	
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <numer_common.h>
#include <numer_grid.h>
#include <numer_complex.h>
#include <numer_mat.h>
#include <numer_eigenfunc.h>
#include <numer_fourier.h>
#include <vector>
#include <algorithm>
#include <execution>




namespace numer{

	namespace qm {

		class HusimiQCalculator {
		private:
			unsigned max_n_;
			unsigned resol_;
			std::vector<double> fcoef_;

		public:
			HusimiQCalculator(unsigned Max_energy_level_, unsigned Resolution_)
				:max_n_(Max_energy_level_ + 1), resol_(Resolution_), fcoef_(max_n_)
			{
				for (unsigned i = 0; i < max_n_; ++i) {
					fcoef_[i] = 1.0 / sqrt(factorial(i));
				}
			}

			template<class AbstractVec>
			mat<double> calculate(AbstractVec&& E_rep_ket_, size_t N_) {
				mat<Complex> a_ket_(resol_, resol_, 0.0);
				const auto cplx_gen = [](double re, double im) { return Complex{ re, im }; };
				const double maxamp = sqrt(max_n_ + 0.5) + sqrt(max_n_) / 2.0;
				const RangeSpec spec{ -maxamp , maxamp , resol_ - 1 };
				BinaryFuncSampler cplx_sampler(cplx_gen, spec, spec);

				using args__ = struct { size_t i__; size_t j__; };
				mat<args__> arg_list = mat<args__>::creat(resol_, resol_, [](size_t i, size_t j) {return args__{ i, j }; });

				const auto calcu_proj = [&](const args__& arg) -> void {
					const Complex a = cplx_sampler(arg.i__, arg.j__).conj();
					Complex an = 1.0, proj = 0.0;
					double nfactor = exp(-a.sqrdAmp() / 2.0);
					for (unsigned k = 0; k < max_n_ && k < N_; ++k) {
						proj += an * E_rep_ket_(k) * fcoef_[k];
						an *= a;
					}
					proj *= nfactor;
					a_ket_[arg.i__][arg.j__] = proj;
				};

				std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), calcu_proj);


				return mat<double>::creat_par(a_ket_, [](const Complex& c) {
					return c.sqrdAmp() / Pi;
					});
			}
		};

	}//namespace qm end

}//namespace numer end