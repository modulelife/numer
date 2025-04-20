//	numer_plot.h	-*- C_++20 -*-
#pragma once
//@brief: plot drawer header
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

#include <cstddef>
#include <cmath>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <numer_grid.h>




namespace numer {

	class Histogram {
	private:
		const size_t hght_;
		const size_t wdth_;
		double min_;
		double max_;
		RGBA bg_clr_{25, 25, 25, 0};
		RGBA line_clr_{ 225, 225, 225, 255 };
		mat<RGBA> rawimg_;

	public:

		Histogram(size_t Height_, size_t Width_, double Min_, double Max_)
			:hght_(Height_), wdth_(Width_), min_(Min_), max_(Max_)
		{
			rawimg_ = std::move(mat<RGBA>(hght_, wdth_, bg_clr_));
		}

		Histogram& clear() {
			rawimg_.set_all_to(bg_clr_);
			return *this;
		}

		Histogram& setBackColor(RGB Clr_) {
			bg_clr_ = attachAlpha(0, Clr_);
			const auto changeBgColor = [&](RGBA& pixel) -> void {
				if (pixel.A > 0) return;
				pixel = bg_clr_;
			};
			rawimg_.modify(changeBgColor);
			return *this;
		}

		Histogram& setLineColor(RGB Clr_) {
			line_clr_ = attachAlpha(line_clr_.A, Clr_);
			return *this;
		}

		Histogram& setLineTrans(uint8_t Alpha_) {
			if (Alpha_ == 0) Alpha_ = 1;
			line_clr_.A = Alpha_;
			return *this;
		}

		template<class AbstractVec>
		Histogram& drawData(AbstractVec&& Vec_real_, size_t N_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ - 1 });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ - 1 });

			size_t xaxis_ypos;
			if (!range.verifyAndIndex(0.0, xaxis_ypos)) {
				if (max_ <= 0.0) xaxis_ypos = hght_ - 1;
				if (min_ >= 0.0) xaxis_ypos = 0;
			}

			size_t val_ypos;
			for (size_t i = 0; i < wdth_; ++i) {
				double val = Vec_real_(idxer(i));
				if (isnan(val)) continue;
				if (!range.verifyAndIndex(val, val_ypos)) {
					val_ypos = val > max_ ? hght_ - 1 : 0;
				}

				size_t s, t;
				if (val_ypos > xaxis_ypos) {
					s = xaxis_ypos;
					t = val_ypos;
				}
				else if (val_ypos < xaxis_ypos) {
					s = val_ypos;
					t = xaxis_ypos;
				}
				else {
					s = xaxis_ypos;
					t = xaxis_ypos;
				}

				for (size_t j = s; j <= t; ++j) {
					rawimg_[hght_ - 1 - j][i] = mixAlpha(rawimg_[hght_ - 1 - j][i], line_clr_);
				}
			}

			return *this;
		}

		template<class AbstractVec, class EntrywiseConverter>
		Histogram& drawData(AbstractVec&& Vec_any_, size_t N_, EntrywiseConverter&& To_real_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ - 1 });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ - 1 });

			size_t xaxis_ypos;
			if (!range.verifyAndIndex(0.0, xaxis_ypos)) {
				if (max_ <= 0.0) xaxis_ypos = hght_ - 1;
				if (min_ >= 0.0) xaxis_ypos = 0;
			}

			size_t val_ypos;
			for (size_t i = 0; i < wdth_; ++i) {
				double val = To_real_(Vec_any_(idxer(i)));
				if (isnan(val)) continue;
				if (!range.verifyAndIndex(val, val_ypos)) {
					val_ypos = val > max_ ? hght_ - 1 : 0;
				}

				size_t s, t;
				if (val_ypos > xaxis_ypos) {
					s = xaxis_ypos;
					t = val_ypos;
				}
				else if(val_ypos < xaxis_ypos) {
					s = val_ypos;
					t = xaxis_ypos;
				}
				else {
					s = xaxis_ypos;
					t = xaxis_ypos;
				}

				for (size_t j = s; j <= t; ++j) {
					rawimg_[hght_ - 1 - j][i] = mixAlpha(rawimg_[hght_ - 1 - j][i], line_clr_);
				}
			}

			return *this;
		}

		template<class AbstractVec, class EntrywiseConverter, typename Ty>
		Histogram& drawData(AbstractVec&& Vec_any_, size_t N_, EntrywiseConverter&& To_real_, const std::function<RGB(Ty)>& Colorize_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ - 1 });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ - 1 });

			size_t xaxis_ypos;
			if (!range.verifyAndIndex(0.0, xaxis_ypos)) {
				if (max_ <= 0.0) xaxis_ypos = hght_ - 1;
				if (min_ >= 0.0) xaxis_ypos = 0;
			}

			size_t val_ypos;
			for (size_t i = 0; i < wdth_; ++i) {
				double val = To_real_(Vec_any_(idxer(i)));
				if (isnan(val)) continue;
				if (!range.verifyAndIndex(val, val_ypos)) {
					val_ypos = val > max_ ? hght_ - 1 : 0;
				}

				size_t s, t;
				if (val_ypos > xaxis_ypos) {
					s = xaxis_ypos;
					t = val_ypos;
				}
				else if (val_ypos < xaxis_ypos) {
					s = val_ypos;
					t = xaxis_ypos;
				}
				else {
					s = xaxis_ypos;
					t = xaxis_ypos;
				}
				
				RGBA line_clr = attachAlpha(line_clr_.A, Colorize_(Vec_any_(idxer(i))));
				for (size_t j = s; j <= t; ++j) {
					rawimg_[hght_ - 1 - j][i] = mixAlpha(rawimg_[hght_ - 1 - j][i], line_clr);
				}
			}

			return *this;
		}

		Histogram& drawHorizLine(double Val_) {
			RangeSampler range(RangeSpec{ min_, max_, hght_ - 1 });
			size_t ypos;
			if (range.verifyAndIndex(Val_, ypos)) {
				for (size_t i = 0; i < wdth_; ++i) {
					rawimg_[hght_ - 1 - ypos][i] = mixAlpha(rawimg_[hght_ - 1 - ypos][i], line_clr_);
				}
			}
			return *this;
		}

		mat<RGB> getImage() const {
			return mat<RGB>::creat(rawimg_, removeAlpha);
		}


	};


	class HeatMapPlot {
	private:
		const size_t hght_;
		const size_t wdth_;

	public:
		HeatMapPlot(size_t Height_, size_t Width_)
			:hght_(Height_), wdth_(Width_)
		{}

		template<class AbstractVec, typename Ty>
		mat<RGB> renderImage(AbstractVec&& Vec_any_, size_t N_, const std::function<RGB(Ty)> Colorize_) {
			
			RangeSampler x_idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ - 1 });
			mat<RGB> plot_img(hght_, wdth_);
			for (size_t i = 0; i < hght_; ++i) {
				for (size_t j = 0; j < wdth_; ++j) {
					plot_img[i][j] = Colorize_(Vec_any_(x_idxer(j)));
				}
			}
			return plot_img;
		}

		template<class AbstractMat, typename Ty>
		mat<RGB> renderImage(AbstractMat&& Mat_any_, size_t Nrow_, size_t Ncol_, const std::function<RGB(Ty)> Colorize_) {

			RangeSampler x_idxer(RangeSpec{ 0, static_cast<double>(Ncol_ - 1), wdth_ - 1 });
			RangeSampler y_idxer(RangeSpec{ 0, static_cast<double>(Nrow_ - 1), hght_ - 1 });
			mat<RGB> plot_img(hght_, wdth_);
			for (size_t i = 0; i < hght_; ++i) {
				for (size_t j = 0; j < wdth_; ++j) {
					plot_img[i][j] = Colorize_(Mat_any_(x_idxer(i), y_idxer(j)));
				}
			}
			return plot_img;
		}

	};


}//namespace numer end




