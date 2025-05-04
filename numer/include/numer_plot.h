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
#include <numer_matrix.h>




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
			:hght_(Height_), wdth_(Width_) {
		}

		template<class AbstractVec, typename Ty>
		mat<RGB> renderImage(AbstractVec&& Vec_any_, size_t N_, const std::function<RGB(Ty)> Colorize_) const {
			
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
		mat<RGB> renderImage(AbstractMat&& Mat_any_, size_t Nrow_, size_t Ncol_, const std::function<RGB(Ty)> Colorize_) const {

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


	class DensityPlot3D {
	private:
		const size_t hght_;
		const size_t wdth_;
		double distance_{ 2.0 };
		double azimuth_{ Pi / 4.0 };
		double elevation_{ Pi / 4.0 };
		double hori_fov_{ Pi / 3.0 };
		double render_dp_{ 1.0 };
		size_t fineness_{ 100 };
		double bright_gain_{ 1.0 };

	public:
		DensityPlot3D(size_t Height_, size_t Width_)
			:hght_(Height_), wdth_(Width_) {
		}

		DensityPlot3D& setCamDistance(double Dist_) {
			distance_ = Dist_;
			return *this;
		}

		//between 0 to 2Pi
		DensityPlot3D& setCamAzimuthAngle(double Azim_) {
			azimuth_ = Azim_;
			return *this;
		}

		//between -Pi/2 to Pi/2
		DensityPlot3D& setCamElevationAngle(double Elev_) {
			elevation_ = Elev_;
			return *this;
		}

		//between 0 to 0.5Pi recommended, must less than Pi
		DensityPlot3D& setHorizontalFOV(double HFOV_) {
			hori_fov_ = HFOV_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		DensityPlot3D& setRenderDepth(double Rdp_) {
			render_dp_ = Rdp_;
			return *this;
		}

		//number of samples along a single ray, aka image quality
		DensityPlot3D& setFineness(size_t Fine_) {
			fineness_ = Fine_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		DensityPlot3D& setBrightnessGain(double Gain_) {
			bright_gain_ = Gain_;
			return *this;
		}

		template<class FieldRelocator, typename Ty>
		mat<RGB> renderImage(
			FieldRelocator&& Field_,
			const std::function<Vec3<double>(Ty)> Grid_colorize_) const
		{
			//convenient defs
			using vec3 = Vec3<double>;
			using vec3t = Vec3t<double>;

			//calculate the basis transform matrix
			constexpr vec3 ex{ 1.0, 0.0, 0.0 }, ey{ 0.0, 1.0, 0.0 }, ez{ 0.0, 0.0, 1.0 };
			constexpr Vec3<vec3t> world_basis{ ex.t(), ey.t(), ez.t() };
			
			vec3 cam_pos{ 
				Spheric_To_Cartes3::x(distance_, Pi / 2.0 - elevation_, azimuth_),
				Spheric_To_Cartes3::y(distance_, Pi / 2.0 - elevation_, azimuth_),
				Spheric_To_Cartes3::z(distance_, Pi / 2.0 - elevation_, azimuth_),
			};

			vec3 ed = -cam_pos / distance_;
			vec3 eh = ez ^ cam_pos;
			eh /= sqrt(eh.t() * eh);
			vec3 ev = eh ^ ed;
			Vec3t<vec3> cam_basis{ ev, eh, ed };

			const auto proj = world_basis * cam_basis;

			//calculate parameters for ray generation
			const double vh_ratio = static_cast<double>(hght_) / static_cast<double>(wdth_);
			const double tan_half_hfov = tan(hori_fov_ / 2.0);
			const double tan_half_vfov = tan_half_hfov * vh_ratio;

			RangeSampler v_slope(RangeSpec{ tan_half_vfov, -tan_half_vfov, hght_ });
			RangeSampler h_slope(RangeSpec{ -tan_half_hfov, tan_half_hfov, wdth_ });

			RangeSampler len_para(RangeSpec{ (render_dp_ >= distance_ ? 0.0 : distance_ - render_dp_), distance_ + render_dp_, fineness_ });

			//render
			mat<RGB> image(hght_, wdth_);

			using args__ = struct { size_t i__; size_t j__; };
			mat<args__> arg_list = mat<args__>::creat(hght_, wdth_, [](size_t i, size_t j) {return args__{ i, j }; });

			const auto frag_shader = [&](const args__& arg) -> void {
				vec3 et{ v_slope(arg.i__), h_slope(arg.j__), 1.0 };
				et /= sqrt(et.t() * et);

				vec3 rgb{ 0.0, 0.0, 0.0 };

				for (size_t i = 0; i < fineness_; ++i) {
					vec3 cam_coord_pos = et * len_para(i);
					vec3 world_coord_pos = proj * cam_coord_pos + cam_pos;

					Ty grid_val = Field_(world_coord_pos[0], world_coord_pos[1], world_coord_pos[2]);
					rgb += Grid_colorize_(grid_val);
				}
				rgb *= len_para.step() * bright_gain_;

				GrayScale to_uint8t(0.0, 1.0);
				RGB pixel{ to_uint8t(rgb[0]), to_uint8t(rgb[1]), to_uint8t(rgb[2]) };
				image[arg.i__][arg.j__] = pixel;
			};

			std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), frag_shader);

			return image;
		}

	};


}//namespace numer end




