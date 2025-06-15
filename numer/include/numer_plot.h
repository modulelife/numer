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
#include <functional>




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
			RangeSampler range(RangeSpec{ min_, max_, hght_ });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ });

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

		template<class AbstractVec>
		Histogram& drawDataLine(AbstractVec&& Vec_real_, size_t N_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ });

			size_t preval_ypos, val_ypos;

			for (size_t i = 0; i < wdth_; ++i) {
				double val = Vec_real_(idxer(i));
				if (isnan(val)) continue;
				if (!range.verifyAndIndex(val, val_ypos)) {
					val_ypos = val > max_ ? hght_ - 1 : 0;
				}
				if (i == 0) preval_ypos = val_ypos;

				size_t s, t;
				if (val_ypos > preval_ypos) {
					s = preval_ypos;
					t = val_ypos;
				}
				else if (val_ypos < preval_ypos) {
					s = val_ypos;
					t = preval_ypos;
				}
				else {
					s = preval_ypos;
					t = preval_ypos;
				}

				for (size_t j = s; j <= t; ++j) {
					rawimg_[hght_ - 1 - j][i] = mixAlpha(rawimg_[hght_ - 1 - j][i], line_clr_);
				}
				preval_ypos = val_ypos;
			}

			return *this;
		}

		template<class AbstractVec, class EntrywiseConverter>
		Histogram& drawData(AbstractVec&& Vec_any_, size_t N_, EntrywiseConverter&& To_real_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ });

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
		Histogram& drawData(AbstractVec&& Vec_any_, size_t N_, EntrywiseConverter&& To_real_, const std::function<RGB(Ty&)> Colorize_) {
			if (N_ == 0) return *this;
			RangeSampler range(RangeSpec{ min_, max_, hght_ });
			RangeSampler idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ });

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
			RangeSampler range(RangeSpec{ min_, max_, hght_ });
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
		mat<RGB> renderImage(AbstractVec&& Vec_any_, size_t N_, const std::function<RGB(Ty&)> Colorize_) const {
			
			RangeSampler x_idxer(RangeSpec{ 0, static_cast<double>(N_ - 1), wdth_ });
			mat<RGB> plot_img(hght_, wdth_);
			for (size_t i = 0; i < hght_; ++i) {
				for (size_t j = 0; j < wdth_; ++j) {
					plot_img[i][j] = Colorize_(Vec_any_(x_idxer(j)));
				}
			}
			return plot_img;
		}

		template<class AbstractMat, typename Ty>
		mat<RGB> renderImage(AbstractMat&& Mat_any_, size_t Nrow_, size_t Ncol_, const std::function<RGB(Ty&)> Colorize_) const {

			RangeSampler x_idxer(RangeSpec{ 0, static_cast<double>(Ncol_ - 1), wdth_ });
			RangeSampler y_idxer(RangeSpec{ 0, static_cast<double>(Nrow_ - 1), hght_ });
			mat<RGB> plot_img(hght_, wdth_);
			for (size_t i = 0; i < hght_; ++i) {
				for (size_t j = 0; j < wdth_; ++j) {
					plot_img[i][j] = Colorize_(Mat_any_(x_idxer(i), y_idxer(j)));
				}
			}
			return plot_img;
		}

	};


	template <class CubeContainer>
	class FieldRelocator {
	private:
		const CubeContainer& base_;
		RangeSampler dim0_rg_, dim1_rg_, dim2_rg_;

	public:
		FieldRelocator(const CubeContainer& Cube_base_, size_t Dim0_size_, size_t Dim1_size_, size_t Dim2_size_)
			: base_(Cube_base_), dim0_rg_(), dim1_rg_(), dim2_rg_()
		{
			constexpr double L_MAX = 1.0;
			const size_t max_size = std::max({ Dim0_size_, Dim1_size_, Dim2_size_ });

			RangeSpec dim0_spec{ -L_MAX, L_MAX, Dim0_size_ }, dim1_spec{ -L_MAX, L_MAX, Dim1_size_ }, dim2_spec{ -L_MAX, L_MAX, Dim2_size_ };


			auto adjust_spec = [max_size](RangeSpec& spec, size_t dim_size) {
				if (max_size != dim_size) {
					double ratio = static_cast<double>(dim_size) / max_size;
					spec.start *= ratio;
					spec.end *= ratio;
				}
				};

			adjust_spec(dim0_spec, Dim0_size_);
			adjust_spec(dim1_spec, Dim1_size_);
			adjust_spec(dim2_spec, Dim2_size_);

			dim0_rg_.setRange(dim0_spec);
			dim1_rg_.setRange(dim1_spec);
			dim2_rg_.setRange(dim2_spec);

		}

		auto operator()(double X_, double Y_, double Z_) {
			using Elem_T = std::decay_t<decltype(base_[0][0][0])>;
			size_t idx0, idx1, idx2;
			if (!dim0_rg_.verifyAndIndex(X_, idx0)) return Elem_T{};
			if (!dim1_rg_.verifyAndIndex(Y_, idx1)) return Elem_T{};
			if (!dim2_rg_.verifyAndIndex(Z_, idx2)) return Elem_T{};
			return base_[idx0][idx1][idx2];
		}
	};


	template <class CubeContainer>
	class FieldSmoothRelocator {
	private:
		const CubeContainer& base_;
		RangeSampler dim0_rg_;
		RangeSampler dim1_rg_;
		RangeSampler dim2_rg_;
		size_t dim0_size_;
		size_t dim1_size_;
		size_t dim2_size_;

	public:
		FieldSmoothRelocator(const CubeContainer& Cube_base_, size_t Dim0_size_, size_t Dim1_size_, size_t Dim2_size_)
			: base_(Cube_base_),
			dim0_size_(Dim0_size_),
			dim1_size_(Dim1_size_),
			dim2_size_(Dim2_size_)
		{
			constexpr double L_MAX = 1.0;
			const size_t max_size = std::max({ Dim0_size_, Dim1_size_, Dim2_size_ });

			RangeSpec dim0_spec{ -L_MAX, L_MAX, Dim0_size_ };
			RangeSpec dim1_spec{ -L_MAX, L_MAX, Dim1_size_ };
			RangeSpec dim2_spec{ -L_MAX, L_MAX, Dim2_size_ };

			auto adjust_spec = [max_size](RangeSpec& spec, size_t dim_size) {
				if (max_size != dim_size) {
					double ratio = static_cast<double>(dim_size) / max_size;
					spec.start *= ratio;
					spec.end *= ratio;
				}
				};

			adjust_spec(dim0_spec, Dim0_size_);
			adjust_spec(dim1_spec, Dim1_size_);
			adjust_spec(dim2_spec, Dim2_size_);

			dim0_rg_.setRange(dim0_spec);
			dim1_rg_.setRange(dim1_spec);
			dim2_rg_.setRange(dim2_spec);
		}

		auto operator()(double X_, double Y_, double Z_) {
			using Elem_T = std::decay_t<decltype(base_[0][0][0])>;

			auto get_continuous_index = [](double coord, const RangeSampler& rg, size_t size)
				-> std::pair<double, bool>
				{
					double normalized = (coord - rg.getStart()) / rg.getDiff();
					if (normalized < 0.0 || normalized > 1.0)
						return { 0.0, false };

					return { normalized * (size - 1), true };
				};

			auto [idx0, valid0] = get_continuous_index(X_, dim0_rg_, dim0_size_);
			auto [idx1, valid1] = get_continuous_index(Y_, dim1_rg_, dim1_size_);
			auto [idx2, valid2] = get_continuous_index(Z_, dim2_rg_, dim2_size_);

			if (!valid0 || !valid1 || !valid2)
				return Elem_T{};

			auto get_parts = [](double index) -> std::pair<size_t, double> {
				size_t i = static_cast<size_t>(std::floor(index));
				double f = index - i;
				return { i, f };
				};

			auto [i0, fx] = get_parts(idx0);
			auto [i1, fy] = get_parts(idx1);
			auto [i2, fz] = get_parts(idx2);

			i0 = (i0 >= dim0_size_ - 1) ? dim0_size_ - 2 : i0;
			i1 = (i1 >= dim1_size_ - 1) ? dim1_size_ - 2 : i1;
			i2 = (i2 >= dim2_size_ - 1) ? dim2_size_ - 2 : i2;

			Elem_T v000 = base_[i0][i1][i2];
			Elem_T v001 = base_[i0][i1][i2 + 1];
			Elem_T v010 = base_[i0][i1 + 1][i2];
			Elem_T v011 = base_[i0][i1 + 1][i2 + 1];
			Elem_T v100 = base_[i0 + 1][i1][i2];
			Elem_T v101 = base_[i0 + 1][i1][i2 + 1];
			Elem_T v110 = base_[i0 + 1][i1 + 1][i2];
			Elem_T v111 = base_[i0 + 1][i1 + 1][i2 + 1];

			Elem_T v00 = v000 * (1 - fz) + v001 * fz;
			Elem_T v01 = v010 * (1 - fz) + v011 * fz;
			Elem_T v10 = v100 * (1 - fz) + v101 * fz;
			Elem_T v11 = v110 * (1 - fz) + v111 * fz;

			Elem_T v0 = v00 * (1 - fy) + v01 * fy;
			Elem_T v1 = v10 * (1 - fy) + v11 * fy;

			return v0 * (1 - fx) + v1 * fx;
		}
	};


	namespace gridclr {

		namespace brt {

			class linear {
			public:
				double operator()(double x) const {
					return x;
				}
			};

			class gradual {
			public:
				double operator()(double x) const {
					return sqrt(x);
				}
			};

			class steep {
			public:
				double operator()(double x) const {
					return x * x;
				}
			};

			class soft {
			public:
				double operator()(double x) const {
					if (x < 0.0) return 0.0;
					if (x < 1.0)return sqrt(x);
					return 1.0;
				}
			};

			class hard {
			public:
				double operator()(double x) const {
					if (x < 0.0) return 0.0;
					if (x < 1.0)return pow(x, 7.0);
					return 1.0;
				}
			};

			class contour {
			public:
				double operator()(double x) const {
					if (x < 0.0) return 0.0;
					if (x < 1.0)return 1.0 * pow(x, 5.0);
					if (x < 2.0)return 1.0 * pow(2.0 - x, 5.0);
					return 0.0;
				}
			};


		}//namespace brt end

		template<class BrightnessStrategy>
		class HeatMap {
		private:
			BrightnessStrategy	brightness_{};
			GrayScale			to_double_{ 0.0, 1.0 };
			double				min_thrs_;
			double				max_thrs_;
			LinearHeatMap		clrmap_;

		public:
			HeatMap()
				: min_thrs_(0.0), max_thrs_(1.0), clrmap_(0.0, 1.0, Color::Glacier()) {
			}
			HeatMap(double Minimum_, double Maximum_)
				: min_thrs_(Minimum_), max_thrs_(Maximum_), clrmap_(Minimum_, Maximum_, Color::Glacier()) {
			}
			HeatMap(double Minimum_, double Maximum_, NormalizedColorMap Color_map_)
				: min_thrs_(Minimum_), max_thrs_(Maximum_), clrmap_(Minimum_, Maximum_, Color_map_) {
			}
			explicit HeatMap(double Minimum_, double Maximum_, LinearHeatMap Heatmap_)
				: min_thrs_(Minimum_), max_thrs_(Maximum_), clrmap_(Heatmap_) {
			}

			Vec3<double> operator()(double x_) const {

				if (x_ < min_thrs_) return Vec3<double>{0.0, 0.0, 0.0};

				double y = (x_ - min_thrs_) / (max_thrs_ - min_thrs_);
				RGB clr = clrmap_(x_);
				Vec3<double> grid_clr{ to_double_(clr.R), to_double_(clr.G), to_double_(clr.B) };
				const double bri = brightness_(y);
				grid_clr *= bri;
				return grid_clr;
			}

			Vec3<double> operator()(Complex c_) const {
				double x_ = c_.sqrdAmp();
				if (x_ < min_thrs_) return Vec3<double>{0.0, 0.0, 0.0};

				double y = (x_ - min_thrs_) / (max_thrs_ - min_thrs_);
				RGB clr = clrmap_(x_);
				Vec3<double> grid_clr{ to_double_(clr.R), to_double_(clr.G), to_double_(clr.B) };
				const double bri = brightness_(y);
				grid_clr *= bri;
				return grid_clr;
			}

		};


		template<class BrightnessStrategy>
		class CRainbow{
		private:
			BrightnessStrategy	brightness_{};
			GrayScale			to_double_{ 0.0, 1.0 };
			double				max_thrs_;
			ComplxRainbowClr	crainbow_;

		public:
			CRainbow()
				: max_thrs_(1.0), crainbow_(1.0) {
			}
			CRainbow(double Maximum_)
				: max_thrs_(Maximum_), crainbow_(sqrt(Maximum_)) {
			}
			explicit CRainbow(double Maximum_, ComplxRainbowClr Crainbow_)
				: max_thrs_(Maximum_), crainbow_(Crainbow_) {
			}


			Vec3<double> operator()(Complex c_) const {
				double x_ = c_.sqrdAmp();

				double y = x_ / max_thrs_;
				RGB clr = crainbow_(c_);
				Vec3<double> grid_clr{ to_double_(clr.R), to_double_(clr.G), to_double_(clr.B) };
				const double bri = brightness_(y);
				grid_clr *= bri;
				return grid_clr;
			}

		};

		template<class BrightnessStrategy>
		class CPhaseClr {
		private:
			BrightnessStrategy	brightness_{};
			GrayScale			to_double_{ 0.0, 1.0 };
			double				max_thrs_;
			ComplxPhaseClr		cphsclr_;

		public:
			CPhaseClr()
				: max_thrs_(1.0), cphsclr_() {
			}
			CPhaseClr(double Maximum_)
				: max_thrs_(Maximum_), cphsclr_() {
			}
			explicit CPhaseClr(double Maximum_, ComplxPhaseClr Cphaseclr_)
				: max_thrs_(Maximum_), cphsclr_(Cphaseclr_) {
			}


			Vec3<double> operator()(Complex c_) const {
				double x_ = c_.sqrdAmp();

				double y = x_ / max_thrs_;
				RGB clr = cphsclr_(c_);
				Vec3<double> grid_clr{ to_double_(clr.R), to_double_(clr.G), to_double_(clr.B) };
				const double bri = brightness_(y);
				grid_clr *= bri;
				return grid_clr;
			}

		};

		class GammaTransp {
		private:
			double min_thrs_;
			double max_thrs_;
			double gm_{ 1.0 };
		public:
			GammaTransp() : min_thrs_(-10.0), max_thrs_(10.0) {}
			GammaTransp(double Minimum_, double Maximum_)
				: min_thrs_(Minimum_), max_thrs_(Maximum_) {
			}
			GammaTransp(double Minimum_, double Maximum_, double Gamma_)
				: min_thrs_(Minimum_), max_thrs_(Maximum_), gm_(Gamma_){
			}

			double operator()(double x_) const {
				if (x_ < min_thrs_) return 0.0;
				if (x_ < max_thrs_)
				{
					double y = (x_ - min_thrs_) / (max_thrs_ - min_thrs_);
					return pow(y, gm_);
				}
				return 1.0;
			}

			double operator()(Complex c_) const {
				double x = c_.sqrdAmp();
				if (x < min_thrs_) return 0.0;
				if (x < max_thrs_)
				{
					double y = (x - min_thrs_) / (max_thrs_ - min_thrs_);
					return pow(y, gm_);
				}
				return 1.0;
			}
		};


	}//namespace gridclr end

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

		template<class RelocatedField, typename Ty>
		mat<RGB> renderImage(
			RelocatedField&& Field_,
			const std::function<Vec3<double>(Ty&)> Grid_colorizer_) const
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
					rgb += Grid_colorizer_(grid_val);
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


	class LuminantPlot3D {
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
		LuminantPlot3D(size_t Height_, size_t Width_)
			:hght_(Height_), wdth_(Width_) {
		}

		LuminantPlot3D& setCamDistance(double Dist_) {
			distance_ = Dist_;
			return *this;
		}

		//between 0 to 2Pi
		LuminantPlot3D& setCamAzimuthAngle(double Azim_) {
			azimuth_ = Azim_;
			return *this;
		}

		//between -Pi/2 to Pi/2
		LuminantPlot3D& setCamElevationAngle(double Elev_) {
			elevation_ = Elev_;
			return *this;
		}

		//between 0 to 0.5Pi recommended, must less than Pi
		LuminantPlot3D& setHorizontalFOV(double HFOV_) {
			hori_fov_ = HFOV_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		LuminantPlot3D& setRenderDepth(double Rdp_) {
			render_dp_ = Rdp_;
			return *this;
		}

		//number of samples along a single ray, aka image quality
		LuminantPlot3D& setFineness(size_t Fine_) {
			fineness_ = Fine_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		LuminantPlot3D& setBrightnessGain(double Gain_) {
			bright_gain_ = Gain_;
			return *this;
		}

		template<class RelocatedField, typename Ty>
		mat<RGB> renderImage(
			RelocatedField&& Field_,
			const std::function<Vec3<double>(Ty&)> Grid_colorizer_) const
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

			RangeSampler len_para(RangeSpec{ (render_dp_ >= distance_ ? (0.01 * distance_) : distance_ - render_dp_), distance_ + render_dp_, fineness_ });

			const double decay_normalizer = distance_ * distance_;
			const auto decay = [=](double d) {
				return decay_normalizer / (d * d);
				};

			//render
			mat<RGB> image(hght_, wdth_);

			using args__ = struct { size_t i__; size_t j__; };
			mat<args__> arg_list = mat<args__>::creat(hght_, wdth_, [](size_t i, size_t j) {return args__{ i, j }; });

			const auto frag_shader = [&](const args__& arg) -> void {
				vec3 et{ v_slope(arg.i__), h_slope(arg.j__), 1.0 };
				et /= sqrt(et.t() * et);

				vec3 rgb{ 0.0, 0.0, 0.0 };

				for (size_t i = 0; i < fineness_; ++i) {
					const double dis = len_para(i);
					vec3 cam_coord_pos = et * dis;
					vec3 world_coord_pos = proj * cam_coord_pos + cam_pos;

					Ty grid_val = Field_(world_coord_pos[0], world_coord_pos[1], world_coord_pos[2]);
					rgb += Grid_colorizer_(grid_val) * decay(dis);
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


	class TranslucentPlot3D {
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
		TranslucentPlot3D(size_t Height_, size_t Width_)
			:hght_(Height_), wdth_(Width_) {
		}

		TranslucentPlot3D& setCamDistance(double Dist_) {
			distance_ = Dist_;
			return *this;
		}

		//between 0 to 2Pi
		TranslucentPlot3D& setCamAzimuthAngle(double Azim_) {
			azimuth_ = Azim_;
			return *this;
		}

		//between -Pi/2 to Pi/2
		TranslucentPlot3D& setCamElevationAngle(double Elev_) {
			elevation_ = Elev_;
			return *this;
		}

		//between 0 to 0.5Pi recommended, must less than Pi
		TranslucentPlot3D& setHorizontalFOV(double HFOV_) {
			hori_fov_ = HFOV_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		TranslucentPlot3D& setRenderDepth(double Rdp_) {
			render_dp_ = Rdp_;
			return *this;
		}

		//number of samples along a single ray, aka image quality
		TranslucentPlot3D& setFineness(size_t Fine_) {
			fineness_ = Fine_;
			return *this;
		}

		//between 1.0 to any, don't set it too large though
		TranslucentPlot3D& setBrightnessGain(double Gain_) {
			bright_gain_ = Gain_;
			return *this;
		}

		template<class RelocatedField, typename Ty>
		mat<RGB> renderImage(
			RelocatedField&& Field_,
			const std::function<Vec3<double>(Ty&)> Grid_colorizer_,
			const std::function<double(Ty&)> Transparency_) const
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

			RangeSampler len_para(RangeSpec{ distance_ + render_dp_, (render_dp_ >= distance_ ? 0.0 : distance_ - render_dp_), fineness_ });

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
					vec3 gird_rgb = Grid_colorizer_(grid_val);
					double transp = Transparency_(grid_val);
					rgb = rgb * (1.0 - transp) + gird_rgb * transp;
				}
				rgb *= bright_gain_;

				GrayScale to_uint8t(0.0, 1.0);
				RGB pixel{ to_uint8t(rgb[0]), to_uint8t(rgb[1]), to_uint8t(rgb[2]) };
				image[arg.i__][arg.j__] = pixel;
				};

			std::for_each(std::execution::par, arg_list.begin(), arg_list.end(), frag_shader);

			return image;
		}

	};


}//namespace numer end




