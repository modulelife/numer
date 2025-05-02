//	numer_matrix.h		-*- C++20 -*-
#pragma once
//@brief: basic linear algebra support for small-size vectors and matrices
//
//@classes:
//	
// 
//
//@description:
//	
//
//@usage:
//
//----------------------------content begins----------------------------

#include <array>
#include <utility>
#include <numer_complex.h>



namespace numer {


    template<class Derived, typename Ty, unsigned Dim>
    class Vec_base {
    protected:
        std::array<Ty, Dim> components_{};

    public:
        using std_array = std::array<Ty, Dim>;

        constexpr Vec_base() : components_{} {}
        constexpr explicit Vec_base(const std_array& Array_)
            : components_(Array_) {
        }
        template<typename... Args>
        constexpr explicit Vec_base(Args... Components_)
            : components_{ static_cast<Ty>(Components_)... } {
        }

        constexpr Ty& operator[](unsigned Idx_) { return components_[Idx_]; }
        constexpr const Ty& operator[](unsigned Idx_) const { return components_[Idx_]; }

        
        constexpr Derived operator-() const {
            Derived result(components_);
            for (unsigned i = 0; i < Dim; ++i) {
                result[i] = -result[i];
            }
            return result;
        }

        constexpr Derived& operator+=(const Derived& Right_) {
            for (unsigned i = 0; i < Dim; ++i) {
                components_[i] += Right_[i];
            }
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& operator-=(const Derived& Right_) {
            for (unsigned i = 0; i < Dim; ++i) {
                components_[i] -= Right_[i];
            }
            return static_cast<Derived&>(*this);
        }

        template<typename Tx>
        constexpr Derived& operator*=(const Tx& Right_) {
            for (unsigned i = 0; i < Dim; ++i) {
                components_[i] *= Right_;
            }
            return static_cast<Derived&>(*this);
        }

        template<typename Tx>
        constexpr Derived& operator/=(const Tx& Right_) {
            for (unsigned i = 0; i < Dim; ++i) {
                components_[i] /= Right_;
            }
            return static_cast<Derived&>(*this);
        }
    };

    
    template<typename Ty, unsigned Dim, bool IsColumn>
    class Vec_impl : public Vec_base<Vec_impl<Ty, Dim, IsColumn>, Ty, Dim> {
    public:
        using Base = Vec_base<Vec_impl, Ty, Dim>;
        using Base::Base;

        constexpr auto transpose() const {
            return Vec_impl<Ty, Dim, !IsColumn>(Base::components_);
        }
    };

    
    template<unsigned Dim, bool IsColumn>
    class Vec_impl<Complex, Dim, IsColumn>
        : public Vec_base<Vec_impl<Complex, Dim, IsColumn>, Complex, Dim> {
    public:
        using Base = Vec_base<Vec_impl, Complex, Dim>;
        using Base::Base;

        constexpr auto conj() const {
            Vec_impl result(Base::components_);
            for (unsigned i = 0; i < Dim; ++i) {
                result[i] = result[i].conj();
            };
            return result;
        }

        constexpr auto transpose() const {
            return Vec_impl<Complex, Dim, !IsColumn>(Base::components_);
        }

        constexpr auto hermiConj() const {
            return this->conj().transpose();
        }
    };

    
    template<typename Ty, unsigned Dim>
    using Vec = Vec_impl<Ty, Dim, true>;// column vector

    template<typename Ty, unsigned Dim>
    using Vect = Vec_impl<Ty, Dim, false>;// row vector

    template<typename Ty>
    using Vec3 = Vec_impl<Ty, 3, true>;// 3d vector




    template<class Derived, typename Ty, unsigned Rows, unsigned Cols>
    class Matrix_base {
    protected:
        std::array<std::array<Ty, Cols>, Rows> entries_{};

    public:
        using std_array2 = std::array<std::array<Ty, Cols>, Rows>;
        using row = std::array<Ty, Cols>;

        constexpr Matrix_base() : entries_{} {}
        constexpr explicit Matrix_base(const std_array2& Array2_)
            : entries_(Array2_) {
        }
        template<typename... RowArgs>
        constexpr Matrix_base(const RowArgs&... Rows_)
            : entries_{ Rows_... } {
        }

        constexpr row& operator[](unsigned Idx_) { return entries_[Idx_]; }
        constexpr const row& operator[](unsigned Idx_) const { return entries_[Idx_]; }

        constexpr Derived operator-() const {
            Derived result(entries_);
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    result[i][j] = -result[i][j];
                }
            }
            return result;
        }

        constexpr Derived& operator+=(const Derived& Right_) {
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    entries_[i][j] += Right_[i][j];
                }
            }
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& operator-=(const Derived& Right_) {
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    entries_[i][j] -= Right_[i][j];
                }
            }
            return static_cast<Derived&>(*this);
        }

        template<typename Tx>
        constexpr Derived& operator*=(const Tx& Right_) {
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    entries_[i][j] *= Right_;
                }
            }
            return static_cast<Derived&>(*this);
        }

        template<typename Tx>
        constexpr Derived& operator/=(const Tx& Right_) {
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    entries_[i][j] /= Right_;
                }
            }
            return static_cast<Derived&>(*this);
        }

    };


    template<typename Ty, unsigned Rows, unsigned Cols>
    class Matrix
        : public Matrix_base<Matrix<Ty, Rows, Cols>, Ty, Rows, Cols> {
    public:
        using Base = Matrix_base<Matrix<Ty, Rows, Cols>, Ty, Rows, Cols>;
        using Base::Base;

        constexpr auto transpose() const {
            Matrix<Ty, Cols, Rows> result{};
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    result[j][i] = Base::entries_[i][j];
                }
            }
            return result;
        }
    };


    template<unsigned Rows, unsigned Cols>
    class Matrix<Complex, Rows, Cols>
        : public Matrix_base<Matrix<Complex, Rows, Cols>, Complex, Rows, Cols> {
    public:
        using Base = Matrix_base<Matrix<Complex, Rows, Cols>, Complex, Rows, Cols>;
        using Base::Base;

        constexpr auto conj() const {
            Matrix<Complex, Rows, Cols> result(Base::entries_);
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    result[i][j] = result[i][j].conj();
                }
            }
            return result;
        }

        constexpr auto transpose() const {
            Matrix<Complex, Cols, Rows> result{};
            for (unsigned i = 0; i < Rows; ++i) {
                for (unsigned j = 0; j < Cols; ++j) {
                    result[j][i] = Base::entries_[i][j];
                }
            }
            return result;
        }

        constexpr auto hermiConj() const {
            return this->conj().transpose();
        }
    };


    template<typename Ty, unsigned Dim>
    using Smatrix = Matrix<Ty, Dim, Dim>; // square matrix



    template<typename Ty, typename Tx, unsigned Dim, bool IsColumn>
    inline constexpr auto operator+(const Vec_impl<Ty, Dim, IsColumn>& Left_, const Vec_impl<Tx, Dim, IsColumn>& Right_);

    template<typename Ty, typename Tx, unsigned Dim, bool IsColumn>
    inline constexpr auto operator-(const Vec_impl<Ty, Dim, IsColumn>& Left_, const Vec_impl<Tx, Dim, IsColumn>& Right_);

    template<typename Ty, typename Tx, unsigned Dim>
    inline constexpr auto operator*(const Vect<Ty, Dim>& Left_, const Vec<Tx, Dim>& Right_);

    template<typename Ty, typename Tx, unsigned Dim>
    inline constexpr auto operator*(const Vec<Ty, Dim>& Left_, const Vect<Tx, Dim>& Right_);

    template<typename Tx, typename Ty, unsigned Dim, bool IsColumn>
    inline constexpr auto operator*(const Tx& Left_, const Vec_impl<Ty, Dim, IsColumn>& Right_);

    template<typename Ty, unsigned Dim, bool IsColumn, typename Tx>
    inline constexpr auto operator*(const Vec_impl<Ty, Dim, IsColumn>& Left_, const Tx& Right_);

    template<typename Ty, unsigned Dim, bool IsColumn, typename Tx>
    inline constexpr auto operator/(const Vec_impl<Ty, Dim, IsColumn>& Left_, const Tx& Right_);

    template<typename Ty, typename Tx>
    inline constexpr auto operator^(const Vec3<Ty>& Left_, const Vec3<Tx>& Right_);//vector product for 3d vectors

    template<typename Tx, typename Ty, unsigned Rows, unsigned Cols>
    inline constexpr auto operator*(const Vect<Tx, Rows>& Left_, const Matrix<Ty, Rows, Cols>& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
    inline constexpr auto operator*(const Matrix<Ty, Rows, Cols>& Left_, const Vec<Tx, Cols>& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
    inline constexpr auto operator+(const Matrix<Ty, Rows, Cols>& Left_, const Matrix<Tx, Rows, Cols>& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
    inline constexpr auto operator-(const Matrix<Ty, Rows, Cols>& Left_, const Matrix<Tx, Rows, Cols>& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Mid, unsigned Cols>
    inline constexpr auto operator*(const Matrix<Ty, Rows, Mid>& Left_, const Matrix<Tx, Mid, Cols>& Right_);

    template<typename Tx, typename Ty, unsigned Rows, unsigned Cols>
    inline constexpr auto operator*(const Tx& Left_, const Matrix<Ty, Rows, Cols>& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
    inline constexpr auto operator*(const Matrix<Ty, Rows, Cols>& Left_, const Tx& Right_);

    template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
    inline constexpr auto operator/(const Matrix<Ty, Rows, Cols>& Left_, const Tx& Right_);

}//namespace numer end


template<typename Ty, typename Tx, unsigned Dim, bool IsColumn>
inline constexpr auto numer::operator+(
    const Vec_impl<Ty, Dim, IsColumn>& Left_,
    const Vec_impl<Tx, Dim, IsColumn>& Right_)
{
    using Tz = decltype(std::declval<Ty>() + std::declval<Tx>());
    Vec_impl<Tz, Dim, IsColumn> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        result[i] = Left_[i] + Right_[i];
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Dim, bool IsColumn>
inline constexpr auto numer::operator-(
    const Vec_impl<Ty, Dim, IsColumn>& Left_,
    const Vec_impl<Tx, Dim, IsColumn>& Right_)
{
    using Tz = decltype(std::declval<Ty>() - std::declval<Tx>());
    Vec_impl<Tz, Dim, IsColumn> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        result[i] = Left_[i] - Right_[i];
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Dim>
inline constexpr auto numer::operator*(
    const Vect<Ty, Dim>& Left_,
    const Vec<Tx, Dim>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Tz result = static_cast<Tz>(0);
    for (unsigned i = 0; i < Dim; ++i) {
        result += Left_[i] * Right_[i];
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Dim>
inline constexpr auto numer::operator*(
    const Vec<Ty, Dim>& Left_,
    const Vect<Tx, Dim>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Smatrix<Tz, Dim> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        for (unsigned j = 0; j < Dim; ++j) {
            result[i][j] = Left_[i] * Right_[j];
        }
    }
    return result;
}

template<typename Tx, typename Ty, unsigned Dim, bool IsColumn>
inline constexpr auto numer::operator*(
    const Tx& Left_,
    const Vec_impl<Ty, Dim, IsColumn>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Vec_impl<Tz, Dim, IsColumn> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        result[i] = Left_ * Right_[i];
    }
    return result;
}

template<typename Ty, unsigned Dim, bool IsColumn, typename Tx>
inline constexpr auto numer::operator*(
    const Vec_impl<Ty, Dim, IsColumn>& Left_,
    const Tx& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Vec_impl<Tz, Dim, IsColumn> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        result[i] = Left_[i] * Right_;
    }
    return result;
}

template<typename Ty, unsigned Dim, bool IsColumn, typename Tx>
inline constexpr auto numer::operator/(
    const Vec_impl<Ty, Dim, IsColumn>& Left_,
    const Tx& Right_)
{
    using Tz = decltype(std::declval<Ty>() / std::declval<Tx>());
    Vec_impl<Tz, Dim, IsColumn> result{};
    for (unsigned i = 0; i < Dim; ++i) {
        result[i] = Left_[i] / Right_;
    }
    return result;
}

template<typename Ty, typename Tx>
inline constexpr auto numer::operator^(
    const Vec3<Ty>& Left_,
    const Vec3<Tx>& Right_)
{
    using Tz = decltype(std::declval<Ty>()* std::declval<Tx>());
    Vec3<Tz> result{};
    result[0] = Left_[1] * Right_[2] - Left_[2] * Right_[1];
    result[1] = Left_[2] * Right_[0] - Left_[0] * Right_[2];
    result[2] = Left_[0] * Right_[1] - Left_[1] * Right_[0];
    return result;
}

template<typename Tx, typename Ty, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator*(
    const Vect<Tx, Rows>& Left_,
    const Matrix<Ty, Rows, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Vect<Tz, Cols> result{};
    for (unsigned j = 0; j < Cols; ++j) {
        Tz component = static_cast<Tz>(0);
        for (unsigned i = 0; i < Rows; ++i) {
            component += Left_[i] * Right_[i][j];
        }
        result[j] = component;
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator*(
    const Matrix<Ty, Rows, Cols>& Left_,
    const Vec<Tx, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>()* std::declval<Tx>());
    Vec<Tz, Rows> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        Tz component = static_cast<Tz>(0);
        for (unsigned j = 0; j < Cols; ++j) {
            component += Left_[i][j] * Right_[j];
        }
        result[i] = component;
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator+(
    const Matrix<Ty, Rows, Cols>& Left_,
    const Matrix<Tx, Rows, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>() + std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result[i][j] = Left_[i][j] + Right_[i][j];
        }
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator-(
    const Matrix<Ty, Rows, Cols>& Left_,
    const Matrix<Tx, Rows, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>() - std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result[i][j] = Left_[i][j] - Right_[i][j];
        }
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Mid, unsigned Cols>
inline constexpr auto numer::operator*(
    const Matrix<Ty, Rows, Mid>& Left_,
    const Matrix<Tx, Mid, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            Tz entry = static_cast<Tz>(0);
            for (unsigned k = 0; k < Mid; ++k) {
                entry += Left_[i][k] * Right_[k][j];
            }
            result[i][j] = entry;
        }
    }
    return result;
}

template<typename Tx, typename Ty, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator*(
    const Tx& Left_,
    const Matrix<Ty, Rows, Cols>& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result[i][j] = Left_ * Right_[i][j];
        }
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator*(
    const Matrix<Ty, Rows, Cols>& Left_,
    const Tx& Right_)
{
    using Tz = decltype(std::declval<Ty>() * std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result[i][j] = Left_[i][j] * Right_;
        }
    }
    return result;
}

template<typename Ty, typename Tx, unsigned Rows, unsigned Cols>
inline constexpr auto numer::operator/(
    const Matrix<Ty, Rows, Cols>& Left_,
    const Tx& Right_)
{
    using Tz = decltype(std::declval<Ty>() / std::declval<Tx>());
    Matrix<Tz, Rows, Cols> result{};
    for (unsigned i = 0; i < Rows; ++i) {
        for (unsigned j = 0; j < Cols; ++j) {
            result[i][j] = Left_[i][j] / Right_;
        }
    }
    return result;
}