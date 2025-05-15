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
#include <cmath>
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

        constexpr auto t() const { return transpose(); }
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

        constexpr auto cc() const { return conj(); }
        constexpr auto t() const { return transpose(); }
        constexpr auto hc() const { return hermiConj(); }
    };

    
    template<typename Ty, unsigned Dim>
    using Vec = Vec_impl<Ty, Dim, true>;// column vector

    template<typename Ty, unsigned Dim>
    using Vect = Vec_impl<Ty, Dim, false>;// row vector

    template<typename Ty>
    using Vec3 = Vec_impl<Ty, 3, true>;// 3d column vector

    template<typename Ty>
    using Vec3t = Vec_impl<Ty, 3, false>;// 3d row vector




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

        constexpr auto t() const { return transpose(); }
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

        constexpr auto cc() const { return conj(); }
        constexpr auto t() const { return transpose(); }
        constexpr auto hc() const { return hermiConj(); }
    };


    template<typename Ty, unsigned Dim>
    using SquareMatrix = Matrix<Ty, Dim, Dim>; // square matrix

    template<typename Ty, unsigned Dim>
    using Carre = Matrix<Ty, Dim, Dim>; // matrice carre



    template<typename Ty, unsigned Dim>
    inline constexpr auto makeIdentityMatrix() {
        SquareMatrix<Ty, Dim> result{};
        Ty zero = static_cast<Ty>(0);
        Ty one = static_cast<Ty>(1);
        for (unsigned i = 0; i < Dim; ++i) {
            for (unsigned j = 0; j < Dim; ++j) {
                if(i == j) result[i][j] = one;
                else result[i][j] = zero;
            }
        }
        return result;
    }

    template<typename Ty, unsigned Rows, unsigned Cols>
    inline void swapRows(Matrix<Ty, Rows, Cols>& M_, unsigned Row1_, unsigned Row2_) {
        for (unsigned i = 0; i < Cols; ++i) {
            std::swap(M_[Row1_][i], M_[Row2_][i]);
        }
    }

    template<typename Ty, unsigned Rows, unsigned Cols>
    inline void scaleRow(Matrix<Ty, Rows, Cols>& M_, const Ty& Multiplier_, unsigned Row_) {
        for (unsigned i = 0; i < Cols; ++i) {
            M_[Row_][i] *=  Multiplier_;
        }
    }

    template<typename Ty, unsigned Rows, unsigned Cols>
    inline void addScaledRowTo(Matrix<Ty, Rows, Cols>& M_, const Ty& Multiplier_, unsigned Src_row_, unsigned Dst_row_) {
        for (unsigned i = 0; i < Cols; ++i) {
            M_[Dst_row_][i] += M_[Src_row_][i] * Multiplier_;
        }
    }

    //you should check reversibility by urself
    template<typename Ty, unsigned Dim>
    inline auto inverse(SquareMatrix<Ty, Dim> A_) {
        SquareMatrix<Ty, Dim> E = makeIdentityMatrix<Ty, Dim>();
        Ty one = static_cast<Ty>(1);

        for (unsigned col = 0; col < Dim; ++col) {
            unsigned pivot_row = col;
            Ty max_val = abs(A_[col][col]);

            for (unsigned i = col + 1; i < Dim; ++i) {
                Ty current = abs(A_[i][col]);
                if (current > max_val) {
                    max_val = current;
                    pivot_row = i;
                }
            }

            if (pivot_row != col) {
                swapRows(A_, col, pivot_row);
                swapRows(E, col, pivot_row);
            }

            const Ty pivot = A_[col][col];
            scaleRow(A_, one / pivot, col);
            scaleRow(E, one / pivot, col);

            for (unsigned i = 0; i < Dim; ++i) {
                if (i != col) {
                    const Ty factor = A_[i][col];
                    addScaledRowTo(A_, -factor, col, i);
                    addScaledRowTo(E, -factor, col, i);
                }
            }
        }

        return E;
    }


    template<typename Ty, unsigned Dim>
    inline constexpr auto trace(SquareMatrix<Ty, Dim> M_) {
        Ty result = static_cast<Ty>(0);
        for (unsigned i = 0; i < Dim; ++i) {
            result += M_[i][i];
        }
        return result;
    }


    template<typename Ty, unsigned Dim>
    inline constexpr auto det(const SquareMatrix<Ty, Dim>& M_) {
        if constexpr (Dim == 1) {
            return M_[0][0];
        }
        else if constexpr (Dim == 2) {
            return M_[0][0] * M_[1][1] - M_[0][1] * M_[1][0];
        }
        else {
            Ty result = static_cast<Ty>(0);
            for (unsigned col = 0; col < Dim; ++col) {
                SquareMatrix<Ty, Dim - 1> sub{};
                for (unsigned i = 1; i < Dim; ++i) {
                    unsigned sub_col = 0;
                    for (unsigned j = 0; j < Dim; ++j) {
                        if (j != col) {
                            if constexpr (Dim - 1 > 0) {
                                sub[i - 1][sub_col] = M_[i][j];
                                sub_col += (sub_col < Dim - 2) ? 1 : 0;
                            }
                        }
                    }
                }
                const Ty sign = (col % 2 == 0) ? Ty(1) : Ty(-1);
                result += sign * M_[0][col] * det(sub);
            }
            return result;
        }
    }



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

    template<typename Ty, unsigned Dim>
    inline auto operator^(const SquareMatrix<Ty, Dim>& M_, const int Expo_);


    //used for the case: t->0, e^At = 1 + At + (At)^2/2 + (At)^3/3! + ... 
    template<unsigned Order, typename Ty, unsigned Dim>
    constexpr auto expm_approx(const SquareMatrix<Ty, Dim>& M_) {
        SquareMatrix<Ty, Dim> I = makeIdentityMatrix<Ty, Dim>();
        SquareMatrix<Ty, Dim> result = I;
        for (unsigned i = Order; i > 0; --i) {
            result = I + M_ * result / static_cast<double>(i);
        }
        return result;
    }



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
    Tz result{};
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
    SquareMatrix<Tz, Dim> result{};
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
        Tz component{};
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
        Tz component{};
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
            Tz entry{};
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

template<typename Ty, unsigned Dim>
inline auto numer::operator^(const SquareMatrix<Ty, Dim>& M_, const int Expo_) {
    if (Expo_ == 0) return makeIdentityMatrix<Ty, Dim>();
    SquareMatrix<Ty, Dim> result = makeIdentityMatrix<Ty, Dim>();
    SquareMatrix<Ty, Dim> base = Expo_ > 0 ? M_ : inverse(M_);
    int n = abs(Expo_);
    while (n > 0) {
        if (n % 2) result = result * base;
        base = base * base;
        n /= 2;
    }
    return result;
}
