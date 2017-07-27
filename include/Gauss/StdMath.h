/*
 * StdMath.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_STDMATH_H
#define GS_STDMATH_H


#include "Vector.h"
#include "Matrix.h"
#include <cmath>
#include <functional>


namespace Gs
{


/* --- Global Functions --- */

#define GS_DECL_STDMATH_FUNC(NAME)                              \
    template <typename T, std::size_t N>                        \
    Vector<T, N> NAME(const Vector<T, N>& x)                    \
    {                                                           \
        Vector<T, N> y { UninitializeTag{} };                   \
        for (std::size_t i = 0; i < N; ++i)                     \
            y[i] = std::NAME(x[i]);                             \
        return y;                                               \
    }                                                           \
    template <typename T, std::size_t Rows, std::size_t Cols>   \
    Matrix<T, Rows, Cols> NAME(const Matrix<T, Rows, Cols>& x)  \
    {                                                           \
        Matrix<T, Rows, Cols> y { UninitializeTag{} };          \
        for (std::size_t i = 0; i < Rows*Cols; ++i)             \
            y[i] = std::NAME(x[i]);                             \
        return y;                                               \
    }

#define GS_DECL_STDMATH_FUNC2(NAME)                                                                 \
    template <typename T, std::size_t N>                                                            \
    Vector<T, N> NAME(const Vector<T, N>& x1, const Vector<T, N>& x2)                               \
    {                                                                                               \
        Vector<T, N> y { UninitializeTag{} };                                                       \
        for (std::size_t i = 0; i < N; ++i)                                                         \
            y[i] = std::NAME(x1[i], x2[i]);                                                         \
        return y;                                                                                   \
    }                                                                                               \
    template <typename T, std::size_t Rows, std::size_t Cols>                                       \
    Matrix<T, Rows, Cols> NAME(const Matrix<T, Rows, Cols>& x1, const Matrix<T, Rows, Cols>& x2)    \
    {                                                                                               \
        Matrix<T, Rows, Cols> y { UninitializeTag{} };                                              \
        for (std::size_t i = 0; i < Rows*Cols; ++i)                                                 \
            y[i] = std::NAME(x1[i], x2[i]);                                                         \
        return y;                                                                                   \
    }


GS_DECL_STDMATH_FUNC( exp )
GS_DECL_STDMATH_FUNC( exp2 )
GS_DECL_STDMATH_FUNC( log )
GS_DECL_STDMATH_FUNC( log10 )
GS_DECL_STDMATH_FUNC( log2 )

GS_DECL_STDMATH_FUNC2( pow )
GS_DECL_STDMATH_FUNC( sqrt )

GS_DECL_STDMATH_FUNC( sin )
GS_DECL_STDMATH_FUNC( cos )
GS_DECL_STDMATH_FUNC( tan )
GS_DECL_STDMATH_FUNC( asin )
GS_DECL_STDMATH_FUNC( acos )
GS_DECL_STDMATH_FUNC( atan )
GS_DECL_STDMATH_FUNC2( atan2 )

GS_DECL_STDMATH_FUNC( sinh )
GS_DECL_STDMATH_FUNC( cosh )
GS_DECL_STDMATH_FUNC( tanh )
GS_DECL_STDMATH_FUNC( asinh )
GS_DECL_STDMATH_FUNC( acosh )
GS_DECL_STDMATH_FUNC( atanh )

GS_DECL_STDMATH_FUNC( ceil )
GS_DECL_STDMATH_FUNC( floor )


#undef GS_DECL_STDMATH_FUNC


} // /namespace Gs


#endif



// ================================================================================
