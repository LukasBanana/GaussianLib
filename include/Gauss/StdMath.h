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

#define GM_DECL_STDMATH_FUNC(NAME)                              \
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

#define GM_DECL_STDMATH_FUNC2(NAME)                                                                 \
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


GM_DECL_STDMATH_FUNC( exp )
GM_DECL_STDMATH_FUNC( exp2 )
GM_DECL_STDMATH_FUNC( log )
GM_DECL_STDMATH_FUNC( log10 )
GM_DECL_STDMATH_FUNC( log2 )

GM_DECL_STDMATH_FUNC2( pow )
GM_DECL_STDMATH_FUNC( sqrt )

GM_DECL_STDMATH_FUNC( sin )
GM_DECL_STDMATH_FUNC( cos )
GM_DECL_STDMATH_FUNC( tan )
GM_DECL_STDMATH_FUNC( asin )
GM_DECL_STDMATH_FUNC( acos )
GM_DECL_STDMATH_FUNC( atan )
GM_DECL_STDMATH_FUNC2( atan2 )

GM_DECL_STDMATH_FUNC( sinh )
GM_DECL_STDMATH_FUNC( cosh )
GM_DECL_STDMATH_FUNC( tanh )
GM_DECL_STDMATH_FUNC( asinh )
GM_DECL_STDMATH_FUNC( acosh )
GM_DECL_STDMATH_FUNC( atanh )

GM_DECL_STDMATH_FUNC( ceil )
GM_DECL_STDMATH_FUNC( floor )


#undef GM_DECL_STDMATH_FUNC


} // /namespace Gs


#endif



// ================================================================================
