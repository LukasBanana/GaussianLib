/*
 * Algebra.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ALGEBRA_H
#define GS_ALGEBRA_H


#include "Macros.h"
#include "Details.h"
#include "Determinant.h"
#include "Inverse.h"
#include "Decl.h"
#include "Real.h"

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>
#include <type_traits>


namespace Gs
{


/* --- Global Functions --- */

//! Returns the value of 1 + 2 + ... + n = n*(n+1)/2.
template <typename T>
T GaussianSum(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSum function only allows integral types");
    return n*(n + T(1))/T(2);
}

//! Returns the value of 1^2 + 2^2 + ... + n^2 = n*(n+1)*(2n+1)/6.
template <typename T>
T GaussianSumSq(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSumSq function only allows integral types");
    return n*(n + T(1))*(n*T(2) + T(1))/T(6);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
T NormalDistribution(const T& x, const T& mean, const T& variance)
{
    return std::exp(-(x - mean)*(x - mean) / (variance + variance)) / std::sqrt(T(2) * T(Gs::pi) * variance);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
T NormalDistribution(const T& x)
{
    return std::exp(-(x*x) / T(2)) / std::sqrt(T(2) * T(Gs::pi));
}

//! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
template <typename T, std::size_t N>
T Angle(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    return std::acos( Dot(lhs, rhs) / (lhs.Length()*rhs.Length()) );
}

//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template <typename T, std::size_t N>
T AngleNorm(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    return std::acos(Dot(lhs, rhs));
}

//! Returns the squared length of the specified vector.
template <typename T, std::size_t N>
T LengthSq(const Vector<T, N>& vec)
{
    return Dot(vec, vec);
}

//! Returns the length (euclidian norm) of the specified vector.
template <typename T, std::size_t N>
T Length(const Vector<T, N>& vec)
{
    return std::sqrt(LengthSq(vec));
}

//! Returns the squared distance between the two vectors 'lhs' and 'rhs'.
template <typename T, std::size_t N>
T DistanceSq(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return LengthSq(result);
}

//! Returns the distance between the two vectors 'lhs' and 'rhs'.
template <typename T, std::size_t N>
T Distance(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return Length(result);
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <typename T, std::size_t N>
T Dot(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    T result = T(0);
    
    for (std::size_t i = 0; i < N; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <typename T>
T Dot(const QuaternionT<T>& lhs, const QuaternionT<T>& rhs)
{
    T result = T(0);
    
    for (std::size_t i = 0; i < QuaternionT<T>::components; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <typename T>
Vector<T, 3> Cross(const Vector<T, 3>& lhs, const Vector<T, 3>& rhs)
{
    return Vector<T, 3>(
        lhs.y*rhs.z - rhs.y*lhs.z,
        rhs.x*lhs.z - lhs.x*rhs.z,
        lhs.x*rhs.y - rhs.x*lhs.y
    );
}

//! Normalizes the specified vector to the unit length of 1.
template <typename T, std::size_t N>
void Normalize(Vector<T, N>& vec)
{
    auto len = LengthSq(vec);
    if (len != T(0) && len != T(1))
    {
        len = T(1) / std::sqrt(len);
        vec *= len;
    }
}

//! Normalizes the specified vector to the unit length of 1.
template <typename T>
void Normalize(QuaternionT<T>& q)
{
    auto len = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
    if (len != T(0) && len != T(1))
    {
        len = T(1) / std::sqrt(len);
        q *= len;
    }
}

//! Resizes the specified vector to the specified length.
template <typename T, std::size_t N>
void Resize(Vector<T, N>& vec, const T& length)
{
    auto len = LengthSq(vec);
    if (len != T(0))
    {
        len = length / std::sqrt(len);
        vec *= len;
    }
}

/**
\brief Computes a linear interpolation between the point 'a' and the point 'b'.
\return Equivalent to: a*(1-t) + b*t
*/
template <typename T, typename I>
void Lerp(T& x, const T& a, const T& b, const I& t)
{
    x = b;
    x -= a;
    x *= t;
    x += a;
}

/**
\brief Computes a linear interpolation between the point 'a' and the point 'b'.
\return Equivalent to: a*(1-t) + b*t
*/
template <typename T, typename I>
T Lerp(const T& a, const T& b, const I& t)
{
    T x = b;
    x -= a;
    x *= t;
    x += a;
    return x;
}

/**
\brief Clamps the value 'x' into the range [minima, maxima].
\return max{ minima, min{ x, maxima } }
*/
template <typename T>
T Clamp(const T& x, const T& minima, const T& maxima)
{
    if (x <= minima)
        return minima;
    if (x >= maxima)
        return maxima;
    return x;
}

/**
\brief Returns the spherical linear interpolation between the two quaternions 'from' and 'to'.
\see QuaternionT::Slerp
*/
template <template <typename> class Quat, typename T>
Quat<T> Slerp(const Quat<T>& from, const Quat<T>& to, const T& t)
{
    Quat<T> q;
    q.Slerp(from, to, t);
    return q;
}

/**
\brief Returns a smooth 'hermite interpolation' in the range [0, 1].
\remarks This hermite interpolation is: 3x^2 - 2x^3.
*/
template <typename T>
T SmoothStep(const T& x)
{
    return x*x * (T(3) - x*T(2));
}

/**
\brief Returns a smooth 'hermite interpolation' in the range [0, 1].
\remarks This hermite interpolation is: 6x^5 - 15x^4 + 10x^3.
*/
template <typename T>
T SmootherStep(const T& x)
{
    return x*x*x * (x*(x*T(6) - T(15)) + T(10));
}


/* --- Global Operators --- */

#ifdef GS_ROW_VECTORS

//! \brief Multiplies the N-dimensional vector with the NxN matrix with.
template <typename T, std::size_t N>
Vector<T, N> operator * (const Vector<T, N>& vec, const Matrix<T, N, N>& mat)
{
    Vector<T, N> result;

    for (std::size_t c = 0; c < N; ++c)
    {
        result[c] = T(0);
        for (std::size_t r = 0; r < N; ++r)
            result[c] += mat(r, c)*vec[r];
    }

    return result;
}

#else

//! \brief Multiplies the NxN matrix with the N-dimensional vector.
template <typename T, std::size_t N>
Vector<T, N> operator * (const Matrix<T, N, N>& mat, const Vector<T, N>& vec)
{
    Vector<T, N> result;

    for (std::size_t r = 0; r < N; ++r)
    {
        result[r] = T(0);
        for (std::size_t c = 0; c < N; ++c)
            result[r] += mat(r, c)*vec[c];
    }

    return result;
}

#endif


} // /namespace Gs


#endif



// ================================================================================
