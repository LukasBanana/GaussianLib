/*
 * Algebra.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ALGEBRA_H__
#define __GS_ALGEBRA_H__


#include "Macros.h"
#include "Details.h"
#include "Determinant.h"
#include "Inverse.h"

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class Vector2T;
template <typename T> class Vector3T;
template <typename T> class Vector4T;


/* --- Global Functions --- */

//! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
T Angle(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos( Dot(lhs, rhs) / (lhs.Length()*rhs.Length()) );
}

//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
T AngleNorm(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos(Dot(lhs, rhs));
}

//! Returns the squared length of the specified vector.
template <template <typename> class Vec, typename T>
T LengthSq(const Vec<T>& vec)
{
    return Dot(vec, vec);
}

//! Returns the length (euclidian norm) of the specified vector.
template <template <typename> class Vec, typename T>
T Length(const Vec<T>& vec)
{
    return std::sqrt(LengthSq(vec));
}

//! Returns the squared distance between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
T DistanceSq(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return LengthSq(result);
}

//! Returns the distance between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
T Distance(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return Length(result);
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
T Dot(const Vec<T>& lhs, const Vec<T>& rhs)
{
    T result = T(0);
    
    for (std::size_t i = 0; i < Vec<T>::components; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T>
Vec<T> Cross(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return Vec<T>(
        lhs.y*rhs.z - rhs.y*lhs.z,
        rhs.x*lhs.z - lhs.x*rhs.z,
        lhs.x*rhs.y - rhs.x*lhs.y
    );
}

//! Normalizes the specified vector to the unit length of 1.
template <template <typename> class Vec, typename T>
void Normalize(Vec<T>& vec)
{
    auto len = LengthSq(vec);
    if (len != Real(0) && len != Real(1))
    {
        len = T(1) / std::sqrt(len);
        vec *= len;
    }
}

/**
\brief Computes a linear interpolation between the point 'a' and the point 'b'.
\return Equivalent to: a*(1-t) + b*t
*/
template <typename T, typename I>
T Lerp(const T& a, const T& b, const I& t)
{
    auto result = b;
    result -= a;
    result *= t;
    result += a;
    return result;
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
Quat<T> Slerp(const Quat<T>& from, Quat<T>& to, const T& t)
{
    Quat<T> q;
    q.Slerp(from, to, t);
    return q;
}


/* --- Global Operators --- */

#ifdef GS_ROW_VECTORS

//! \brief Multiplies the N-dimensional vector with the NxN matrix with.
template <template <typename> class Vec, typename T, std::size_t N>
Vec<T> operator * (const Vec<T>& vec, const Matrix<T, N, N>& mat)
{
    static_assert(
        Vec<T>::components == N,
        "only N-dimensional vectors can be multiplied with NxN matrices"
    );

    Vec<T> result;

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
template <template <typename> class Vec, typename T, std::size_t N>
Vec<T> operator * (const Matrix<T, N, N>& mat, const Vec<T>& vec)
{
    static_assert(
        Vec<T>::components == N,
        "only NxN matrices can be multiplied with N-dimensional vectors"
    );

    Vec<T> result;

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
