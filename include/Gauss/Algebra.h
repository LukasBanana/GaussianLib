/*
 * Algebra.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ALGEBRA_H
#define GS_ALGEBRA_H


#include <Gauss/Macros.h>
#include <Gauss/Details.h>
#include <Gauss/Determinant.h>
#include <Gauss/Inverse.h>
#include <Gauss/Decl.h>
#include <Gauss/Real.h>
#include <Gauss/Tags.h>
#include <Gauss/ScalarType.h>

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <limits>


namespace Gs
{


/* --- Global Functions --- */

//! Returns the value of 1 + 2 + ... + n = n*(n+1)/2.
template <typename T>
GS_NODISCARD
T GaussianSum(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSum function only allows integral types");
    return n*(n + T(1))/T(2);
}

//! Returns the value of 1^2 + 2^2 + ... + n^2 = n*(n+1)*(2n+1)/6.
template <typename T>
GS_NODISCARD
T GaussianSumSq(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSumSq function only allows integral types");
    return n*(n + T(1))*(n*T(2) + T(1))/T(6);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
GS_NODISCARD
T NormalDistribution(const T& x, const T& mean, const T& variance)
{
    return std::exp(-(x - mean)*(x - mean) / (variance + variance)) / std::sqrt(T(2) * T(Gs::pi) * variance);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
GS_NODISCARD
T NormalDistribution(const T& x)
{
    return std::exp(-(x*x) / T(2)) / std::sqrt(T(2) * T(Gs::pi));
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar Dot(const TVector& lhs, const TVector& rhs)
{
    TScalar result = TScalar(0);

    for (std::size_t i = 0; i < VectorType<TVector>::elements; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <typename TVector>
void Cross(TVector& result, const TVector& lhs, const TVector& rhs)
{
    static_assert(VectorType<TVector>::elements == 3, "Vector type must have exactly three components");
    result[0] = lhs[1]*rhs[2] - rhs[1]*lhs[2];
    result[1] = rhs[0]*lhs[2] - lhs[0]*rhs[2];
    result[2] = lhs[0]*rhs[1] - rhs[0]*lhs[1];
}

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <typename TVector>
GS_NODISCARD
TVector Cross(const TVector& lhs, const TVector& rhs)
{
    TVector result{ UninitializeTag{} };
    Cross<TVector>(result, lhs, rhs);
    return result;
}

//! Returns the squared length of the specified vector.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar LengthSq(const TVector& vec)
{
    return Dot<TVector, TScalar>(vec, vec);
}

//! Returns the length (euclidian norm) of the specified vector.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar Length(const TVector& vec)
{
    return std::sqrt(LengthSq<TVector, TScalar>(vec));
}

//! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar Angle(const TVector& lhs, const TVector& rhs)
{
    return std::acos( Dot<TVector, TScalar>(lhs, rhs) / (Length<TVector, TScalar>(lhs) * Length<TVector, TScalar>(rhs)) );
}

//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar AngleNorm(const TVector& lhs, const TVector& rhs)
{
    return std::acos(Dot<TVector, TScalar>(lhs, rhs));
}

//! Returns the squared distance between the two vectors 'lhs' and 'rhs'.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar DistanceSq(const TVector& lhs, const TVector& rhs)
{
    TVector result = rhs;
    result -= lhs;
    return LengthSq<TVector, TScalar>(result);
}

//! Returns the distance between the two vectors 'lhs' and 'rhs'.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TScalar Distance(const TVector& lhs, const TVector& rhs)
{
    TVector result = rhs;
    result -= lhs;
    return Length<TVector, TScalar>(result);
}

//! Returns the reflected vector of the incident vector for the specified surface normal.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TVector Reflect(const TVector& incident, const TVector& normal)
{
    /* Compute reflection as: I - N x Dot(N, I) x 2 */
    TVector v = normal;
    v *= (Dot<TVector, TScalar>(normal, incident) * TScalar(-2));
    v += incident;
    return v;
}

//! Normalizes the specified vector to the unit length of 1.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
void Normalize(TVector& outVec, const TVector& inVec)
{
    TScalar len = LengthSq<TVector, TScalar>(inVec);
    if (len != TScalar(0) && len != TScalar(1))
    {
        len = TScalar(1) / std::sqrt(len);
        outVec = inVec * len;
    }
    else
        outVec = inVec;
}

//! Normalizes the specified vector to the unit length of 1.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TVector Normalize(const TVector& vec)
{
    TVector result{ UninitializeTag{} };
    Normalize<TVector, TScalar>(result, vec);
    return result;
}

//! Resizes the specified vector to the specified length.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
void Resize(TVector& outVec, const TVector& inVec, const TScalar& length)
{
    TScalar len = LengthSq<TVector, TScalar>(inVec);
    if (len != TScalar(0))
    {
        len = length / std::sqrt(len);
        outVec = inVec * len;
    }
    else
        outVec = inVec;
}

//! Resizes the specified vector to the specified length.
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
GS_NODISCARD
TVector Resize(const TVector& vec, const TScalar& length)
{
    TVector result{ UninitializeTag{} };
    Resize<TVector, TScalar>(result, vec, length);
    return result;
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
    /* Return (b - a) * t + a */
    T x = b;
    x -= a;
    x *= t;
    x += a;
    return x;
}

/**
\brief Mixes the two values by their scalings.
\return Equivalent to: v0*scale0 + v1*scale1
*/
template <typename T, typename I>
T Mix(const T& v0, const T& v1, const I& scale0, const I& scale1)
{
    return v0*scale0 + v1*scale1;
}

/**
\brief Clamps the input value 'x' into the range [0, 1].
\return max{ 0, min{ x, 1 } }
*/
template <typename T>
T Saturate(const T& x)
{
    return std::max<T>(T(0), std::min<T>(x, T(1)));
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
template <typename TVector, typename TScalar = typename VectorType<TVector>::ScalarType>
TVector Slerp(const TVector& from, const TVector& to, const TScalar& t)
{
    TScalar omega, cosom, sinom;
    TScalar scale0, scale1;

    /* Calculate cosine */
    cosom = Dot<TVector, TScalar>(from, to);

    /* Adjust signs (if necessary) */
    if (cosom < TScalar(0))
    {
        cosom = -cosom;
        scale1 = TScalar(-1);
    }
    else
        scale1 = TScalar(1);

    /* Calculate coefficients */
    if ((TScalar(1) - cosom) > std::numeric_limits<TScalar>::epsilon())
    {
        /* Standard case (slerp) */
        omega = std::acos(cosom);
        sinom = std::sin(omega);
        scale0 = std::sin((TScalar(1) - t) * omega) / sinom;
        scale1 *= std::sin(t * omega) / sinom;
    }
    else
    {
        /* 'from' and 'to' quaternions are very close, so we can do a linear interpolation */
        scale0 = TScalar(1) - t;
        scale1 *= t;
    }

    /* Calculate final values */
    return Mix(from, to, scale0, scale1);
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

//! Returns the reciprocal of the specified scalar value.
template <typename T>
T Rcp(const T& x)
{
	return T(1) / x;
}

//! Returns the per-component reciprocal of the specified N-dimensional vector.
template <typename T, std::size_t N>
Vector<T, N> Rcp(const Vector<T, N>& vec)
{
	Vector<T, N> vecRcp{ UninitializeTag{} };

	for (std::size_t i = 0; i < N; ++i)
		vecRcp[i] = T(1) / vec[i];

	return vecRcp;
}

//! Returns the per-component reciprocal of the specified NxM-dimensional matrix.
template <typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Rcp(const Matrix<T, N, M>& mat)
{
	Matrix<T, N, M> matRcp{ UninitializeTag{} };

	for (std::size_t i = 0; i < N*M; ++i)
		matRcp[i] = T(1) / mat[i];

	return matRcp;
}

//! Rescales the specified value 't' from the first range [lower0, upper0] into the second range [lower1, upper1].
template <typename T, typename I>
T Rescale(const T& t, const I& lower0, const I& upper0, const I& lower1, const I& upper1)
{
    /* Return (((t - lower0) / (upper0 - lower0)) * (upper1 - lower1) + lower1) */
    T x = t;
    x -= T(lower0);
    x /= (upper0 - lower0);
    x *= (upper1 - lower1);
    x += T(lower1);
    return x;
}


/* --- Global Operators --- */

/**
\brief Multiplies the N-dimensional row-vector with the NxM matrix.
\remarks This is equivalent to: Transpose(rhs) * lhs.
*/
template <typename T, std::size_t Rows, std::size_t Cols>
Vector<T, Cols> operator * (const Vector<T, Rows>& lhs, const Matrix<T, Rows, Cols>& rhs)
{
    Vector<T, Cols> result;

    for (std::size_t c = 0; c < Cols; ++c)
    {
        result[c] = T(0);
        for (std::size_t r = 0; r < Rows; ++r)
            result[c] += rhs(r, c)*lhs[r];
    }

    return result;
}

/**
\brief Multiplies the NxM matrix with the M-dimensional column-vector.
\remarks This is equivalent to: rhs * Transpose(lhs).
*/
template <typename T, std::size_t Rows, std::size_t Cols>
Vector<T, Rows> operator * (const Matrix<T, Rows, Cols>& lhs, const Vector<T, Cols>& rhs)
{
    Vector<T, Rows> result;

    for (std::size_t r = 0; r < Rows; ++r)
    {
        result[r] = T(0);
        for (std::size_t c = 0; c < Cols; ++c)
            result[r] += lhs(r, c)*rhs[c];
    }

    return result;
}


} // /namespace Gs


#endif



// ================================================================================
