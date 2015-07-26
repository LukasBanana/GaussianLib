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

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class SparseMatrix4T;

template <typename T, std::size_t Rows, std::size_t Cols> class Matrix;


/* --- Global Functions --- */

//! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> T Angle(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos( Dot(lhs, rhs) / (lhs.Length()*rhs.Length()) );
}

//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> T AngleNorm(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos(Dot(lhs, rhs));
}

//! Returns the squared length of the specified vector.
template <template <typename> class Vec, typename T> T LengthSq(const Vec<T>& vec)
{
    return Dot(vec, vec);
}

//! Returns the length (euclidian norm) of the specified vector.
template <template <typename> class Vec, typename T> T Length(const Vec<T>& vec)
{
    return std::sqrt(LengthSq(vec));
}

//! Returns the squared distance between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> T DistanceSq(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return LengthSq(result);
}

//! Returns the distance between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> T Distance(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return Length(result);
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> T Dot(const Vec<T>& lhs, const Vec<T>& rhs)
{
    T result = T(0);
    
    for (std::size_t i = 0; i < Vec<T>::components; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <template <typename> class Vec, typename T> Vec<T> Cross(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return Vec<T>(
        lhs.y*rhs.z - rhs.y*lhs.z,
        rhs.x*lhs.z - lhs.x*rhs.z,
        lhs.x*rhs.y - rhs.x*lhs.y
    );
}

//! Normalizes the specified vector to the unit length of 1.
template <template <typename> class Vec, typename T> void Normalize(Vec<T>& vec)
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
template <typename T, typename I> T Lerp(const T& a, const T& b, const I& t)
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
template <typename T> T Clamp(const T& x, const T& minima, const T& maxima)
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

/**
\brief Computes a free rotation around an axis and stores the result into the matrix 'm'.
\tparam M Specifies the matrix type. This should be Matrix3, Matrix4, or SparseMatrix4.
This type must implement the following interface:
\code
static const std::size_t rows;    // >= 3
static const std::size_t columns; // >= 3
\endcode
\tparam Vec Specifies the vector type. This should be Vector3, or Vector4.
\tparam T Specifies the data type. This should should be float or double.
\param[out] m Specifies the resulting matrix.
\param[in] axis Specifies the rotation axis. This must be normalized!
\param[in] angle Specifies the rotation angle (in radians).
*/
template <class M, template <typename> class V, typename T>
void MakeFreeRotation(M& mat, const V<T>& axis, const T& angle)
{
    static_assert(
        M::rows >= 3 && M::columns >= 3,
        "free rotation can only be computed for matrices with at least 3 rows and 3 column"
    );

    /* Setup rotation values */
    const T  c  = std::cos(angle);
    const T  s  = std::sin(angle);
    const T  cc = T(1) - c;

    const T& x  = axis.x;
    const T& y  = axis.y;
    const T& z  = axis.z;

    /* Perform matrix rotation */
    mat(0, 0) = x*x*cc + c;
    mat(0, 1) = y*x*cc + z*s;
    mat(0, 2) = x*z*cc - y*s;

    mat(1, 0) = x*y*cc - z*s;
    mat(1, 1) = y*y*cc + c;
    mat(1, 2) = y*z*cc + x*s;

    mat(2, 0) = x*z*cc + y*s;
    mat(2, 1) = y*z*cc - x*s;
    mat(2, 2) = z*z*cc + c;
}

template <class M, template <typename> class Q, typename T>
void MatrixToQuaternion(M mat, Q<T>& quaternion)
{
    static_assert(
        M::rows >= 3 && M::columns >= 3,
        "matrices can only be converted to quaternions, when they have at least 3 rows and 3 column"
    );

    /* Make sure the matrix has an identitiy scaling */
    //mat.SetScale({ 1, 1, 1 });

    /* Only get the trace of the 3x3 upper left matrix */
    const T trace = mat(0, 0) + mat(1, 1) + mat(2, 2) + T(1);
    
    if (trace > T(0))
    {
        const T s = T(2) * std::sqrt(trace);
        quaternion.x = (mat(1, 2) - mat(2, 1)) / s;
        quaternion.y = (mat(2, 0) - mat(0, 2)) / s;
        quaternion.z = (mat(0, 1) - mat(1, 0)) / s;
        quaternion.w = T(0.25) * s;
    }
    else
    {
        if (mat(0, 0) > mat(1, 1) && mat(0, 0) > mat(2, 2))
        {
            const T s = T(2) * std::sqrt(T(1) + mat(0, 0) - mat(1, 1) - mat(2, 2));
            quaternion.x = T(0.25) * s;
            quaternion.y = (mat(1, 0) + mat(0, 1) ) / s;
            quaternion.z = (mat(0, 2) + mat(2, 0) ) / s;
            quaternion.w = (mat(1, 2) - mat(2, 1) ) / s;
        }
        else if (mat(1, 1) > mat(2, 2))
        {
            const T s = T(2) * std::sqrt(T(1) + mat(1, 1) - mat(0, 0) - mat(2, 2));
            quaternion.x = (mat(1, 0) + mat(0, 1) ) / s;
            quaternion.y = T(0.25) * s;
            quaternion.z = (mat(2, 1) + mat(1, 2) ) / s;
            quaternion.w = (mat(2, 0) - mat(0, 2) ) / s;
        }
        else
        {
            const T s = T(2) * std::sqrt(T(1) + mat(2, 2) - mat(0, 0) - mat(1, 1));
            quaternion.x = (mat(2, 0) + mat(0, 2) ) / s;
            quaternion.y = (mat(2, 1) + mat(1, 2) ) / s;
            quaternion.z = T(0.25) * s;
            quaternion.w = (mat(0, 1) - mat(1, 0) ) / s;
        }
    }

    Normalize(quaternion);
}


/* --- Determinant Functions --- */

/**
\brief Computes the determinant of an arbitrary NxN matrix.
\tparam M Specifies the matrix type. This should be "Matrix".
\tparam T Specifies the data type. This should be float or double.
\tparam Rows Specifies the rows of the matrix.
\tparam Cols Specifies the columns of the matrix.
\remarks The template arguments 'Rows' and 'Cols' must be equal, otherwise a compile time error will occur,
since a determinant is only defined for squared matrices.
\param[in] m Specifies the squared matrix for which the determinant is to be computed.
*/
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
T Determinant(const M<T, Rows, Cols>& m)
{
    static_assert(Rows == Cols, "determinants can only be computed for squared matrices");
    using Helper = Details::MatrixHelper<M, T, Rows, Cols>;
    return Helper::OrderedDeterminant(Helper::MatrixToArray(m), Rows);
}

//! Computes the determinant of the specified 1x1 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 1, 1>& m)
{
    return m(0, 0);
}

//! Computes the determinant of the specified 2x2 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 2, 2>& m)
{
    return m(0, 0)*m(1, 1) - m(1, 0)*m(0, 1);
}

//! Computes the determinant of the specified 3x3 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 3, 3>& m)
{
    return
        ( m(0, 0) * m(1, 1) * m(2, 2) ) + ( m(1, 0) * m(2, 1) * m(0, 2) ) + ( m(2, 0) * m(0, 1) * m(1, 2) ) -
        ( m(0, 2) * m(1, 1) * m(2, 0) ) - ( m(0, 2) * m(2, 1) * m(0, 0) ) - ( m(2, 2) * m(0, 1) * m(1, 0) );
}

//! Computes the determinant of the specified 4x4 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 4, 4>& m)
{
    return
        ( m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) ) * ( m(2, 2) * m(3, 3) - m(3, 2) * m(2, 3) ) -
        ( m(0, 0) * m(2, 1) - m(2, 0) * m(0, 1) ) * ( m(1, 2) * m(3, 3) - m(3, 2) * m(1, 3) ) +
        ( m(0, 0) * m(3, 1) - m(3, 0) * m(0, 1) ) * ( m(1, 2) * m(2, 3) - m(2, 2) * m(1, 3) ) +
        ( m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1) ) * ( m(0, 2) * m(3, 3) - m(3, 2) * m(0, 3) ) -
        ( m(1, 0) * m(3, 1) - m(3, 0) * m(1, 1) ) * ( m(0, 2) * m(2, 3) - m(2, 2) * m(0, 3) ) +
        ( m(2, 0) * m(3, 1) - m(3, 0) * m(2, 1) ) * ( m(0, 2) * m(1, 3) - m(1, 2) * m(0, 3) );
}

//! Computes the determinant of the specified sparse 4x4 matrix 'm'.
template <typename T> T Determinant(const SparseMatrix4T<T>& m)
{
    return
        ( m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) ) * ( m(2, 2) ) -
        ( m(0, 0) * m(2, 1) - m(2, 0) * m(0, 1) ) * ( m(1, 2) ) +
        ( m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1) ) * ( m(0, 2) );
}


/* --- "Inverse" Functions --- */

//! Computes the inverse of the specified matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
bool Inverse(M<T, Rows, Cols>& inv, const M<T, Rows, Cols>& m)
{
    static_assert(Rows == Cols, "inverses can only be computed for squared matrices");
    //using Helper = Details::MatrixHelper<M, T, Rows, Cols>;
    //return Helper::OrderedInverse(inv, Helper::MatrixToArray(m), Rows);
    return false;//!!!
}

//! Computes the inverse of the specified 2x2 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
bool Inverse(M<T, 2, 2>& inv, const M<T, 2, 2>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv(0, 0) = d * (  m(1, 1) );
    inv(0, 1) = d * ( -m(0, 1) );
    inv(1, 0) = d * ( -m(1, 0) );
    inv(1, 1) = d * (  m(0, 0) );

    return true;
}

//! Computes the inverse of the specified 3x3 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
bool Inverse(M<T, 3, 3>& inv, const M<T, 3, 3>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv(0, 0) = d * ( m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2) );
    inv(1, 0) = d * ( m(2, 0) * m(1, 2) - m(1, 0) * m(2, 2) );
    inv(2, 0) = d * ( m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1) );
    inv(0, 1) = d * ( m(1, 1) * m(0, 2) - m(0, 1) * m(2, 2) );
    inv(1, 1) = d * ( m(0, 0) * m(2, 2) - m(2, 0) * m(0, 2) );
    inv(2, 1) = d * ( m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) );
    inv(0, 2) = d * ( m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2) );
    inv(1, 2) = d * ( m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2) );
    inv(2, 2) = d * ( m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) );

    return true;
}

//! Computes the inverse of the specified 4x4 matrix 'm'.
template <template <typename, std::size_t, std::size_t> class M, typename T>
bool Inverse(M<T, 4, 4>& inv, const M<T, 4, 4>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv(0, 0) = d * ( m(1, 1) * (m(2, 2) * m(3, 3) - m(3, 2) * m(2, 3)) + m(2, 1) * (m(3, 2) * m(1, 3) - m(1, 2) * m(3, 3)) + m(3, 1) * (m(1, 2) * m(2, 3) - m(2, 2) * m(1, 3)) );
    inv(1, 0) = d * ( m(1, 2) * (m(2, 0) * m(3, 3) - m(3, 0) * m(2, 3)) + m(2, 2) * (m(3, 0) * m(1, 3) - m(1, 0) * m(3, 3)) + m(3, 2) * (m(1, 0) * m(2, 3) - m(2, 0) * m(1, 3)) );
    inv(2, 0) = d * ( m(1, 3) * (m(2, 0) * m(3, 1) - m(3, 0) * m(2, 1)) + m(2, 3) * (m(3, 0) * m(1, 1) - m(1, 0) * m(3, 1)) + m(3, 3) * (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) );
    inv(3, 0) = d * ( m(1, 0) * (m(3, 1) * m(2, 2) - m(2, 1) * m(3, 2)) + m(2, 0) * (m(1, 1) * m(3, 2) - m(3, 1) * m(1, 2)) + m(3, 0) * (m(2, 1) * m(1, 2) - m(1, 1) * m(2, 2)) );
    inv(0, 1) = d * ( m(2, 1) * (m(0, 2) * m(3, 3) - m(3, 2) * m(0, 3)) + m(3, 1) * (m(2, 2) * m(0, 3) - m(0, 2) * m(2, 3)) + m(0, 1) * (m(3, 2) * m(2, 3) - m(2, 2) * m(3, 3)) );
    inv(1, 1) = d * ( m(2, 2) * (m(0, 0) * m(3, 3) - m(3, 0) * m(0, 3)) + m(3, 2) * (m(2, 0) * m(0, 3) - m(0, 0) * m(2, 3)) + m(0, 2) * (m(3, 0) * m(2, 3) - m(2, 0) * m(3, 3)) );
    inv(2, 1) = d * ( m(2, 3) * (m(0, 0) * m(3, 1) - m(3, 0) * m(0, 1)) + m(3, 3) * (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) + m(0, 3) * (m(3, 0) * m(2, 1) - m(2, 0) * m(3, 1)) );
    inv(3, 1) = d * ( m(2, 0) * (m(3, 1) * m(0, 2) - m(0, 1) * m(3, 2)) + m(3, 0) * (m(0, 1) * m(2, 2) - m(2, 1) * m(0, 2)) + m(0, 0) * (m(2, 1) * m(3, 2) - m(3, 1) * m(2, 2)) );
    inv(0, 2) = d * ( m(3, 1) * (m(0, 2) * m(1, 3) - m(1, 2) * m(0, 3)) + m(0, 1) * (m(1, 2) * m(3, 3) - m(3, 2) * m(1, 3)) + m(1, 1) * (m(3, 2) * m(0, 3) - m(0, 2) * m(3, 3)) );
    inv(1, 2) = d * ( m(3, 2) * (m(0, 0) * m(1, 3) - m(1, 0) * m(0, 3)) + m(0, 2) * (m(1, 0) * m(3, 3) - m(3, 0) * m(1, 3)) + m(1, 2) * (m(3, 0) * m(0, 3) - m(0, 0) * m(3, 3)) );
    inv(2, 2) = d * ( m(3, 3) * (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) + m(0, 3) * (m(1, 0) * m(3, 1) - m(3, 0) * m(1, 1)) + m(1, 3) * (m(3, 0) * m(0, 1) - m(0, 0) * m(3, 1)) );
    inv(3, 2) = d * ( m(3, 0) * (m(1, 1) * m(0, 2) - m(0, 1) * m(1, 2)) + m(0, 0) * (m(3, 1) * m(1, 2) - m(1, 1) * m(3, 2)) + m(1, 0) * (m(0, 1) * m(3, 2) - m(3, 1) * m(0, 2)) );
    inv(0, 3) = d * ( m(0, 1) * (m(2, 2) * m(1, 3) - m(1, 2) * m(2, 3)) + m(1, 1) * (m(0, 2) * m(2, 3) - m(2, 2) * m(0, 3)) + m(2, 1) * (m(1, 2) * m(0, 3) - m(0, 2) * m(1, 3)) );
    inv(1, 3) = d * ( m(0, 2) * (m(2, 0) * m(1, 3) - m(1, 0) * m(2, 3)) + m(1, 2) * (m(0, 0) * m(2, 3) - m(2, 0) * m(0, 3)) + m(2, 2) * (m(1, 0) * m(0, 3) - m(0, 0) * m(1, 3)) );
    inv(2, 3) = d * ( m(0, 3) * (m(2, 0) * m(1, 1) - m(1, 0) * m(2, 1)) + m(1, 3) * (m(0, 0) * m(2, 1) - m(2, 0) * m(0, 1)) + m(2, 3) * (m(1, 0) * m(0, 1) - m(0, 0) * m(1, 1)) );
    inv(3, 3) = d * ( m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) + m(1, 0) * (m(2, 1) * m(0, 2) - m(0, 1) * m(2, 2)) + m(2, 0) * (m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2)) );

    return true;
}

//! Computes the inverse of the specified sparse 4x4 matrix 'm'.
template <typename T> bool Inverse(SparseMatrix4T<T>& inv, const SparseMatrix4T<T>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv(0, 0) = d * ( m(1, 1) * m(2, 2) + m(2, 1) * ( -m(1, 2) ) );
    inv(1, 0) = d * ( m(1, 2) * m(2, 0) + m(2, 2) * ( -m(1, 0) ) );
    inv(2, 0) = d * ( m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1) );
  /*inv(3, 0) = 0;*/
    inv(0, 1) = d * ( m(2, 1) * m(0, 2) + m(0, 1) * ( -m(2, 2) ) );
    inv(1, 1) = d * ( m(2, 2) * m(0, 0) + m(0, 2) * ( -m(2, 0) ) );
    inv(2, 1) = d * ( m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1) );
  /*inv(3, 1) = 0;*/
    inv(0, 2) = d * ( m(0, 1) * m(1, 2) + m(1, 1) * ( -m(0, 2) ) );
    inv(1, 2) = d * ( m(0, 2) * m(1, 0) + m(1, 2) * ( -m(0, 0) ) );
    inv(2, 2) = d * ( m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) );
  /*inv(3, 2) = 0;*/
    inv(0, 3) = d * ( m(0, 1) * (m(2, 2) * m(1, 3) - m(1, 2) * m(2, 3)) + m(1, 1) * (m(0, 2) * m(2, 3) - m(2, 2) * m(0, 3)) + m(2, 1) * (m(1, 2) * m(0, 3) - m(0, 2) * m(1, 3)) );
    inv(1, 3) = d * ( m(0, 2) * (m(2, 0) * m(1, 3) - m(1, 0) * m(2, 3)) + m(1, 2) * (m(0, 0) * m(2, 3) - m(2, 0) * m(0, 3)) + m(2, 2) * (m(1, 0) * m(0, 3) - m(0, 0) * m(1, 3)) );
    inv(2, 3) = d * ( m(0, 3) * (m(2, 0) * m(1, 1) - m(1, 0) * m(2, 1)) + m(1, 3) * (m(0, 0) * m(2, 1) - m(2, 0) * m(0, 1)) + m(2, 3) * (m(1, 0) * m(0, 1) - m(0, 0) * m(1, 1)) );
  /*inv(3, 3) = 1;*/

    return true;
}


/* --- Global Operators --- */

//! \brief Multiplies the NxN matrix with the N-dimensional vector.
template <template <typename> class Vec, template <typename, std::size_t, std::size_t> class Mat, typename T, std::size_t N>
Vec<T> operator * (const Mat<T, N, N>& mat, const Vec<T>& vec)
{
    static_assert(Vec<T>::components == N, __GS_FILE_LINE__ "function only allows multiplication of an NxN matrix with an N-dimensional vector");

    Vec<T> result;

    for (std::size_t r = 0; r < N; ++r)
    {
        result[r] = 0;
        for (std::size_t c = 0; c < N; ++c)
            result[r] += mat(r, c)*vec[c];
    }

    return result;
}


} // /namespace Gs


#endif



// ================================================================================
