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
#include "SparseMatrix4.h"

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>


namespace Gs
{


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
\return min{}
*/
template <typename T> T Clamp(const T& x, const T& minima, const T& maxima)
{
    if (x <= minima)
        return minima;
    if (x >= maxima)
        return maxima;
    return x;
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
