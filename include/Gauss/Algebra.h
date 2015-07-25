/*
 * Algebra.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ALGEBRA_H__
#define __GS_ALGEBRA_H__


#include "Macros.h"

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

// Forward declaration
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
T Determinant(const M<T, Rows, Cols>&);

namespace Details
{

template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
class DeterminantHelper
{

    protected:
        
        friend T Gs::Determinant<M, T, Cols, Rows>(const M<T, Rows, Cols>&);

        static T OrderedDeterminant(const std::vector<T>& mat, std::size_t order)
        {
            if (order == 1)
                return mat[0];

            std::vector<T> minor((order - 1)*(order - 1));

            T det = T(0);

            for (std::size_t i = 0; i < order; ++i)
            {
                GetMinorMatrix(mat, minor, i, order);
                if (i % 2 == 1)
                    det -= mat[i] * OrderedDeterminant(minor, order - 1);
                else
                    det += mat[i] * OrderedDeterminant(minor, order - 1);
            }
    
            return det;
        }

    private:

        static void GetMinorMatrix(const std::vector<T>& mat, std::vector<T>& minor, std::size_t column, std::size_t order)
        {
            for (std::size_t r = 1; r < order; ++r)
            {
                for (std::size_t c = 0, i = 0; c < order; ++c)
                {
                    if (c != column)
                    {
                        minor[(r - 1)*(order - 1) + i] = mat[r*order + c];
                        ++i;
                    }
                }
            }
        }

};

} // /namespace Details

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
    std::vector<T> mat(Rows*Cols);
    for (std::size_t i = 0; i < Rows*Cols; ++i)
        mat[i] = m[i];
    return Details::DeterminantHelper<M, T, Rows, Cols>::OrderedDeterminant(mat, Rows);
}

//! Computes the determinant of a 1x1 matrix.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 1, 1>& m)
{
    return m(0, 0);
}

//! Computes the determinant of a 2x2 matrix.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 2, 2>& m)
{
    return m(0, 0)*m(1, 1) - m(1, 0)*m(0, 1);
}

//! Computes the determinant of a 3x3 matrix.
template <template <typename, std::size_t, std::size_t> class M, typename T>
T Determinant(const M<T, 3, 3>& m)
{
    return
        ( m(0, 0) * m(1, 1) * m(2, 2) ) + ( m(1, 0) * m(2, 1) * m(0, 2) ) + ( m(2, 0) * m(0, 1) * m(1, 2) ) -
        ( m(0, 2) * m(1, 1) * m(2, 0) ) - ( m(0, 2) * m(2, 1) * m(0, 0) ) - ( m(2, 2) * m(0, 1) * m(1, 0) );
}

//! Computes the determinant of a 4x4 matrix.
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
