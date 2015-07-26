/*
 * Determinant.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_DETERMINANT_H__
#define __GS_DETERMINANT_H__


#include "Details.h"


namespace Gs
{


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


} // /namespace Gs


#endif



// ================================================================================
