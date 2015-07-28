/*
 * Translate.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_TRANSLATE_H__
#define __GS_TRANSLATE_H__


#include "Decl.h"


namespace Gs
{


namespace Details
{

#define __GS_ASSERT_MxN_MATRIX__(info, T, m, n)             \
    static_assert(                                          \
        T::rows >= m && T::columns >= n,                    \
        info " requires at least a " #m "x" #n " matrix"    \
    )

template <typename M, typename T>
void Translate(M& m, const Vector3T<T>& v)
{
    #ifdef GS_ROW_VECTORS
    __GS_ASSERT_MxN_MATRIX__("translation with row vectors", M, 4, 3);
    m(3, 0) += ( m(0, 0)*v.x + m(1, 0)*v.y + m(2, 0)*v.z );
    m(3, 1) += ( m(0, 1)*v.x + m(1, 1)*v.y + m(2, 1)*v.z );
    m(3, 2) += ( m(0, 2)*v.x + m(1, 2)*v.y + m(2, 2)*v.z );
    #else
    __GS_ASSERT_MxN_MATRIX__("translation with column vectors", M, 3, 4);
    m(0, 3) += ( m(0, 0)*v.x + m(0, 1)*v.y + m(0, 2)*v.z );
    m(1, 3) += ( m(1, 0)*v.x + m(1, 1)*v.y + m(1, 2)*v.z );
    m(2, 3) += ( m(2, 0)*v.x + m(2, 1)*v.y + m(2, 2)*v.z );
    #endif
}


} // /namespace Details


//! Translates the specified matrix 'm' by the vector 'v'.
template <typename T>
void Translate(Matrix<T, 4, 4>& m, const Vector3T<T>& v)
{
    /* Translate x, y, z */
    Details::Translate(m, v);

    /* Also translate w */
    #ifdef GS_ROW_VECTORS
    m(3, 3) += ( m(0, 3)*v.x + m(1, 3)*v.y + m(2, 3)*v.z );
    #else
    m(3, 3) += ( m(3, 0)*v.x + m(3, 1)*v.y + m(3, 2)*v.z );
    #endif
}

//! Computes the inverse of the specified sparse 4x4 matrix 'm'.
template <typename T>
void Translate(SparseMatrix4T<T>& m, const Vector3T<T>& v)
{
    Details::Translate(m, v);
}


} // /namespace Gs


#endif



// ================================================================================
