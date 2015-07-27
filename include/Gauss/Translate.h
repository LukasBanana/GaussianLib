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


/* --- Details --- */

namespace Details
{


template <typename M, typename T> void Translate(M& m, const Vector3T<T>& v)
{
    static_assert(
        M::rows >= 3 && M::columns >= 4,
        "matrices can only be translated, if they have at least 3 rows and 4 columns"
    );
    m(0, 3) += ( m(0, 0)*v.x + m(0, 1)*v.y + m(0, 2)*v.z );
    m(1, 3) += ( m(1, 0)*v.x + m(1, 1)*v.y + m(1, 2)*v.z );
    m(2, 3) += ( m(2, 0)*v.x + m(2, 1)*v.y + m(2, 2)*v.z );
}


} // /namespace Details


/* --- Global Functions --- */

//! Translates the specified matrix 'm' by the vector 'v'.
template <typename T> void Translate(Matrix<T, 4, 4>& m, const Vector3T<T>& v)
{
    Details::Translate(m, v);
}

//! Computes the inverse of the specified sparse 4x4 matrix 'm'.
template <typename T> void Translate(SparseMatrix4T<T>& m, const Vector3T<T>& v)
{
    Details::Translate(m, v);
}


} // /namespace Gs


#endif



// ================================================================================
