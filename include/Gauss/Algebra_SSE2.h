/*
 * Algebra_SSE2.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ALGEBRA_SSE2_H
#define GS_ALGEBRA_SSE2_H


#include "Algebra.h"

#include <emmintrin.h>


namespace Gs
{


/* --- Global Functions --- */

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <>
inline double Dot(const Vector<double, 4>& lhs, const Vector<double, 4>& rhs)
{
    __m128d v[2];
    v[0] = _mm_mul_pd(lhs.m128[0], rhs.m128[0]);
    v[0] = _mm_hadd_pd(v[0], v[0]);
    v[1] = _mm_mul_pd(lhs.m128[1], rhs.m128[1]);
    v[1] = _mm_hadd_pd(v[1], v[1]);
    return _mm_cvtsd_f64(v[0]) + _mm_cvtsd_f64(v[1]);
}


} // /namespace Gs


#endif



// ================================================================================
