/*
 * Algebra_SSE.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ALGEBRA_SSE_H
#define GS_ALGEBRA_SSE_H


#include "Algebra.h"

#include <pmmintrin.h>


namespace Gs
{


/* --- Global Functions --- */

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <>
float Dot(const Vector<float, 4>& lhs, const Vector<float, 4>& rhs)
{
    union
    {
        float  f[4];
        __m128 v;
    };
    v = _mm_mul_ps(lhs.m128, rhs.m128);
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
    return f[0];
}

//! Returns the per-component reciprocal of the specified N-dimensional vector.
template <>
Vector<float, 4> Rcp(const Vector<float, 4>& vec)
{
    return _mm_rcp_ps(vec.m128);
}


} // /namespace Gs


#endif



// ================================================================================
