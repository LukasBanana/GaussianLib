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

namespace Details
{


/* --- Internal functions --- */

inline float DotProductM128f(__m128 lhs, __m128 rhs)
{
    __m128 v;
    v = _mm_mul_ps(lhs, rhs);
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
    return _mm_cvtss_f32(v);
}


} // /namespace Details


/* --- Global Functions --- */

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template <>
inline float Dot(const Vector<float, 4>& lhs, const Vector<float, 4>& rhs)
{
    return Details::DotProductM128f(lhs.m128, rhs.m128);
}

//! Returns the dot or rather scalar product between the two quaternions 'lhs' and 'rhs'.
template <>
inline float Dot(const QuaternionT<float>& lhs, const QuaternionT<float>& rhs)
{
    return Details::DotProductM128f(lhs.m128, rhs.m128);
}

//! Returns the per-component reciprocal of the specified N-dimensional vector.
template <>
inline Vector<float, 4> Rcp(const Vector<float, 4>& vec)
{
    return _mm_rcp_ps(vec.m128);
}

/**
\brief Mixes the two values by their scalings.
\return Equivalent to: v0*scale0 + v1*scale1
*/
template <>
inline QuaternionT<float> Mix(const QuaternionT<float>& v0, const QuaternionT<float>& v1, const float& scale0, const float& scale1)
{
    return QuaternionT<float>
    {
        _mm_add_ps(
            _mm_mul_ps(_mm_set_ps1(scale0), v0.m128), 
            _mm_mul_ps(_mm_set_ps1(scale1), v1.m128)
        )
    };
}


} // /namespace Gs


#endif



// ================================================================================
