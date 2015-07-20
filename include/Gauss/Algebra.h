/*
 * Algebra.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ALGEBRA_H__
#define __GS_ALGEBRA_H__


#include "Vector3.h"


namespace Gs
{


/* --- Global Functions --- */

//! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
template < template <typename> class Vec, typename T > T Angle(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos( Dot(lhs, rhs) / (lhs.Length()*rhs.Length()) );
}

//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template < template <typename> class Vec, typename T > T AngleNorm(const Vec<T>& lhs, const Vec<T>& rhs)
{
    return std::acos(Dot(lhs, rhs));
}

//! Returns the squared length of the specified vector.
template < template <typename> class Vec, typename T > T LengthSq(const Vec<T>& vec)
{
    return Dot(vec, vec);
}

//! Returns the length (euclidian norm) of the specified vector.
template < template <typename> class Vec, typename T > T Length(const Vec<T>& vec)
{
    return std::sqrt(LengthSq(vec));
}

//! Returns the squared distance between the two vectors 'lhs' and 'rhs'.
template < template <typename> class Vec, typename T > T DistanceSq(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return LengthSq(result);
}

//! Returns the distance between the two vectors 'lhs' and 'rhs'.
template < template <typename> class Vec, typename T > T Distance(const Vec<T>& lhs, const Vec<T>& rhs)
{
    auto result = rhs;
    result -= lhs;
    return Length(result);
}

//! Returns the dot or rather scalar product between the two vectors 'lhs' and 'rhs'.
template < template <typename> class Vec, typename T > T Dot(const Vec<T>& lhs, const Vec<T>& rhs)
{
    T result = T(0);
    
    for (size_t i = 0; i < Vec<T>::components; ++i)
        result += lhs[i]*rhs[i];

    return result;
}

//! Normalizes the specified vector to the unit length of 1.
template < template <typename> class Vec, typename T > void Normalize(Vec<T>& vec)
{
    auto len = vec.LengthSq();
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


} // /namespace Gs


#endif



// ================================================================================
