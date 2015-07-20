/*
 * Equal.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_EQUAL_H__
#define __GS_EQUAL_H__


#include <cmath>


namespace Gs
{


/* --- Global Functions --- */

template <typename T> inline bool IsEqual(const T& lhs, const T& rhs)
{
    return lhs == rhs;
}

template <> inline bool IsEqual<float>(const float& lhs, const float& rhs)
{
    return std::abs(lhs - rhs) < epsilon32;
}

template <> inline bool IsEqual<double>(const double& lhs, const double& rhs)
{
    return std::abs(lhs - rhs) < epsilon64;
}

template <typename T> inline bool IsZero(const T& value)
{
    return value == T(0);
}

template <> inline bool IsZero<float>(const float& value)
{
    return std::abs(value) < epsilon32;
}

template <> inline bool IsZero<double>(const double& value)
{
    return std::abs(value) < epsilon64;
}

template <typename T> inline bool IsOne(const T& value)
{
    return IsEqual(value, T(1));
}


} // /namespace Gs


#endif



// ================================================================================
