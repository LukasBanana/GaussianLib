/*
 * Equals.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_EQUALS_H__
#define __GS_EQUALS_H__


#include "Real.h"
#include "Epsilon.h"

#include <cmath>


namespace Gs
{


/* --- Global Functions --- */

template <typename T>
bool Equals(const T& lhs, const T& rhs)
{
    return lhs == rhs;
}

template <>
bool Equals<float>(const float& lhs, const float& rhs)
{
    return std::abs(lhs - rhs) <= Epsilon<float>();
}

template <>
bool Equals<double>(const double& lhs, const double& rhs)
{
    return std::abs(lhs - rhs) <= Epsilon<double>();
}

template <typename T, std::size_t N>
bool Equals(const Vector<T, N>& lhs, const Vector<T, N>& rhs)
{
    for (std::size_t i = 0; i < N; ++i)
    {
        if (!Equals(lhs[i], rhs[i]))
            return false;
    }
    return true;
}


} // /namespace Gs


#endif



// ================================================================================
