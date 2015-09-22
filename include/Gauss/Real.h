/*
 * Real.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_REAL_H__
#define __GS_REAL_H__


#include <cmath>


namespace Gs
{


/* --- Constants --- */

#ifdef GS_HIGH_PRECISION_FLOAT

using Real = double;

#else

using Real = float;

#endif

static const Real pi = Real(3.14159265358979323846);

//! Converts radian to degree.
template <typename T>
T Rad2Deg(const T& x)
{
    return x*T(180)/T(pi);
}

//! Converts degree to radian.
template <typename T>
T Deg2Rad(const T& x)
{
    return x*T(pi)/T(180);
}


} // /namespace Gs


#endif



// ================================================================================
