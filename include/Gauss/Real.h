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

static const float epsilon32 = 0.00001f;
static const double epsilon64 = 0.00000001;


#ifdef GD_HIGH_PRECISION_FLOAT

using Real = double;
static const Real epsilon = epsilon64;

#else

using Real = float;
static const Real epsilon = epsilon32;

#endif


} // /namespace Gs


#endif



// ================================================================================
