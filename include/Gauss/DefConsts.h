/*
 * DefConsts.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_DEF_CONSTS_H__
#define __GS_DEF_CONSTS_H__


#include "Epsilon.h"


namespace Gs
{


template <>
const float Epsilon<float>::value = 1.0e-6f;//0.000001f;

template <>
const double Epsilon<double>::value = 1.0e-8;//0.00000001;


} // /namespace Gs


#endif



// ================================================================================
