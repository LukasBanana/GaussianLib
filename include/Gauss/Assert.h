/*
 * Assert.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ASSERT_H__
#define __GS_ASSERT_H__


#include "Config.h"

#include <cassert>


#ifdef GS_ENABLE_ASSERT
#   define GS_ASSERT(expr) assert((expr))
#else
#   define GS_ASSERT(expr)
#endif


#endif



// ================================================================================
