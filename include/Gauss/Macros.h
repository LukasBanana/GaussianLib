/*
 * Macros.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_MACROS_H
#define GS_MACROS_H


#include "Config.h"


#define __GS_TOSTRING_SUB__(x)  #x
#define __GS_TOSTRING__(x)      __GS_TOSTRING_SUB__(x)
#define __GS_FILE_LINE__        __FILE__ " (" __GS_TOSTRING__(__LINE__) "): "

#ifdef GS_ROW_VECTORS

#define __GS_ASSERT_MxN_MATRIX__(info, T, n, m)             \
    static_assert(                                          \
        T::rows >= m && T::columns >= n,                    \
        info " requires at least a " #m "x" #n " matrix"    \
    )

#else

#define __GS_ASSERT_MxN_MATRIX__(info, T, m, n)             \
    static_assert(                                          \
        T::rows >= m && T::columns >= n,                    \
        info " requires at least a " #m "x" #n " matrix"    \
    )

#endif


#endif



// ================================================================================
