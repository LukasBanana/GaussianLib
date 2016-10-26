/*
 * Assert.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ASSERT_H
#define GS_ASSERT_H


#include "Config.h"

#include <cassert>


#define _GS_TOSTRING_PRIMARY_(x) #x
#define _GS_TOSTRING_(x) _GS_TOSTRING_PRIMARY_(x)

#ifdef GS_ENABLE_ASSERT
#   ifdef GS_ASSERT_EXCEPTION
#       include <exception>
#       define GS_ASSERT(expr)                                  \
            if (!(expr))                                        \
            {                                                   \
                throw std::runtime_error(                       \
                    "assertion failed: (" #expr "), file "      \
                    __FILE__ ", line " _GS_TOSTRING_(__LINE__)  \
                );                                              \
            }
#   else
#       define GS_ASSERT(expr) assert((expr))
#   endif
#else
#   define GS_ASSERT(expr)
#endif


#endif



// ================================================================================
