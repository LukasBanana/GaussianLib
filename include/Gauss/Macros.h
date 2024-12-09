/*
 * Macros.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_MACROS_H
#define GS_MACROS_H


#include <Gauss/Config.h>


#define GS_TOSTRING_PRIMARY(x)  #x
#define GS_TOSTRING(x)          GS_TOSTRING_PRIMARY(x)
#define GS_FILE_LINE            __FILE__ " (" GS_TOSTRING(__LINE__) "): "

#ifdef GS_ROW_VECTORS

#define GS_ASSERT_MxN_MATRIX(info, T, n, m)                 \
    static_assert(                                          \
        T::rows >= m && T::columns >= n,                    \
        info " requires at least a " #m "x" #n " matrix"    \
    )

#else

#define GS_ASSERT_MxN_MATRIX(info, T, m, n)                 \
    static_assert(                                          \
        T::rows >= m && T::columns >= n,                    \
        info " requires at least a " #m "x" #n " matrix"    \
    )

#endif

#if __cplusplus >= 201703L // C++17
#   define GS_NODISCARD [[nodiscard]]
#elif defined __GNUC__ || defined __clang__ // GNU/Clang extensions
#   define GS_NODISCARD __attribute__((warn_unused_result))
#elif _MSC_VER >= 1700 // MSVC extensions
#   define GS_NODISCARD _Check_return_
#else
#   define GS_NODISCARD
#endif


#endif



// ================================================================================
