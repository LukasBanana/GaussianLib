/*
 * Swizzle.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_SWIZZLE_H__
#define __GS_SWIZZLE_H__


#include "Config.h"


#ifdef GS_ENABLE_SWIZZLE_OPERATOR

#define __GS_SWIZZLE_VECTOR_OP_SUB__(S, V, OP)                                  \
    template <typename T> V<T> operator OP (const S<T>& lhs, const S<T>& rhs)   \
    {                                                                           \
        return V<T>(lhs) + V<T>(rhs);                                           \
    }                                                                           \
    template <typename T> V<T> operator OP (const V<T>& lhs, const S<T>& rhs)   \
    {                                                                           \
        return lhs + V<T>(rhs);                                                 \
    }                                                                           \
    template <typename T> V<T> operator OP (const S<T>& lhs, const V<T>& rhs)   \
    {                                                                           \
        return V<T>(lhs) + rhs;                                                 \
    }

#define __GS_SWIZZLE_VECTOR_OP__(N, OP) \
    __GS_SWIZZLE_VECTOR_OP_SUB__(SwizzleRef##N, Vector##N##T, OP)

#define __GS_SWIZZLE_VECTOR_OP_ALL__(N) \
    __GS_SWIZZLE_VECTOR_OP__(N, +)      \
    __GS_SWIZZLE_VECTOR_OP__(N, -)      \
    __GS_SWIZZLE_VECTOR_OP__(N, *)      \
    __GS_SWIZZLE_VECTOR_OP__(N, /)

#define __GS_SWIZZLE_INTERFACE__(N, OP_DECL)                    \
    OP_DECL(=);                                                 \
    OP_DECL(+=);                                                \
    OP_DECL(-=);                                                \
    OP_DECL(*=);                                                \
    OP_DECL(/=);                                                \
    SwizzleRef##N<T>& operator = (const Vector##N##T<T>& rhs);  \
    operator Vector##N##T<T> () const;


#endif


#endif



// ================================================================================
