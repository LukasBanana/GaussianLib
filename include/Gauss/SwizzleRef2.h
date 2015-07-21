/*
 * SwizzleRef2.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_SWIZZLE_REF2_H__
#define __GS_SWIZZLE_REF2_H__


#include "Swizzle.h"


namespace Gs
{


template <typename T> class Vector2T;

//! Helper class to support the 'swizzle'-lile operator with 2 components.
template <typename T> class SwizzleRef2
{
            
    public:

        SwizzleRef2(T& x, T& y) :
            x_( x ),
            y_( y )
        {
        }

        #define __GS_SWIZZLE_OP__(OP)                               \
            SwizzleRef2<T>& operator OP (const SwizzleRef2<T>& rhs) \
            {                                                       \
                x_ OP rhs.x_;                                       \
                y_ OP rhs.y_;                                       \
                return *this;                                       \
            }

        __GS_SWIZZLE_INTERFACE__(2, __GS_SWIZZLE_OP__);

        #undef __GS_SWIZZLE_OP__

    private:

        T& x_, &y_;

};


#define __GS_SWIZZLE_REF2__(v0, v1)             \
    SwizzleRef2<T> v0##v1()                     \
    {                                           \
        return SwizzleRef2<T>(v0, v1);          \
    }                                           \
    SwizzleRef2<const T> v0##v1() const         \
    {                                           \
        return SwizzleRef2<const T>(v0, v1);    \
    }


} // /namespace Gs


#endif



// ================================================================================
