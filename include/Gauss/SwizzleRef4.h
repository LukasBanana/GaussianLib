/*
 * SwizzleRef4.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_SWIZZLE_REF4_H__
#define __GS_SWIZZLE_REF4_H__


#include "Swizzle.h"


namespace Gs
{


template <typename T> class Vector4T;

//! Helper class to support the 'swizzle'-lile operator with 4 components.
template <typename T> class SwizzleRef4
{
            
    public:

        SwizzleRef4(T& x, T& y, T&z, T& w) :
            x_( x ),
            y_( y ),
            z_( z ),
            w_( w )
        {
        }

        #define __GS_SWIZZLE_OP__(OP)                               \
            SwizzleRef4<T>& operator OP (const SwizzleRef4<T>& rhs) \
            {                                                       \
                x_ OP rhs.x_;                                       \
                y_ OP rhs.y_;                                       \
                z_ OP rhs.z_;                                       \
                w_ OP rhs.w_;                                       \
                return *this;                                       \
            }

        __GS_SWIZZLE_INTERFACE__(4, __GS_SWIZZLE_OP__);

        #undef __GS_SWIZZLE_OP__

    private:

        T& x_, &y_, &z_, &w_;

};


#define __GS_SWIZZLE_REF4__(v0, v1, v2, v3)             \
    SwizzleRef4<T> v0##v1##v2##v3()                     \
    {                                                   \
        return SwizzleRef4<T>(v0, v1, v2, v3);          \
    }                                                   \
    SwizzleRef4<const T> v0##v1##v2##v3() const         \
    {                                                   \
        return SwizzleRef4<const T>(v0, v1, v2, v3);    \
    }


} // /namespace Gs


#endif



// ================================================================================
