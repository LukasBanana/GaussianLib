/*
 * SwizzleRef3.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_SWIZZLE_REF3_H__
#define __GS_SWIZZLE_REF3_H__


#include "Swizzle.h"


namespace Gs
{


template <typename T> class Vector3T;

//! Helper class to support the 'swizzle'-lile operator with 3 components.
template <typename T> class SwizzleRef3
{
            
    public:

        SwizzleRef3(T& x, T& y, T&z) :
            x_( x ),
            y_( y ),
            z_( z )
        {
        }

        #define __GS_SWIZZLE_OP__(OP)                               \
            SwizzleRef3<T>& operator OP (const SwizzleRef3<T>& rhs) \
            {                                                       \
                x_ OP rhs.x_;                                       \
                y_ OP rhs.y_;                                       \
                z_ OP rhs.z_;                                       \
                return *this;                                       \
            }

        __GS_SWIZZLE_INTERFACE__(3, __GS_SWIZZLE_OP__);

        #undef __GS_SWIZZLE_OP__

    private:

        T& x_, &y_, &z_;

};


#define __GS_SWIZZLE_REF3__(v0, v1, v2)             \
    SwizzleRef3<T> v0##v1##v2()                     \
    {                                               \
        return SwizzleRef3<T>(v0, v1, v2);          \
    }                                               \
    SwizzleRef3<const T> v0##v1##v2() const         \
    {                                               \
        return SwizzleRef3<const T>(v0, v1, v2);    \
    }


} // /namespace Gs


#endif



// ================================================================================
