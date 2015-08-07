/*
 * Vector3.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_VECTOR3_H__
#define __GS_VECTOR3_H__


#include "Vector.h"
#include "Algebra.h"
#include "Swizzle.h"

#include <cmath>


namespace Gs
{


/**
Base 3D vector class with components: x, y, and z.
\tparam T Specifies the data type of the vector components.
This should be a primitive data type such as float, double, int etc.
*/
template <typename T>
class Vector<T, 3>
{
    
    public:
        
        //! Specifies the number of vector components.
        static const std::size_t components = 3;

        #ifndef GS_DISABLE_AUTO_INIT
        Vector() :
            x( T(0) ),
            y( T(0) ),
            z( T(0) )
        {
        }
        #else
        Vector() = default;
        #endif

        Vector(const Vector<T, 3>& rhs) :
            x( rhs.x ),
            y( rhs.y ),
            z( rhs.z )
        {
        }

        explicit Vector(const T& scalar) :
            x( scalar ),
            y( scalar ),
            z( scalar )
        {
        }

        Vector(const T& x, const T& y, const T& z) :
            x( x ),
            y( y ),
            z( z )
        {
        }

        Vector(UninitializeTag)
        {
            // do nothing
        }

        Vector<T, 3>& operator += (const Vector<T, 3>& rhs)
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Vector<T, 3>& operator -= (const Vector<T, 3>& rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        Vector<T, 3>& operator *= (const Vector<T, 3>& rhs)
        {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            return *this;
        }

        Vector<T, 3>& operator /= (const Vector<T, 3>& rhs)
        {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            return *this;
        }

        Vector<T, 3>& operator *= (const T& rhs)
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }

        Vector<T, 3>& operator /= (const T& rhs)
        {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }

        /**
        \brief Returns the specified vector component.
        \param[in] component Specifies the vector component index. This must be 0, 1, or 2.
        */
        T& operator [] (std::size_t component)
        {
            GS_ASSERT(component < (Vector<T, 3>::components));
            return *((&x) + component);
        }

        /**
        \brief Returns the specified vector component.
        \param[in] component Specifies the vector component index. This must be 0, 1, or 2.
        */
        const T& operator [] (std::size_t component) const
        {
            GS_ASSERT(component < (Vector<T, 3>::components));
            return *((&x) + component);
        }

        //! Returns the squared length of this vector.
        T LengthSq() const
        {
            return Gs::LengthSq(*this);
        }

        //! Returns the length of this vector.
        T Length() const
        {
            return Gs::Length(*this);
        }

        /**
        Normalizes this vector to the unit length of 1.
        \see Normalized
        \see Length
        */
        void Normalize()
        {
            Gs::Normalize(*this);
        }

        /**
        Returns a normalized instance of this vector.
        \see Normalize
        */
        Vector<T, 3> Normalized() const
        {
            auto vec = *this;
            vec.Normalize();
            return vec;
        }

        /**
        Resizes this vector to the specified length.
        \see Normalize
        \see Length
        */
        void Resize(const T& length)
        {
            Gs::Resize(*this, length);
        }

        /**
        Returns a type casted instance of this vector.
        \tparam C Specifies the static cast type.
        */
        template <typename C>
        Vector<C, 3> Cast() const
        {
            return Vector<C, 3>(
                static_cast<C>(x),
                static_cast<C>(y),
                static_cast<C>(z)
            );
        }

        //! Returns a pointer to the first element of this vector.
        T* Ptr()
        {
            return &x;
        }

        //! Returns a constant pointer to the first element of this vector.
        const T* Ptr() const
        {
            return &x;
        }

        #ifdef GS_ENABLE_SWIZZLE_OPERATOR
        #   include "SwizzleVec2Op2.h"
        #   include "SwizzleVec2Op3.h"
        #   include "SwizzleVec2Op4.h"
        #   include "SwizzleVec3Op2.h"
        #   include "SwizzleVec3Op3.h"
        #   include "SwizzleVec3Op4.h"
        #endif

        T x, y, z;

};


/* --- Type Alias --- */

template <typename T> using Vector3T = Vector<T, 3>;

using Vector3   = Vector3T<Real>;
using Vector3f  = Vector3T<float>;
using Vector3d  = Vector3T<double>;
using Vector3i  = Vector3T<int>;
using Vector3ui = Vector3T<unsigned int>;
using Vector3b  = Vector3T<char>;
using Vector3ub = Vector3T<unsigned char>;


} // /namespace Gs


#endif



// ================================================================================
