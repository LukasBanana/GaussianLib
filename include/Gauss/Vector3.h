/*
 * Vector3.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_VECTOR3_H__
#define __GS_VECTOR3_H__


#include "Real.h"
#include "Assert.h"

#include <cmath>


namespace Gs
{


/**
Base 3D vector class with components: x, y, and z.
\tparam T Specifies the data type of the vector components.
This should be a primitive data type such as float, double, int etc.
*/
template <typename T> class Vector3T
{
    
    public:
        
        //! Specifies the number of vector components.
        static const size_t components = 3;

        Vector3T() :
            x{ T(0) },
            y{ T(0) },
            z{ T(0) }
        {
        }
        Vector3T(const Vector3T<T>& rhs) :
            x{ rhs.x },
            y{ rhs.y },
            z{ rhs.z }
        {
        }
        explicit Vector3T(const T& scalar) :
            x{ scalar },
            y{ scalar },
            z{ scalar }
        {
        }
        Vector3T(const T& x, const T& y, const T& z) :
            x{ x },
            y{ y },
            z{ z }
        {
        }

        Vector3T<T>& operator += (const Vector3T<T>& rhs)
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Vector3T<T>& operator -= (const Vector3T<T>& rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        Vector3T<T>& operator *= (const Vector3T<T>& rhs)
        {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            return *this;
        }

        Vector3T<T>& operator /= (const Vector3T<T>& rhs)
        {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            return *this;
        }

        Vector3T<T>& operator *= (const T& rhs)
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }

        Vector3T<T>& operator /= (const T& rhs)
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
        T& operator [] (size_t component)
        {
            GS_ASSERT(component < Vector3T<T>::components);
            return *((&x) + component);
        }

        /**
        \brief Returns the specified vector component.
        \param[in] component Specifies the vector component index. This must be 0, 1, or 2.
        */
        const T& operator [] (size_t component) const
        {
            GS_ASSERT(component < Vector3T<T>::components);
            return *((&x) + component);
        }

        //! Returns the squared length of this vector.
        T LengthSq() const;
        //! Returns the length of this vector.
        T Length() const;

        /**
        Normalizes the vector to the unit length of 1.
        \see Normalized
        \see Length
        */
        void Normalize();
        /**
        Returns a normalized instance of this vector.
        \see Normalize
        */
        Vector3T<T> Normalized() const;

        /**
        Returns a type casted instance of this vector.
        \tparam C Specifies the static cast type.
        */
        template <typename C> Vector3T<C> Cast() const
        {
            return Vector3T<C>(
                static_cast<C>(x),
                static_cast<C>(y),
                static_cast<C>(z)
            );
        }

        T x, y, z;

};


/* --- Global Functions --- */

//! Returns the cross or rather vector product between the two vectors 'lhs' and 'rhs'.
template <typename T> Vector3T<T> Cross(const Vector3T<T>& lhs, const Vector3T<T>& rhs)
{
    return Vector3T<T>(
        lhs.y*rhs.z - rhs.y*lhs.z,
        rhs.x*lhs.z - lhs.x*rhs.z,
        lhs.x*rhs.y - rhs.x*lhs.y
    );
}


/* --- Global Operators --- */

template <typename T> Vector3T<T> operator + (const Vector3T<T>& lhs, const Vector3T<T>& rhs)
{
    auto result = lhs;
    result += rhs;
    return result;
}

template <typename T> Vector3T<T> operator - (const Vector3T<T>& lhs, const Vector3T<T>& rhs)
{
    auto result = lhs;
    result -= rhs;
    return result;
}

template <typename T> Vector3T<T> operator * (const Vector3T<T>& lhs, const Vector3T<T>& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T> Vector3T<T> operator / (const Vector3T<T>& lhs, const Vector3T<T>& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T> Vector3T<T> operator * (const Vector3T<T>& lhs, const T& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T> Vector3T<T> operator * (const T& lhs, const Vector3T<T>& rhs)
{
    auto result = rhs;
    result *= lhs;
    return result;
}

template <typename T> Vector3T<T> operator / (const Vector3T<T>& lhs, const T& rhs)
{
    auto result = lhs;
    result /= rhs;
    return result;
}


/* --- Class Functions --- */

template <typename T> T Vector3T<T>::LengthSq() const
{
    return Gs::LengthSq(*this);
}

template <typename T> T Vector3T<T>::Length() const
{
    return Gs::Length(*this);
}

template <typename T> void Vector3T<T>::Normalize()
{
    Gs::Normalize(*this);
}

template <typename T> Vector3T<T> Vector3T<T>::Normalized() const
{
    auto vec = *this;
    vec.Normalize();
    return vec;
}


/* --- Type Alias --- */

using Vector3 = Vector3T<Real>;
using Vector3f = Vector3T<float>;
using Vector3d = Vector3T<double>;
using Vector3i = Vector3T<int>;


} // /namespace Gs


#endif



// ================================================================================
