/*
 * Vector4_SSE2.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_VECTOR4_SSE2_H
#define GS_VECTOR4_SSE2_H


#include "Vector4.h"

#include <emmintrin.h>


namespace Gs
{


//! Template specialization with SSE support for 4D single-precision floating-point vectors.
template <>
class alignas(32) Vector<double, 4>
{
    
    public:
        
        using T = double;

        //! Specifies the number of vector components.
        static const std::size_t components = 4;

        #ifndef GS_DISABLE_AUTO_INIT
        Vector() :
            m128 { _mm_setzero_pd(), _mm_setzero_pd() }
        {
        }
        #else
        Vector() = default;
        #endif
        
        Vector(__m128d rhs0, __m128d rhs1) :
            m128 { rhs0, rhs1 }
        {
        }

        Vector(const Vector<T, 4>& rhs) :
            m128 { rhs.m128[0], rhs.m128[1] }
        {
        }

        explicit Vector(const Vector<T, 2>& xy, const Vector<T, 2>& zw) :
            x { xy.x },
            y { xy.y },
            z { zw.x },
            w { zw.y }
        {
        }

        explicit Vector(const Vector<T, 2>& xy, const T& z, const T& w) :
            x { xy.x },
            y { xy.y },
            z { z    },
            w { w    }
        {
        }

        explicit Vector(const Vector<T, 3>& xyz, const T& w) :
            x { xyz.x },
            y { xyz.y },
            z { xyz.z },
            w { w     }
        {
        }

        explicit Vector(const T& scalar) :
            x { scalar },
            y { scalar },
            z { scalar },
            w { scalar }
        {
        }

        Vector(const T& x, const T& y, const T& z, const T& w) :
            m128 { { x, y }, { z, w } }
        {
        }

        Vector(UninitializeTag)
        {
            // do nothing
        }

        Vector<T, 4>& operator += (const Vector<T, 4>& rhs)
        {
            m128[0] = _mm_add_pd(m128[0], rhs.m128[0]);
            m128[1] = _mm_add_pd(m128[1], rhs.m128[1]);
            return *this;
        }

        Vector<T, 4>& operator -= (const Vector<T, 4>& rhs)
        {
            m128[0] = _mm_sub_pd(m128[0], rhs.m128[0]);
            m128[1] = _mm_sub_pd(m128[1], rhs.m128[1]);
            return *this;
        }

        Vector<T, 4>& operator *= (const Vector<T, 4>& rhs)
        {
            m128[0] = _mm_mul_pd(m128[0], rhs.m128[0]);
            m128[1] = _mm_mul_pd(m128[1], rhs.m128[1]);
            return *this;
        }

        Vector<T, 4>& operator /= (const Vector<T, 4>& rhs)
        {
            m128[0] = _mm_div_pd(m128[0], rhs.m128[0]);
            m128[1] = _mm_div_pd(m128[1], rhs.m128[1]);
            return *this;
        }

        Vector<T, 4>& operator *= (const T rhs)
        {
            m128[0] = _mm_mul_pd(m128[0], _mm_set1_pd(rhs));
            m128[1] = _mm_mul_pd(m128[1], _mm_set1_pd(rhs));
            return *this;
        }

        Vector<T, 4>& operator /= (const T rhs)
        {
            m128[0] = _mm_div_pd(m128[0], _mm_set1_pd(rhs));
            m128[1] = _mm_div_pd(m128[1], _mm_set1_pd(rhs));
            return *this;
        }

        Vector<T, 4> operator - () const
        {
            return Vector<T, 4>(-x, -y, -z, -w);
        }

        /**
        \brief Returns the specified vector component.
        \param[in] component Specifies the vector component index. This must be 0, 1, 2, or 3.
        */
        T& operator [] (std::size_t component)
        {
            GS_ASSERT(component < (Vector<T, 4>::components));
            return *((&x) + component);
        }

        /**
        \brief Returns the specified vector component.
        \param[in] component Specifies the vector component index. This must be 0, 1, 2, or 3.
        */
        const T& operator [] (std::size_t component) const
        {
            GS_ASSERT(component < (Vector<T, 4>::components));
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
        Vector<T, 4> Normalized() const
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
        void Resize(T length)
        {
            Gs::Resize(*this, length);
        }

        /**
        Returns a type casted instance of this vector.
        \tparam C Specifies the static cast type.
        */
        template <typename C>
        Vector<C, 4> Cast() const
        {
            return Vector<C, 4>(
                static_cast<C>(x),
                static_cast<C>(y),
                static_cast<C>(z),
                static_cast<C>(w)
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
        #   include "SwizzleVec4Op2.h"
        #   include "SwizzleVec4Op3.h"
        #   include "SwizzleVec4Op4.h"
        #endif
        
        union
        {
            struct
            {
                T x, y, z, w;
            };
            __m128d m128[2];
        };

};


/* --- Global Operators --- */

template <>
inline Vector<double, 4> operator + (const Vector<double, 4>& lhs, const Vector<double, 4>& rhs)
{
    return Vector<double, 4>(_mm_add_pd(lhs.m128[0], rhs.m128[0]), _mm_add_pd(lhs.m128[1], rhs.m128[1]));
}

template <>
inline Vector<double, 4> operator - (const Vector<double, 4>& lhs, const Vector<double, 4>& rhs)
{
    return Vector<double, 4>(_mm_sub_pd(lhs.m128[0], rhs.m128[0]), _mm_sub_pd(lhs.m128[1], rhs.m128[1]));
}

template <>
inline Vector<double, 4> operator * (const Vector<double, 4>& lhs, const Vector<double, 4>& rhs)
{
    return Vector<double, 4>(_mm_mul_pd(lhs.m128[0], rhs.m128[0]), _mm_mul_pd(lhs.m128[1], rhs.m128[1]));
}

template <>
inline Vector<double, 4> operator / (const Vector<double, 4>& lhs, const Vector<double, 4>& rhs)
{
    return Vector<double, 4>(_mm_div_pd(lhs.m128[0], rhs.m128[0]), _mm_div_pd(lhs.m128[1], rhs.m128[1]));
}

template <>
inline Vector<double, 4> operator * (const Vector<double, 4>& lhs, const double& rhs)
{
    const auto rhs128 = _mm_set1_pd(rhs);
    return Vector<double, 4>(_mm_mul_pd(lhs.m128[0], rhs128), _mm_mul_pd(lhs.m128[1], rhs128));
}

template <>
inline Vector<double, 4> operator * (const double& lhs, const Vector<double, 4>& rhs)
{
    const auto lhs128 = _mm_set1_pd(lhs);
    return Vector<double, 4>(_mm_mul_pd(lhs128, rhs.m128[0]), _mm_mul_pd(lhs128, rhs.m128[1]));
}

template <>
inline Vector<double, 4> operator / (const Vector<double, 4>& lhs, const double& rhs)
{
    const auto rhs128 = _mm_set1_pd(rhs);
    return Vector<double, 4>(_mm_div_pd(lhs.m128[0], rhs128), _mm_div_pd(lhs.m128[1], rhs128));
}

template <>
inline Vector<double, 4> operator / (const double& lhs, const Vector<double, 4>& rhs)
{
    const auto lhs128 = _mm_set1_pd(lhs);
    return Vector<double, 4>(_mm_div_pd(lhs128, rhs.m128[0]), _mm_mul_pd(lhs128, rhs.m128[1]));
}


} // /namespace Gs


#endif



// ================================================================================
