/*
 * Quaternion_SSE.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_QUATERNION_SSE_H
#define GS_QUATERNION_SSE_H


#include "Quaternion.h"
#include "Vector4_SSE.h"

#include <xmmintrin.h>


namespace Gs
{


//! Template specialization with SSE support for 4D single-precision floating-point quaternions.
template <>
class alignas(alignof(__m128)) QuaternionT<float>
{
    
        using T = float;

    public:
        
        //! Specifies the typename of the scalar components.
        using ScalarType = T;

        //! Specifies the number of quaternion components. This is just for the internal template interface.
        static const std::size_t components = 4;

        #ifndef GS_DISABLE_AUTO_INIT
        QuaternionT() :
            m128 ( _mm_setzero_ps() )
        {
        }
        #else
        QuaternionT() = default;
        #endif

        QuaternionT(__m128 rhs) :
            m128 ( rhs )
        {
        }

        QuaternionT(const QuaternionT<T>& rhs) :
            m128 ( rhs.m128 )
        {
        }

        QuaternionT(const T& x, const T& y, const T& z, const T& w) :
            m128 ( _mm_set_ps(w, z, y, x) )
        {
        }

        template <template <typename, std::size_t, std::size_t> class M, std::size_t Rows, std::size_t Cols>
        explicit QuaternionT(const M<T, Rows, Cols>& matrix)
        {
            Gs::MatrixToQuaternion(*this, matrix);
        }

        explicit QuaternionT(UninitializeTag)
        {
            // do nothing
        }

        QuaternionT<T>& operator += (const QuaternionT<T>& rhs)
        {
            m128 = _mm_add_ps(m128, rhs.m128);
            return *this;
        }

        QuaternionT<T>& operator -= (const QuaternionT<T>& rhs)
        {
            m128 = _mm_sub_ps(m128, rhs.m128);
            return *this;
        }

        QuaternionT<T>& operator *= (const QuaternionT<T>& rhs)
        {
            *this = (*this * rhs);
            return *this;
        }

        QuaternionT<T>& operator *= (const T& rhs)
        {
            m128 = _mm_mul_ps(m128, _mm_set_ps1(rhs));
            return *this;
        }

        /**
        \brief Returns the specified quaternion component.
        \param[in] component Specifies the quaternion component index. This must be 0, 1, 2, or 3.
        */
        T& operator [] (std::size_t component)
        {
            GS_ASSERT(component < QuaternionT<T>::components);
            return *((&x) + component);
        }

        /**
        \brief Returns the specified quaternion component.
        \param[in] component Specifies the quaternion component index. This must be 0, 1, 2, or 3.
        */
        const T& operator [] (std::size_t component) const
        {
            GS_ASSERT(component < QuaternionT<T>::components);
            return *((&x) + component);
        }

        /**
        \brief Normalizes the quaternion to the unit length of 1.
        \see Normalized
        \see Length
        */
        void Normalize()
        {
            Gs::Normalize(*this);
        }

        /**
        \brief Returns a normalized instance of this quaternion.
        \see Normalize
        */
        QuaternionT<T> Normalized() const
        {
            auto quat = *this;
            quat.Normalize();
            return quat;
        }

        //! Sets this quaternion to the identity quaternion (x = y = z = 0, w = 1).
        void LoadIdentity()
        {
            x = y = z = T(0);
            w = T(1);
        }

        //! Makes this quaternion to its inverse.
        void MakeInverse()
        {
            x = -x;
            y = -y;
            z = -z;
        }

        //! Returns the inverse of this quaternion.
        QuaternionT<T> Inverse() const
        {
            return QuaternionT<T>{ -x, -y, -z, w };
        }

        /**
        \brief Computes a spherical linear interpolation between the two quaternions and stores the result into this quaternion.
        \param[in] t Specifies the interpolation factor. This should be in the range [0.0, 1.0].
        */
        void Slerp(const QuaternionT<T>& from, const QuaternionT<T>& to, const T& t)
        {
            *this = Gs::Slerp(from, to, t);
        }

        //! Sets the quaternion to an euler rotation with the specified angles (in radian).
        void SetEulerAngles(const Vector<T, 3>& angles)
        {
            const T cr = std::cos(angles.x/T(2));
            const T cp = std::cos(angles.y/T(2));
            const T cy = std::cos(angles.z/T(2));

            const T sr = std::sin(angles.x/T(2));
            const T sp = std::sin(angles.y/T(2));
            const T sy = std::sin(angles.z/T(2));

            const T cpcy = cp * cy;
            const T spsy = sp * sy;
            const T cpsy = cp * sy;
            const T spcy = sp * cy;

            x = sr * cpcy - cr * spsy;
            y = cr * spcy + sr * cpsy;
            z = cr * cpsy - sr * spcy;
            w = cr * cpcy + sr * spsy;

            Normalize();
        }

        void GetEulerAngles(Vector<T, 3>& angles) const
        {
            Vector<T, 4> sq { _mm_mul_ps(m128, m128) };
            angles.x = std::atan2(2.0f * (y*z + x*w), -sq.x - sq.y + sq.z + sq.w);
            angles.y = std::asin(Clamp(2.0f * (y*w - x*z), -1.0f, 1.0f));
            angles.z = std::atan2(2.0f * (x*y + z*w), sq.x - sq.y - sq.z + sq.w);
        }

        /**
        \brief Sets the rotation of this quaternion by the specified euler axis.
        \param[in] aixs Specifies the aixs. This must be normalized!
        \param[in] angle Specifies the rotation angle (in radians).
        */
        void SetAngleAxis(const Vector<T, 3>& axis, const T& angle)
        {
            const T halfAngle   = angle * 0.5f;
            const T sine        = std::sin(halfAngle);

            m128 = _mm_set_ps(
                std::cos(halfAngle),
                sine * axis.z,
                sine * axis.y,
                sine * axis.x
            );
        }

        void GetAngleAxis(Vector<T, 3>& axis, T& angle) const
        {
            const T scale = std::sqrt(x*x + y*y + z*z);

            if ( ( std::abs(scale) <= std::numeric_limits<T>::epsilon() ) || w > 1.0f || w < -1.0f )
            {
                axis.x  = 0.0f;
                axis.y  = 1.0f;
                axis.z  = 0.0f;
                angle   = 0.0f;
            }
            else
            {
                const T invScale = 1.0f / scale;
                axis.x  = x * invScale;
                axis.y  = y * invScale;
                axis.z  = z * invScale;
                angle   = 2.0f * std::acos(w);
            }
        }

        Matrix3T<T> ToMatrix3() const
        {
            Matrix3T<T> result(UninitializeTag{});
            Gs::QuaternionToMatrix(result, *this);
            return result;
        }

        Matrix3T<T> ToMatrix3Transposed() const
        {
            Matrix3T<T> result(UninitializeTag{});
            Gs::QuaternionToMatrixTransposed(result, *this);
            return result;
        }

        /**
        Returns a type casted instance of this quaternion.
        \tparam C Specifies the static cast type.
        */
        template <typename C>
        QuaternionT<C> Cast() const
        {
            return QuaternionT<C>(
                static_cast<C>(x),
                static_cast<C>(y),
                static_cast<C>(z),
                static_cast<C>(w)
            );
        }

        //! Returns a pointer to the first element of this quaternion.
        T* Ptr()
        {
            return &x;
        }

        //! Returns a constant pointer to the first element of this quaternion.
        const T* Ptr() const
        {
            return &x;
        }

        //! Returns a new quaternion, rotated with the specified euler angles.
        static QuaternionT<T> EulerAngles(const Vector<T, 3>& angles)
        {
            QuaternionT<T> result;
            result.SetEulerAngles(angles);
            return result;
        }

        //! Returns a new quaternion, rotated with the specified angle axis.
        static QuaternionT<T> AngleAxis(const Vector<T, 3>& axis, const T& angle)
        {
            QuaternionT<T> result;
            result.SetAngleAxis(axis, angle);
            return result;
        }

        union
        {
            __m128 m128;
            struct
            {
                T x, y, z, w;
            };
        };

};


/* --- Global Operators --- */

template <>
inline QuaternionT<float> operator + (const QuaternionT<float>& lhs, const QuaternionT<float>& rhs)
{
    return QuaternionT<float> { _mm_add_ps(lhs.m128, rhs.m128) };
}

template <>
inline QuaternionT<float> operator - (const QuaternionT<float>& lhs, const QuaternionT<float>& rhs)
{
    return QuaternionT<float> { _mm_sub_ps(lhs.m128, rhs.m128) };
}

template <>
inline QuaternionT<float> operator * (const QuaternionT<float>& lhs, const float& rhs)
{
    return QuaternionT<float> { _mm_mul_ps(lhs.m128, _mm_set_ps1(rhs)) };
}

template <>
inline QuaternionT<float> operator * (const float& lhs, const QuaternionT<float>& rhs)
{
    return QuaternionT<float> { _mm_mul_ps(_mm_set_ps1(lhs), rhs.m128) };
}


} // /namespace Gs


#endif



// ================================================================================
