/*
 * AffineMatrix4.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_AFFINE_MATRIX4_H__
#define __GS_AFFINE_MATRIX4_H__


#include "Real.h"
#include "Assert.h"
#include "Macros.h"
#include "Matrix.h"
#include "Tags.h"

#include <cmath>
#include <cstring>
#include <algorithm>


namespace Gs
{


#ifdef GS_ROW_MAJOR_STORAGE
#   define __GS_FOREACH_ROW_COL__(r, c)                                     \
        for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)     \
        for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
#else
#   define __GS_FOREACH_ROW_COL__(r, c)                                     \
        for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)  \
        for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
#endif

/**
This is an affine 4x4 matrix for affine transformations,
i.e. it can only contain translations, scaling, rotations and shearing.
It only stores a 3x4 matrix where the 4th row is always implicitly (0, 0, 0, 1).
\tparam T Specifies the data type of the matrix components.
This should be a primitive data type such as float, double, int etc.
\remarks The macro GS_ROW_MAJOR_STORAGE can be defined, to use row-major storage layout.
By default column-major storage layout is used.
The macro GS_ROW_VECTORS can be defined, to use row vectors. By default column vectors are used.
Here is an example, how an affine 4x4 matrix is laid-out with column- and row vectors:
\code
// Affine 4x4 matrix with column vectors:
// / x1 y1 z1 w1 \
// | x2 y2 z2 w2 |
// | x3 y3 z3 w3 |
// \  0  0  0  1 /

// Affine 4x4 matrix with row vectors:
// / x1 x2 x3 0 \
// | y1 y2 y3 0 |
// | z1 z2 z3 0 |
// \ w1 w2 w3 1 /

// In both cases, (w1, w2, w3, 1) stores the position in an affine transformation.
\endcode
*/
template <typename T>
class AffineMatrix4T
{
    
    public:
        
        static const std::size_t rows           = 4;
        static const std::size_t columns        = 4;
        static const std::size_t elements       = AffineMatrix4T<T>::rows*AffineMatrix4T<T>::columns;

        #ifdef GS_ROW_VECTORS
        static const std::size_t rowsSparse     = 4;
        static const std::size_t columnsSparse  = 3;
        #else
        static const std::size_t rowsSparse     = 3;
        static const std::size_t columnsSparse  = 4;
        #endif

        static const std::size_t elementsSparse = AffineMatrix4T<T>::rowsSparse*AffineMatrix4T<T>::columnsSparse;

        using ThisType = AffineMatrix4T<T>;

        //! Transposed matrix type, i.e. AffineMatrix4T<T> becomes Matrix<T, 4, 4>.
        using TransposedType = Matrix<T, AffineMatrix4T<T>::rows, AffineMatrix4T<T>::columns>;

        class Initializer
        {
            
            public:
                
                Initializer(ThisType& matrix) :
                    matrix_ ( matrix ),
                    element_( 0      )
                {
                }

                Initializer& operator , (const T& nextValue)
                {
                    matrix_(element_ / AffineMatrix4T<T>::columnsSparse, element_ % AffineMatrix4T<T>::columnsSparse) = nextValue;
                    ++element_;
                    return *this;
                }

            private:

                ThisType&   matrix_;
                std::size_t element_;

        };

        AffineMatrix4T()
        {
            #ifndef GS_ENABLE_AUTO_INIT
            Reset();
            #endif
        }

        AffineMatrix4T(const ThisType& rhs)
        {
            *this = rhs;
        }

        #ifdef GS_ROW_VECTORS

        AffineMatrix4T(
            const T& m11, const T& m12, const T& m13,
            const T& m21, const T& m22, const T& m23,
            const T& m31, const T& m32, const T& m33,
            const T& m41, const T& m42, const T& m43)
        {
            (*this)(0, 0) = m11; (*this)(0, 1) = m12; (*this)(0, 2) = m13;
            (*this)(1, 0) = m21; (*this)(1, 1) = m22; (*this)(1, 2) = m23;
            (*this)(2, 0) = m31; (*this)(2, 1) = m32; (*this)(2, 2) = m33;
            (*this)(3, 0) = m41; (*this)(3, 1) = m42; (*this)(3, 2) = m43;
        }

        #else

        AffineMatrix4T(
            const T& m11, const T& m12, const T& m13, const T& m14,
            const T& m21, const T& m22, const T& m23, const T& m24,
            const T& m31, const T& m32, const T& m33, const T& m34)
        {
            (*this)(0, 0) = m11; (*this)(0, 1) = m12; (*this)(0, 2) = m13; (*this)(0, 3) = m14;
            (*this)(1, 0) = m21; (*this)(1, 1) = m22; (*this)(1, 2) = m23; (*this)(1, 3) = m24;
            (*this)(2, 0) = m31; (*this)(2, 1) = m32; (*this)(2, 2) = m33; (*this)(2, 3) = m34;
        }

        #endif

        AffineMatrix4T(UninitializeTag)
        {
            // do nothing
        }

        /**
        \brief Returns a reference to the element at the specified entry.
        \param[in] row Specifies the row index. This must be in the range [0, 2], or [0, 3] if GS_ROW_VECTORS is defined.
        \param[in] col Specifies the column index. This must be in the range [0, 3], or [0, 2] if GS_ROW_VECTORS is defined.
        */
        T& operator () (std::size_t row, std::size_t col)
        {
            GS_ASSERT(row < AffineMatrix4T<T>::rowsSparse);
            GS_ASSERT(col < AffineMatrix4T<T>::columnsSparse);
            #ifdef GS_ROW_MAJOR_STORAGE
            return m_[row*AffineMatrix4T<T>::columnsSparse + col];
            #else
            return m_[col*AffineMatrix4T<T>::rowsSparse + row];
            #endif
        }

        /**
        \brief Returns a constant reference to the element at the specified entry.
        \param[in] row Specifies the row index. This must be in the range [0, 2], or [0, 3] if GS_ROW_VECTORS is defined.
        \param[in] col Specifies the column index. This must be in the range [0, 3], or [0, 2] if GS_ROW_VECTORS is defined.
        */
        const T& operator () (std::size_t row, std::size_t col) const
        {
            GS_ASSERT(row < AffineMatrix4T<T>::rowsSparse);
            GS_ASSERT(col < AffineMatrix4T<T>::columnsSparse);
            #ifdef GS_ROW_MAJOR_STORAGE
            return m_[row*AffineMatrix4T<T>::columnsSparse + col];
            #else
            return m_[col*AffineMatrix4T<T>::rowsSparse + row];
            #endif
        }

        T& operator [] (std::size_t element)
        {
            GS_ASSERT(element < AffineMatrix4T<T>::elementsSparse);
            return m_[element];
        }

        const T& operator [] (std::size_t element) const
        {
            GS_ASSERT(element < AffineMatrix4T<T>::elementsSparse);
            return m_[element];
        }

        ThisType& operator *= (const ThisType& rhs)
        {
            *this = (*this * rhs);
            return *this;
        }

        ThisType& operator = (const ThisType& rhs)
        {
            for (std::size_t i = 0; i < ThisType::elementsSparse; ++i)
                m_[i] = rhs.m_[i];
            return *this;
        }

        #ifdef GS_ROW_VECTORS

        T& At(std::size_t col, std::size_t row)
        {
            return (*this)(row, col);
        }

        const T& At(std::size_t col, std::size_t row) const
        {
            return (*this)(row, col);
        }

        #else

        T& At(std::size_t row, std::size_t col)
        {
            return (*this)(row, col);
        }

        const T& At(std::size_t row, std::size_t col) const
        {
            return (*this)(row, col);
        }

        #endif

        void Reset()
        {
            for (std::size_t i = 0; i < ThisType::elementsSparse; ++i)
                m_[i] = T(0);
        }

        void LoadIdentity()
        {
            __GS_FOREACH_ROW_COL__(r, c)
            {
                (*this)(r, c) = (r == c ? T(1) : T(0));
            }
        }

        static ThisType Identity()
        {
            ThisType result;
            result.LoadIdentity();
            return result;
        }

        TransposedType Transposed() const
        {
            TransposedType result;

            __GS_FOREACH_ROW_COL__(r, c)
            {
                result(c, r) = (*this)(r, c);
            }

            for (std::size_t i = 0; i < ThisType::columnsSparse; ++i)
                result(i, ThisType::rowsSparse);

            return result;
        }

        T Determinant() const
        {
            return Gs::Determinant(*this);
        }

        //! Returns the trace of this matrix: M(0, 0) + M(1, 1) + M(2, 2) + 1.
        T Trace() const
        {
            return (*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2) + T(1);
        }

        AffineMatrix4T<T> Inverse() const
        {
            AffineMatrix4T<T> inv{ *this };
            inv.MakeInverse();
            return inv;
        }

        bool MakeInverse()
        {
            AffineMatrix4T<T> in{ *this };
            return Gs::Inverse(*this, in);
        }

        //! Returns a pointer to the first element of this matrix.
        T* Ptr()
        {
            return &(m_[0]);
        }

        //! Returns a constant pointer to the first element of this matrix.
        const T* Ptr() const
        {
            return &(m_[0]);
        }

        /* --- Extended functions for affine transformations --- */

        void SetPosition(const Vector3T<T>& position)
        {
            At(0, 3) = position.x;
            At(1, 3) = position.y;
            At(2, 3) = position.z;
        }

        Vector3T<T> GetPosition() const
        {
            return Vector3T<T>(At(0, 3), At(1, 3), At(2, 3));
        }

    private:
        
        T m_[ThisType::elementsSparse];

};

#undef __GS_FOREACH_ROW_COL__


/* --- Global Operators --- */

template <typename T> AffineMatrix4T<T> operator * (const AffineMatrix4T<T>& lhs, const AffineMatrix4T<T>& rhs)
{
    AffineMatrix4T<T> result(UninitializeTag{});

    #ifdef GS_ROW_VECTORS

    for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
    {
        for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
        {
            /* Only accumulate with 'rowsSparse' here! */
            result(r, c) = T(0);
            for (std::size_t i = 0; i < AffineMatrix4T<T>::columnsSparse; ++i)
                result(r, c) += lhs(r, i)*rhs(i, c);
        }

        /* Accumulate the rest of the current column of 'lhs' and the implicit 1 of 'rhs' */
        result(AffineMatrix4T<T>::rowsSparse - 1, c) += lhs(AffineMatrix4T<T>::rowsSparse - 1, c);
    }

    #else

    for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
    {
        for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
        {
            /* Only accumulate with 'rowsSparse' here! */
            result(r, c) = T(0);
            for (std::size_t i = 0; i < AffineMatrix4T<T>::rowsSparse; ++i)
                result(r, c) += lhs(r, i)*rhs(i, c);
        }

        /* Accumulate the rest of the current row of 'lhs' and the implicit 1 of 'rhs' */
        result(r, AffineMatrix4T<T>::columnsSparse - 1) += lhs(r, AffineMatrix4T<T>::columnsSparse - 1);
    }

    #endif

    return result;
}

template <typename T, typename I>
typename AffineMatrix4T<T>::Initializer operator << (AffineMatrix4T<T>& matrix, const I& firstValue)
{
    typename AffineMatrix4T<T>::Initializer initializer(matrix);
    initializer , static_cast<T>(firstValue);
    return initializer;
}


/* --- Type Alias --- */

using AffineMatrix4     = AffineMatrix4T<Real>;
using AffineMatrix4f    = AffineMatrix4T<float>;
using AffineMatrix4d    = AffineMatrix4T<double>;
using AffineMatrix4i    = AffineMatrix4T<int>;


} // /namespace Gs


#endif



// ================================================================================
