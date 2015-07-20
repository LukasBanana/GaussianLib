/*
 * Matrix.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_MATRIX_H__
#define __GS_MATRIX_H__


#include "Real.h"
#include "Assert.h"

#include <cmath>
#include <cstring>


namespace Gs
{


#define __GS_ASSERT_NxN_MATRIX__ \
    static_assert(Rows == Cols, __FUNCTION__ " can only be used with NxN matrices")

#ifdef GS_MATRIX_COLUMN_MAJOR
#   define __GS_FOREACH_ROW_COL__(r, c)     \
        for (size_t r = 0; r < Rows; ++r)   \
        for (size_t c = 0; c < Rows; ++c)
#else
#   define __GS_FOREACH_ROW_COL__(r, c)     \
        for (size_t c = 0; c < Rows; ++c)   \
        for (size_t r = 0; r < Rows; ++r)
#endif


/**
\brief Base matrix class.
\tparam T Specifies the data type of the matrix components.
This should be a primitive data type such as float, double, int etc.
\remarks The macro GS_MATRIX_COLUMN_MAJOR can be defined, to use column-major matrices.
By default row-major matrices are used.
*/
template <typename T, size_t Rows, size_t Cols> class Matrix
{
    
    public:
        
        static const size_t rows        = Rows;
        static const size_t columns     = Cols;

        //! Transposed matrix type, i.e. NxM becomes MxN.
        using ThisType          = Matrix<T, Rows, Cols>;
        using TransposedType    = Matrix<T, Cols, Rows>;

        class Initializer
        {
            
            public:
                
                Initializer(ThisType& matrix) :
                    matrix  { matrix },
                    element { 0      }
                {
                }

                Initializer& operator , (const T& nextValue)
                {
                    matrix(element / Cols, element % Rows) = nextValue;
                    ++element;
                    return *this;
                }

            private:

                ThisType&   matrix;
                size_t      element;

        };

        Matrix()
        {
            Reset();
        }
        Matrix(const ThisType& rhs)
        {
            *this = rhs;
        }

        T& operator () (size_t row, size_t col)
        {
            GS_ASSERT(row < Rows);
            GS_ASSERT(col < Cols);
            #ifdef GS_MATRIX_COLUMN_MAJOR
            return m_[col*Rows + row];
            #else
            return m_[row*Cols + col];
            #endif
        }

        const T& operator () (size_t row, size_t col) const
        {
            GS_ASSERT(row < Rows);
            GS_ASSERT(col < Cols);
            #ifdef GS_MATRIX_COLUMN_MAJOR
            return m_[col*Rows + row];
            #else
            return m_[row*Cols + col];
            #endif
        }

        T& operator [] (size_t element)
        {
            GS_ASSERT(element < Rows*Cols);
            return m_[element];
        }

        const T& operator [] (size_t element) const
        {
            GS_ASSERT(element < Rows*Cols);
            return m_[element];
        }

        ThisType operator *= (const ThisType& rhs)
        {
            __GS_ASSERT_NxN_MATRIX__;
            *this = (*this * rhs);
            return *this;
        }

        ThisType& operator = (const ThisType& rhs)
        {
            for (size_t i = 0; i < Rows*Cols; ++i)
                m_[i] = rhs.m_[i];
            return *this;
        }

        void Reset()
        {
            for (size_t i = 0; i < Rows*Cols; ++i)
                m_[i] = T(0);
        }

        void LoadIdentitiy()
        {
            __GS_ASSERT_NxN_MATRIX__;
            __GS_FOREACH_ROW_COL__(r, c)
            {
                (*this)(r, c) = (r == c ? T(1) : T(0));
            }
        }

        TransposedType Transposed() const
        {
            TransposedType result;

            __GS_FOREACH_ROW_COL__(r, c)
            {
                result(r, c) = (*this)(c, r);
            }

            return result;
        }

        void Transpose()
        {
            __GS_ASSERT_NxN_MATRIX__;

            for (size_t i = 0; i + 1 < Cols; ++i)
            {
                for (size_t j = 1; j + i < Cols; ++j)
                {
                    std::swap(
                        m_[i*(Cols + 1) + j],
                        m_[(j + i)*Cols + i]
                    );
                }
            }
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

    private:
        
        T m_[Rows*Cols];

};

#undef __GS_ASSERT_NxN_MATRIX__
#undef __GS_FOREACH_ROW_COL__


/* --- Global Operators --- */

template <typename T, size_t Rows, size_t Cols>
Matrix<T, Cols, Cols> operator * (const Matrix<T, Rows, Cols>& lhs, const Matrix<T, Cols, Rows>& rhs)
{
    Matrix<T, Cols, Cols> result;

    for (size_t r = 0; r < Cols; ++r)
    {
        for (size_t c = 0; c < Cols; ++c)
        {
            result(r, c) = T(0);
            for (size_t i = 0; i < Cols; ++i)
                result(r, c) += lhs(r, i)*rhs(i, c);
        }
    }

    return result;
}

template <typename T, typename I, size_t Rows, size_t Cols>
typename Matrix<T, Rows, Cols>::Initializer operator << (Matrix<T, Rows, Cols>& matrix, const I& firstValue)
{
    typename Matrix<T, Rows, Cols>::Initializer initializer(matrix);
    initializer , firstValue;
    return initializer;
}


/* --- Type Alias --- */

#define __GS_DEF_MATRIX_TYPES_NxN__(n)                          \
    template <typename T> using Matrix##n##T = Matrix<T, n, n>; \
    using Matrix##n = Matrix##n##T<Real>;                       \
    using Matrix##n##f = Matrix##n##T<float>;                   \
    using Matrix##n##d = Matrix##n##T<double>;                  \
    using Matrix##n##i = Matrix##n##T<int>

__GS_DEF_MATRIX_TYPES_NxN__(2);
__GS_DEF_MATRIX_TYPES_NxN__(3);
__GS_DEF_MATRIX_TYPES_NxN__(4);

#undef __GS_DEF_MATRIX_TYPES_NxN__


} // /namespace Gs


#endif



// ================================================================================
