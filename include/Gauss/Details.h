/*
 * Details.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_DETAILS_H__
#define __GS_DETAILS_H__


#include <vector>


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class SparseMatrix4T;
template <typename T, std::size_t Rows, std::size_t Cols> class Matrix;

// Determinant
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
T Determinant(const M<T, Rows, Cols>&);

// Inverse
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
bool Inverse(M<T, Rows, Cols>&, const M<T, Rows, Cols>&);


namespace Details
{


//! Internal class for implementation details.
template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
class MatrixHelper
{

        MatrixHelper() = delete;

    protected:
        
        friend T Gs::Determinant<M, T, Cols, Rows>(const M<T, Rows, Cols>&);
        friend bool Gs::Inverse<M, T, Cols, Rows>(M<T, Rows, Cols>&, const M<T, Rows, Cols>&);

        static std::vector<T> MatrixToArray(const M<T, Rows, Cols>& mat)
        {
            std::vector<T> vec(Rows*Cols);

            for (std::size_t r = 0, i = 0; r < Rows; ++r)
            {
                for (std::size_t c = 0; c < Cols; ++c)
                    vec[i++] = mat(r, c);
            }

            return vec;
        }

        static T OrderedDeterminant(const std::vector<T>& mat, std::size_t order)
        {
            if (order == 1)
                return mat[0];

            std::vector<T> minor((order - 1)*(order - 1));

            T det = T(0);

            for (std::size_t i = 0; i < order; ++i)
            {
                GetMinorMatrix(mat, minor, 0, i, order);
                if (i % 2 == 1)
                    det -= mat[i] * OrderedDeterminant(minor, order - 1);
                else
                    det += mat[i] * OrderedDeterminant(minor, order - 1);
            }
    
            return det;
        }

        static void OrderedInverse(const std::vector<T>& mat, std::size_t order)
        {
            T det = OrderedDeterminant(mat, order);

            //todo...
        }

    private:

        static void GetMinorMatrix(
            const std::vector<T>& mat, std::vector<T>& minor, std::size_t row, std::size_t column, std::size_t order)
        {
            for (std::size_t r = 1, i = 0; r < order; ++r)
            {
                if (r != row)
                {
                    for (std::size_t c = 0, j = 0; c < order; ++c)
                    {
                        if (c != column)
                        {
                            minor[i*(order - 1) + j] = mat[r*order + c];
                            ++j;
                        }
                    }
                    ++i;
                }
            }
        }

};


} // /namespace Details


} // /namespace Gs


#endif



// ================================================================================
