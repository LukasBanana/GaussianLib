/*
 * OStream.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_OSTREAM_H__
#define __GS_OSTREAM_H__


#include "Vector2.h"
#include "Vector3.h"
#include "Vector4.h"
#include "Matrix.h"
#include "SparseMatrix4.h"

#include <iostream>


namespace Gs
{


template < template <typename> class Vec, typename T> std::ostream& operator << (std::ostream& stream, const Vec<T>& vec)
{
    stream << "( ";

    for (size_t i = 0; i < Vec<T>::components; ++i)
        stream << vec[i] << (i + 1 < Vec<T>::components ? " | " : " )");

    return stream;
}

template <typename T, size_t Rows, size_t Cols> std::ostream& operator << (std::ostream& stream, const Matrix<T, Rows, Cols>& mat)
{
    for (size_t r = 0; r < Rows; ++r)
    {
        stream << (r == 0 ? '/' : r + 1 == Rows ? '\\' : '|') << ' ';
        for (size_t c = 0; c < Cols; ++c)
            stream << mat(r, c) << (c + 1 < Cols ? '\t' : ' ');
        stream << (r == 0 ? '\\' : r + 1 == Rows ? '/' : '|')  << std::endl;
    }

    return stream;
}

template <typename T> std::ostream& operator << (std::ostream& stream, const SparseMatrix4T<T>& mat)
{
    for (size_t r = 0; r < SparseMatrix4T<T>::rowsSparse; ++r)
    {
        stream << (r == 0 ? '/' : '|') << ' ';
        for (size_t c = 0; c < SparseMatrix4T<T>::columnsSparse; ++c)
            stream << mat(r, c) << (c + 1 < SparseMatrix4T<T>::columnsSparse ? '\t' : ' ');
        stream << (r == 0 ? '\\' : '|')  << std::endl;
    }
    
    stream << "\\ 0\t0\t0\t1 /" << std::endl;

    return stream;
}


} // /namespace Gs


#endif



// ================================================================================
