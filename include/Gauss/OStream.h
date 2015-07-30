/*
 * OStream.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_OSTREAM_H__
#define __GS_OSTREAM_H__


#include "Decl.h"

#include <iostream>
#include <algorithm>
#include <array>
#include <string>
#include <sstream>


namespace Gs
{


namespace Details
{


template <typename T>
std::size_t Length(const T& value)
{
    /*
    Use std::stringstream instead of std::to_string to guarantee
    to have the same output of the length as in the final writing
    */
    std::stringstream stream;
    stream << value;
    return stream.str().size();
}

template <template <typename> class Vec, typename T>
std::ostream& ShiftVectorOStream(std::ostream& stream, const Vec<T>& vec)
{
    stream << "( ";

    for (std::size_t i = 0; i < Vec<T>::components; ++i)
        stream << vec[i] << (i + 1 < Vec<T>::components ? " | " : " )");

    return stream;
}

template <class M, typename T>
std::ostream& ShiftAffineMatrixOStream(std::ostream& stream, const M& mat)
{
    /* Determine longest elements for each row */
    std::array<std::size_t, M::columnsSparse> lengths;

    for (std::size_t c = 0; c < M::columnsSparse; ++c)
    {
        lengths[c] = 0;
        for (std::size_t r = 0; r < M::rowsSparse; ++r)
            lengths[c] = std::max(lengths[c], Details::Length(mat(r, c)));
    }

    #ifdef GS_ROW_VECTORS

    /* Write each row */
    for (std::size_t r = 0; r < M::rowsSparse; ++r)
    {
        stream << (r == 0 ? '/' : r + 1 == M::rowsSparse ? '\\' : '|');

        for (std::size_t c = 0; c < M::columnsSparse; ++c)
        {
            stream << std::string(lengths[c] + 1u - Details::Length(mat(r, c)), ' ');
            stream << mat(r, c) << ' ';
        }

        /* Write implicit column */
        stream << ' ' << (r + 1 == M::rowsSparse ? '1' : '0');
        stream << ' ' << (r == 0 ? '\\' : r + 1 == M::rowsSparse ? '/' : '|')  << std::endl;
    }

    #else

    /* Write each row */
    for (std::size_t r = 0; r < M::rowsSparse; ++r)
    {
        stream << (r == 0 ? '/' : '|');

        for (std::size_t c = 0; c < M::columnsSparse; ++c)
        {
            stream << std::string(lengths[c] + 1u - Details::Length(mat(r, c)), ' ');
            stream << mat(r, c) << ' ';
        }

        stream << (r == 0 ? '\\' : '|')  << std::endl;
    }

    /* Write implicit row */
    stream << '\\';
    
    for (std::size_t c = 0; c < M::columnsSparse; ++c)
    {
        stream << std::string(lengths[c], ' ');
        stream << (c + 1 == M::columnsSparse ? '1' : '0') << ' ';
    }

    stream << '/' << std::endl;

    #endif

    return stream;
}


} // /namespace Details


template <typename T>
std::ostream& operator << (std::ostream& stream, const Vector2T<T>& vec)
{
    return Details::ShiftVectorOStream(stream, vec);
}

template <typename T>
std::ostream& operator << (std::ostream& stream, const Vector3T<T>& vec)
{
    return Details::ShiftVectorOStream(stream, vec);
}

template <typename T>
std::ostream& operator << (std::ostream& stream, const Vector4T<T>& vec)
{
    return Details::ShiftVectorOStream(stream, vec);
}

template <typename T>
std::ostream& operator << (std::ostream& stream, const QuaternionT<T>& vec)
{
    return Details::ShiftVectorOStream(stream, vec);
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::ostream& operator << (std::ostream& stream, const Matrix<T, Rows, Cols>& mat)
{
    /* Determine longest elements for each row */
    std::array<std::size_t, Matrix<T, Rows, Cols>::columns> lengths;

    for (std::size_t c = 0; c < Matrix<T, Rows, Cols>::columns; ++c)
    {
        lengths[c] = 0;
        for (std::size_t r = 0; r < Matrix<T, Rows, Cols>::rows; ++r)
            lengths[c] = std::max(lengths[c], Details::Length(mat(r, c)));
    }

    /* Write each row */
    for (std::size_t r = 0; r < Rows; ++r)
    {
        stream << (r == 0 ? '/' : r + 1 == Rows ? '\\' : '|');
        
        for (std::size_t c = 0; c < Cols; ++c)
        {
            stream << std::string(lengths[c] + 1u - Details::Length(mat(r, c)), ' ');
            stream << mat(r, c) << ' ';
        }

        stream << (r == 0 ? '\\' : r + 1 == Rows ? '/' : '|')  << std::endl;
    }

    return stream;
}

template <typename T>
std::ostream& operator << (std::ostream& stream, const AffineMatrix4T<T>& mat)
{
    return Details::ShiftAffineMatrixOStream<AffineMatrix4T<T>, T>(stream, mat);
}

template <typename T>
std::ostream& operator << (std::ostream& stream, const AffineMatrix3T<T>& mat)
{
    return Details::ShiftAffineMatrixOStream<AffineMatrix3T<T>, T>(stream, mat);
}


} // /namespace Gs


#endif



// ================================================================================
