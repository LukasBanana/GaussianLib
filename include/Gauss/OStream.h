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
#include "AffineMatrix4.h"

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


} // /namespace Details


template < template <typename> class Vec, typename T>
std::ostream& operator << (std::ostream& stream, const Vec<T>& vec)
{
    stream << "( ";

    for (std::size_t i = 0; i < Vec<T>::components; ++i)
        stream << vec[i] << (i + 1 < Vec<T>::components ? " | " : " )");

    return stream;
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::ostream& operator << (std::ostream& stream, const Matrix<T, Rows, Cols>& mat)
{
    /* Determine longest elements for each row */
    std::array<std::size_t, AffineMatrix4T<T>::columns> lengths;

    for (std::size_t c = 0; c < AffineMatrix4T<T>::columns; ++c)
    {
        lengths[c] = 0;
        for (std::size_t r = 0; r < AffineMatrix4T<T>::rows; ++r)
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
    /* Determine longest elements for each row */
    std::array<std::size_t, AffineMatrix4T<T>::columnsSparse> lengths;

    for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
    {
        lengths[c] = 0;
        for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
            lengths[c] = std::max(lengths[c], Details::Length(mat(r, c)));
    }

    #ifdef GS_ROW_VECTORS

    /* Write each row */
    for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
    {
        stream << (r == 0 ? '/' : r + 1 == AffineMatrix4T<T>::rowsSparse ? '\\' : '|');

        for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
        {
            stream << std::string(lengths[c] + 1u - Details::Length(mat(r, c)), ' ');
            stream << mat(r, c) << ' ';
        }

        /* Write implicit column */
        stream << ' ' << (r + 1 == AffineMatrix4T<T>::rowsSparse ? '1' : '0');
        stream << ' ' << (r == 0 ? '\\' : r + 1 == AffineMatrix4T<T>::rowsSparse ? '/' : '|')  << std::endl;
    }

    #else

    /* Write each row */
    for (std::size_t r = 0; r < AffineMatrix4T<T>::rowsSparse; ++r)
    {
        stream << (r == 0 ? '/' : '|');

        for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
        {
            stream << std::string(lengths[c] + 1u - Details::Length(mat(r, c)), ' ');
            stream << mat(r, c) << ' ';
        }

        stream << (r == 0 ? '\\' : '|')  << std::endl;
    }

    /* Write implicit row */
    stream << '\\';
    
    for (std::size_t c = 0; c < AffineMatrix4T<T>::columnsSparse; ++c)
    {
        stream << std::string(lengths[c], ' ');
        stream << (c + 1 == AffineMatrix4T<T>::columnsSparse ? '1' : '0') << ' ';
    }

    stream << '/' << std::endl;

    #endif

    return stream;
}


} // /namespace Gs


#endif



// ================================================================================
