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

#include <iostream>


template < template <typename> class Vec, typename T> std::ostream& operator << (std::ostream& stream, const Vec<T>& v)
{
    stream << "( ";

    for (size_t i = 0; i < Vec<T>::components; ++i)
        stream << v[i] << (i + 1 < Vec<T>::components ? " | " : " )");

    return stream;
}


#endif



// ================================================================================
