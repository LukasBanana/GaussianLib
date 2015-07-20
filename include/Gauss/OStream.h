/*
 * OStream.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_OSTREAM_H__
#define __GS_OSTREAM_H__


#include "Vector3.h"

#include <iostream>


template <typename T> std::ostream& operator << (std::ostream& stream, const Gs::Vector3T<T>& v)
{
    stream << "( " << v.x << " | " << v.y << " | " << v.z << " )";
    return stream;
}


#endif



// ================================================================================
