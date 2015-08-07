/*
 * Decl.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_DECL_H__
#define __GS_DECL_H__


#include <cstddef>


namespace Gs
{


/* --- Forward Declarations --- */

// matrices
template <typename T> class AffineMatrix3T;
template <typename T> class AffineMatrix4T;

template <typename T> class ProjectionMatrix4T;

template <typename T, std::size_t Rows, std::size_t Cols> class Matrix;

// vectors
template <typename T, std::size_t N> class Vector;

// quaternions
template <typename T> class QuaternionT;


} // /namespace Gs


#endif



// ================================================================================
