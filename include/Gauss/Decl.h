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

template <typename T> class SparseMatrix4T;
template <typename T, std::size_t Rows, std::size_t Cols> class Matrix;

template <typename T> class Vector2T;
template <typename T> class Vector3T;
template <typename T> class Vector4T;

template <typename T> class QuaternionT;


} // /namespace Gs


#endif



// ================================================================================
