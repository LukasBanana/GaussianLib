/*
 * TransformVector.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_TRANSFORM_VECTOR_H__
#define __GS_TRANSFORM_VECTOR_H__


#include "Macros.h"


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class Vector2T;
template <typename T> class Vector3T;
template <typename T> class Vector4T;


template <typename M, typename T>
Vector2T<T> TransformVector(const M& mat, const Vector2T<T>& vec)
{
    __GS_ASSERT_MxN_MATRIX__("2D vector transformation by matrix", M, 2, 2);
    #ifdef GS_ROW_VECTORS
    return Vector2T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0) + mat(2, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1) + mat(2, 1)
    );
    #else
    return Vector2T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + mat(0, 2),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + mat(1, 2)
    );
    #endif
}

template <typename M, typename T>
Vector3T<T> TransformVector(const M& mat, const Vector3T<T>& vec)
{
    #ifdef GS_ROW_VECTORS
    __GS_ASSERT_MxN_MATRIX__("3D row vector transformation by matrix", M, 4, 3);
    return Vector3T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0) + vec.z*mat(2, 0) + mat(3, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1) + vec.z*mat(2, 1) + mat(3, 1),
        vec.x*mat(0, 2) + vec.y*mat(1, 2) + vec.z*mat(2, 2) + mat(3, 2)
    );
    #else
    __GS_ASSERT_MxN_MATRIX__("3D column vector transformation by matrix", M, 3, 4);
    return Vector3T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + mat(2, 3)
    );
    #endif
}

template <typename M, typename T>
Vector4T<T> TransformVector(const M& mat, const Vector4T<T>& vec)
{
    #ifdef GS_ROW_VECTORS
    __GS_ASSERT_MxN_MATRIX__("4D row vector transformation by matrix", M, 4, 4);
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0) + vec.z*mat(2, 0) + vec.w*mat(3, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1) + vec.z*mat(2, 1) + vec.w*mat(3, 1),
        vec.x*mat(0, 2) + vec.y*mat(1, 2) + vec.z*mat(2, 2) + vec.w*mat(3, 2),
        vec.x*mat(0, 3) + vec.y*mat(1, 3) + vec.z*mat(2, 3) + vec.w*mat(3, 3)
    );
    #else
    __GS_ASSERT_MxN_MATRIX__("4D column vector transformation by matrix", M, 4, 4);
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + vec.w*mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + vec.w*mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + vec.w*mat(2, 3),
        vec.x*mat(3, 0) + vec.y*mat(3, 1) + vec.z*mat(3, 2) + vec.w*mat(3, 3)
    );
    #endif
}

template <typename T>
Vector4T<T> TransformVector(const SparseMatrix4T<T>& mat, const Vector4T<T>& vec)
{
    #ifdef GS_ROW_VECTORS
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0) + vec.z*mat(2, 0) + vec.w*mat(3, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1) + vec.z*mat(2, 1) + vec.w*mat(3, 1),
        vec.x*mat(0, 2) + vec.y*mat(1, 2) + vec.z*mat(2, 2) + vec.w*mat(3, 2),
        vec.w
    );
    #else
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + vec.w*mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + vec.w*mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + vec.w*mat(2, 3),
        vec.w
    );
    #endif
}


} // /namespace Gs


#endif



// ================================================================================
