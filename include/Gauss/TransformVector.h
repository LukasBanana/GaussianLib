/*
 * TransformVector.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_TRANSFORM_VECTOR_H__
#define __GS_TRANSFORM_VECTOR_H__


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class Vector2T;
template <typename T> class Vector3T;
template <typename T> class Vector4T;


template <class M, typename T>
Vector2T<T> TransformVector(const M& mat, const Vector2T<T>& vec)
{
    static_assert(
        M::rows >= 2 && M::columns >= 3,
        "2D vector can only be transformed with matrix which has at least 2 rows and 3 columns"
    );
    return Vector2T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + mat(0, 2),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + mat(1, 2)
    );
}

template <class M, typename T>
Vector3T<T> TransformVector(const M& mat, const Vector3T<T>& vec)
{
    static_assert(
        M::rows >= 3 && M::columns >= 4,
        "3D vector can only be transformed with matrix which has at least 3 rows and 4 columns"
    );
    return Vector3T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + mat(2, 3)
    );
}

template <template <typename, std::size_t, std::size_t> class M, typename T, std::size_t Rows, std::size_t Cols>
Vector4T<T> TransformVector(const M<T, Rows, Cols>& mat, const Vector4T<T>& vec)
{
    static_assert(
        M<T, Rows, Cols>::rows >= 4 && M<T, Rows, Cols>::columns >= 4,
        "4D vector can only be transformed with matrix which has at least 4 rows and 4 columns"
    );
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + vec.w*mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + vec.w*mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + vec.w*mat(2, 3),
        vec.x*mat(3, 0) + vec.y*mat(3, 1) + vec.z*mat(3, 2) + vec.w*mat(3, 3)
    );
}

template <typename T>
Vector4T<T> TransformVector(const SparseMatrix4T<T>& mat, const Vector4T<T>& vec)
{
    return Vector4T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2) + vec.w*mat(0, 3),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2) + vec.w*mat(1, 3),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2) + vec.w*mat(2, 3),
        vec.w
    );
}


} // /namespace Gs


#endif



// ================================================================================
