/*
 * RotateVector.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ROTATE_VECTOR_H__
#define __GS_ROTATE_VECTOR_H__


namespace Gs
{


/* --- Forward Declarations --- */

template <typename T> class Vector2T;
template <typename T> class Vector3T;
template <typename T> class Vector4T;


template <class M, typename T>
Vector2T<T> RotateVector(const M& mat, const Vector2T<T>& vec)
{
    static_assert(
        M::rows >= 2 && M::columns >= 2,
        "2D vector can only be rotated with matrix which has at least 2 rows and 2 columns"
    );
    return Vector2T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1),
        vec.x*mat(1, 0) + vec.y*mat(1, 1)
    );
}

template <class M, typename T>
Vector3T<T> RotateVector(const M& mat, const Vector3T<T>& vec)
{
    static_assert(
        M::rows >= 3 && M::columns >= 3,
        "3D vector can only be rotated with matrix which has at least 3 rows and 3 columns"
    );
    return Vector3T<T>(
        vec.x*mat(0, 0) + vec.y*mat(0, 1) + vec.z*mat(0, 2),
        vec.x*mat(1, 0) + vec.y*mat(1, 1) + vec.z*mat(1, 2),
        vec.x*mat(2, 0) + vec.y*mat(2, 1) + vec.z*mat(2, 2)
    );
}

template <class M, typename T>
Vector2T<T> RotateVectorInverse(const M& mat, const Vector2T<T>& vec)
{
    static_assert(
        M::rows >= 2 && M::columns >= 2,
        "2D vector can only be rotated with matrix which has at least 2 rows and 2 columns"
    );
    return Vector2T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1)
    );
}

template <class M, typename T>
Vector3T<T> RotateVectorInverse(const M& mat, const Vector3T<T>& vec)
{
    static_assert(
        M::rows >= 3 && M::columns >= 3,
        "3D vector can only be rotated with matrix which has at least 3 rows and 3 columns"
    );
    return Vector3T<T>(
        vec.x*mat(0, 0) + vec.y*mat(1, 0) + vec.z*mat(2, 0),
        vec.x*mat(0, 1) + vec.y*mat(1, 1) + vec.z*mat(2, 1),
        vec.x*mat(0, 2) + vec.y*mat(1, 2) + vec.z*mat(2, 2)
    );
}


} // /namespace Gs


#endif



// ================================================================================
