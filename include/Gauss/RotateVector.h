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
    __GS_ASSERT_MxN_MATRIX__("2D vector rotation by matrix", M, 2, 2);
    return Vector2T<T>(
        vec.x*mat.At(0, 0) + vec.y*mat.At(0, 1),
        vec.x*mat.At(1, 0) + vec.y*mat.At(1, 1)
    );
}

template <class M, typename T>
Vector3T<T> RotateVector(const M& mat, const Vector3T<T>& vec)
{
    __GS_ASSERT_MxN_MATRIX__("2D vector rotation by matrix", M, 3, 3);
    return Vector3T<T>(
        vec.x*mat.At(0, 0) + vec.y*mat.At(0, 1) + vec.z*mat.At(0, 2),
        vec.x*mat.At(1, 0) + vec.y*mat.At(1, 1) + vec.z*mat.At(1, 2),
        vec.x*mat.At(2, 0) + vec.y*mat.At(2, 1) + vec.z*mat.At(2, 2)
    );
}

template <class M, typename T>
Vector2T<T> RotateVectorInverse(const M& mat, const Vector2T<T>& vec)
{
    __GS_ASSERT_MxN_MATRIX__("2D vector inverse rotation by matrix", M, 2, 2);
    return Vector2T<T>(
        vec.x*mat.At(0, 0) + vec.y*mat.At(1, 0),
        vec.x*mat.At(0, 1) + vec.y*mat.At(1, 1)
    );
}

template <class M, typename T>
Vector3T<T> RotateVectorInverse(const M& mat, const Vector3T<T>& vec)
{
    __GS_ASSERT_MxN_MATRIX__("2D vector inverse rotation by matrix", M, 3, 3);
    return Vector3T<T>(
        vec.x*mat.At(0, 0) + vec.y*mat.At(1, 0) + vec.z*mat.At(2, 0),
        vec.x*mat.At(0, 1) + vec.y*mat.At(1, 1) + vec.z*mat.At(2, 1),
        vec.x*mat.At(0, 2) + vec.y*mat.At(1, 2) + vec.z*mat.At(2, 2)
    );
}


} // /namespace Gs


#endif



// ================================================================================
