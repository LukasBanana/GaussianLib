/*
 * Rotate.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_ROTATE_H
#define GS_ROTATE_H


#include <Gauss/Decl.h>
#include <Gauss/Macros.h>
#include <Gauss/Vector3.h>
#include <Gauss/Translate.h>


namespace Gs
{


namespace Details
{


template <class M, typename T>
void FreeRotation(M& mat, const Vector3T<T>& axis, const T& angle)
{
    GS_ASSERT_MxN_MATRIX("free rotation", M, 3, 3);

    /* Setup rotation values */
    const T  c  = std::cos(angle);
    const T  s  = std::sin(angle);
    const T  cc = T(1) - c;

    const T& x  = axis.x;
    const T& y  = axis.y;
    const T& z  = axis.z;

    /* Perform matrix rotation */
    mat.At(0, 0) = x*x*cc + c;
    mat.At(1, 0) = x*y*cc - z*s;
    mat.At(2, 0) = x*z*cc + y*s;

    mat.At(0, 1) = y*x*cc + z*s;
    mat.At(1, 1) = y*y*cc + c;
    mat.At(2, 1) = y*z*cc - x*s;

    mat.At(0, 2) = x*z*cc - y*s;
    mat.At(1, 2) = y*z*cc + x*s;
    mat.At(2, 2) = z*z*cc + c;
}


} // /namespace Details


/**
\brief Computes a free rotation around an axis and stores the result into the matrix 'm'.
\tparam M Specifies the matrix type. This should be Matrix3, Matrix4, or AffineMatrix4.
This type must implement the following interface:
\code
static const std::size_t rows;    // >= 3
static const std::size_t columns; // >= 3
\endcode
\tparam T Specifies the data type. This should should be float or double.
\param[in,out] mat Specifies the matrix which is to be rotated.
\param[in] axis Specifies the rotation axis. This must be normalized!
\param[in] angle Specifies the rotation angle (in radians).
*/
template <typename T>
void FreeRotation(Matrix<T, 4, 4>& mat, const Vector3T<T>& axis, const T& angle)
{
    Details::FreeRotation(mat, axis, angle);
}

template <typename T>
void FreeRotation(Matrix<T, 3, 3>& mat, const Vector3T<T>& axis, const T& angle)
{
    Details::FreeRotation(mat, axis, angle);
}

template <typename T>
void FreeRotation(AffineMatrix4T<T>& mat, const Vector3T<T>& axis, const T& angle)
{
    Details::FreeRotation(mat, axis, angle);
}

/**
\brief Rotates the matrix 'mat' around the specified axis.
\param[in,out] mat Specifies the matrix that is to be modified.
\param[in] axis Specifies the rotation axis. This must be normalized!
\param[in] angle Specifies the rotation angle (in radians).
\see FreeRotation
*/
template <class M, typename T>
void RotateFree(M& mat, const Vector3T<T>& axis, const T& angle)
{
    auto rotation = M::Identity();
    FreeRotation(rotation, axis, angle);
    mat *= rotation;
}

/**
\brief Rotates the matrix 'mat' around the specified axis and a pivot. The pivot offsets the center of the matrix.
\param[in,out] mat Specifies the matrix that is to be modified.
\param[in] axis Specifies the rotation axis. This must be normalized!
\param[in] angle Specifies the rotation angle (in radians).
\param[in] pivot Specifies the pivot to offset the center of rotation.
\see FreeRotation
*/
template <class M, typename T>
void RotateFree(M& mat, const Vector3T<T>& axis, const T& angle, const Vector3T<T>& pivot)
{
    /* Calculate offset between matrix center and pivot */
    Matrix<T, 3, 3> rotation;
    FreeRotation(rotation, axis, angle);
    const Vector3T<T> offset = rotation * pivot;

    /* Move and then rotate around the pivot */
    Translate(mat, pivot - offset);
    RotateFree(mat, axis, angle);
}

//! Rotates the matrix at the X-axis with the specified angle (in radians).
template <class M, typename T>
void RotateX(M& mat, const T& angle)
{
    const T c = std::cos(angle);
    const T s = std::sin(angle);

    /* Temporaries */
    const T m01 = mat.At(0, 1);
    const T m11 = mat.At(1, 1);
    const T m21 = mat.At(2, 1);

    /* Rotation */
    mat.At(0, 1) = m01*c + mat.At(0, 2)*s;
    mat.At(1, 1) = m11*c + mat.At(1, 2)*s;
    mat.At(2, 1) = m21*c + mat.At(2, 2)*s;

    mat.At(0, 2) = mat.At(0, 2)*c - m01*s;
    mat.At(1, 2) = mat.At(1, 2)*c - m11*s;
    mat.At(2, 2) = mat.At(2, 2)*c - m21*s;
}

//! Rotates the matrix at the Y-axis with the specified angle (in radians).
template <class M, typename T>
void RotateY(M& mat, const T& angle)
{
    const T c = std::cos(angle);
    const T s = std::sin(angle);

    /* Temporaries */
    const T m00 = mat.At(0, 0);
    const T m10 = mat.At(1, 0);
    const T m20 = mat.At(2, 0);

    /* Rotation */
    mat.At(0, 0) = m00*c + mat.At(0, 2)*s;
    mat.At(1, 0) = m10*c + mat.At(1, 2)*s;
    mat.At(2, 0) = m20*c + mat.At(2, 2)*s;

    mat.At(0, 2) = mat.At(0, 2)*c - m00*s;
    mat.At(1, 2) = mat.At(1, 2)*c - m10*s;
    mat.At(2, 2) = mat.At(2, 2)*c - m20*s;
}

//! Rotates the matrix at the Z-axis with the specified angle (in radians).
template <class M, typename T>
void RotateZ(M& mat, const T& angle)
{
    const T c = std::cos(angle);
    const T s = std::sin(angle);

    /* Temporaries */
    const T m00 = mat.At(0, 0);
    const T m10 = mat.At(1, 0);
    const T m20 = mat.At(2, 0);

    /* Rotation */
    mat.At(0, 0) = m00*c + mat.At(0, 1)*s;
    mat.At(1, 0) = m10*c + mat.At(1, 1)*s;
    mat.At(2, 0) = m20*c + mat.At(2, 1)*s;

    mat.At(0, 1) = mat.At(0, 1)*c - m00*s;
    mat.At(1, 1) = mat.At(1, 1)*c - m10*s;
    mat.At(2, 1) = mat.At(2, 1)*c - m20*s;
}


} // /namespace Gs


#endif



// ================================================================================
