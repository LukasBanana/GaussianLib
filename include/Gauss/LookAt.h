/*
 * LookAt.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_LOOKAT_H
#define GS_LOOKAT_H


#include <Gauss/Decl.h>
#include <Gauss/Macros.h>
#include <Gauss/Matrix.h>
#include <Gauss/AffineMatrix4.h>
#include <Gauss/Vector3.h>
#include <Gauss/Algebra.h>


namespace Gs
{


namespace Details
{


template <class M, typename T>
void LookAt(M& m, const Vector3T<T>& position, const Vector3T<T>& target, const Vector3T<T>& upVector)
{
    GS_ASSERT_MxN_MATRIX("Gs::LookAt()", M, 3, 4);

    Vector3T<T> zAxis = (target - position).Normalized();
    Vector3T<T> xAxis = Cross(upVector.Normalized(), zAxis).Normalized();
    Vector3T<T> yAxis = Cross(zAxis, xAxis);

    m.At(0, 0) = xAxis.x;
    m.At(1, 0) = xAxis.y;
    m.At(2, 0) = xAxis.z;

    m.At(0, 1) = yAxis.x;
    m.At(1, 1) = yAxis.y;
    m.At(2, 1) = yAxis.z;

    m.At(0, 2) = zAxis.x;
    m.At(1, 2) = zAxis.y;
    m.At(2, 2) = zAxis.z;

    m.At(0, 3) = position.x;
    m.At(1, 3) = position.y;
    m.At(2, 3) = position.z;
}


} // /namespace Details


//! Constructs a 4x4 matrix that is position and orientated to look at a specific target location.
template <typename T>
void LookAt(Matrix<T, 4, 4>& m, const Vector3T<T>& position, const Vector3T<T>& target, const Vector3T<T>& upVector)
{
    Details::LookAt(m, position, target, upVector);
    m.At(3, 0) = 0.0f;
    m.At(3, 1) = 0.0f;
    m.At(3, 2) = 0.0f;
    m.At(3, 3) = 1.0f;
}

//! Constructs an affine 4x4 matrix that is position and orientated to look at a specific target location.
template <typename T>
void LookAt(AffineMatrix4T<T>& m, const Vector3T<T>& position, const Vector3T<T>& target, const Vector3T<T>& upVector)
{
    Details::LookAt(m, position, target, upVector);
}



} // /namespace Gs


#endif



// ================================================================================
