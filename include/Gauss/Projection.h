/*
 * Projection.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_PROJECTION_H
#define GS_PROJECTION_H


#include <Gauss/Decl.h>
#include <Gauss/Macros.h>
#include <Gauss/Matrix.h>
#include <Gauss/ProjectionMatrix4.h>


namespace Gs
{


template <class M, typename T>
void Orthogonal(M& mat, const T& left, const T& right, const T& top, const T& bottom, const T& near, const T& far, int flags = 0)
{
    GS_ASSERT_MxN_MATRIX("Gs::Orthogonal()", M, 4, 4);

    bool rightHanded    = (( flags & ProjectionFlags::RightHanded ) != 0);
    bool unitCube       = (( flags & ProjectionFlags::UnitCube    ) != 0);

    mat.LoadIdentity();

    mat.At(0, 0) = T(2)/(right - left);
    mat.At(1, 1) = T(2)/(top - bottom);
    mat.At(2, 2) = (unitCube ? T(2)/(far - near) : T(1)/(far - near));

    mat.At(0, 3) = -(right + left) / (right - left);
    mat.At(1, 3) = -(top + bottom) / (top - bottom);
    mat.At(2, 3) = (unitCube ? -(far + near)/(far - near) : -near/(far - near));

    if (rightHanded)
        mat.At(2, 2) = -mat.At(2, 2);
}


} // /namespace Gs


#endif



// ================================================================================
