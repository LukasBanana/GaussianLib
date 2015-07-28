/*
 * Projection.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_PROJECTION_H__
#define __GS_PROJECTION_H__

// <<< extension header >>>


#include "Decl.h"
#include "Macros.h"
#include "Tags.h"

#include <cmath>


namespace Gs
{


/**
Builds a 4x4 left-handed (LH) perspective projection matrix.
\param[in] aspect Specifies the aspect ratio. This can should be either resolution.x/resolution.y or resolution.y/resolution.x.
\param[in] near Specifies the near clipping plane. Should be greater than 0.0.
\param[in] far Specifies the far clipping plane. Should be greater than 0.0 and greater than 'near'.
\param[in] fov Specifies the field-of-view (FOV) (in radians).
\param[in] fovVertical Specifies whether the field-of-view is vertical (true) of horizontal (false). By default true.
\remarks This function builds a left-handed projection matrix, as it is used in Direct3D, i.e. the z values can be
transformed to normalized-device-coordinates (NDC) in the range [0.0, 1.0]. Moreover, positive z values point into the screen.
*/
template <typename T>
Matrix<T, 4, 4> BuildPerspectiveProjectionLH(const T& aspect, const T& near, const T& far, const T& fov, bool fovVertical = true)
{
    Matrix<T, 4, 4> m(UninitializeTag{});

    T h, w;

    if (fovVertical)
    {
        h = T(1) / std::tan(fov / T(2));
        w = h / aspect;
    }
    else
    {
        w = T(1) / std::tan(fov / T(2));
        h = w / aspect;
    }

    m.At(0, 0) = w;
    m.At(1, 0) = T(0);
    m.At(2, 0) = T(0);
    m.At(3, 0) = T(0);

    m.At(0, 1) = T(0);
    m.At(1, 1) = h;
    m.At(2, 1) = T(0);
    m.At(3, 1) = T(0);

    m.At(0, 2) = T(0);
    m.At(1, 2) = T(0);
    m.At(2, 2) = far/(far - near);
    m.At(3, 2) = T(1);

    m.At(0, 3) = T(0);
    m.At(1, 3) = T(0);
    m.At(2, 3) = T(0);
    m.At(3, 3) = -(far*near)/(far - near);

    return m;
}

/**
Builds a 4x4 right-handed (RH) perspective projection matrix.
\param[in] aspect Specifies the aspect ratio. This can should be either resolution.x/resolution.y or resolution.y/resolution.x.
\param[in] near Specifies the near clipping plane. Should be greater than 0.0.
\param[in] far Specifies the far clipping plane. Should be greater than 0.0 and greater than 'near'.
\param[in] fov Specifies the field-of-view (FOV) (in radians).
\param[in] fovVertical Specifies whether the field-of-view is vertical (true) of horizontal (false). By default true.
\remarks This function builds a right-handed projection matrix, as it is used in OpenGL, i.e. the z values can be
transformed to normalized-device-coordinates (NDC) in the range [-1.0, 1.0]. Moreover, positive z values point out of the screen.
*/
template <typename T>
Matrix<T, 4, 4> BuildPerspectiveProjectionRH(const T& aspect, const T& near, const T& far, const T& fov, bool fovVertical = true)
{
    Matrix<T, 4, 4> m(UninitializeTag{});

    T h, w;

    if (fovVertical)
    {
        h = T(1) / std::tan(fov / T(2));
        w = h / aspect;
    }
    else
    {
        w = T(1) / std::tan(fov / T(2));
        h = w / aspect;
    }

    m.At(0, 0) = w;
    m.At(1, 0) = T(0);
    m.At(2, 0) = T(0);
    m.At(3, 0) = T(0);

    m.At(0, 1) = T(0);
    m.At(1, 1) = h;
    m.At(2, 1) = T(0);
    m.At(3, 1) = T(0);

    m.At(0, 2) = T(0);
    m.At(1, 2) = T(0);
    m.At(2, 2) = -(far + near)/(far - near);
    m.At(3, 2) = T(-1);

    m.At(0, 3) = T(0);
    m.At(1, 3) = T(0);
    m.At(2, 3) = T(0);
    m.At(3, 3) = -(T(2)*far*near)/(far - near);

    return m;
}

/**
Builds a 4x4 left-handed (LH) orthogonal projection matrix.
\param[in] width Specifies the view width.
\param[in] height Specifies the view height.
\param[in] near Specifies the near clipping plane. Should be greater than 0.0.
\param[in] far Specifies the far clipping plane. Should be greater than 0.0 and greater than 'near'.
\remarks This function builds a right-handed projection matrix, as it is used in OpenGL, i.e. the z values can be
transformed to normalized-device-coordinates (NDC) in the range [-1.0, 1.0]. Moreover, positive z values point out of the screen.
*/
template <typename T>
Matrix<T, 4, 4> BuildOrthogonalProjectionLH(const T& width, const T& height, const T& near, const T& far)
{
    Matrix<T, 4, 4> m(UninitializeTag{});

    m.At(0, 0) = T(2)/width;
    m.At(1, 0) = T(0);
    m.At(2, 0) = T(0);
    m.At(3, 0) = T(0);

    m.At(0, 1) = T(0);
    m.At(1, 1) = T(2)/height;
    m.At(2, 1) = T(0);
    m.At(3, 1) = T(0);

    m.At(0, 2) = T(0);
    m.At(1, 2) = T(0);
    m.At(2, 2) = T(1)/(far - near);
    m.At(3, 2) = T(0);

    m.At(0, 3) = T(0);
    m.At(1, 3) = T(0);
    m.At(2, 3) = -near/(far - near);
    m.At(3, 3) = T(1);

    return m;
}

/**
Builds a 4x4 right-handed (RH) orthogonal projection matrix.
\param[in] width Specifies the view width.
\param[in] height Specifies the view height.
\param[in] near Specifies the near clipping plane. Should be greater than 0.0.
\param[in] far Specifies the far clipping plane. Should be greater than 0.0 and greater than 'near'.
\remarks This function builds a right-handed projection matrix, as it is used in OpenGL, i.e. the z values can be
transformed to normalized-device-coordinates (NDC) in the range [-1.0, 1.0]. Moreover, positive z values point out of the screen.
*/
template <typename T>
Matrix<T, 4, 4> BuildOrthogonalProjectionRH(const T& width, const T& height, const T& near, const T& far)
{
    Matrix<T, 4, 4> m(UninitializeTag{});

    m.At(0, 0) = T(2)/width;
    m.At(1, 0) = T(0);
    m.At(2, 0) = T(0);
    m.At(3, 0) = T(0);

    m.At(0, 1) = T(0);
    m.At(1, 1) = T(2)/height;
    m.At(2, 1) = T(0);
    m.At(3, 1) = T(0);

    m.At(0, 2) = T(0);
    m.At(1, 2) = T(0);
    m.At(2, 2) = -T(2)/(far - near);
    m.At(3, 2) = T(0);

    m.At(0, 3) = T(0);
    m.At(1, 3) = T(0);
    m.At(2, 3) = -(far + near)/(far - near);
    m.At(3, 3) = T(1);

    return m;
}


} // /namespace Gs


#endif



// ================================================================================
