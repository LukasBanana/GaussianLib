/*
 * Rotate.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_ROTATE_H__
#define __GS_ROTATE_H__


#include "Decl.h"
#include "Macros.h"


namespace Gs
{


/**
\brief Computes a free rotation around an axis and stores the result into the matrix 'm'.
\tparam M Specifies the matrix type. This should be Matrix3, Matrix4, or SparseMatrix4.
This type must implement the following interface:
\code
static const std::size_t rows;    // >= 3
static const std::size_t columns; // >= 3
\endcode
\tparam Vec Specifies the vector type. This should be Vector3, or Vector4.
\tparam T Specifies the data type. This should should be float or double.
\param[out] m Specifies the resulting matrix.
\param[in] axis Specifies the rotation axis. This must be normalized!
\param[in] angle Specifies the rotation angle (in radians).
*/
template <class M, template <typename> class V, typename T>
void FreeRotation(M& mat, const V<T>& axis, const T& angle)
{
    __GS_ASSERT_MxN_MATRIX__("free rotation", M, 3, 3);

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


} // /namespace Gs


#endif



// ================================================================================
