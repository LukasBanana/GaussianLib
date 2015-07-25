/*
 * GLSLTypes.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_GLSL_TYPES_H__
#define __GS_GLSL_TYPES_H__


#include <Gauss/Gauss.h>


namespace Gs
{


/* --- Type Alias --- */

using vec2      = Vector2f;
using vec3      = Vector3f;
using vec4      = Vector4f;

using dvec2     = Vector2d;
using dvec3     = Vector3d;
using dvec4     = Vector4d;

using ivec2     = Vector2i;
using ivec3     = Vector3i;
using ivec4     = Vector4i;

using uvec2     = Vector2ui;
using uvec3     = Vector3ui;
using uvec4     = Vector4ui;

using bvec2     = Vector2T<bool>;
using bvec3     = Vector3T<bool>;
using bvec4     = Vector4T<bool>;

using mat2      = Matrix2f;
using mat2x3    = Matrix<float, 2, 3>;
using mat2x4    = Matrix<float, 2, 4>;

using mat3x2    = Matrix<float, 3, 2>;
using mat3      = Matrix3f;
using mat3x4    = Matrix<float, 3, 4>;

using mat4x2    = Matrix<float, 4, 2>;
using mat4x3    = Matrix<float, 4, 3>;
using mat4      = Matrix4f;

using dmat2     = Matrix2d;
using dmat2x3   = Matrix<double, 2, 3>;
using dmat2x4   = Matrix<double, 2, 4>;

using dmat3x2   = Matrix<double, 3, 2>;
using dmat3     = Matrix3d;
using dmat3x4   = Matrix<double, 3, 4>;

using dmat4x2   = Matrix<double, 4, 2>;
using dmat4x3   = Matrix<double, 4, 3>;
using dmat4     = Matrix4d;


} // /namespace Gs


#endif



// ================================================================================
