/*
 * Config.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_CONFIG_H
#define GS_CONFIG_H


#ifdef _DEBUG
#   define GS_ENABLE_ASSERT 1
#endif

#ifndef GS_ENABLE_ASSERT
#define GS_ENABLE_ASSERT 0
#endif

//! Uses double precision floating-points for real types.
#ifndef GS_REAL_DOUBLE
#define GS_REAL_DOUBLE 0
#endif

//! Uses a runtime exception instead of an assert.
#ifndef GS_ASSERT_EXCEPTION
#define GS_ASSERT_EXCEPTION 0
#endif

//! Enables the swizzle operator in the vector classes. If undefined, swizzle operator is disabled (default).
#ifndef GS_ENABLE_SWIZZLE_OPERATOR
#define GS_ENABLE_SWIZZLE_OPERATOR 0
#endif

//! Enables the inverse matrix operator, i.e. allows "A^-1" expressions as shortcut for "A.Inverse()".
#ifndef GS_ENABLE_INVERSE_OPERATOR
#define GS_ENABLE_INVERSE_OPERATOR 0
#endif

//! Disables automatic data initialization. If undefined, automatic initialization is enabled (default).
#ifndef GS_DISABLE_AUTO_INIT
#define GS_DISABLE_AUTO_INIT 0
#endif

//! Enables row-major storage. If undefined, column-major storage is used (default).
#ifndef GS_ROW_MAJOR_STORAGE
#define GS_ROW_MAJOR_STORAGE 0
#endif

//! Enables row vectors. If undefined, column vectors are used (default).
#ifndef GS_ROW_VECTORS
#define GS_ROW_VECTORS 0
#endif


#endif



// ================================================================================
