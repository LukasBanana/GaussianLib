/*
 * Config.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_CONFIG_H__
#define __GS_CONFIG_H__


#ifdef _DEBUG
#   define GS_ENABLE_ASSERT
#endif

//! Enables the swizzle operator in the vector classes. If undefined, swizzle operator is disabled (default).
#define GS_ENABLE_SWIZZLE_OPERATOR

//! Disables automatic data initialization. If undefined, automatic initialization is enabled (default). 
//#define GS_DISABLE_AUTO_INIT

//! Enables column-major storage. If undefined, column-major storage is used (default).
//#define GS_COLUMN_MAJOR_STORAGE

//! Enables row vectors. If undefined, column vectors are used (default).
//#define GS_ROW_VECTORS //!UNUSED!


#endif



// ================================================================================
