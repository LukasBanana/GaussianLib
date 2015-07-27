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

//! Enables the swizzle operator in the vector classes.
#define GS_ENABLE_SWIZZLE_OPERATOR

//! Disables automatic data initialization.
//#define GS_DISABLE_AUTO_INIT

//! Enables column-major storage. If disabled, column-major storage is used.
//#define GS_COLUMN_MAJOR_STORAGE

//! Enables column vectors. If disabled, row vectors are used.
#define GS_COLUMN_VECTORS


#endif



// ================================================================================
