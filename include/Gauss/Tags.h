/*
 * Tags.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_TAGS_H__
#define __GS_TAGS_H__


namespace Gs
{


/**
\brief Common uninitialize tag.
\remarks This can be used to explicitly construct an uninitialized vector, matrix or quaternion:
\code
Gs::Matrix4 m(Gs::UninitializeTag{});
// ...
//m(0, 1) = ...
\endcode
*/
struct UninitializeTag {};


} // /namespace Gs


#endif



// ================================================================================
