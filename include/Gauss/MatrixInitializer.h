/*
 * MatrixInitializer.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_MATRIX_INITIALIZER_H__
#define __GS_MATRIX_INITIALIZER_H__


#include <cstddef>


namespace Gs
{


//! Common matrix initializer.
template <class M, typename T, std::size_t Cols>
class MatrixInitializer
{
            
    public:

        MatrixInitializer(M& matrix) :
            matrix_ ( matrix ),
            element_( 0      )
        {
        }

        MatrixInitializer<M, T, Cols>& operator , (const T& nextValue)
        {
            matrix_(element_ / Cols, element_ % Cols) = nextValue;
            ++element_;
            return *this;
        }

    private:

        M&          matrix_;
        std::size_t element_;

};


} // /namespace Gs


#endif



// ================================================================================
