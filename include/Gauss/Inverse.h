/*
 * Inverse.h
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef __GS_INVERSE_H__
#define __GS_INVERSE_H__


#include "Decl.h"
#include "Details.h"
#include "Determinant.h"


namespace Gs
{


//! Computes the inverse of the specified matrix 'm'.
template <typename T, std::size_t N>
bool Inverse(Matrix<T, N, N>& inv, const Matrix<T, N, N>& m)
{
    //using Helper = Details::MatrixHelper<M, T, Rows, Cols>;
    //return Helper::OrderedInverse(inv, Helper::MatrixToArray(m), Rows);
    return false;//!!!
}

//! Computes the inverse of the specified 2x2 matrix 'm'.
template <typename T>
bool Inverse(Matrix<T, 2, 2>& inv, const Matrix<T, 2, 2>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv.At(0, 0) = d * (  m.At(1, 1) );
    inv.At(0, 1) = d * ( -m.At(0, 1) );

    inv.At(1, 0) = d * ( -m.At(1, 0) );
    inv.At(1, 1) = d * (  m.At(0, 0) );

    return true;
}

//! Computes the inverse of the specified 3x3 matrix 'm'.
template <typename T>
bool Inverse(Matrix<T, 3, 3>& inv, const Matrix<T, 3, 3>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv.At(0, 0) = d * ( m.At(1, 1) * m.At(2, 2) - m.At(2, 1) * m.At(1, 2) );
    inv.At(1, 0) = d * ( m.At(1, 2) * m.At(2, 0) - m.At(1, 0) * m.At(2, 2) );
    inv.At(2, 0) = d * ( m.At(1, 0) * m.At(2, 1) - m.At(2, 0) * m.At(1, 1) );

    inv.At(0, 1) = d * ( m.At(0, 2) * m.At(2, 1) - m.At(0, 1) * m.At(2, 2) );
    inv.At(1, 1) = d * ( m.At(0, 0) * m.At(2, 2) - m.At(0, 2) * m.At(2, 0) );
    inv.At(2, 1) = d * ( m.At(2, 0) * m.At(0, 1) - m.At(0, 0) * m.At(2, 1) );

    inv.At(0, 2) = d * ( m.At(0, 1) * m.At(1, 2) - m.At(0, 2) * m.At(1, 1) );
    inv.At(1, 2) = d * ( m.At(1, 0) * m.At(0, 2) - m.At(0, 0) * m.At(1, 2) );
    inv.At(2, 2) = d * ( m.At(0, 0) * m.At(1, 1) - m.At(1, 0) * m.At(0, 1) );

    return true;
}

//! Computes the inverse of the specified 4x4 matrix 'm'.
template <typename T>
bool Inverse(Matrix<T, 4, 4>& inv, const Matrix<T, 4, 4>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    //INCOMPLETE
    /* Compute inverse matrix */
    inv.At(0, 0) = d * ( m.At(1, 1) * (m.At(2, 2) * m.At(3, 3) - m.At(3, 2) * m.At(2, 3)) + m.At(2, 1) * (m.At(3, 2) * m.At(1, 3) - m.At(1, 2) * m.At(3, 3)) + m.At(3, 1) * (m.At(1, 2) * m.At(2, 3) - m.At(2, 2) * m.At(1, 3)) );
    inv.At(1, 0) = d * ( m.At(1, 2) * (m.At(2, 0) * m.At(3, 3) - m.At(3, 0) * m.At(2, 3)) + m.At(2, 2) * (m.At(3, 0) * m.At(1, 3) - m.At(1, 0) * m.At(3, 3)) + m.At(3, 2) * (m.At(1, 0) * m.At(2, 3) - m.At(2, 0) * m.At(1, 3)) );
    inv.At(2, 0) = d * ( m.At(1, 3) * (m.At(2, 0) * m.At(3, 1) - m.At(3, 0) * m.At(2, 1)) + m.At(2, 3) * (m.At(3, 0) * m.At(1, 1) - m.At(1, 0) * m.At(3, 1)) + m.At(3, 3) * (m.At(1, 0) * m.At(2, 1) - m.At(2, 0) * m.At(1, 1)) );
    inv.At(3, 0) = d * ( m.At(1, 0) * (m.At(3, 1) * m.At(2, 2) - m.At(2, 1) * m.At(3, 2)) + m.At(2, 0) * (m.At(1, 1) * m.At(3, 2) - m.At(3, 1) * m.At(1, 2)) + m.At(3, 0) * (m.At(2, 1) * m.At(1, 2) - m.At(1, 1) * m.At(2, 2)) );
    
    inv.At(0, 1) = d * ( m.At(2, 1) * (m.At(0, 2) * m.At(3, 3) - m.At(3, 2) * m.At(0, 3)) + m.At(3, 1) * (m.At(2, 2) * m.At(0, 3) - m.At(0, 2) * m.At(2, 3)) + m.At(0, 1) * (m.At(3, 2) * m.At(2, 3) - m.At(2, 2) * m.At(3, 3)) );
    inv.At(1, 1) = d * ( m.At(2, 2) * (m.At(0, 0) * m.At(3, 3) - m.At(3, 0) * m.At(0, 3)) + m.At(3, 2) * (m.At(2, 0) * m.At(0, 3) - m.At(0, 0) * m.At(2, 3)) + m.At(0, 2) * (m.At(3, 0) * m.At(2, 3) - m.At(2, 0) * m.At(3, 3)) );
    inv.At(2, 1) = d * ( m.At(2, 3) * (m.At(0, 0) * m.At(3, 1) - m.At(3, 0) * m.At(0, 1)) + m.At(3, 3) * (m.At(2, 0) * m.At(0, 1) - m.At(0, 0) * m.At(2, 1)) + m.At(0, 3) * (m.At(3, 0) * m.At(2, 1) - m.At(2, 0) * m.At(3, 1)) );
    inv.At(3, 1) = d * ( m.At(2, 0) * (m.At(3, 1) * m.At(0, 2) - m.At(0, 1) * m.At(3, 2)) + m.At(3, 0) * (m.At(0, 1) * m.At(2, 2) - m.At(2, 1) * m.At(0, 2)) + m.At(0, 0) * (m.At(2, 1) * m.At(3, 2) - m.At(3, 1) * m.At(2, 2)) );
    
    inv.At(0, 2) = d * ( m.At(3, 1) * (m.At(0, 2) * m.At(1, 3) - m.At(1, 2) * m.At(0, 3)) + m.At(0, 1) * (m.At(1, 2) * m.At(3, 3) - m.At(3, 2) * m.At(1, 3)) + m.At(1, 1) * (m.At(3, 2) * m.At(0, 3) - m.At(0, 2) * m.At(3, 3)) );
    inv.At(1, 2) = d * ( m.At(3, 2) * (m.At(0, 0) * m.At(1, 3) - m.At(1, 0) * m.At(0, 3)) + m.At(0, 2) * (m.At(1, 0) * m.At(3, 3) - m.At(3, 0) * m.At(1, 3)) + m.At(1, 2) * (m.At(3, 0) * m.At(0, 3) - m.At(0, 0) * m.At(3, 3)) );
    inv.At(2, 2) = d * ( m.At(3, 3) * (m.At(0, 0) * m.At(1, 1) - m.At(1, 0) * m.At(0, 1)) + m.At(0, 3) * (m.At(1, 0) * m.At(3, 1) - m.At(3, 0) * m.At(1, 1)) + m.At(1, 3) * (m.At(3, 0) * m.At(0, 1) - m.At(0, 0) * m.At(3, 1)) );
    inv.At(3, 2) = d * ( m.At(3, 0) * (m.At(1, 1) * m.At(0, 2) - m.At(0, 1) * m.At(1, 2)) + m.At(0, 0) * (m.At(3, 1) * m.At(1, 2) - m.At(1, 1) * m.At(3, 2)) + m.At(1, 0) * (m.At(0, 1) * m.At(3, 2) - m.At(3, 1) * m.At(0, 2)) );

    inv.At(0, 3) = d * ( m.At(0, 1) * (m.At(2, 2) * m.At(1, 3) - m.At(1, 2) * m.At(2, 3)) + m.At(1, 1) * (m.At(0, 2) * m.At(2, 3) - m.At(2, 2) * m.At(0, 3)) + m.At(2, 1) * (m.At(1, 2) * m.At(0, 3) - m.At(0, 2) * m.At(1, 3)) );
    inv.At(1, 3) = d * ( m.At(0, 2) * (m.At(2, 0) * m.At(1, 3) - m.At(1, 0) * m.At(2, 3)) + m.At(1, 2) * (m.At(0, 0) * m.At(2, 3) - m.At(2, 0) * m.At(0, 3)) + m.At(2, 2) * (m.At(1, 0) * m.At(0, 3) - m.At(0, 0) * m.At(1, 3)) );
    inv.At(2, 3) = d * ( m.At(0, 3) * (m.At(2, 0) * m.At(1, 1) - m.At(1, 0) * m.At(2, 1)) + m.At(1, 3) * (m.At(0, 0) * m.At(2, 1) - m.At(2, 0) * m.At(0, 1)) + m.At(2, 3) * (m.At(1, 0) * m.At(0, 1) - m.At(0, 0) * m.At(1, 1)) );
    inv.At(3, 3) = d * ( m.At(0, 0) * (m.At(1, 1) * m.At(2, 2) - m.At(2, 1) * m.At(1, 2)) + m.At(1, 0) * (m.At(2, 1) * m.At(0, 2) - m.At(0, 1) * m.At(2, 2)) + m.At(2, 0) * (m.At(0, 1) * m.At(1, 2) - m.At(1, 1) * m.At(0, 2)) );

    return true;
}

//! Computes the inverse of the specified affine 3x3 matrix 'm'.
template <typename T> bool Inverse(AffineMatrix3T<T>& inv, const AffineMatrix3T<T>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv.At(0, 0) = d * (  m.At(1, 1) );
    inv.At(1, 0) = d * ( -m.At(1, 0) );
  /*inv.At(2, 0) = 0*/;

    inv.At(0, 1) = d * ( -m.At(0, 1) );
    inv.At(1, 1) = d * (  m.At(0, 0) );
  /*inv.At(2, 1) = 0*/;

    inv.At(0, 2) = d * ( m.At(0, 1) * m.At(1, 2) - m.At(0, 2) * m.At(1, 1) );
    inv.At(1, 2) = d * ( m.At(1, 0) * m.At(0, 2) - m.At(0, 0) * m.At(1, 2) );
  /*inv.At(2, 2) = 1*/;

    return true;
}

//! Computes the inverse of the specified affine 4x4 matrix 'm'.
template <typename T> bool Inverse(AffineMatrix4T<T>& inv, const AffineMatrix4T<T>& m)
{
    /* Compute inverse determinant */
    T d = Determinant(m);

    if (d == T(0))
        return false;

    d = T(1) / d;

    /* Compute inverse matrix */
    inv.At(0, 0) = d * ( m.At(1, 1) * m.At(2, 2) + m.At(2, 1) * ( -m.At(1, 2) ) );
    inv.At(1, 0) = d * ( m.At(1, 2) * m.At(2, 0) + m.At(2, 2) * ( -m.At(1, 0) ) );
    inv.At(2, 0) = d * ( m.At(1, 0) * m.At(2, 1) - m.At(2, 0) * m.At(1, 1) );
  /*inv.At(3, 0) = 0;*/
    inv.At(0, 1) = d * ( m.At(2, 1) * m.At(0, 2) + m.At(0, 1) * ( -m.At(2, 2) ) );
    inv.At(1, 1) = d * ( m.At(2, 2) * m.At(0, 0) + m.At(0, 2) * ( -m.At(2, 0) ) );
    inv.At(2, 1) = d * ( m.At(2, 0) * m.At(0, 1) - m.At(0, 0) * m.At(2, 1) );
  /*inv.At(3, 1) = 0;*/
    inv.At(0, 2) = d * ( m.At(0, 1) * m.At(1, 2) + m.At(1, 1) * ( -m.At(0, 2) ) );
    inv.At(1, 2) = d * ( m.At(0, 2) * m.At(1, 0) + m.At(1, 2) * ( -m.At(0, 0) ) );
    inv.At(2, 2) = d * ( m.At(0, 0) * m.At(1, 1) - m.At(1, 0) * m.At(0, 1) );
  /*inv.At(3, 2) = 0;*/
    inv.At(0, 3) = d * ( m.At(0, 1) * (m.At(2, 2) * m.At(1, 3) - m.At(1, 2) * m.At(2, 3)) + m.At(1, 1) * (m.At(0, 2) * m.At(2, 3) - m.At(2, 2) * m.At(0, 3)) + m.At(2, 1) * (m.At(1, 2) * m.At(0, 3) - m.At(0, 2) * m.At(1, 3)) );
    inv.At(1, 3) = d * ( m.At(0, 2) * (m.At(2, 0) * m.At(1, 3) - m.At(1, 0) * m.At(2, 3)) + m.At(1, 2) * (m.At(0, 0) * m.At(2, 3) - m.At(2, 0) * m.At(0, 3)) + m.At(2, 2) * (m.At(1, 0) * m.At(0, 3) - m.At(0, 0) * m.At(1, 3)) );
    inv.At(2, 3) = d * ( m.At(0, 3) * (m.At(2, 0) * m.At(1, 1) - m.At(1, 0) * m.At(2, 1)) + m.At(1, 3) * (m.At(0, 0) * m.At(2, 1) - m.At(2, 0) * m.At(0, 1)) + m.At(2, 3) * (m.At(1, 0) * m.At(0, 1) - m.At(0, 0) * m.At(1, 1)) );
  /*inv.At(3, 3) = 1;*/

    return true;
}


} // /namespace Gs


#endif



// ================================================================================
