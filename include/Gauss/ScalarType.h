/*
 * Typename.h
 *
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#ifndef GS_TYPENAME_H
#define GS_TYPENAME_H


#include <Gauss/Decl.h>
#include <cstddef>


namespace Gs
{


//! Provides the scalar type of the scalar, vector, or matrix type specified by 'T'.
template <typename T>
struct VectorType
{
    using ScalarType = T;
    static constexpr std::size_t elements = T::elements;
};

template <typename T, std::size_t N>
struct VectorType< T[N] >
{
    using ScalarType = T;
    static constexpr std::size_t elements = N;
};

template <typename T, std::size_t N>
struct VectorType< Vector<T, N> >
{
    using ScalarType = T;
    static constexpr std::size_t elements = N;
};

template <typename T, std::size_t Rows, std::size_t Cols>
struct VectorType< Matrix<T, Rows, Cols> >
{
    using ScalarType = T;
    static constexpr std::size_t elements = Matrix<T, Rows, Cols>::elements;
};

template <typename T>
struct VectorType<AffineMatrix3T<T>>
{
    using ScalarType = T;
    static constexpr std::size_t elements = AffineMatrix3T<T>::elements;
};

template <typename T>
struct VectorType< AffineMatrix4T<T> >
{
    using ScalarType = T;
    static constexpr std::size_t elements = AffineMatrix4T<T>::elements;
};

template <typename T>
struct VectorType< ProjectionMatrix4T<T> >
{
    using ScalarType = T;
    static constexpr std::size_t elements = ProjectionMatrix4T<T>::elements;
};

template <typename T>
struct VectorType< QuaternionT<T> >
{
    using ScalarType = T;
    static constexpr std::size_t elements = QuaternionT<T>::components;
};

template <typename T>
struct ScalarType
{
    using Type = typename VectorType<T>::ScalarType;
};


} // /namespace Gs


#endif



// ================================================================================
