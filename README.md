GaussianLib - A basic linear algebra C++ library for 2D and 3D applications
===========================================================================

License
-------

[3-Clause BSD License](https://github.com/LukasBanana/GaussianLib/blob/master/LICENSE.txt)

Getting Started:
----------------
[Getting Started with GaussianLib.pdf](https://github.com/LukasBanana/GaussianLib/blob/master2/docu/GettingStarted/Getting%20Started%20with%20GaussianLib.pdf)

Example
-------

```cpp
// Optional macro to switch between column- or row major matrix storage layout:
// #define GS_ROW_MAJOR_STORAGE

// Optional macro to switch between column- or row vectors:
// #define GS_ROW_VECTORS

#include <Gauss/Gauss.h>     // include gaussian lib main header
#include <Gauss/DefConsts.h> // implement constant definitions
#include <iostream>

static const Gs::Real pi = Gs::Real(3.141592654);

int main()
{
    // Initialize 4 dimensional vectors a and b
    Gs::Vector4 a(1, 2, 3, 4), b(-12, 0.5f, 0, 1);
    const Gs::Vector2 c(42, 19);
    
    // Simple vector access
    c.x = a.x;
    a.y = c[1]; // equivalent to a.y = c.y
    a[2] += a.z // equivalent to a.z += a.z
    
    // 'Swizzle operator' like functionality
    Gs::Vector3 d = a.xyw()*2.0f + b.zxy() - c.xxy();

    // Declare 3x4 matrix A and 4x3 matrix B
    Gs::Matrix<double, 3, 4> A;
    Gs::Matrix<double, 4, 3> B;
    
    // Initialize 3x4 matrix A
    A << 1, 2, 0, -12,
         0, 0, 1, 4,
         0, 1, 0, 5;
         
    // Initialize 4x3 matrix B with the transposed matrix of A
    B = A.Transposed();
    
    /*
                        | / b11 b12 b13 \
                        | | b21 b22 b23 |
              x         | | b31 b32 b33 |
                        | \ b41 b42 b43 /
    --------------------|------------------
    / a11 a12 a13 a14 \ | / c11 c12 c13 \
    | a21 a22 a23 a24 | | | c21 c22 c23 |
    \ a31 a32 a33 a34 / | \ c31 c32 c33 /
    */
    Gs::Matrix3<double> C = A * B;
    
    // Invert matrix C
    C.MakeInverse();
    
    // Declare affine 4x4 matrix (only stores 3x4 elements,
    // or 4x3 elements whether GS_ROW_VECTORS is defined or not).
    // This requires less storage and most functions (such as "Inverse")
    // are much faster than with a common 4x4 matrix.
    Gs::AffineMatrix4 D = Gs::AffineMatrix4::Identity();
    
    // Set some transformations for the affine matrix
    D.SetPosition(Gs::Vector3(1, -2, 5));
    D.RotateX(pi*0.5);
    D.RotateZ(-pi*0.25);
    D.Scale(Gs::Vector3(1, 2, 3));
    
    // Declare quaternions
    Gs::Quaternion q0, q1;
    q0 = Gs::Quaternion::EulerAngles(Gs::Vector3(pi*-0.25, pi*0.8, 0));
    q1.SetEulerAngles(Gs::Vector3(pi*0.5, 0, 0));
    
    // Spherical-linear-interpolation (Slerp) with two quaternions
    Gs::Quaternion p = Slerp(q0, q1, 0.5);
    
    // Print matrices to standard output
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "C = " << std::endl << C << std::endl;
    
    // Print vectors to standard output
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a*b = " << a*b << std::endl;
    std::cout << "a . b = " << Dot(a, b) << std::endl;
    std::cout << "a x b = " << Cross(a, b) << std::endl;
    std::cout << "|a| = " << a.Length() << std::endl;
    std::cout << "a / |a| = " << a.Normalized() << std::endl;
    
    return 0;
}
```


