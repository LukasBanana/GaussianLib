GaussianLib - A basic linear algebra C++ library for 2D and 3D applications
===========================================================================

License
-------

[3-Clause BSD License](https://github.com/LukasBanana/GaussianLib/blob/master/LICENSE.txt)

Status
------

**Alpha**

Example
-------

```cpp
#include <Gauss/Gauss.h>
#include <iostream>

int main()
{
    // Initialize 4 dimensional vectors a and b
    Gs::Vector4 a(1, 2, 3, 4), b(-12, 0.5f, 0, 1);
    
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
    C.Invert();
    
    // Print matrices to standard output
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "C = " << std::endl << C << std::endl;
    
    // Print vectors to standard output
    std::cout << "a = " << a << ", b = " << b << ", a*b = " << a*b << std::endl;
    std::cout << "a . b = " << Dot(a, b) << "a x b = " << Cross(a, b) << std::endl;
    std::cout << "|| a || = " << a.Length() << ", a / || a || = " << a.Normalized() << std::endl;
    
    return 0;
}
```


