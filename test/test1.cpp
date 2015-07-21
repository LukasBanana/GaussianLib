/*
 * test1.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#include <Gauss/Gauss.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>


using namespace Gs;


int main()
{
    std::cout << "GaussianLib Test 1" << std::endl;
    std::cout << "==================" << std::endl;

    // --- vector tests ---

    const Vector4 a(1, 2, 3, 4), b(-4, 0, 2);

    Vector4 c = a.zzzw() + a.xyxy() - b.yxzw();
    c.yxzw() = a;
    Vector2 d = c.xx(), e = c.xw() + b.yz();
    d.yx() += c.zw();

    // --- matrix tests ---

    Matrix4 A = Matrix4::Identity(), B = Matrix4::Identity();

    A << 1, -5, 0, 12.5f,
         0,  1, 6, -7.8f;
    B << 1, 0, 0, 0,
         0, 0, 1, 0,
         0, 1, 0, 0,
         0, 0, 0, 1;

    Matrix<Real, 3, 4> C;
    C << 1, 0, 0, 12,
         0, 1, 0, -4,
         0, 0, 1, 5;
    
    Matrix<Real, 4, 1> D;
    D << 4,
         2,
         0,
         1;

    auto E = C * D;


    // --- output ---

    std::cout << "C = " << std::endl << C << std::endl;
    std::cout << "D = " << std::endl << D << std::endl;
    std::cout << "E = " << std::endl << E << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "a / || a || = " << a.Normalized() << std::endl;
    std::cout << "|| a - b || = " << Distance(a, b) << std::endl;
    std::cout << "a ANGLE b = " << Angle(a, b) << std::endl;
    std::cout << "a.x = " << a[0] << ", a.y = " << a[1] << ", a.z = " << a[2] << std::endl;
    std::cout << "a DOT b = " << Dot(a, b) << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "A^T = " << std::endl << A.Transposed() << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "A * B = " << std::endl << (A * B) << std::endl;
    std::cout << "A * a = " << std::endl << (A * a) << std::endl;

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

