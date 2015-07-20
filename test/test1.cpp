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

    Vector4 a(1, 2, 3), b(-4, 0, 2);

    Matrix4 A, B;
    A.LoadIdentitiy();
    B.LoadIdentitiy();

    A << 1.0f, -5.0f, 0.0f, 12.5f,
         0.0f,  1.0f, 6.0f, -7.8f;
    B << 1, 0, 0, 0,
         0, 0, 1, 0,
         0, 1, 0, 0,
         0, 0, 0, 1;

    Matrix<Real, 4, 3> C;
    Matrix<Real, 3, 4> D;
    auto E = C * D;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "|| a || = " << a.Normalized() << std::endl;
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

