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
    A(0, 2) = 5;
    A.Transpose();

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

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

