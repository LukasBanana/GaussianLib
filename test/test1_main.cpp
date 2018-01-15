/*
 * test1_main.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#include "test1.h"

#include <iostream>
#include <cstdlib>

using namespace Gs;


void sseVector4Test1_2()
{
    Vector4f v0;
    Vector4f v1(3.0f, -0.5f, 0.0f, 1.0f);
    Vector4f v2(0.5f, 1.0f, 2.5f, 1.0f);

    std::cout << "v0          = " << v0 << std::endl;
    std::cout << "v1          = " << v1 << std::endl;
    std::cout << "v2          = " << v2 << std::endl;
    std::cout << "v1 + v2     = " << v1 + v2 << std::endl;
    std::cout << "v1 * v2     = " << v1 * v2 << std::endl;
    std::cout << "v1 * -5     = " << v1 * -5.0f << std::endl;
    std::cout << "Dot(v1, v2) = " << Dot(v1, v2) << std::endl;
    std::cout << "|v1|        = " << v1.Length() << std::endl;
    std::cout << "v1 / |v1|   = " << v1.Normalized() << std::endl;
}

void sseVector4Test2_2()
{
    Vector4d v0;
    Vector4d v1(3.0, -0.5, 0.0, 1.0);
    Vector4d v2(0.5, 1.0, 2.5, 1.0);

    std::cout << "v0          = " << v0 << std::endl;
    std::cout << "v1          = " << v1 << std::endl;
    std::cout << "v2          = " << v2 << std::endl;
    std::cout << "v1 + v2     = " << v1 + v2 << std::endl;
    std::cout << "v1 * v2     = " << v1 * v2 << std::endl;
    std::cout << "v1 * -5     = " << v1 * -5.0 << std::endl;
    std::cout << "Dot(v1, v2) = " << Dot(v1, v2) << std::endl;
    std::cout << "|v1|        = " << v1.Length() << std::endl;
    std::cout << "v1 / |v1|   = " << v1.Normalized() << std::endl;
}

int main()
{
    std::cout << "GaussianLib Test 1" << std::endl;
    std::cout << "==================" << std::endl;

    try
    {
        /*commonTest1();
        affineMatrixTest1();
        affineMatrixTest2();
        quaternionTest1();
        quaternionTest2();
        matrixVectorTest1();
        complexTest1();
        projectionTest1();
        equalsTest1();
        vectorTest1();
        epsilonTest1();
        sphericalTest1();
        crossProductTest1();
        rotateVectorTest1();
        sortingTest1();
        flipTest1();
        rotateMatrixTest1();
		rcpTest1();
        stdMathTest1();
        vector3Test1();*/
        sseVector4Test1();
        sseVector4Test1_2();
        sseVector4Test2();
        sseVector4Test2_2();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

