/*
 * test1_main.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#include "test1.h"

#include <iostream>
#include <cstdlib>


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
        sseVector4Test2();
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

