/*
 * test1.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

//#define GS_MATRIX_COLUMN_MAJOR

#include <Gauss/Gauss.h>
#include <Gauss/HLSLTypes.h>
#include <Gauss/GLSLTypes.h>

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

    #ifdef GS_ENABLE_SWIZZLE_OPERATOR
    Vector4 c = a.zzzw() + a.xyxy() - b.yxzw();
    c.yxzw() = a;
    Vector2 d = c.xx(), e = c.xw() + b.yz();
    d.yx() += c.zw();
    #endif

    // --- quaternion tests ---

    Quaternion q0;
    q0.SetEulerAngles(Vector3(3.14f*0.5f, 0.0f, 0.0f));

    Vector3 eulerAngles;
    q0.GetEulerAngles(eulerAngles);

    Vector3 v0 = q0 * Vector3(1, 2, 3);
    Quaternion q1 = q0 * 3.0f;

    // --- matrix tests ---

    Matrix4 A = Matrix4::Identity(), B = Matrix4::Identity();

    A << 1, -5, 0, 12.5f,
         0,  1, 6, -7.8f;
    B << 1, 0, 0, 0,
         0, 0, 1, 0,
         0, 1, 0, 0,
         0, 0, 0, 1;

    Matrix<Real, 3, 4> C(UninitializeTag{});
    C << 1, 0, 0, 12,
         0, 1, 0, -4,
         0, 0, 1, 5;
    
    Matrix<Real, 4, 1> D;
    D << 4,
         2,
         0,
         1;

    auto E = C * D;

    Matrix2 m2x2 = Matrix2::Identity();
    Matrix3 m3x3 = Matrix3::Identity();

    Matrix<float, 6, 6> hugeMatrix; hugeMatrix.LoadIdentity();

    // --- sparse matrix tests ---

    SparseMatrix4 As(
        1, 0, 0, 4,
        0, 1, 0, -2,
        0, 0, 1, 5
    );
    SparseMatrix4 Bs(
        1, 0, 0, 6,
        0, 0, 1, 0,
        0, 1, 0, 0
    );

    As = Bs * As;

    // --- output ---

    #if 0

    std::cout << "As = " << std::endl << As << std::endl;
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

    #endif

    std::cout << "m2x2 = " << std::endl << m2x2 << std::endl;
    std::cout << "Inverse(m2x2) = " << std::endl << m2x2.Inverse() << std::endl;
    std::cout << "Determinant(m2x2) = " << m2x2.Determinant() << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "Inverse(A) = " << std::endl << A.Inverse() << std::endl;
    std::cout << "hugeMatrix = " << std::endl << hugeMatrix << std::endl;
    std::cout << "Determinant(hugeMatrix) = " << hugeMatrix.Determinant() << std::endl;

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

