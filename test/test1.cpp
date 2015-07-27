/*
 * test1.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

//#define GS_COLUMN_MAJOR_STORAGE

#include <Gauss/Gauss.h>
#include <Gauss/HLSLTypes.h>
#include <Gauss/GLSLTypes.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <complex>


using namespace Gs;


static const Real pi = Real(3.141592654);

static void commonTest1()
{
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

    #if 1

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

    #if 1

    std::cout << "m2x2 = " << std::endl << m2x2 << std::endl;
    std::cout << "Inverse(m2x2) = " << std::endl << m2x2.Inverse() << std::endl;
    std::cout << "Determinant(m2x2) = " << m2x2.Determinant() << std::endl;

    #endif
    
    #if 1

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "Inverse(A) = " << std::endl << A.Inverse() << std::endl;
    std::cout << "A*Inverse(A) = " << std::endl << A*A.Inverse() << std::endl;
    std::cout << "Inverse(A)*A = " << std::endl << A.Inverse()*A << std::endl;

    #endif

    #if 1

    std::cout << "hugeMatrix = " << std::endl << hugeMatrix << std::endl;
    std::cout << "Determinant(hugeMatrix) = " << hugeMatrix.Determinant() << std::endl;

    #endif
}

static void sparesMatrixTest1()
{
    Matrix4 A;
    A << 1, 0, -2, 6,
         0, 8, 0, -4,
         0, 1, 2, 0,
         0, 0, 0, 1;

    SparseMatrix4 B;
    B << 1, 0, -2, 6,
         0, 8, 0, -4,
         0, 1, 2, 0;

    auto B2 = B;
    B2.RotateFree(Vector4(1, 1, 1).Normalized(), pi*0.5f);
    //B.MakeInverse();

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "Inv(A) = " << std::endl << A.Inverse() << std::endl;
    std::cout << "A*Inv(A) = " << std::endl << A*A.Inverse() << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "Inv(B) = " << std::endl << B.Inverse() << std::endl;
    std::cout << "B*Inv(B) = " << std::endl << B*B.Inverse() << std::endl;
    std::cout << "| A | = " << A.Determinant() << std::endl;
    std::cout << "| B | = " << B.Determinant() << std::endl;
    std::cout << "Trace A = " << A.Trace() << std::endl;
    std::cout << "Trace B = " << B.Trace() << std::endl;
    std::cout << std::endl << "B1 = " << std::endl << B << std::endl;
    std::cout << std::endl << "B2 = " << std::endl << B2 << std::endl;
    std::cout << std::endl << "Lerp(B1, B2, 0.5) = " << std::endl << Lerp(B, B2, 0.5f) << std::endl;
    std::cout << std::endl << "5 * I3 = " << std::endl << 5.0f * Matrix3::Identity() << std::endl;
}

static void quaternionTest1()
{
    Quaternion q0, q1;

    q0.SetEulerAngles(Vector3(pi*0.5f, 0, 0));
    q1.SetEulerAngles(Vector3(pi*1.0f, 0, 0));
    
    std::cout << "q0 = " << q0 << std::endl;
    std::cout << "q1 = " << q1 << std::endl;

    /*for (int i = 0; i <= 10; ++i)
    {
        auto t = static_cast<Real>(i) / 10;
        std::cout << "Slerp(" << t << ") = " << Slerp(q0, q1, t) << std::endl;
    }*/

    Matrix3 m = Matrix3::Identity();
    m.RotateFree(Vector3(1, 0, 1), pi*0.5f);

    std::cout << "m = " << std::endl << m << std::endl;
    std::cout << "Quaterion(m) = " << Quaternion(m) << std::endl;
}

static void matrixVectorTest1()
{
    Matrix4 A;
    SparseMatrix4 B;

    A << 1, 0, 0, 4,
         0, 1, 0, 2,
         0, 0, 1, -5,
         0, 0, 0, 1;
    
    B << 1, 0, 0, 4,
         0, 1, 0, 2,
         0, 0, 1, -5;

    auto a = TransformVector(A, Vector4(0, 0, 0, 1));
    auto b = TransformVector(B, Vector4(0, 0, 0, 1));

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
}

static void complexTest1()
{
    using complex = std::complex<Real>;

    Vector4T<complex> a(complex(0, 0), complex(0, 1), complex(4, -2));
    Matrix4T<complex> A = Matrix4T<complex>::Identity();

    auto b = A * a;

}

int main()
{
    std::cout << "GaussianLib Test 1" << std::endl;
    std::cout << "==================" << std::endl;

    //commonTest1();
    //sparesMatrixTest1();
    //quaternionTest1();
    //matrixVectorTest1();
    complexTest1();

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

