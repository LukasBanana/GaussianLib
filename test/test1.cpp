/*
 * test1.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

//#define GS_ROW_MAJOR_STORAGE
#define GS_ENABLE_SWIZZLE_OPERATOR
#define GS_HIGH_PRECISION_FLOAT
//#define GS_ROW_VECTORS

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
    Vector4 c = a.zzzw()*Real(2) + a.xyxy() - b.yxzw();
    Vector2 d = c.xx(), e = c.xw() + b.yz();
    #endif

    // --- quaternion tests ---

    Quaternion q0;
    q0.SetEulerAngles(Vector3(3.14f*0.5f, 0.0f, 0.0f));

    Vector3 eulerAngles;
    q0.GetEulerAngles(eulerAngles);

    Vector3 v0 = q0 * Vector3(1, 2, 3);
    Quaternion q1 = q0 * Real(3);

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
    hugeMatrix.Inverse();

    // --- affine matrix tests ---

    AffineMatrix4 As(
        1, 0, 0, 4,
        0, 1, 0, -2,
        0, 0, 1, 5
    );
    AffineMatrix4 Bs(
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

    #ifdef GS_ROW_VECTORS
    std::cout << "a * A = " << std::endl << (a * A) << std::endl;
    #else
    std::cout << "A * a = " << std::endl << (A * a) << std::endl;
    #endif

    #endif

    #if 0

    std::cout << "m2x2 = " << std::endl << m2x2 << std::endl;
    std::cout << "Inverse(m2x2) = " << std::endl << m2x2.Inverse() << std::endl;
    std::cout << "Determinant(m2x2) = " << m2x2.Determinant() << std::endl;

    #endif
    
    #if 0

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "Inverse(A) = " << std::endl << A.Inverse() << std::endl;
    std::cout << "A*Inverse(A) = " << std::endl << A*A.Inverse() << std::endl;
    std::cout << "Inverse(A)*A = " << std::endl << A.Inverse()*A << std::endl;

    #endif

    #if 0

    std::cout << "hugeMatrix = " << std::endl << hugeMatrix << std::endl;
    std::cout << "Determinant(hugeMatrix) = " << hugeMatrix.Determinant() << std::endl;

    #endif
}

static void affineMatrixTest1()
{
    Matrix4 A;
    A << 1, 0, -2, 6,
         0, 8, 0, -4,
         0, 1, 2, 0,
         0, 0, 0, 1;

    AffineMatrix4 B;
    B << 1, 0, -2, 6,
         0, 8, 0, -4,
         0, 1, 2, 0;

    auto A2 = A;
    RotateFree(A2, Vector3(1, 1, 1).Normalized(), pi*0.5f);
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
    std::cout << std::endl << "A1 = " << std::endl << A << std::endl;
    std::cout << std::endl << "A2 = " << std::endl << A2 << std::endl;
    std::cout << std::endl << "Lerp(A1, A2, 0.5) = " << std::endl << Lerp(A, A2, 0.5f) << std::endl;
    std::cout << std::endl << "5 * I3 = " << std::endl << Real(5) * Matrix3::Identity() << std::endl;
}

static void affineMatrixTest2()
{
    AffineMatrix3 A;

    A << 1, 0, -2,
         0, 8, 3;
    
    AffineMatrix4 B;
    B << 1, 7, 9, -6,
         2, -4, 0, 1,
         3, 1, 6, -2;

    std::cout << "AffineMatrix3:" << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "A^-1 = " << std::endl << A.Inverse() << std::endl;
    std::cout << "A*A^-1 = " << std::endl << A*A.Inverse() << std::endl;

    std::cout << "AffineMatrix4:" << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "B^-1 = " << std::endl << B.Inverse() << std::endl;
    std::cout << "B*B^-1 = " << std::endl << B*B.Inverse() << std::endl;
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
    RotateFree(m, Vector3(1, 0, 1), pi*0.5f);

    std::cout << "m = " << std::endl << m << std::endl;
    std::cout << "Quaterion(m) = " << Quaternion(m) << std::endl;
}

static void matrixVectorTest1()
{
    Matrix4 A = Matrix4::Identity();
    AffineMatrix4 B = AffineMatrix4::Identity();

    Translate(A, Vector3(4, 2, -5.0f/3.0f));
    Scale(B, Vector3(1, 2, 3));
    Translate(B, Vector3(4, 2, -5.0f/3.0f));

    #if 1
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "A^-1 = " << std::endl << A.Inverse() << std::endl;
    std::cout << "A*A^-1 = " << std::endl << A*A.Inverse() << std::endl;
    #endif

    auto a = TransformVector(A, Vector4(0, 0, 0, 1));
    auto b = TransformVector(B, Vector4(0, 0, 0, 1));

    //Translate(A, a.xyz()*2.0f);
    //Translate(B, a.yxz());

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
}

static void complexTest1()
{
    using complex = std::complex<Real>;

    Vector4T<complex> a(complex(0, 0), complex(0, 1), complex(4, -2));
    Matrix4T<complex> A = Matrix4T<complex>::Identity();

    #ifdef GS_ROW_VECTORS
    auto b = a * A;
    #else
    auto b = A * a;
    #endif
}

static void projectionTest1()
{
    const Real w = 800, h = 600, near = 1.0f, far = 100.0f, fov = 74.0f*pi/180.0f;

    auto P = ProjectionMatrix4::Perspective(w/h, near, far, fov);
    auto Q = ProjectionMatrix4::Orthogonal(w, h, near, far);
    auto R = ProjectionMatrix4::Planar(w, h);

    Vector4 a(50, 0, 0, 1);

    std::cout << "Perspective Projection P = " << std::endl << P << std::endl;
    std::cout << "Orthogonal  Projection Q = " << std::endl << Q << std::endl;
    std::cout << "Planar      Projection R = " << std::endl << R << std::endl;
    std::cout << "P*P^-1 = " << std::endl << P*P.Inverse() << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "Project(R, a) = " << (R * a).xy() << std::endl;
}

static void equalsTest1()
{
    Vector3 a(1, 2, 3), b(4, 5, 6), c(1, 2, 3);
    Real x = 1, y = 2, z = 1;

    auto YesNo = [](bool b)
    {
        return b ? "Yes" : "No";
    };

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "a equals b ? " << YesNo(Equals(a, b)) << std::endl;
    std::cout << "a equals c ? " << YesNo(Equals(a, c)) << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl;
    std::cout << "x equals y ? " << YesNo(Equals(x, y)) << std::endl;
    std::cout << "x equals z ? " << YesNo(Equals(x, z)) << std::endl;
}

int main()
{
    std::cout << "GaussianLib Test 1" << std::endl;
    std::cout << "==================" << std::endl;

    //commonTest1();
    //affineMatrixTest1();
    //affineMatrixTest2();
    //quaternionTest1();
    //matrixVectorTest1();
    //complexTest1();
    //projectionTest1();
    equalsTest1();

    #ifdef _WIN32
    system("pause");
    #endif

    return 0;
}

