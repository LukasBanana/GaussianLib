/*
 * test1.cpp
 * 
 * This file is part of the "GaussianLib" project (Copyright (c) 2015 by Lukas Hermanns)
 * See "LICENSE.txt" for license information.
 */

#include "test1.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <complex>
#include <ctime>

#ifdef _WIN32
#   define NOMINMAX
#   include <Windows.h>
#   undef near
#   undef far
#endif


#ifdef _MSC_VER
#pragma warning (disable:4101) // "unreferenced local variable"
#endif

using namespace Gs;


#ifdef _WIN32

Timer::Timer()
{
    QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&startTime_));
}

Timer::~Timer()
{
    __int64 endTime, frequency;
    QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&endTime));
    QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&frequency));
    printf("%0.1f ms\n", static_cast<float>(endTime - startTime_) / static_cast<float>(frequency) * 1000.0f);
}

#endif


void commonTest1()
{
    // --- vector tests ---

    const Vector4 a(1, 2, 3, 4), b(-4, 0, 2, 1);

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

    Matrix<Real, 3, 4> C { UninitializeTag{} };
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

    std::cout << "A * a = " << (A * a) << std::endl;
    std::cout << "a * A = " << (a * A) << std::endl;
    std::cout << "C * b = " << (C * b) << std::endl;
    std::cout << "b * C^T = " << (b * C.Transposed()) << std::endl;
    std::cout << "Vector3(a) * C = " << (Vector3(a) * C) << std::endl;

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

void affineMatrixTest1()
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
    std::cout << "Inv(A) = " << std::endl << (A^-1) << std::endl;
    std::cout << "A*Inv(A) = " << std::endl << A*(A^-1) << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
    std::cout << "B^T = " << std::endl << B.Transposed() << std::endl;
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

void affineMatrixTest2()
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

void quaternionTest1()
{
    Quaternion q0, q1;

    q0.SetEulerAngles(Vector3(pi*0.5f, 0, 0));
    q1.SetEulerAngles(Vector3(pi*1.0f, 0, 0));
    
    std::cout << "q0 = " << q0 << std::endl;
    std::cout << "q1 = " << q1 << std::endl;

    for (int i = 0; i <= 10; ++i)
    {
        auto t = static_cast<Real>(i) / 10;
        std::cout << "Slerp(" << t << ") = " << Slerp(q0, q1, t) << std::endl;
    }

    Matrix3 m = Matrix3::Identity();
    RotateFree(m, Vector3(1, 0, 1), pi*0.5f);

    std::cout << "m = " << std::endl << m << std::endl;
    std::cout << "Quaterion(m) = " << Quaternion(m) << std::endl;
}

void quaternionTest2()
{
    Quaternion q0, q1;

    q0.SetEulerAngles(Vector3(pi*0.5f, pi*0.25f, 0));
    q1.SetEulerAngles(Vector3(pi*1.0f, 0, 0));
    
    std::cout << "q0       = " << q0 << std::endl;
    std::cout << "q1       = " << q1 << std::endl;

    std::cout << "q0 * q0  = " << q0*q0 << std::endl;
    std::cout << "Dot(q0, q0)  = " << Dot(q0, q0) << std::endl;

    q0 *= q0;
    std::cout << "q0 *= q0 = " << q0 << std::endl;
}

void matrixVectorTest1()
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

void complexTest1()
{
    using complex = std::complex<Real>;

    Vector4T<complex> a(complex(0, 0), complex(0, 1), complex(4, -2), complex());
    Matrix4T<complex> A = Matrix4T<complex>::Identity();

    auto b = A * a;
}

void projectionTest1()
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
    std::cout << "Project(R, a) = ";
    std::cout << (R * a).xy() << std::endl;
}

void equalsTest1()
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

void vectorTest1()
{
    auto x = Vector<Real, 10>();
    auto A = Matrix<Real, 10, 10>::Identity();
    
    x[0] = 12;
    x[1] = 5;
    x[2] = 3;
    x[3] = -9;
    x[4] = 16;
    x[5] = -4;

    A(0, 1) = 4;
    A(5, 3) = 18;
    A(2, 7) = -5;
    A(4, 4) = 6;
    A(8, 8) = -2;

    std::cout << "x = " << x << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "A*x = " << A*x << std::endl;
    std::cout << "Trace(A) = " << A.Trace() << std::endl;
}

void epsilonTest1()
{
    auto eps1 = Epsilon<float>();
    auto eps2 = Epsilon<double>();
    float x = 0.1f;

    std::cout << "epsilon<float> = " << eps1 << std::endl;
    std::cout << "epsilon<double> = " << eps2 << std::endl;
    std::cout << x << ((x < eps1) ? " < " : " >= ") << eps1 << std::endl;
    std::cout << "pi = " << pi << std::endl;
}

void sphericalTest1()
{
    auto ShowCoords = [](const Vector3& v)
    {
        std::cout << "cartesian coordinate a = " << v << std::endl;
        std::cout << "spherical coordinate b = " << Spherical(v) << std::endl;
        std::cout << "cartesian coordinate b = " << Vector3(Spherical(v)) << std::endl << std::endl;
    };

    ShowCoords({ 1, 2, 3 });
    ShowCoords({ -4, 0, 8 });
    ShowCoords({ 0, 0, 0 });
    ShowCoords({ 0.01f, -0.08f, 0.0005f });
    ShowCoords({ -3948, 4933, -239382 });
}

void crossProductTest1()
{
    Vector3 a(1, 0, 0);
    Vector3 b(0, 1, 0);
    b.Normalize();

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "Cross(a, b) = " << Cross(a, b).Normalized() << std::endl;
}

void rotateVectorTest1()
{
    Vector3 vec(1, 2, 3);
    Vector3 axis(0, 1, 0);

    Real angle = Deg2Rad(45.0f);
    auto x = RotateVectorAroundAxis(vec, axis, angle);

    std::cout << "vec = " << vec << ", axis = " << axis << ", angle = " << angle << ", x = " << x << std::endl;
}

void sortingTest1()
{
    Vector3 a(1, 2, 3);
    Vector3 b(2, 3, 4);

    Matrix3 A, B;
    A.LoadIdentity();
    B.LoadIdentity();
    A(2, 1) = 5;

    std::cout << std::boolalpha;
    std::cout << "a = " << a << ", b = " << b << ", a < b: " << Compare(a, b) << std::endl;
    std::cout << "A = " << std::endl << A << "B = " << std::endl << B << "A < B: " << Compare(A, B) << std::endl;
}

void flipTest1()
{
    Matrix4 A;
    A <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16;

    AffineMatrix4 B;
    #ifdef GS_ROW_VECTORS
    B <<  1,  2,  3,
          4,  5,  6, 
          7,  8,  9,
         10, 11, 12;
    #else
    B <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12;
    #endif

    Gs::FlipAxis(A, 1);
    Gs::FlipAxis(B, 1);

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;
}

void rotateMatrixTest1()
{
    AffineMatrix4 A;
    A << 1, 0, 0,   7,
         0, 1, 0,   3,
         0, 0, 1, -12;

    A.RotateX(45*Gs::pi/180);

    std::cout << "A = " << std::endl << A << std::endl;
}

void rcpTest1()
{
	const Real x = Gs::pi;
	const Vector3 v(1, 2, 3);
	const Matrix4 A;

	std::cout << "Rcp(" << x << ") = " << Rcp(x) << std::endl;
	std::cout << "Rcp(" << v << ") = " << Rcp(v) << std::endl;
	std::cout << "Rcp(" << std::endl << A << ") = " << std::endl << Gs::Rcp(A) << std::endl;
}

void typeTest1()
{
    typename ScalarType<Vector3f>::Type a;
    typename ScalarType<Vector2T<double>>::Type b;
    typename ScalarType<Vector<int, 5>>::Type c;
    typename ScalarType<Matrix<int, 3, 2>>::Type d;
    typename ScalarType<Matrix4d>::Type e;
    typename ScalarType<AffineMatrix4f>::Type f;
}

void stdMathTest1()
{
    Vector3f a(0.5f, 1.0f, 1.5f), b(1, 2, 3);
    Matrix2f A; A << 1, 2, 3, 4;

    std::cout << "sin(" << a << ") = " << sin(a) << std::endl;
    std::cout << "acos(" << a << ") = " << acos(a) << std::endl;
    std::cout << "floor(" << a << ") = " << floor(a) << std::endl;
    std::cout << "ceil(" << a << ") = " << ceil(a) << std::endl;
    std::cout << "sin(" << std::endl << A << ") = " << std::endl << sin(A) << std::endl;
}

void sseVector4Test1()
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

void sseVector4Test2()
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

void vector3Test1()
{
    Vector3f v(2.5f, 3.4f, 9.5f);

    std::cout << "v                               = " << v << std::endl;
    std::cout << "1 / v                           = " << 1.0f / v << std::endl;
    std::cout << "v / v.x                         = " << v / v.x << std::endl;
    
    v /= v.x;
    Vector2f v2D(v);
    Vector4f v4D(v2D, v2D);

    std::cout << "v /= v.x                        = " << v << std::endl;
    std::cout << "Vector2(v)                      = " << v2D << std::endl;
    std::cout << "Vector4(Vector2(v), Vector2(v)) = " << v4D << std::endl;
}


// Quick simple RNG based on Xorhash
class RandomNumberGenerator
{

    private:

        std::uint32_t state_ = 0;

    public:

        RandomNumberGenerator()
        {
            seed(47);
        }

        void seed(std::uint32_t seed)
        {
            // Thomas Wang's integer hash, as reported by Bob Jenkins
            seed = (seed ^ 61) ^ (seed >> 16);
            seed *= 9;
            seed = seed ^ (seed >> 4);
            seed *= 0x27d4eb2d;
            seed = seed ^ (seed >> 15);
            state_ = seed;
        }

        std::uint32_t randUInt()
        {
            // Xorshift algorithm from George Marsaglia's paper
            state_ ^= (state_ << 13);
            state_ ^= (state_ >> 17);
            state_ ^= (state_ << 5);
            return state_;
        }

        std::int32_t randInt()
        {
            return static_cast<std::int32_t>(randUInt());
        }

        float randFloat()
        {
            return static_cast<float>(randUInt()) * (1.0f / 4294967296.0f);
        }

        float randFloat(float min, float max)
        {
            return min + (min - max)*randFloat();
        }

};

void performanceTest1()
{
    #ifdef _WIN32

    std::cout << "performance test:" << std::endl;

    RandomNumberGenerator rng;
    rng.seed(static_cast<std::uint32_t>(time(0)));

    for (int i = 0; i < 10; ++i)
        std::cout << "random number: " << rng.randFloat() << std::endl;

    #else

    std::cout << "performance test only available on Win32" << std::endl;

    #endif
}

