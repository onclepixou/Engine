#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "tennin.h"
#include "geometry.h"



//----------------------------------//
typedef struct Matrix4x4 Matrix4x4;
//----------------------------------//


//------------------------------------------------------------------------------------//
// Matrix4x4 //
struct Matrix4x4{
    Float m[4][4];

    // class constructor
    Matrix4x4();
    Matrix4x4(Float mat[4][4]);
    Matrix4x4(Float t00, Float t01, Float t02, Float t03,
              Float t10, Float t11, Float t12, Float t13,
              Float t20, Float t21, Float t22, Float t23,
              Float t30, Float t31, Float t32, Float t33);

    // operator overload
    bool operator==(const Matrix4x4& mat)const;
    bool operator!=(const Matrix4x4& mat)const;

    //  class method
    std::string print();
    bool isIdentity();
    static Matrix4x4 Mul(const Matrix4x4& m1, const Matrix4x4& m2);


};

// Matrix4x4 related function
    std::ostream& operator<<(std::ostream& o, Matrix4x4 m);
    Matrix4x4 Transpose(const Matrix4x4& mat);
    Matrix4x4 Inverse(const Matrix4x4& m);
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Transform //

class Transform{
    public: 

        // class constructor
        Transform();
        Transform(const Float mat[4][4]);
        Transform(const Matrix4x4& m);
        Transform(const Matrix4x4& m, const Matrix4x4& mInv);

        // operator overloading
        bool operator==(const Transform& t)const;
        bool operator!=(const Transform& t)const;
        Transform operator*(const Transform& t2)const;

        template <typename T>
        Point3<T> operator()(const Point3<T>& p) const;

        template <typename T>
        Point3<T> operator()(const Point3<T>& p, Vector3<T>* absError) const;

        template <typename T>
        Vector3<T> operator()(const Vector3<T>& v) const;

        template <typename T>
        Normal3<T> operator()(const Normal3<T>& n) const;

        Ray operator()(const Ray& r) const;
        RayDifferential operator()(const RayDifferential& rd) const;  
        Bounds3f operator()(const Bounds3f& b)const;


        // class method
        std::string print();
        bool isIdentity();
        bool hasScale()const;
        Matrix4x4 getM()const;
        Matrix4x4 getInv()const;
        bool swapsHandedness()const;

        // friend function
        friend Transform Inverse(const Transform& t);
        friend Transform Transpose(const Transform& t);


    private:
        Matrix4x4 m;
        Matrix4x4 mInv;
};

// Transform related function
//------------------------------------------------------------------------------------//

    template <typename T>
    Transform Translate(Vector3<T> v);

    template <typename T>
    Transform Scale(Vector3<T> v);

    Transform RotateX(Float theta);
    Transform RotateY(Float theta);
    Transform RotateZ(Float theta);
    Transform Rotate (Float theta,const Vect3f axis);

    Transform LookAt(const Point3f pos, const Point3f look, const Vect3f up);

 

    std::ostream& operator<<(std::ostream& o, Transform t);


#endif