#ifndef TRANSFORM_H
#define TRANSFORM_H
#include "tennin.h"
#include "geometry.h"

#define Float_AS_DOUBLE
//#define Float_AS_FLOAT

#ifdef Float_AS_DOUBLE
    typedef double Float;
#endif

#ifdef Float_AS_FLOAT
    typedef float Float;
#endif

typedef struct Matrix4x4 Matrix4x4;

struct Matrix4x4{
    Float m[4][4];

    Matrix4x4();

    Matrix4x4(Float mat[4][4]);

    Matrix4x4(Float t00, Float t01, Float t02, Float t03,
              Float t10, Float t11, Float t12, Float t13,
              Float t20, Float t21, Float t22, Float t23,
              Float t30, Float t31, Float t32, Float t33);


    bool operator==(const Matrix4x4& mat)const;

    bool operator!=(const Matrix4x4& mat)const;


    std::string print();
    bool isIdentity();


    static Matrix4x4 Mul(const Matrix4x4& m1, const Matrix4x4& m2){
        Matrix4x4 r;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] +
                            m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] +
                            m1.m[i][3] * m2.m[3][j] ;
            }
        }
        return r;
    }


};

    std::ostream& operator<<(std::ostream& o, Matrix4x4 m);

class Transform{
    public: 
        Transform();
        Transform(const Float mat[4][4]);
        Transform(const Matrix4x4& m);
        Transform(const Matrix4x4& m, const Matrix4x4& mInv);

        std::string print();

        bool operator==(const Transform& t)const;
        bool operator!=(const Transform& t)const;
        bool isIdentity();


        friend Transform Inverse(const Transform& t);
        friend Transform Transpose(const Transform& t);


    private:
        Matrix4x4 m;
        Matrix4x4 mInv;
};



    template <typename T>
    Transform Translate(Vector3<T>& v){
        Matrix4x4 m(1,0,0,v.x,
                    0,1,0,v.y,
                    0,0,1,v.z,
                    0,0,0,1);

        Matrix4x4 mInv(1,0,0, -v.x,
                       0,1,0, -v.y,
                       0,0,1, -v.z,
                       0,0,0, 1);

        return Transform(m, mInv);
    }

    template <typename T>
    Transform Scale(Vector3<T>& v){
        Matrix4x4 m(v.x,0,  0,0,
                    0,v.y,  0,0,
                    0,  0,v.z,0,
                    0,  0,  0,1);

        Matrix4x4 mInv(1/v.x,0,0,0,
                       0,1/v.y,0,0,
                       0,0,1/v.z,0,
                       0,0,0,1);

        return Transform(m, mInv);
    }


    std::ostream& operator<<(std::ostream& o, Transform t);


#endif