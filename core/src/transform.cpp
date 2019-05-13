
#include "transform.h"


// Matrix4x4
    //----------------------------------------------//
    // Matrix4x4 struct constructor

        Matrix4x4::Matrix4x4(){
            m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
            m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
            m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
            m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
        }

        Matrix4x4::Matrix4x4(Float mat[4][4]){
            m[0][0] = mat[0][0]; m[0][1] = mat[0][1]; m[0][2] = mat[0][2]; m[0][3] = mat[0][3];
            m[1][0] = mat[1][0]; m[1][1] = mat[1][1]; m[1][2] = mat[1][2]; m[1][3] = mat[1][3];
            m[2][0] = mat[2][0]; m[2][1] = mat[2][1]; m[2][2] = mat[2][2]; m[2][3] = mat[2][3];
            m[3][0] = mat[3][0]; m[3][1] = mat[3][1]; m[3][2] = mat[3][2]; m[3][3] = mat[3][3];
        }

        Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03,
                             Float t10, Float t11, Float t12, Float t13,
                             Float t20, Float t21, Float t22, Float t23,
                             Float t30, Float t31, Float t32, Float t33){

            m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
            m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
            m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
            m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;      

        }


    // Matrix4x4 class operator overload

        bool Matrix4x4::operator==(const Matrix4x4& mat)const{
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(m[i][j] != mat.m[i][j])
                        return false;
                }
            }
            return true;
        }


        bool Matrix4x4::operator!=(const Matrix4x4& mat)const{
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(m[i][j] != mat.m[i][j])
                        return true;
                }
            }
            return false;
        }

 
    // Matrix4x4 class method
        
        std::string Matrix4x4::print(){
            char name[10000];
            sprintf(name,"|%.2F %.2F %.2F %.2F|\n"
                         "|%.2F %.2F %.2F %.2F|\n" 
                         "|%.2F %.2F %.2F %.2F|\n"
                         "|%.2F %.2F %.2F %.2F|\n",  m[0][0], m[0][1], m[0][2], m[0][3],
                                                     m[1][0], m[1][1], m[1][2], m[1][3],
                                                     m[2][0], m[2][1], m[2][2], m[2][3],
                                                     m[3][0], m[3][1], m[3][2], m[3][3]);
            return std::string(name);
        }

        bool Matrix4x4::isIdentity(){
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(((i == j) && (m[i][j] != 1)))
                        return false;
                    if(((i != j) && (m[i][j] != 0)))
                        return false;
                }
            }
            return true;
        }

        // static method of Matrix4x4
        Matrix4x4 Matrix4x4::Mul(const Matrix4x4& m1, const Matrix4x4& m2){
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

    // Matrix4x4 related functions

        Matrix4x4 Transpose(const Matrix4x4& mat){
            return Matrix4x4(mat.m[0][0], mat.m[1][0], mat.m[2][0], mat.m[3][0],
                             mat.m[0][1], mat.m[1][1], mat.m[2][1], mat.m[3][1],
                             mat.m[0][2], mat.m[1][2], mat.m[2][2], mat.m[3][2],
                             mat.m[0][3], mat.m[1][4], mat.m[2][3], mat.m[3][3]);
        }

        Matrix4x4 Inverse(const Matrix4x4& m){
            int index_row[4], index_col[4];
            int index_pivot[4] = {0, 0, 0, 0};
            Float minv[4][4];

            std::memcpy(minv, m.m, 4 * 4 * sizeof(Float));
            for(int i = 0 ; i < 4; i++){
                int irow = 0, icol = 0;
                Float big = 0.0f;

                //pivot selection
                for(int j = 0; j < 4; j++){
                    if(index_pivot[j] != 1){
                        for(int k = 0; k < 4; k++){
                            if(index_pivot[k] == 0){
                                if(std::abs(minv[j][k]) >= big){
                                    big = Float(std::abs(minv[j][k]));
                                    irow = j;
                                    icol = k;
                                }
                            }
                                else if(index_pivot[k] > 1){
                                    std::cout << "Singular matrix : aborting" << std::endl;
                                    std::exit(0); //provisoire
                                }
                        }
                    }
                }

                ++index_pivot[icol];
                // swap row irow and col icol
                if(irow != icol){
                    for(int k = 0; k < 4; ++k){
                            std::swap(minv[irow][k], minv[icol][k]);
                    }
                }


                index_row[i] = irow;
                index_col[i] = icol;
                if(minv[icol][icol] == 0.f){
                    std::cout << "Singular matrix : aborting" << std::endl;
                    std::exit(0); //provisoire
                }

                // scale row i_col according to pivot
                Float pivinv = 1. / minv[icol][icol];
                minv[icol][icol] = 1.;
                for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;


                // Linear operation to zero column
                for (int j = 0; j < 4; j++) {
                    if (j != icol) {
                        Float save = minv[j][icol];
                        minv[j][icol] = 0;
                        for (int k = 0; k < 4; k++){ 
                            minv[j][k] -= minv[icol][k] * save;
                        }
                    }
                }


            }

            for(int j = 3; j >= 0; j--) {
                if (index_row[j] != index_col[j]) {
                    for (int k = 0; k < 4; k++){
                        std::swap(minv[k][index_row[j]], minv[k][index_col[j]]);
                    }
                }
            }

            return Matrix4x4(minv);
        }


        std::ostream& operator<<(std::ostream& o, Matrix4x4 m){
            o << m.print();
            return o;
        }

//----------------------------------------------//


// Transform
    //----------------------------------------------//
    // Transform struct constructor

        Transform::Transform(): m(Matrix4x4()), mInv(Matrix4x4()){};

        Transform::Transform(const Matrix4x4& m) : m(m), mInv(Inverse(m)){}

        Transform::Transform(const Float mat[4][4]){
            m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                          mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                          mat[2][0], mat[2][1], mat[0][2], mat[0][3],
                          mat[3][0], mat[3][1], mat[0][2], mat[0][3]);

            mInv = Inverse(m);
        }

        Transform::Transform(const Matrix4x4& m, const Matrix4x4& mInv) : m(m), mInv(mInv){};


    // Transform operator overloading
        bool Transform::operator==(const Transform& t)const{
            if( (m == t.m) && (mInv == t.mInv) ){
                return true;
            }
            return false;
        }

        bool Transform::operator!=(const Transform& t)const{
            if( (m == t.m) && (mInv == t.mInv) ){
                return false;
            }
            return true;
        }

        Transform Transform::operator*(const Transform& t2)const{
            return Transform(Matrix4x4::Mul(m, t2.m), Matrix4x4::Mul(t2.mInv, mInv));
        }

        template <typename T>
        Point3<T> Transform::operator()(const Point3<T> &p) const{
            T x = p.x; 
            T y = p.y; 
            T z = p.z;

            T xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
            T yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
            T zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
            T wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];

            if(wp == 1) 
                return Point3<T>(xp, yp, zp);
            else{
                assert(wp != 0);
                return (Point3<T>(xp, yp, zp)/wp);
            }
        }

        template <typename T>
        Point3<T> Transform::operator()(const Point3<T>& p, Vector3<T>* absError) const{
            T x = p.x; 
            T y = p.y; 
            T z = p.z;

            T xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
            T yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
            T zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
            T wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];

            T xAbsSum = std::abs(m.m[0][0] * x ) + std::abs(m.m[0][1] * y ) + std::abs(m.m[0][2] * z ) + std::abs(m.m[0][3]);
            T yAbsSum = std::abs(m.m[1][0] * x ) + std::abs(m.m[1][1] * y ) + std::abs(m.m[1][2] * z ) + std::abs(m.m[1][3]);
            T zAbsSum = std::abs(m.m[2][0] * x ) + std::abs(m.m[2][1] * y ) + std::abs(m.m[2][2] * z ) + std::abs(m.m[2][3]);

            *absError = gamma(3) * Vector3<T>(xAbsSum, yAbsSum, zAbsSum);

            if(wp == 1) 
                return Point3<T>(xp, yp, zp);
            else {
                assert(wp != 0);
                return (Point3<T>(xp, yp, zp)/wp);  
            }          
        }

        template <typename T>
        Vector3<T> Transform::operator()(const Vector3<T> &v) const{
            T x = v.x; 
            T y = v.y; 
            T z = v.z;
            return Vector3<T>(m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z,
                              m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z,
                              m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z);
        }


        template <typename T>
        Normal3<T> Transform::operator()(const Normal3<T> &n) const{
            T x = n.x; 
            T y = n.y; 
            T z = n.z;
            return Normal3<T>(mInv.m[0][0]*x + mInv.m[0][1]*y + mInv.m[0][2]*z,
                              mInv.m[1][0]*x + mInv.m[1][1]*y + mInv.m[1][2]*z,
                              mInv.m[2][0]*x + mInv.m[2][1]*y + mInv.m[2][2]*z);            
        }


        Ray Transform::operator()(const Ray& r) const{
            Vect3f oError;
            Point3f o = (*this)(r.o, &oError);
            Vect3f d = (*this)(r.d);

            Float lengthSquared = d.lengthSquared();
            Float tMax = r.tMax;

            if(lengthSquared > 0){
                Float dt = Dot(Abs(d), oError)/ lengthSquared;
                o += (d * dt);
                tMax -= dt;
            }

            return Ray(o, d, tMax, r.time, r.medium);
        }


        RayDifferential Transform::operator()(const RayDifferential& r) const{

            Ray tr = (*this)(Ray(r));
            RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
            ret.hasDifferentials = r.hasDifferentials;
            ret.rxOrigin = (*this)(r.rxOrigin) ;
            ret.ryOrigin = (*this)(r.ryOrigin) ; 
            ret.rxDirection = (*this)(r.rxDirection) ; 
            ret.ryDirection = (*this)(r.ryDirection) ; 

            return ret;
        } 

        Bounds3f Transform::operator()(const Bounds3f& b)const{
            const Transform& M = (*this);
            Bounds3f ret(M(b.corner(0)));
            for(int i = 1; i < 8; i++){
                ret = Union(ret, M(b.corner(i)));
            }
            return ret;
        }


    // Transform class method
        bool Transform::isIdentity(){
            return (m.isIdentity() && mInv.isIdentity());
        }

        std::string Transform::print(){
            std::string txt = "Transform :\n----\nMatrix\n" + m.print() + "\nInverse\n" + mInv.print() + "\n" + "----" ;
            return txt;
        }

        bool Transform::hasScale()const{
            Float la2 = (*this)(Vect3f(1, 0, 0)).lengthSquared();
            Float lb2 = (*this)(Vect3f(0, 1, 0)).lengthSquared();
            Float lc2 = (*this)(Vect3f(0, 0, 1)).lengthSquared();

            return ( (!floatEqual(la2, 1.0f, 0.001)) || (!floatEqual(lb2, 1.0f, 0.001)) || (!floatEqual(lc2, 1.0f, 0.001) ));
        }

        Matrix4x4 Transform::getM()const{return m;}
        Matrix4x4 Transform::getInv()const{return mInv;}

        bool Transform::swapsHandedness()const{
            Float det = 
                m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
                m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
                m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
            
            return (det < 0);
        }


 // Transform related function

        Transform Inverse(const Transform& t){
            return Transform(t.mInv, t.m);
        }

        Transform Transpose(const Transform& t){
            return Transform(Transpose(t.m), Transpose(t.mInv));
        }

        std::ostream& operator<<(std::ostream& o, Transform t){
            o << t.print();
            return o;
        }


        template <typename T>
        Transform Translate(Vector3<T> v){
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
        Transform Scale(Vector3<T> v){
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


        Transform RotateX(Float theta){
            Float ct = std::cos(Radians(theta));
            Float st = std::sin(Radians(theta));

            Matrix4x4 m(  1,   0,    0,    0,
                          0,   ct,  -st,   0,
                          0,   st,   ct,   0,
                          0,    0,    0,   1);
            
            return Transform(m, Transpose(m));
        }


        Transform RotateY(Float theta){
            Float ct = std::cos(Radians(theta));
            Float st = std::sin(Radians(theta));

            Matrix4x4 m( ct,   0,   st,   0,
                          0,   1,    0,   0,
                        -st,   0,   ct,   0,
                          0,   0,    0,   1);
            
            return Transform(m, Transpose(m));
        }


        Transform RotateZ(Float theta){
            Float ct = std::cos(Radians(theta));
            Float st = std::sin(Radians(theta));

            Matrix4x4 m( ct,  -st,   0,   0,
                         st,   ct,   0,   0,
                          0,   0,    1,   0,
                          0,   0,    0,   1);
            
            return Transform(m, Transpose(m));
        }

        // Rodrigue's formula
        Transform Rotate (Float theta,const Vect3f axis){
            Vect3f a = Normalize(axis);
            Float ct = std::cos(Radians(theta));
            Float st = std::sin(Radians(theta));
            Matrix4x4 m;

            Float wx = a.x;
            Float wy = a.y;
            Float wz = a.z;

            m.m[0][0] = ct + (wx * wx * (1 - ct));         m.m[0][1] = (wx * wy * (1 - ct)) - (wz * st);   m.m[0][2] = (wy * st) + (wx * wz * (1 - ct)); 
            m.m[1][0] = (wz * st) + (wx * wy * (1 - ct));  m.m[1][1] = ct + (wy * wy * (1 - ct));          m.m[1][2] = (-wx * st) + (wy * wz * (1 - ct));
            m.m[2][0] = (-wy * st) + (wx * wz * (1-ct));   m.m[2][1] = (wx * st) + (wy * wz * (1 - ct));   m.m[2][2] = ct + (wz * wz * (1 - ct));    

            return Transform(m, Transpose(m));
        }


        Transform LookAt(const Point3f pos, const Point3f look, const Vect3f up){
            Matrix4x4 cameraToWorld;

            cameraToWorld.m[0][3] = pos.x;
            cameraToWorld.m[1][3] = pos.y;
            cameraToWorld.m[2][3] = pos.z;
            cameraToWorld.m[3][3] = 1;

            Vect3f dir = Normalize(look - pos);
            Vect3f left = Normalize(Cross(Normalize(up), dir));
            Vect3f newUp = Cross(dir, left);


            cameraToWorld.m[0][0] = left.x;
            cameraToWorld.m[1][0] = left.y;
            cameraToWorld.m[2][0] = left.z;
            cameraToWorld.m[3][0] = 0;

            cameraToWorld.m[0][1] = newUp.x;
            cameraToWorld.m[1][1] = newUp.y;
            cameraToWorld.m[2][1] = newUp.z;
            cameraToWorld.m[3][1] = 0;

            cameraToWorld.m[0][2] = dir.x;
            cameraToWorld.m[1][2] = dir.y;
            cameraToWorld.m[2][2] = dir.z;
            cameraToWorld.m[3][2] = 0;


            return Transform(Inverse(cameraToWorld), cameraToWorld);

        }

        
//----------------------------------------------//

template Transform Translate(Vector3<int> v);
template Transform Translate(Vector3<Float> v);
template Transform Scale(Vector3<int> v);
template Transform Scale(Vector3<Float> v);


