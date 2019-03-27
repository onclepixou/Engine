#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <assert.h>


#include "transform.h"

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

    Transform::Transform(const Matrix4x4& m) : m(m), mInv(Inverse(m)){}

    Transform::Transform(const Float mat[4][4]){
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                      mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                      mat[2][0], mat[2][1], mat[0][2], mat[0][3],
                      mat[3][0], mat[3][1], mat[0][2], mat[0][3]);

        mInv = Inverse(m);
    }


int main(){
    Float Mat1[4][4] = {2,2,2,2,2,2,2,2,0,2,2,2,2,2,2,2};
    Float Mat2[4][4] = {4,5,6,1,2,3,2,3,1,7,4,3,4,8,5,3};

    Matrix4x4 m1(Mat1) ;
    Matrix4x4 m2(Mat2) ;
    Matrix4x4 m3 = Inverse(m2);

    std::cout << Matrix4x4::Mul(m1, m2) << std::endl;

    std::cout << (m1 == m2) << std::endl;

    std::cout << "inverse matrix" << std::endl;
    std::cout << m3 << std::endl;
}