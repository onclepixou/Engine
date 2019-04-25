#include "inc/geometry.h"
#include "inc/transform.h"

int main(){


    Transform A = Translate(Vect3f(1.0, 2.0, 3.0));
    Transform B = Rotate(90, Vect3f(1.0, 0, 0));


    Transform C = A * B;

    std::cout << C << std::endl;
    std::cout << Matrix4x4::Mul(A.getM(), B.getM()) << std::endl;
    std::cout << Matrix4x4::Mul(B.getInv(), A.getInv()) << std::endl;

    return EXIT_SUCCESS;


}