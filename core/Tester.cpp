#include "inc/geometry.h"
#include "inc/transform.h"

int main(){

    Point3i p(-2, 7, -5);
    Point3i p2 = Abs(p);
    std::cout << p2 << std::endl;
    return EXIT_SUCCESS;


}