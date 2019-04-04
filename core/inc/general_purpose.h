#include <cmath>
#include <ostream>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <assert.h>
#include <string>
#include <cstring>    


#define Float_AS_DOUBLE
//#define Float_AS_FLOAT

#ifdef Float_AS_DOUBLE
    typedef double Float;
#endif

#ifdef Float_AS_FLOAT
    typedef float Float;
#endif

static constexpr Float MaxFloat = std::numeric_limits<Float>::max();
static constexpr Float Infinity = std::numeric_limits<Float>::infinity();



inline Float Lerp(Float t, Float v1, Float v2){
    return (((1 - t) * v1) + (t * v2));
}