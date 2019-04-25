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

#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)

inline Float gamma(int n) {
    return ((n * MachineEpsilon)/(1 - (n * MachineEpsilon)));
}



inline Float Lerp(Float t, Float v1, Float v2){
    return (((1 - t) * v1) + (t * v2));
}

// check that test is in [ref - precision, ref + precision]
inline bool floatEqual(Float test, Float ref, Float precision){
    test = std::abs(test); ref = std:: abs(ref); precision = std::abs(precision);
    Float min = ref - precision;
    Float max = ref + precision;
    if((min <= test) && ( test <= max))
        return true;
    return false;
}

inline Float Radians(Float deg){
    return ((deg * M_PI) / 180.0);
}