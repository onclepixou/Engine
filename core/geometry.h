#include <cmath>
#include <ostream>
#include <limits>

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

// class declaration //
template <typename T> class Point2;
template <typename T> class Point3;
template <typename T> class Vector2;
template <typename T> class Vector3;
template <typename T> class Normal3;
template <typename T> class Bounds2;
template <typename T> class Bounds3;
class Ray;
class Medium;
//----------------------------------//

typedef Vector2<int> Vect2i;
typedef Vector2<Float> Vect2f;
typedef Vector3<int> Vect3i;
typedef Vector3<Float> Vect3f;

typedef Point2<int> Point2i;
typedef Point2<Float> Point2f;
typedef Point3<int> Point3i;
typedef Point3<Float> Point3f;

typedef Normal3<int> Normal3i;
typedef Normal3<Float> Normal3f;

typedef Bounds2<int> Bounds2i;
typedef Bounds2<Float> Bounds2f;
typedef Bounds3<int> Bounds3i;
typedef Bounds3<Float> Bounds3f;


template <typename T> class Vector3{

    public:
        T x;
        T y;
        T z;

        Vector3() : x(0), y(0), z(0){};
        Vector3(T a, T b, T c) : x(a), y(b), z(c){assert(!hasNan());};
        Vector3(const Vector2<T>& v) : x(v.x), y(v.y), z(0) {assert(!hasNan());};
        explicit Vector3(const Normal3<T> &n): x(n.x), y(n.y), z(n.z){assert(!hasNan());};
        template <typename U> explicit operator Vector3<U>() const {return Vector3<U>(x, y, z);}
        Vector3<T> operator+(const Vector3<T>& rhs) const{ return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z); };
        Vector3<T>& operator+=(const Vector3<T>& rhs) { x += rhs.x; y += rhs.y; z+= rhs.z;return *this;};
        Vector3<T> operator-(const Vector3<T>& rhs) const{ return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z); };
        Vector3<T>& operator-=(const Vector3<T>& rhs) { x -= rhs.x; y -= rhs.y; z-= rhs.z;return *this;};
        Vector3<T> operator*(T s) const {return Vector3<T>( s * x, s * y, s * z);};
        Vector3<T>& operator*=(T s){ x = s * x; y = s * y; z = s * z; return *this;};
        Vector3<T> operator/(T i)const{assert(i != 0); T inv = 1/i; return Vector3<T>(x * inv, y * inv, z * inv);}
        Vector3<T>& operator/=(T i){assert(i != 0); T inv = 1/i; x *= inv; y *= inv; z *= inv; return *this;}
        Vector3<T> operator-(){return Vector3<T>(-x, -y, -z);};
        T  operator[](int i);
        T& operator[](int i) const;

        bool hasNan()const;       
        std::string tostring();
        Vector3<T> abs(){return Vector3<T>(std::abs(x), std::abs(y), std::abs(z));}
        Float lenghtSquared()const{return x * x + y * y + z * z;};
        Float lenght()const{return std::sqrt(lenghtSquared());};
        bool isNormalized(){return ((lenght() - 1 ) <= 1e-6);}

};

template <typename T> class Vector2{

    public: 
        T x;
        T y;

        Vector2() : x(0), y(0){};
        Vector2(T a, T b) : x(a), y(b){assert(!hasNan());};
        Vector2(const Vector3<T>& v) : x(v.x), y(v.y){assert(!hasNan());};
        template <typename U> explicit operator Vector2<U>() const {return Vector2<U>(x, y);}
        Vector2<T> operator+(const Vector2<T>& rhs) const{ return Vector2<T>(x + rhs.x, y + rhs.y); };
        Vector2<T>& operator+=(const Vector2<T>& rhs) { x += rhs.x; y += rhs.y;return *this;};
        Vector2<T> operator-(const Vector2<T>& rhs) const{ return Vector2<T>(x - rhs.x, y - rhs.y); };
        Vector2<T>& operator-=(const Vector2<T>& rhs) { x -= rhs.x; y -= rhs.y;return *this;};
        Vector2<T> operator*(T s) const {return Vector2<T>( s * x, s * y);};
        Vector2<T>& operator*=(T s){  x = s * x; y = s * y; return *this;};
        Vector2<T> operator/(T i)const{assert(i != 0); T inv = 1/i; return Vector2<T>(x * inv, y * inv);};
        Vector2<T>& operator/=(T i){assert(i != 0); T inv = 1/i; x *= inv; y *= inv; return *this;};
        Vector2<T> operator-(){return Vector2<T>(-x, -y);};
        T operator[](int i)const;
        T& operator[](int i);

        bool hasNan() const;
        std::string tostring();
        Vector2<T> abs(){return Vector2<T>(std::abs(x), std::abs(y));};
        Float lenghtSquared()const{return x * x + y * y;};
        Float lenght()const{return std::sqrt(lenghtSquared());};
        bool isNormalized(){return ((lenght() - 1 ) <= 1e-6);}


};



template <typename T> class Point3{

    public:
        T x;
        T y;
        T z;

        Point3():x(0), y(0), z(0){};
        Point3(T a, T b, T c): x(a), y(b), z(c){ assert(!hasNan());};
        Point3(const Point2<T>& p):x(p.x), y(p.y), z(0){ assert(!hasNan());};
        template <typename U> explicit Point3(const Point3<U> p): x((T)p.x), y((T)p.y), z((T)p.z){assert(!hasNan());};
        Point3<T> operator+(const Point3<T>& p)const{return Point3<T>((p.x + x), (p.y + y), (p.z + z));};
        Point3<T>& operator+=(const Point3<T>& p){x += p.x; y += p.y; z += p.z; return (*this); };
        Point3<T> operator+(const Vector3<T>& v)const{return Point3<T>((v.x + x), (v.y + y), (v.z + z));};
        Point3<T>& operator+=(const Vector3<T>& v){x += v.x; y += v.y; z += v.z; return (*this); };
        Vector3<T> operator-(const Point3<T>& p)const{return Vector3<T>((x-p.x), (y-p.y), (z-p.z));}
        Point3<T> operator-(const Vector3<T>& v)const{return Point3<T>((x-v.x), (y-v.y), (z-v.z));};
        Point3<T>& operator-=(const Vector3<T>& v){x -= v.x; y -= v.y; z -= v.z; return (*this);}
        Point3<T> operator*(const T& s) const {return Point3<T>( s * x, s * y, s * z);};
        Point3<T>& operator*=(const T& s){ x = s * x; y = s * y; z = s * z; return *this;};
        T  operator[](int i)const;
        T& operator[](int i);

        bool hasNan()const;
        std::string tostring();
        
         
};

template <typename T> class Point2{

    public:
        T x;
        T y;

        Point2():x(0), y(0){};
        Point2(T a, T b): x(a), y(b){ assert(!hasNan());};
        explicit Point2(const Point3<T>& p): x(p.x), y(p.y){assert(!hasNan());};
        template <typename U> explicit Point2(const Point2<U> p): x((T)p.x), y((T)p.y){assert(!hasNan());};
        Point2<T> operator+(const Point2<T>& p) const{return Point3<T>((p.x + x), (p.y + y));};
        Point2<T>& operator+=(const Point2<T>& p){x += p.x; y += p.y; return *this;};
        Point2<T> operator+(const Vector2<T>& v)const{return Point2<T>((v.x + x), (v.y + y));};
        Point2<T>& operator+=(const Vector2<T>& v){x += v.x; y += v.y; return (*this); };
        Vector2<T> operator-(const Point2<T>& p)const{return Vector2<T>((x-p.x), (y-p.y));};
        Point2<T> operator-(const Vector2<T>& v)const{return Point2<T>((x-v.x), (y-v.y));};
        Point2<T>& operator-=(const Vector2<T>& v){x -= v.x; y -= v.y; return (*this);};
        Point2<T> operator*(const T& s) const {return Point2<T>( s * x, s * y);};
        Point2<T>& operator*=(const T& s){ x = s * x; y = s * y; return *this;};
        T  operator[](int i)const;
        T& operator[](int i);

        bool hasNan()const;
        std::string tostring();
         
};


template <typename T> class Normal3{

    public:
        T x;
        T y;
        T z;

        Normal3() : x(0), y(0), z(0){};
        Normal3(T a, T b, T c) : x(a), y(b), z(c){assert(!hasNan());};
        explicit Normal3(const Vector3<T>& v) : x(v.x), y(v.y), z(v.z){assert(!hasNan());};
        template <typename U> explicit operator Normal3<U>() const {return Normal3<U>(x, y, z);}
        Normal3<T> operator+(const Normal3<T>& rhs) const{ return Normal3<T>(x + rhs.x, y + rhs.y, z + rhs.z); };
        Normal3<T>& operator+=(const Normal3<T>& rhs) { x += rhs.x; y += rhs.y; z+= rhs.z;return *this;};
        Normal3<T> operator-(const Normal3<T>& rhs) const{ return Normal3<T>(x - rhs.x, y - rhs.y, z - rhs.z); };
        Normal3<T>& operator-=(const Normal3<T>& rhs) { x -= rhs.x; y -= rhs.y; z-= rhs.z;return *this;};
        Normal3<T> operator*(T s) const {return Normal3<T>( s * x, s * y, s * z);};
        Normal3<T>& operator*=(T s){ x = s * x; y = s * y; z = s * z; return *this;};
        Normal3<T> operator/(T i)const{assert(i != 0); T inv = 1/i; return Normal3<T>(x * inv, y * inv, z * inv);}
        Normal3<T>& operator/=(T i){assert(i != 0); T inv = 1/i; x *= inv; y *= inv; z *= inv; return *this;}
        Normal3<T> operator-(){return Normal3<T>(-x, -y, -z);};
        T  operator[](int i);
        T& operator[](int i) const;
        
        bool hasNan()const;       
        std::string tostring();
        Normal3<T> abs(){return Normal3<T>(std::abs(x), std::abs(y), std::abs(z));}
        Float lenghtSquared()const{return x * x + y * y + z * z;};
        Float lenght()const{return std::sqrt(lenghtSquared());};
        bool isNormalized(){return ((lenght() - 1 ) <= 1e-6);}

};

class Ray{
    public: 
        Point3<Float> o;
        Vector3<Float> d;
        mutable Float tMax;
        Float time;
        const Medium *medium;

        Ray():tMax(Infinity), time(0.0f), medium(nullptr){};
        Ray(const Point3<Float>& o, const Vector3<Float>& d, Float tMax = Infinity,
            Float time = 0.0f, const Medium* medium = nullptr) : o(o), d(d), tMax(tMax), time(time), medium(medium){};
        
        Point3<Float> operator()(Float t){return o + (d * t);}

};

class RayDifferential : public Ray{
    public: 
        bool hasDifferentials;
        Point3<Float> rxOrigin, ryOrigin;
        Vector3<Float> rxDirection, ryDirection;

        RayDifferential(){hasDifferentials = false;};
        RayDifferential(const Point3<Float>& o, const Vector3<Float>& d, 
                        Float tMax = Infinity, Float time = 0.0f, const Medium* medium = nullptr)
                        : Ray(o, d, tMax, time, medium){hasDifferentials = false;};
        RayDifferential(const Ray& r) : Ray(r){hasDifferentials = false; };

        void scaleDifferentials(Float s);
};


template <typename T> class Bounds2{
    public: 
        Point2<T> pMin;
        Point2<T> pMax;
};

template <typename T> class Bounds3{

    public: 
        Point3<T> pMin;
        Point3<T> pMax;

        Bounds3(){
            T minNum = std::numeric_limits<T>::lowest();
            T maxNum = std::numeric_limits<T>::max();
            pMin = Point3<T>(maxNum, maxNum, maxNum);
            pMax = Point3<T>(minNum, minNum, minNum);
        }
        Bounds3(const Point3<T>& p) : pMin(p), pMax(p){};
        Bounds3(const Point3<T>& p1, const Point3<T>& p2)
         : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
           pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)){};

        Point3<T> corner(int i) const;
        Vector3<T> diagonal() const;
        T surfaceArea()const;
        T volume()const;
        int maximumExtent()const;
        Point3<T> lerp(const Point3<Float>& t) const;
        Vector3<T> offset(const Point3<T>& p) const;
        void boundingSphere(Point3<T>* center, Float* radius)const;


        const Point3<T>& operator[](int i)const;
        Point3<T>& operator[](int i);


};


