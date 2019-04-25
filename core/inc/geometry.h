#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "tennin.h"

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


//------------------------------------------------------------------------------------//
// VECTOR 3 //

template <typename T> class Vector3{

    public:
        T x;
        T y;
        T z;

        // Constructors
        Vector3();
        Vector3(T a, T b, T c);
        Vector3(const Vector2<T>& v);
        explicit Vector3(const Normal3<T> &n);
        template <typename U> explicit operator Vector3<U>() const;

        // Operator Overloading
        T  operator[](int i)const;
        T& operator[](int i);
        bool operator==(const Vector3<T>& v);
        bool operator!=(const Vector3<T>& v);
        Vector3<T> operator+(const Vector3<T>& rhs) const;
        Vector3<T>& operator+=(const Vector3<T>& rhs);
        Vector3<T> operator-(const Vector3<T>& rhs) const;
        Vector3<T>& operator-=(const Vector3<T>& rhs);
        Vector3<T> operator*(T s) const;
        Vector3<T>& operator*=(T s);
        Vector3<T> operator/(T i)const;
        Vector3<T>& operator/=(T i);
        Vector3<T> operator-();

        // Public methods
        bool hasNan()const;
        std::string tostring();
        Vector3<T> abs();
        Float lengthSquared()const;
        Float length()const;
        bool isNormalized();

};

// Vector3 related function
    template <typename T> 
    inline Vector3<T> operator*(T s, Vector3<T> v){
        return Vector3<T>(s * v.x, s * v.y, s * v.z);
    }

    template <typename T>
    inline T Dot(const Vector3<T>& a, const Vector3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}

    template <typename T>
    inline Vector3<T> Abs(const Vector3<T>& a){return Vector3<T>(std::abs(a.x), std::abs(a.y), std::abs(a.z));}


    template <typename T>
    inline T AbsDot(const Vector3<T>& a, const Vector3<T>& b){return std::abs(dot(a, b));}

    template <typename T>
    inline Vector3<T> Cross(const Vector3<T>& a, const Vector3<T>& b){

        Float ax = a.x, ay = a.y, az = a.z;
        Float bx = a.x, by = b.y, bz = b.z;
        return Vector3<T>((ay * bz - az * by), (az * bx - ax * bz), (ax * by - ay * bx));
    }

    template <typename T>
    inline Vector3<T> Normalize(const Vector3<T>& v){return (v/v.length());}

    template <typename T> 
    inline T MinComponent(const Vector3<T>& v){
        return std::min(std::min(v.x, v.y), v.z);
    }

    template <typename T> 
    inline T MaxComponent(const Vector3<T>& v){
        return std::max(std::max(v.x, v.y), v.z);
    }

    template <typename T> 
    inline int MinDimension(const Vector3<T>& v){
        return (v.x < v.y) ? ((v.y < v.z) ? 1 : 2) : ((v.x < v.z) ? 0 : 2);
    }

    template <typename T> 
    inline int MaxDimension(const Vector3<T>& v){
        return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
    }

    template <typename T>
    inline Vector3<T> Min(const Vector3<T>& v1, const Vector3<T>& v2){
        return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    template <typename T>
    inline Vector3<T> Max(const Vector3<T>& v1, const Vector3<T>& v2){
        return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    template <typename T>
    inline Vector3<T> Permute(const Vector3<T>& v, int x, int y, int z){
        return Vector3<T>(v[x], v[y], v[z]);
    }

    template <typename T>
    inline void CoordinateSystem(const Vector3<T>& v1, Vector3<T>* v2, Vector3<T>* v3){
        if(std::abs(v1.x) > std::abs(v1.y))
            *v2 = (Vector3<T>(-v1.z, 0, v1.x) / std::sqrt((v1.x * v1.x) + (v1.z * v1.z)) );
        else
            *v2 = (Vector3<T>(0, v1.z, -v1.y) / std::sqrt((v1.y * v1.y) + (v1.z * v1.z)) );
        *v3 = cross(v1, *v2);

        return;
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Vector3<T>& v){
        os << v.tostring();
        return os;
    }

//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// VECTOR 2 //

template <typename T> class Vector2{

    public: 
        T x;
        T y;

        // Constructors
        Vector2();
        Vector2(T a, T b);
        Vector2(const Vector3<T>& v);
        template <typename U> explicit operator Vector2<U>() const;

        // Operator Overloading
        T operator[](int i)const;
        T& operator[](int i);
        bool operator==(const Vector2<T>& v);
        bool operator!=(const Vector2<T>& v);

        Vector2<T> operator+(const Vector2<T>& rhs)const;
        Vector2<T>& operator+=(const Vector2<T>& rhs);
        Vector2<T> operator-(const Vector2<T>& rhs) const;
        Vector2<T>& operator-=(const Vector2<T>& rhs);
        Vector2<T> operator*(T s) const;
        Vector2<T>& operator*=(T s);
        Vector2<T> operator/(T i)const;
        Vector2<T>& operator/=(T i);
        Vector2<T> operator-();

        // Public methods
        bool hasNan() const;
        std::string tostring();
        Vector2<T> abs();
        Float lengthSquared()const;
        Float length()const;
        bool isNormalized()const;


};

// Vector2 related function
    template <typename T> 
    inline Vector2<T> operator*(T s, Vector2<T> v){
        return Vector2<T>(s * v.x, s * v.y);
    }

    template <typename T>
    inline T Dot(const Vector2<T>& a, const Vector2<T>& b){return a.x * b.x + a.y * b.y;}

    template <typename T>
    inline Vector2<T> Abs(const Vector2<T>& a){return Vector2<T>(std::abs(a.x), std::abs(a.y));}

    template <typename T>
    inline T AbsDot(const Vector2<T>& a, const Vector2<T>& b){return std::abs(dot(a, b));}

    template <typename T>
    inline Vector2<T> Cross(const Vector2<T>& a, const Vector2<T>& b){

        Float ax = a.x, ay = a.y, az = a.z;
        Float bx = a.x, by = b.y, bz = b.z;
        return Vector2<T>((ax * by - ay * bx), (ay * bx - ax * by));
    }

    template <typename T>
    inline Vector2<T> Normalize(const Vector2<T>& v){return (v/v.length());}

    template <typename T> 
    inline T MinComponent(const Vector2<T>& v){
        return std::min(v.x, v.y);
    }

    template <typename T> 
    inline T MaxComponent(const Vector2<T>& v){
        return std::max(v.x, v.y);
    }

    template <typename T> 
    inline int MinDimension(const Vector2<T>& v){
        return ( (v.x > v.y )? 1 : 0 );
    }

    template <typename T> 
    inline int MaxDimension(const Vector2<T>& v){
        return ( (v.x > v.y )? 0 : 1 );
    }

    template <typename T>
    inline Vector2<T> Min(const Vector2<T>& v1, const Vector2<T>& v2){
        return Vector2<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y));
    }

    template <typename T>
    inline Vector2<T> Max(const Vector2<T>& v1, const Vector2<T>& v2){
        return Vector2<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y));
    }

    template <typename T>
    inline Vector2<T> Permute(const Vector2<T>& v, int x, int y){
        return Vector2<T>(v[x], v[y]);
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Vector2<T>& v){
        os << v.tostring();
        return os;
    }
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Point 3 //

template <typename T> class Point3{

    public:
        T x;
        T y;
        T z;

        // Constructors
        Point3();
        Point3(T a, T b, T c);
        Point3(const Point2<T>& p);
        template <typename U> explicit Point3(const Point3<U> p);

        // Operator Overloading
        T  operator[](int i)const;
        T& operator[](int i);
        bool operator==(const Point3<T>& p);
        bool operator!=(const Point3<T>& p);
        Point3<T> operator+(const Point3<T>& p)const;
        Point3<T>& operator+=(const Point3<T>& p);
        Point3<T> operator+(const Vector3<T>& v)const;
        Point3<T>& operator+=(const Vector3<T>& v);
        Vector3<T> operator-(const Point3<T>& p)const;
        Point3<T> operator-(const Vector3<T>& v)const;
        Point3<T>& operator-=(const Vector3<T>& v);
        Point3<T> operator*(const T& s) const;
        Point3<T>& operator*=(const T& s);
        Point3<T> operator/(const T& n)const;
        Point3<T>& operator/=(const T& n);

        // Public methods
        bool hasNan()const;
        std::string tostring();
        
         
};

// Point3 related function

    template <typename T>
    inline Float Distance(const Point3<T>& p1, const Point3<T>& p2){
        return (p1 - p2).length();
    }

    template <typename T>
    inline Float DistanceSquared(const Point3<T>& p1, const Point3<T>& p2){
        return (p1 - p2).lengthSquared();
    }

    template <typename T>
    inline Point3<T> operator*(T s, const Point3<T> p){
        return Point3<T>(s * p.x, s * p.y, s * p.z); 
    } 

    template <typename T>
    inline Point3<T> Lerp(Float t, const Point3<T>& p0, const Point3<T>& p1){
        return (1-t) * p0 + p1;
    }

    template <typename T>
    inline Point3<T> Min(const Point3<T>& p1, const Point3<T>& p2){
        return Point3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
    }

    template <typename T>
    inline Point3<T> Max(const Point3<T>& p1, const Point3<T>& p2){
        return Point3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
    }

    template <typename T>
    inline Point3<T> Floor(const Point3<T>& p){
        return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
    }

    template <typename T>
    inline Point3<T> Ceil(const Point3<T>& p){
        return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
    }

    template <typename T>
    inline Point3<T> Abs(const Point3<T>& p){
        return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
    }

    template <typename T>
    inline Point3<T> Permute(const Point3<T>& p, int x, int y, int z){
        return Point3<T>(p[x], p[y], p[z]);
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Point3<T>& p){
        os << p.tostring();
        return os;
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Point3<T> p){
        os << p.tostring();
        return os;
    }
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Point 2 //

template <typename T> class Point2{

    public:
        T x;
        T y;

        // Constructors
        Point2();
        Point2(T a, T b);
        explicit Point2(const Point3<T>& p);
        template <typename U> explicit Point2(const Point2<U> p);

        // Operator Overloading
        T  operator[](int i)const;
        T& operator[](int i);
        bool operator==(const Point3<T>& p);
        bool operator!=(const Point3<T>& p);
        Point2<T> operator+(const Point2<T>& p) const;
        Point2<T>& operator+=(const Point2<T>& p);
        Point2<T> operator+(const Vector2<T>& v)const;
        Point2<T>& operator+=(const Vector2<T>& v);
        Vector2<T> operator-(const Point2<T>& p)const;
        Point2<T> operator-(const Vector2<T>& v)const;
        Point2<T>& operator-=(const Vector2<T>& v);
        Point2<T> operator*(const T& s) const;
        Point2<T>& operator*=(const T& s);
        Point2<T> operator/(const T& n)const;
        Point2<T>& operator/=(const T& n);


        // Public methods
        bool hasNan()const;
        std::string tostring();
         
};

// Point2 related function
    template <typename T>
    inline Float Distance(const Point2<T>& p1, const Point2<T>& p2){
        return (p1 - p2).length();
    }

    template <typename T>
    inline Float DistanceSquared(const Point2<T>& p1, const Point2<T>& p2){
        return (p1 - p2).lengthSquared();
    }

    template <typename T>
    inline Point2<T> operator*(T s, const Point2<T> p){
        return Point2<T>(s * p.x, s * p.y); 
    }

    template <typename T>
    inline Point2<T> Lerp(Float t, const Point2<T>& p0, const Point2<T>& p1){
        return (1-t) * p0 + p1;
    }

    template <typename T>
    inline Point2<T> Min(const Point2<T>& p1, const Point2<T>& p2){
        return Point2<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
    }

    template <typename T>
    inline Point2<T> Max(const Point2<T>& p1, const Point2<T>& p2){
        return Point2<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
    }

    template <typename T>
    inline Point2<T> Floor(const Point2<T>& p){
        return Point2<T>(std::floor(p.x), std::floor(p.y));
    }

    template <typename T>
    inline Point2<T> Ceil(const Point2<T>& p){
        return Point2<T>(std::ceil(p.x), std::ceil(p.y));
    }

    template <typename T>
    inline Point2<T> Abs(const Point2<T>& p){
        return Point2<T>(std::abs(p.x), std::abs(p.y));
    }

    template <typename T>
    inline Point2<T> Permute(const Point2<T>& v, int x, int y){
        return Point2<T>(v[x], v[y]);
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Point2<T>& p){
        os << p.tostring();
        return os;
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Point2<T> p){
        os << p.tostring();
        return os;
    }
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Normal3

template <typename T> class Normal3{

    public:
        T x;
        T y;
        T z;

        // Constructors
        Normal3();
        Normal3(T a, T b, T c);
        explicit Normal3(const Vector3<T>& v);
        template <typename U> explicit operator Normal3<U>() const;

        // Operator Overloading
        T  operator[](int i)const;
        T& operator[](int i);
        bool operator==(const Normal3<T>& v);
        bool operator!=(const Normal3<T>& v);
        Normal3<T> operator+(const Normal3<T>& rhs) const;
        Normal3<T>& operator+=(const Normal3<T>& rhs);
        Normal3<T> operator-(const Normal3<T>& rhs) const;
        Normal3<T>& operator-=(const Normal3<T>& rhs);
        Normal3<T> operator*(T s) const;
        Normal3<T>& operator*=(T s);
        Normal3<T> operator/(T i)const;
        Normal3<T>& operator/=(T i);
        Normal3<T> operator-();

        
        bool hasNan()const;       
        std::string tostring();
        Normal3<T> abs();
        Float lengthSquared()const;
        Float length()const;
        bool isNormalized();

};

// Normal3 related functions

    template <typename T> 
    inline Normal3<T> operator*(T s, Normal3<T> v){
        return Normal3<T>(s * v.x, s * v.y, s * v.z);
    }


    template <typename T>
    inline Normal3<T> Normalize(const Normal3<T>& v){return (v/v.length());}

    template <typename T>
    inline T Dot(const Normal3<T>& a, const Normal3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}

    template <typename T>
    inline T AbsDot(const Normal3<T>& a, const Normal3<T>& b){return std::abs(dot(a, b));}

    template <typename T>
    inline T Dot(const Normal3<T>& a, const Vector3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}

    template <typename T>
    inline T AbsDot(const Normal3<T>& a, const Vector3<T>& b){return std::abs(dot(a, b));}

    template <typename T>
    inline T Dot(const Vector3<T>& a, const Normal3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}

    template <typename T>
    inline T AbsDot(const Vector3<T>& a, const Normal3<T>& b){return std::abs(dot(a, b));}

    template <typename T> 
    inline T MinComponent(const Normal3<T>& v){
        return std::min(std::min(v.x, v.y), v.z);
    }

    template <typename T> 
    inline T MaxComponent(const Normal3<T>& v){
        return std::max(std::max(v.x, v.y), v.z);
    }

    template <typename T> 
    inline int MinDimension(const Normal3<T>& v){
        return (v.x < v.y) ? ((v.y < v.z) ? 1 : 2) : ((v.x < v.z) ? 0 : 2);
    }

    template <typename T> 
    inline int MaxDimension(const Normal3<T>& v){
        return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
    }

    template <typename T>
    inline Normal3<T> Min(const Normal3<T>& v1, const Normal3<T>& v2){
        return Normal3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    template <typename T>
    inline Normal3<T> Max(const Normal3<T>& v1, const Normal3<T>& v2){
        return Normal3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    template <typename T>
    inline Normal3<T> Permute(const Normal3<T>& v, int x, int y, int z){
        return Normal3<T>(v[x], v[y], v[z]);
    }


    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Normal3<T>& v){
        os << v.tostring();
        return os;
    }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, Normal3<T> v){
        os << v.tostring();
        return os;
    }

    template <typename T>
    inline Normal3<T> Faceforward(const Normal3<T>& n, const Vector3<T>& v){
        return ((dot(n,v) < 0.0f) ? -n : n);
    }

    template <typename T>
    inline void CoordinateSystem(const Normal3<T>& v1, Normal3<T>* v2, Normal3<T>* v3){
        if(std::abs(v1.x) > std::abs(v1.y))
            *v2 = (Vector3<T>(-v1.z, 0, v1.x) / std::sqrt((v1.x * v1.x) + (v1.z * v1.z)) );
        else
            *v2 = (Vector3<T>(0, v1.z, -v1.y) / std::sqrt((v1.y * v1.y) + (v1.z * v1.z)) );
        *v3 = cross(v1, *v2);

        return;
    }

//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Ray
    class Ray{
        public: 
            Point3<Float> o;
            Vector3<Float> d;
            mutable Float tMax;
            Float time;
            const Medium *medium;

            // class constructor

            Ray();
            Ray(const Point3<Float>& o, const Vector3<Float>& d, Float tMax = Infinity,
                Float time = 0.0f, const Medium* medium = nullptr);

            // operator overloading
            Point3<Float> operator()(Float t);

    };
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// RayDifferential
class RayDifferential : public Ray{
    public: 
        bool hasDifferentials;
        Point3<Float> rxOrigin, ryOrigin;
        Vector3<Float> rxDirection, ryDirection;

        // class constructor

        RayDifferential();
        RayDifferential(const Point3<Float>& o, const Vector3<Float>& d, 
                        Float tMax = Infinity, Float time = 0.0f, const Medium* medium = nullptr);
        RayDifferential(const Ray& r);

        void scaleDifferentials(Float s);
};
//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// Bound2

template <typename T> class Bounds2{
    public: 
        Point2<T> pMin;
        Point2<T> pMax;

        // class constructor

        Bounds2();
        Bounds2(const Point2<T>& p);
        Bounds2(const Point2<T>& p1, const Point2<T>& p2);

        // operator overloading

        const Point2<T>& operator[](int i)const;
        Point2<T>& operator[](int i);
        bool operator==(const Bounds2<T>& b);
        bool operator!=(const Bounds2<T>& b);

        // class method

        Point2<T> corner(int i) const;
        Vector2<T> diagonal() const;
        T surfaceArea()const;
        int maximumExtent()const;
        Point2<T> lerp(const Point2<Float>& t) const;
        Vector2<T> offset(const Point2<T>& p) const;
        void boundingCircle(Point2<T>* center, Float* radius)const;


};

// Bounds 2 related functions

        template <typename T>
        Bounds2<T> Union(const Bounds2<T>& b, const Point2<T>& p);
    
        template <typename T> 
        Bounds2<T> Union(const Bounds2<T>& b1, const Bounds2<T>& b2);
    
        template <typename T> 
        Bounds2<T> Intersect(const Bounds2<T>& b1, const Bounds2<T>& b2);
    
        template <typename T>
        bool Overlaps(const Bounds2<T>& b1, const Bounds2<T>& b2);
    
        template <typename T>
        bool Inside(const Point2<T>& p, const Bounds2<T>& b);
    
        template <typename T>
        bool InsideExclusive(const Point2<T>& p, const Bounds2<T>& b);
    
        template <typename T, typename U>
        inline Bounds2<T> Expand(const Bounds2<T>& b, U delta){
            return Bounds2<T>(b.pMin - Vector2<U>(delta, delta),
                              b.pMax + Vector2<U>(delta, delta));
        }

//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//

template <typename T> class Bounds3{

    public: 
        Point3<T> pMin;
        Point3<T> pMax;

        // Class constructor

        Bounds3();
        Bounds3(const Point3<T>& p);
        Bounds3(const Point3<T>& p1, const Point3<T>& p2);

        // Operator overloading

        const Point3<T>& operator[](int i)const;
        Point3<T>& operator[](int i);
        bool operator==(const Bounds3<T>& b);
        bool operator!=(const Bounds3<T>& b);

        // Class method 

        Point3<T> corner(int i) const;
        Vector3<T> diagonal() const;
        T surfaceArea()const;
        T volume()const;
        int maximumExtent()const;
        Point3<T> lerp(const Point3<Float>& t) const;
        Vector3<T> offset(const Point3<T>& p) const;
        void boundingSphere(Point3<T>* center, Float* radius)const;

};


// Bounds3 related function //
        template <typename T>
        Bounds3<T> Union(const Bounds3<T>& b, const Point3<T>& p);

        template <typename T> 
        Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2);

        template <typename T> 
        Bounds3<T> Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2);

        template <typename T>
        bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2);

        template <typename T>
        bool Inside(const Point3<T>& p, const Bounds3<T>& b);

        template <typename T>
        bool InsideExclusive(const Point3<T>& p, const Bounds3<T>& b);

        template <typename T, typename U>
        inline Bounds3<T> Expand(const Bounds3<T>& b, U delta){
            return Bounds3<T>(b.pMin - Vector3<U>(delta, delta, delta),
                              b.pMax + Vector3<U>(delta, delta, delta));
        }


//------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------//

class boundIterator : public std::forward_iterator_tag{

    public:

        // class constructor

        boundIterator(const Bounds2i& b, const Point2i& p);

        //operator overloading

        boundIterator operator++();
        boundIterator operator++(int);
        bool operator==(const boundIterator& bi);
        bool operator!=(const boundIterator& bi);
        Point2i operator*() const ;

    private:
        void advance();

        const Bounds2i* b;
        Point2i pt;
};

//------------------------------------------------------------------------------------//


#endif
