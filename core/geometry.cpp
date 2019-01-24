#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>

#include "geometry.h"



// Global inline //
//----------------------------------------------//
inline Float Lerp(Float t, Float v1, Float v2){
    return (((1 - t) * v1) + (t * v2));
}
//----------------------------------------------//

//----------------------------------------------//
// Vector2 class method //

template <class T> 
T Vector2<T>::operator[](int i)const{
    assert( (i == 0) || (i == 1) );
    if(i == 0){return x;}
    if(i == 1){return y;}
}

template <class T> 
T& Vector2<T>::operator[](int i){
    assert( (i == 0) || (i == 1) );
    if(i == 0){return x;}
    if(i == 1){return y;}
}

template <class T>
bool Vector2<T>::hasNan()const{
    return (std::isnan(x) || std::isnan(y));
}

template <>
std::string Vector2<Float>::tostring(){
    char text[50];
    sprintf(text, "[%.3f, %.3f]", x,y);
    return std::string(text);
}


template <>
std::string Vector2<int>::tostring(){
    char text[50];
    sprintf(text, "[%d, %d]", x,y);
    return std::string(text);
}

// Vector2 related function //

template <typename T> 
inline Vector2<T> operator*(T s, Vector2<T> v){
    return Vector2<T>(s * v.x, s * v.y);
}

template <typename T>
inline T Dot(const Vector2<T>& a, const Vector2<T>& b){return a.x * b.x + a.y * b.y;}

template <typename T>
inline T AbsDot(const Vector2<T>& a, const Vector2<T>& b){return std::abs(dot(a, b));}

template <typename T>
inline Vector2<T> Cross(const Vector2<T>& a, const Vector2<T>& b){

    Float ax = a.x, ay = a.y, az = a.z;
    Float bx = a.x, by = b.y, bz = b.z;
    return Vector2<T>((ax * by - ay * bx), (ay * bx - ax * by));
}

template <typename T>
inline Vector2<T> Normalize(const Vector2<T>& v){return (v/v.lenght());}

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

//----------------------------------------------//


//----------------------------------------------//
// Vector3 class method //

template <class T>
bool Vector3<T>::hasNan()const{
    return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
}


template <class T>
T Vector3<T>::operator[](int i){
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <class T>
T& Vector3<T>::operator[](int i)const{
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <>
std::string Vector3<Float>::tostring(){
    char text[50];
    sprintf(text, "[%.3f, %.3f, %.3f]", x, y, z);
    return std::string(text);
}


template <>
std::string Vector3<int>::tostring(){
    char text[50];
    sprintf(text, "[%d, %d, %d]", x, y, z);
    return std::string(text);
}

// Vector3 related function //

template <typename T> 
inline Vector3<T> operator*(T s, Vector3<T> v){
    return Vector3<T>(s * v.x, s * v.y, s * v.z);
}

template <typename T>
inline T Dot(const Vector3<T>& a, const Vector3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}

template <typename T>
inline T AbsDot(const Vector3<T>& a, const Vector3<T>& b){return std::abs(dot(a, b));}

template <typename T>
inline Vector3<T> Cross(const Vector3<T>& a, const Vector3<T>& b){

    Float ax = a.x, ay = a.y, az = a.z;
    Float bx = a.x, by = b.y, bz = b.z;
    return Vector3<T>((ay * bz - az * by), (az * bx - ax * bz), (ax * by - ay * bx));
}

template <typename T>
inline Vector3<T> Normalize(const Vector3<T>& v){return (v/v.lenght());}

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

//----------------------------------------------//



//----------------------------------------------//
// Point2 class method //

template <class T>
bool Point2<T>::hasNan()const{
    return (std::isnan(x) || std::isnan(y)) ;
}

template <class T>
T  Point2<T>::operator[](int i)const{
    assert(( i == 0 ) || ( i == 1 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
}

template <class T>
T& Point2<T>::operator[](int i){
    assert(( i == 0 ) || ( i == 1 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
}

template <>
std::string Point2<Float>::tostring(){
    char text[50];
    sprintf(text, "(%.3f, %.3f)", x,y);
    return std::string(text);
}


template <>
std::string Point2<int>::tostring(){
    char text[50];
    sprintf(text, "(%d, %d)", x,y);
    return std::string(text);
}

// Point2 related function //


template <typename T>
inline Float Distance(const Point2<T>& p1, const Point2<T>& p2){
    return (p1 - p2).lenght();
}

template <typename T>
inline Float DistanceSquared(const Point2<T>& p1, const Point2<T>& p2){
    return (p1 - p2).lenghtSquared();
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

//----------------------------------------------//



//----------------------------------------------//
// Point3 class method //
template <class T>
bool Point3<T>::hasNan()const{
    return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
}

template <class T>
T  Point3<T>::operator[](int i)const{
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <class T>
T& Point3<T>::operator[](int i){
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <>
std::string Point3<Float>::tostring(){
    char text[50];
    sprintf(text, "(%.3f, %.3f, %.3f)", x,y,z);
    return std::string(text);
}


template <>
std::string Point3<int>::tostring(){
    char text[50];
    sprintf(text, "(%d, %d, %d)", x,y,z);
    return std::string(text);
}

// Point3 related function //

template <typename T>
inline Float Distance(const Point3<T>& p1, const Point3<T>& p2){
    return (p1 - p2).lenght();
}

template <typename T>
inline Float DistanceSquared(const Point3<T>& p1, const Point3<T>& p2){
    return (p1 - p2).lenghtSquared();
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
//----------------------------------------------//


//----------------------------------------------//
// Normal3 class method //

template <class T>
bool Normal3<T>::hasNan()const{
    return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
}
template <class T>
T Normal3<T>::operator[](int i){
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <class T>
T& Normal3<T>::operator[](int i)const{
    assert(( i >= 0 ) && ( i <= 2 ));
    if(i == 0){return x;}
    if(i == 1){return y;}
    if(i == 2){return z;}
}

template <>
std::string Normal3<Float>::tostring(){
    char text[50];
    sprintf(text, "|%.3f, %.3f, %.3f|", x, y, z);
    return std::string(text);
}

template <>
std::string Normal3<int>::tostring(){
    char text[50];
    sprintf(text, "|%d, %d, %d|", x, y, z);
    return std::string(text);
}

// Normal3 related function //
//----------------------------------------------//

template <typename T> 
inline Normal3<T> operator*(T s, Normal3<T> v){
    return Normal3<T>(s * v.x, s * v.y, s * v.z);
}


template <typename T>
inline Normal3<T> Normalize(const Normal3<T>& v){return (v/v.lenght());}

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
//----------------------------------------------//


//----------------------------------------------//
// RayDifferential class method //
    void RayDifferential::scaleDifferentials(Float s){
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    };
// RayDifferential related function //
//----------------------------------------------//


//----------------------------------------------//
// Bounds3 class method //
template <typename T>
const Point3<T>& Bounds3<T>::operator[](int i)const{
    assert((i == 0) || (i == 1));
    return ((i==0)?pMin:pMax);
}

template <typename T>
Point3<T>& Bounds3<T>::operator[](int i){
    assert((i == 0) || (i == 1));
    return ((i==0)?pMin:pMax);
}

template <typename T>
Point3<T> Bounds3<T>::corner(int i) const{
    return Point3<T>((*this)[(i & 1)].x,
                     (*this)[(i & 2)? 1 : 0].y,
                     (*this)[(i & 4)? 1 : 0].z);
}

template <typename T>
Vector3<T> Bounds3<T>::diagonal()const{
    return (pMax - pMin);
}

template <typename T>
T Bounds3<T>::surfaceArea()const{
    Vector3<T> v = diagonal();
    return (v.x * v.y + v.x * v.z + v.y * v.z);
}

template <typename T>
T Bounds3<T>::volume()const{
    Vector3<T> v = diagonal();
    return (v.x * v.y * v.z);
}

template <typename T>
int Bounds3<T>::maximumExtent()const{
    Vector3<T> v = diagonal();
    if((v.x > v.y) && (v.x > v.z)){return 0;}
    return (v.y > v.z) ? 1 : 2;
}


// Bounds3 related function //
template <typename T>
Bounds3<T> Union(const Bounds3<T>& b, const Point3<T>& p){
    return Bounds3<T>(Point3<T>(std::min(b.pMin.x, p.x), 
                                std::min(b.pMin.y, p.y),
                                std::min(b.pMin.z, p.z)),
                      Point3<T>(std::max(b.pMax.x, p.x),
                                std::max(b.pMax.y, p.y),
                                std::max(b.pMax.z, p.z)));
}

template <typename T> 
Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2){
    return Bounds3<T>(Point3<T>(std::min(b1.pMin.x, b2.pMin.x), 
                                std::min(b1.pMin.y, b2.pMin.y),
                                std::min(b1.pMin.z, b2.pMin.z)),
                      Point3<T>(std::max(b1.pMax.x, b2.pMax.x),
                                std::max(b1.pMax.y, b2.pMax.y),
                                std::max(b1.pMax.z, b2.pMax.z)));    
}

template <typename T> 
Bounds3<T> Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2){
    return Bounds3<T>(Point3<T>(std::max(b1.pMin.x, b2.pMin.x), 
                                std::max(b1.pMin.y, b2.pMin.y),
                                std::max(b1.pMin.z, b2.pMin.z)),
                      Point3<T>(std::min(b1.pMax.x, b2.pMax.x),
                                std::min(b1.pMax.y, b2.pMax.y),
                                std::min(b1.pMax.z, b2.pMax.z)));    
}

template <typename T>
bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2){
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
} 

template <typename T>
bool Inside(const Point3<T>& p, const Bounds3<T>& b){
    bool x = ( (p.x >= b.pMin.x) && (p.x <= b.pMax.x) );
    bool y = ( (p.y >= b.pMin.y) && (p.y <= b.pMax.y) );
    bool z = ( (p.z >= b.pMin.z) && (p.z <= b.pMax.z) );

    return (x && y && z);
}

template <typename T>
bool InsideExclusive(const Point3<T>& p, const Bounds3<T>& b){
    bool x = ( (p.x >= b.pMin.x) && (p.x < b.pMax.x) );
    bool y = ( (p.y >= b.pMin.y) && (p.y < b.pMax.y) );
    bool z = ( (p.z >= b.pMin.z) && (p.z < b.pMax.z) );

    return (x && y && z);
}

template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T>& b, U delta){
    return Bounds3<T>(b.pMin - Vector3<U>(delta, delta, delta),
                      b.pMax + Vector3<U>(delta, delta, delta));
}

template <typename T>
Point3<T> Bounds3<T>::lerp(const Point3<Float>& t) const{
    return Point3<T>(Lerp(t.x, pMin.x, pMax.x), 
                     Lerp(t.y, pMin.y, pMax.y),
                     Lerp(t.z, pMin.y, pMax.z));
}

template <typename T>
Vector3<T> Bounds3<T>::offset(const Point3<T>& p) const{
    Vector3<T> o = p - pMin;
    if(pMax.x > pMin.x){ o.x /= (pMax.x - pMin.x); }
    if(pMax.y > pMin.y){ o.y /= (pMax.y - pMin.y); }
    if(pMax.z > pMin.z){ o.z /= (pMax.z - pMin.z); }
    return o;
}

template <typename T>
void Bounds3<T>::boundingSphere(Point3<T>* center, Float* radius)const{
    *center = pMax - pMin;
    *radius =  Inside(*center, *this) ? Distance(*center, pMax) : 0;
    return;
}
//----------------------------------------------//

int main(int argc, char** argv){

    
    Normal3f n = Normal3f(5.0, 2.3, 25);
    Vect3f v = Vect3f(n);
    Point3f p1 = Point3f(0.7, 0.1, 0.2);
    Point3f p2 = Point3f(1.7, 10, 0.3);
    Ray r = Ray(p1, v);
    RayDifferential rd = RayDifferential(r);
    std::cout << rd.o << rd.d << std::endl;

    Bounds3f b1 = Bounds3f(p1, p2);
    std::cout << b1.maximumExtent() << std::endl;

    std::cout<<"hello"<<std::endl;
    //std::cout << p2 << std::endl;

    return EXIT_SUCCESS;
}
