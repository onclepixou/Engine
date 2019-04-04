#include "geometry.h"


//VECTOR2
    //----------------------------------------------//
    // Vector2 class constructor
    template <class T>
    Vector2<T>::Vector2() : x(0), y(0){};

    template <class T>
    Vector2<T>::Vector2(T a, T b) : x(a), y(b){
        assert(!hasNan());
    }

    template <class T>
    Vector2<T>::Vector2(const Vector3<T>& v) : x(v.x), y(v.y){
        assert(!hasNan());
    };

    template <class T>
    template <typename U> 
    Vector2<T>::operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }



    // Vector2 class operator overload //

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

    template <typename T>
    bool Vector2<T>::operator==(const Vector2<T>& v){
        if(( x == v.x ) && ( y == v.y ))
            return true;
        return false;
    }

    template <typename T>
    bool Vector2<T>::operator!=(const Vector2<T>& v){
        if(( x != v.x ) || ( y != v.y ))
            return true;
        return false;
    }

    template <class T>
    Vector2<T> Vector2<T>::operator+(const Vector2<T>& rhs) const{ 
        return Vector2<T>(x + rhs.x, y + rhs.y); 
    } 

    template <class T>
    Vector2<T>& Vector2<T>::operator+=(const Vector2<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    template <class T>
    Vector2<T> Vector2<T>::operator-(const Vector2<T>& rhs) const{ 
        return Vector2<T>(x - rhs.x, y - rhs.y); 
    }

    template <class T>
    Vector2<T>& Vector2<T>::operator-=(const Vector2<T>& rhs) { 
        x -= rhs.x; 
        y -= rhs.y;
        return *this;
    }

    template <class T>
    Vector2<T> Vector2<T>::operator*(T s) const {
        return Vector2<T>( s * x, s * y);
    }

    template <class T>
    Vector2<T>& Vector2<T>::operator*=(T s){  
        x = s * x; 
        y = s * y; 
        return *this;
    }

    template <class T>
    Vector2<T> Vector2<T>::operator/(T i)const{
        assert(i != 0); 
        T inv = 1/i; 
        return Vector2<T>(x * inv, y * inv);
    }

    template <class T>
    Vector2<T>& Vector2<T>::operator/=(T i){
        assert(i != 0); 
        T inv = 1/i; 
        x *= inv; 
        y *= inv; 
        return *this;
    }

    template <class T>
    Vector2<T> Vector2<T>::operator-(){
        return Vector2<T>(-x, -y);
    }

    // Vector2 class method //

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

    template <class T>
    Vector2<T> Vector2<T>::abs(){
        return Vector2<T>(std::abs(x), std::abs(y));
    }

    template <class T>
    Float Vector2<T>::lenghtSquared()const{
        return x * x + y * y;
    }

    template <class T>
    Float Vector2<T>::lenght()const{
        return std::sqrt(lenghtSquared());
    }

    template <class T>
    bool Vector2<T>::isNormalized()const{
        return ((lenght() - 1 ) <= 1e-6);
    }

//----------------------------------------------//

//VECTOR3

    // Vector3 class constructor //

    template <typename T>
    Vector3<T>::Vector3() : x(0), y(0), z(0){};

    template <typename T>
    Vector3<T>::Vector3(T a, T b, T c) : x(a), y(b), z(c){assert(!hasNan());};

    template <typename T>
    Vector3<T>::Vector3(const Vector2<T>& v) : x(v.x), y(v.y), z(0) {assert(!hasNan());};

    template <typename T>
    Vector3<T>::Vector3(const Normal3<T> &n): x(n.x), y(n.y), z(n.z){assert(!hasNan());};

    template <typename T>
    template <typename U> 
    Vector3<T>::operator Vector3<U>() const {return Vector3<U>(x, y, z);}

    // Vector3 class operator overload //

    template <typename T>
    T Vector3<T>::operator[](int i)const{
        assert(( i >= 0 ) && ( i <= 2 ));
        if(i == 0){return x;}
        if(i == 1){return y;}
        if(i == 2){return z;}
    }

    template <typename T>
    T& Vector3<T>::operator[](int i){
        assert(( i >= 0 ) && ( i <= 2 ));
        if(i == 0){return x;}
        if(i == 1){return y;}
        if(i == 2){return z;}
    }

    template <typename T>
    bool Vector3<T>::operator==(const Vector3<T>& v){
        if(( x == v.x ) && ( y == v.y ) && ( z == v.z))
            return true;
        return false;
    }

    template <typename T>
    bool Vector3<T>::operator!=(const Vector3<T>& v){
        if(( x != v.x ) || ( y != v.y ) || ( z != v.z))
            return true;
        return false;
    }

    template <typename T>
    Vector3<T> Vector3<T>::operator+(const Vector3<T>& rhs) const{ 
        return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z); 
    }

    template <typename T>
    Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& rhs) {
         x += rhs.x; y += rhs.y; z+= rhs.z;return *this;
    }

    template <typename T>
    Vector3<T> Vector3<T>::operator-(const Vector3<T>& rhs) const{ 
        return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z); 
    }

    template <typename T>
    Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& rhs) {
         x -= rhs.x; y -= rhs.y; z-= rhs.z;return *this;
    }

    template <typename T>
    Vector3<T> Vector3<T>::operator*(T s) const {
        return Vector3<T>( s * x, s * y, s * z);
    }

    template <typename T>
    Vector3<T>& Vector3<T>::operator*=(T s){
        x = s * x; y = s * y;
        z = s * z; 
        return *this;
    }

    template <typename T>
    Vector3<T> Vector3<T>::operator/(T i)const{
        assert(i != 0);
        T inv = 1/i; 
        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    template <typename T>
    Vector3<T>& Vector3<T>::operator/=(T i){
        assert(i != 0); 
        T inv = 1/i; 
        x *= inv; 
        y *= inv; 
        z *= inv; 
        return *this;
    }

    template <typename T>
    Vector3<T> Vector3<T>::operator-(){
        return Vector3<T>(-x, -y, -z);
    }

    // Vector3 class method //

    template <typename T>
    bool Vector3<T>::hasNan()const{
        return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
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

    template <typename T>
    Vector3<T> Vector3<T>::abs(){
        return Vector3<T>(std::abs(x), std::abs(y), std::abs(z));
    }

    template <typename T>
    Float Vector3<T>::lenghtSquared()const{
        return x * x + y * y + z * z;
    }

    template <typename T>
    Float Vector3<T>::lenght()const{
        return std::sqrt(lenghtSquared());
    }

    template <typename T>
    bool Vector3<T>::isNormalized(){
        return ((lenght() - 1 ) <= 1e-6);
    }
//----------------------------------------------//

//POINT2
    //----------------------------------------------//

    // Point2 class constructor

    template <class T>
    Point2<T>::Point2():x(0), y(0){};

    template <class T>
    Point2<T>::Point2(T a, T b): x(a), y(b){ assert(!hasNan());};

    template <class T>
    Point2<T>::Point2(const Point3<T>& p): x(p.x), y(p.y){
        assert(!hasNan());
    }

    template <class T>
    template <typename U> 
    Point2<T>::Point2(const Point2<U> p): x((T)p.x), y((T)p.y){
        assert(!hasNan());
    }

    // Point2 class operator overloading



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

    template <class T>
    bool Point2<T>::operator==(const Point3<T>& p){
        return ((x == p.x) && ( y == p.y ));
    }

    template <class T>
    bool Point2<T>::operator!=(const Point3<T>& p){
        return ((x != p.x) || ( y != p.y ));
    }

    template <class T>
    Point2<T> Point2<T>::operator+(const Point2<T>& p) const{
        return Point2<T>((p.x + x), (p.y + y));
    }

    template <class T>
    Point2<T>& Point2<T>::operator+=(const Point2<T>& p){
        x += p.x; 
        y += p.y; 
        return *this;
    }

    template <class T>
    Point2<T> Point2<T>::operator+(const Vector2<T>& v)const{
        return Point2<T>((v.x + x), (v.y + y));
    }

    template <class T>
    Point2<T>& Point2<T>::operator+=(const Vector2<T>& v){
        x += v.x; 
        y += v.y; 
        return (*this); 
    }

    template <class T>
    Vector2<T> Point2<T>::operator-(const Point2<T>& p)const{
        return Vector2<T>((x-p.x), (y-p.y));
    }

    template <class T>
    Point2<T> Point2<T>::operator-(const Vector2<T>& v)const{
        return Point2<T>((x-v.x), (y-v.y));
    } 

    template <class T>
    Point2<T>& Point2<T>::operator-=(const Vector2<T>& v){
        x -= v.x; 
        y -= v.y; 
        return (*this);
    }

    template <class T>
    Point2<T> Point2<T>::operator*(const T& s) const {
        return Point2<T>( s * x, s * y);
    } 

    template <class T>
    Point2<T>& Point2<T>::operator*=(const T& s){ 
        x = s * x; 
        y = s * y; 
        return *this;
    }

    template <class T>
    Point2<T> Point2<T>::operator/(const T& n)const{
        return Point2<T>(x/n, y/n);
    }
    
    template <class T>
    Point2<T>& Point2<T>::operator/=(const T& n){
        x/=n;
        y/=n;
        return *this;
    }


    // Point2 class method //

    template <class T>
    bool Point2<T>::hasNan()const{
        return (std::isnan(x) || std::isnan(y)) ;
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

//----------------------------------------------//

//POINT3
    //----------------------------------------------//
    // Point3 class constructor //

    template <typename T>
    Point3<T>::Point3():x(0), y(0), z(0){};

    template <typename T>
    Point3<T>::Point3(T a, T b, T c): x(a), y(b), z(c){ 
        assert(!hasNan());
    }

    template <typename T>
    Point3<T>::Point3(const Point2<T>& p):x(p.x), y(p.y), z(0){
        assert(!hasNan());
    }

    template <typename T>
    template <typename U> 
    Point3<T>::Point3(const Point3<U> p): x((T)p.x), y((T)p.y), z((T)p.z){
        assert(!hasNan());
    }



    // Point3 class operator overload //


    template <class T>
    Point3<T> Point3<T>::operator+(const Point3<T>& p)const{
        return Point3<T>((p.x + x), (p.y + y), (p.z + z));
    }

    template <class T>
    Point3<T>& Point3<T>::operator+=(const Point3<T>& p){
        x += p.x; 
        y += p.y; 
        z += p.z; 
        return (*this); 
    }

    template <class T>
    Point3<T> Point3<T>::operator+(const Vector3<T>& v)const{
        return Point3<T>((v.x + x), (v.y + y), (v.z + z));
    } 


    template <class T>
    Point3<T>& Point3<T>::operator+=(const Vector3<T>& v){
        x += v.x; 
        y += v.y; 
        z += v.z; 
        return (*this); 
    }

    template <class T>
    Vector3<T> Point3<T>::operator-(const Point3<T>& p)const{
        return Vector3<T>((x-p.x), (y-p.y), (z-p.z));
    }

    template <class T>
    Point3<T> Point3<T>::operator-(const Vector3<T>& v)const{
        return Point3<T>((x-v.x), (y-v.y), (z-v.z));
    }

    template <class T>
    Point3<T>& Point3<T>::operator-=(const Vector3<T>& v){
        x -= v.x; 
        y -= v.y; 
        z -= v.z; 
        return (*this);
    }

    template <class T>
    Point3<T> Point3<T>::operator*(const T& s) const {
        return Point3<T>( s * x, s * y, s * z);
    }

    template <class T>
    Point3<T>& Point3<T>::operator*=(const T& s){ 
        x = s * x; 
        y = s * y; 
        z = s * z; 
        return *this;
    }

    template <class T>
    bool Point3<T>::operator==(const Point3<T>& p){
        return ((x == p.x) && ( y == p.y ) && (z == p.z));
    }

    template <class T>
    bool Point3<T>::operator!=(const Point3<T>& p){
        return ((x != p.x) || ( y != p.y ) || (z != p.z));
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

    template <class T>
    Point3<T> Point3<T>::operator/(const T& n)const{
        return Point3<T>(x/n, y/n, z/n);
    }
    
    template <class T>
    Point3<T>& Point3<T>::operator/=(const T& n){
        x/=n;
        y/=n;
        z/=n;
        return *this;
    }


    // Point3 class method //

    template <class T>
    bool Point3<T>::hasNan()const{
        return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
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
//----------------------------------------------//

//NORMAL3
    //----------------------------------------------//

    // Normal3 class constructor //

    template <class T>
    Normal3<T>::Normal3() : x(0), y(0), z(0){};

    template <class T>
    Normal3<T>::Normal3(T a, T b, T c) : x(a), y(b), z(c){
        assert(!hasNan());
    }

    template <class T>
    Normal3<T>::Normal3(const Vector3<T>& v) : x(v.x), y(v.y), z(v.z){
        assert(!hasNan());
    }

    template <class T>
    template <typename U> 
    Normal3<T>::operator Normal3<U>() const {
        return Normal3<U>(x, y, z);
    }



    // Normal3 class operator overload //

    template <class T>
    T Normal3<T>::operator[](int i)const{
        assert(( i >= 0 ) && ( i <= 2 ));
        if(i == 0){return x;}
        if(i == 1){return y;}
        if(i == 2){return z;}
    }

    template <class T>
    T& Normal3<T>::operator[](int i){
        assert(( i >= 0 ) && ( i <= 2 ));
        if(i == 0){return x;}
        if(i == 1){return y;}
        if(i == 2){return z;}
    }

    template <class T>
    bool Normal3<T>::operator==(const Normal3<T>& n){
        if(( x == n.x ) && ( y == n.y ) && ( z == n.z))
            return true;
        return false;
    }

    template <class T>
    bool Normal3<T>::operator!=(const Normal3<T>& n){
        if(( x != n.x ) || ( y != n.y ) || ( z != n.z))
            return true;
        return false;
    }

    template <class T>
    Normal3<T> Normal3<T>::operator+(const Normal3<T>& rhs) const{ 
        return Normal3<T>(x + rhs.x, y + rhs.y, z + rhs.z); 
    }

    template <class T>
    Normal3<T>& Normal3<T>::operator+=(const Normal3<T>& rhs) {
         x += rhs.x; 
         y += rhs.y; 
         z+= rhs.z;
         return *this;
    }

    template <class T>
    Normal3<T> Normal3<T>::operator-(const Normal3<T>& rhs) const{ 
        return Normal3<T>(x - rhs.x, y - rhs.y, z - rhs.z); 
    }

    template <class T>
    Normal3<T>& Normal3<T>::operator-=(const Normal3<T>& rhs) { 
        x -= rhs.x; 
        y -= rhs.y; 
        z-= rhs.z;
        return *this;
    }


    template <class T>
    Normal3<T> Normal3<T>::operator*(T s) const {
        return Normal3<T>( s * x, s * y, s * z);
    }

    template <class T>
    Normal3<T>& Normal3<T>::operator*=(T s){ 
        x = s * x; 
        y = s * y; 
        z = s * z; 
        return *this;
    }

    template <class T>
    Normal3<T> Normal3<T>::operator/(T i)const{
        assert(i != 0); 
        T inv = 1/i; 
        return Normal3<T>(x * inv, y * inv, z * inv);
    }


    template <class T>
    Normal3<T>& Normal3<T>::operator/=(T i){
        assert(i != 0); 
        T inv = 1/i; 
        x *= inv; 
        y *= inv; 
        z *= inv; 
        return *this;
    }

    template <class T>
    Normal3<T> Normal3<T>::operator-(){
        return Normal3<T>(-x, -y, -z);
    }

    // Normal3 class method //

    template <class T>
    bool Normal3<T>::hasNan()const{
        return (std::isnan(x) || std::isnan(y) || std::isnan(z) ) ;
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

    template <class T>
    Normal3<T> Normal3<T>::abs(){
        return Normal3<T>(std::abs(x), std::abs(y), std::abs(z));
    }

    template <class T>
    Float Normal3<T>::lenghtSquared()const{
        return x * x + y * y + z * z;
    }

    template <class T>
    Float Normal3<T>::lenght()const{
        return std::sqrt(lenghtSquared());
    }

    template <class T>
    bool Normal3<T>::isNormalized(){
        return ((lenght() - 1 ) <= 1e-6);
    }


//----------------------------------------------//

//RAY
    //----------------------------------------------//

    // RAY class constructor //

    Ray::Ray():tMax(Infinity), time(0.0f), medium(nullptr){};
    Ray::Ray(const Point3<Float>& o, const Vector3<Float>& d, Float tMax,
            Float time, const Medium* medium) : o(o), d(d), tMax(tMax), time(time), medium(medium){};

    // RAY class operator overloading

    Point3<Float> Ray::operator()(Float t){return o + (d * t);}



//----------------------------------------------//

//RAYDIF
    //----------------------------------------------//

    // RayDifferential class constructor

    RayDifferential::RayDifferential(){hasDifferentials = false;};

    RayDifferential::RayDifferential(const Point3<Float>& o, const Vector3<Float>& d, 
                        Float tMax , Float time , const Medium* medium )
                        : Ray(o, d, tMax, time, medium){hasDifferentials = false;};

//    RayDifferential::RayDifferential(const Point3<Float>& o, const Vector3<Float>& d, 
//                        Float tMax = Infinity, Float time = 0.0f, const Medium* medium = nullptr)
//                        : Ray(o, d, tMax, time, medium){hasDifferentials = false;};

    RayDifferential::RayDifferential(const Ray& r) : Ray(r){hasDifferentials = false; };


    // RayDifferential class method //

    void RayDifferential::scaleDifferentials(Float s){
            rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    };

    // RayDifferential related function //
//----------------------------------------------//

//BOUNDS2
    //----------------------------------------------//

    // Bounds2 class constructor //

        template <class T>
        Bounds2<T>::Bounds2(){
            T minNum = std::numeric_limits<T>::lowest();
            T maxNum = std::numeric_limits<T>::max();
            pMin = Point2<T>(maxNum, maxNum);
            pMax = Point2<T>(minNum, minNum);
        }

        template <class T>
        Bounds2<T>::Bounds2(const Point2<T>& p) : pMin(p), pMax(p){};

        template <class T>
        Bounds2<T>::Bounds2(const Point2<T>& p1, const Point2<T>& p2): 
            pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y)),
            pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y)){};

    // Bounds2 operator overloading //

        template <typename T>
        const Point2<T>& Bounds2<T>::operator[](int i)const{
            assert((i == 0) || (i == 1));
            return ((i==0)?pMin:pMax);
        }

        template <typename T>
        Point2<T>& Bounds2<T>::operator[](int i){
            assert((i == 0) || (i == 1));
            return ((i==0)?pMin:pMax);
        }

        template <typename T>
        bool Bounds2<T>::operator==(const Bounds2<T>& b){
            return (( pMin == b.pMin ) && ( pMax == b.pMax));
        }

        template <typename T>
        bool Bounds2<T>::operator!=(const Bounds2<T>& b){
            return (( pMin != b.pMin ) || ( pMax != b.pMax));
        }

    // Bounds2 class method

        template <typename T>
        Point2<T> Bounds2<T>::corner(int i) const{
            assert((i >= 0) && (i <= 3));
            return Point2<T>((*this)[(i & 1)].x,
                             (*this)[(i & 2)? 1 : 0].y);
        }

        template <typename T>
        Vector2<T> Bounds2<T>::diagonal()const{
            return (pMax - pMin);
        }

        template <typename T>
        T Bounds2<T>::surfaceArea()const{
            Vector2<T> v = diagonal();
            return (v.x * v.y );
        }

        template <typename T>
        int Bounds2<T>::maximumExtent()const{
            Vector2<T> v = diagonal();
            return ((v.x > v.y)?0:1);
        }

        template <typename T>
        Point2<T> Bounds2<T>::lerp(const Point2<Float>& t) const{
            return Point2<T>(Lerp(t.x, pMin.x, pMax.x), 
                             Lerp(t.y, pMin.y, pMax.y));
        }

        template <typename T>
        Vector2<T> Bounds2<T>::offset(const Point2<T>& p) const{
            Vector2<T> o = p - pMin;
            if(pMax.x > pMin.x){ o.x /= (pMax.x - pMin.x); }
            if(pMax.y > pMin.y){ o.y /= (pMax.y - pMin.y); }
            return o;
        }

        template <typename T>
        void Bounds2<T>::boundingCircle(Point2<T>* center, Float* radius)const{
            *center = (pMax + pMin)/2;
            *radius =  Inside(*center, *this) ? Distance(*center, pMax) : 0;
            return;    
        }

    // Bounds2 related function //
    //----------------------------------------------//

        template <typename T>
        Bounds2<T> Union(const Bounds2<T>& b, const Point2<T>& p){
            return Bounds2<T>(Point2<T>(std::min(b.pMin.x, p.x), 
                                        std::min(b.pMin.y, p.y)),
                              Point2<T>(std::max(b.pMax.x, p.x),
                                        std::max(b.pMax.y, p.y)));
        }

        template <typename T> 
        Bounds2<T> Union(const Bounds2<T>& b1, const Bounds2<T>& b2){
            return Bounds2<T>(Point2<T>(std::min(b1.pMin.x, b2.pMin.x), 
                                        std::min(b1.pMin.y, b2.pMin.y)),
                              Point2<T>(std::max(b1.pMax.x, b2.pMax.x),
                                        std::max(b1.pMax.y, b2.pMax.y)));    
        }

        template <typename T> 
        Bounds2<T> Intersect(const Bounds2<T>& b1, const Bounds2<T>& b2){
            return Bounds3<T>(Point3<T>(std::max(b1.pMin.x, b2.pMin.x), 
                                        std::max(b1.pMin.y, b2.pMin.y)),
                              Point3<T>(std::min(b1.pMax.x, b2.pMax.x),
                                        std::min(b1.pMax.y, b2.pMax.y)));  
        }

        template <typename T>
        bool Overlaps(const Bounds2<T>& b1, const Bounds2<T>& b2){
            bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
            bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
            return (x && y );
        } 

        template <typename T>
        bool Inside(const Point2<T>& p, const Bounds2<T>& b){
            bool x = ( (p.x >= b.pMin.x) && (p.x <= b.pMax.x) );
            bool y = ( (p.y >= b.pMin.y) && (p.y <= b.pMax.y) );
            return (x && y );
        }

        template <typename T>
        bool InsideExclusive(const Point2<T>& p, const Bounds2<T>& b){
            bool x = ( (p.x >= b.pMin.x) && (p.x < b.pMax.x) );
            bool y = ( (p.y >= b.pMin.y) && (p.y < b.pMax.y) );
            return (x && y );
        }

//----------------------------------------------//

//BOUNDS3
    // Bounds3 class method //
    //----------------------------------------------//

        template <typename T>
        Bounds3<T>::Bounds3(){
                T minNum = std::numeric_limits<T>::lowest();
                T maxNum = std::numeric_limits<T>::max();
                pMin = Point3<T>(maxNum, maxNum, maxNum);
                pMax = Point3<T>(minNum, minNum, minNum);
        }  

        template <typename T>
        Bounds3<T>::Bounds3(const Point3<T>& p) : pMin(p), pMax(p){};

        template <typename T>
        Bounds3<T>::Bounds3(const Point3<T>& p1, const Point3<T>& p2)
             : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
               pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)){};  


    // Bounds3 operator overloading//

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
        bool Bounds3<T>::operator==(const Bounds3<T>& b){return (( pMin == b.pMin ) && ( pMax == b.pMax));}

        template <typename T>
        bool Bounds3<T>::operator!=(const Bounds3<T>& b){return (( pMin != b.pMin ) || ( pMax != b.pMax));}


        template <typename T>
        Point3<T> Bounds3<T>::corner(int i) const{
            assert((i >= 0) && (i <= 7));
            return Point3<T>((*this)[(i & 1)].x,
                             (*this)[(i & 2)? 1 : 0].y,
                             (*this)[(i & 4)? 1 : 0].z);
        }

    // Bounds3 public method

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
            *center = (pMax + pMin)/2;
            *radius =  Inside(*center, *this) ? Distance(*center, pMax) : 0;
            return;
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
//----------------------------------------------//

//BOUND_IT

    // BoundIterator class constructor //

    boundIterator::boundIterator(const Bounds2i& b, const Point2i& p) : b(&b), pt(p){};

    // BoundIterator class method //

        void boundIterator::advance(){
            ++pt.x;
            if(pt.x == b->pMax.x){
                pt.x = b->pMin.x;
                ++pt.y;
            }
        }


    // BoundIterator operator overload //

        boundIterator boundIterator::operator++(){
            advance();
            return *this;
        }

        boundIterator boundIterator::operator++(int){
            boundIterator old = *this;
            advance();
            return old;
        }

        bool boundIterator::operator==(const boundIterator& bi){
            return ((pt == bi.pt) && (b == bi.b));
        }

        bool boundIterator::operator!=(const boundIterator& bi){
            return ((pt != bi.pt) || (b != bi.b));
        }

        Point2i boundIterator::operator*() const { return pt; }

    // BoundIterator related function //

    inline boundIterator begin(const Bounds2i& b){
        return boundIterator(b, b.pMin);
    }

    inline boundIterator end(const Bounds2i& b){
        Point2i end  = Point2i(b.pMin.x, b.pMax.y);
        if((b.pMin.x >= b.pMax.x)||(b.pMin.y >= b.pMax.y)){
            end = b.pMin;
        }
        return boundIterator(b, end);
    }
//----------------------------------------------//


template class Vector3<int>;
template class Vector3<Float>;
template class Vector2<int>;
template class Vector2<Float>;

template class Point3<int>;
template class Point3<Float>;
template class Point2<int>;
template class Point2<Float>;

template class Normal3<int>;
template class Normal3<Float>;

template class Bounds3<int>;
template class Bounds3<Float>;
template class Bounds2<int>;
template class Bounds2<Float>;



