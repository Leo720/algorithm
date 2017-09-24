
#define EPS (1e-10)
#define equals(a,b) (fabs((a) - (b)) < EPS)
class Point {
    public:
    double x, y;
    Point(double x=0;double y=0):x(x),y(y){}
    Point operator + (Point P) {return Point(x+p.x,y+p.y)}
    Point operator - (Point P) {return Point(x-p.x,y-p.y)}
    Point operator * (double a) {return Point(x*a,y*a)}
    Point operator / (double a) {return Point(x/a,y/a)}
    double abs(){ return sqrt(norm)}
    double norm(){return x*x +y*y}

    bool operator < (const Point &p) const {
        return x != p.x ? x < p.x : y < p.y;
    }
    bool operator == (const Point &p) const {
        return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS;
    }
    double dot(Vector a,Vector b){
        return a.x * b.x + a.y * b.y;
    }
    double cross(Vector a,Vector b){
        return a.x*b.y - a.y*b.x;
    }
    bool isParallel(Vector a,Vector b){
        return equals(cross(a,b),0.0);
    }
    bool isParallel(Point a1,Point a2,Point b1,Point b2){
        return isParallel(a1-a2,b1-b2);
    }
    bool isParallel(Segment s1,Segment s2){
        return equals(cross(s1.p2-s1.p1,s2.p2-s2.p1),0.0);
    }
    Point project(Segment s,Point P){
        Vector base = s.p2 - s.p1;
        double r = dot(p - s.p1,base)/norm(base);
        return sp1+base*r;
    }
    Point reflect(Segment s,Point p){
        return p + (project(s,p) - p)*2.0;
    }
    double getDistance(Point a,Point b){
        return abs(a-b);
    }
    double getDistanceLP(Line l, Point p){
        return abs(cross(l.p2 - l.p1,p-l.p1)/abs(l.p2-l.p1));
    }
    double getDistanceSP(Segment s,Point p){
        if (dot(s.p2 - s.p1,p - s.p1)<0.0) return abs(p-s.p1);
        if (dot(s.p1 - s.p2,p - s.p2)<0.0) return abs(p-sp2);
        return getDistanceLP(s,p);
    }
    
}
typedef Point Vector;

struct Segment {
    Point p1,p2;
};
typedef Segment Line;
class Circle {
    Point c;
    double r;
    Circle(Point c = Point(), double r= 0.0:c(c),r(r)){}
};