
#define EPS (1e-10)
#define equals(a,b) (fabs((a) - (b)) < EPS)
static const int COUNTER_CLOCKWISE = 1; //反时钟方向
static const int CLOCKWISE = -1; //时钟方向
static const int ONLINE_BACK = 2; //顺直线方向
static const int ONLINE_FRONT = -2; //反直线方向
static const int ON_SEGMENT = 0; //线段中间
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
    //点乘
    double dot(Vector a,Vector b){
        return a.x * b.x + a.y * b.y;
    }
    //X乘
    double cross(Vector a,Vector b){
        return a.x*b.y - a.y*b.x;
    }
    //是否垂直
    bool isOrthogonal(Vector a,Vector b){
        return equals(dot(a,b),0);
    }
    bool isOrthogonal(Point a1,Point a2,Point b1,Point b2){
        return isOrthogonal(a1-a2,b1-b2);
    }
    bool isOrthogonal(Segment s1,Segment s2){
        return equals(dot(s1.p2-s1.p1,s2.p2-s2.p1),0.0);
    }
    //是否平行
    bool isParallel(Vector a,Vector b){
        return equals(cross(a,b),0.0);
    }
    bool isParallel(Point a1,Point a2,Point b1,Point b2){
        return isParallel(a1-a2,b1-b2);
    }
    bool isParallel(Segment s1,Segment s2){
        return equals(cross(s1.p2-s1.p1,s2.p2-s2.p1),0.0);
    }
    //投影点
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
    //点到线段垂直距离
    double getDistanceLP(Line l, Point p){
        return abs(cross(l.p2 - l.p1,p-l.p1)/abs(l.p2-l.p1));
    }
    //点到线段最短距离
    double getDistanceSP(Segment s,Point p){
        if (dot(s.p2 - s.p1,p - s.p1)<0.0) return abs(p-s.p1);
        if (dot(s.p1 - s.p2,p - s.p2)<0.0) return abs(p-sp2);
        return getDistanceLP(s,p);
    }
    //线段之间最短距离
    double getDistance(Segment s1,Segment s2){
        if (intersect(s1,s2)) return 0.0;
        return min(min(getDistanceSP(s1,s2.p1),getDistanceSP(s1,s2.p2),
                   min(getDistanceSP(s2,s1.p1),getDistanceSP(s2,s1.p2))));
    }
    //判断点跟直线位置关系
    int ccw(Point p0,Point p1,Point p2){
        Vector a = p1 - p0;
        Vector b = p2 - p0;
        if( cross(a,b) > 0) return COUNTER_CLOCKWISE；
        if( cross(a,b) < 0) return CLOCKWISE;
        if( dot(a, b) < 0) return ONLINE_BACK;
        if( a.norm() < b.norm()) return ONLINE_FRONT;
        return ON_SEGMENT;
    }
    bool intersect(Point p1,Point p2,Point p3,Point p4){
        return ( ccw(p1, p2 ,p3)* ccw(p1,p2,p4) <= 0 && 
                 ccw(p3, p4 ,p1)* ccw(p3,p4,p2) <= 0);
    }
    bool intersect(Segment s1,Segment s2){
        return intersect(s1.p1,s1.p2,s2.p1,s2.p2);
    }
    Point getCrossPoint(Segment s1,Segment s2){
        Vector base = s2.p2 - s2.p1;
        double d1 = abs(cross(base, s1.p1 - s2.p1));
        double d2 = abs(cross(base, s1.p2 - s2.p1));
        double t = d1/(d1 + d2);
        return s1.p1 + (s1.p2 - s1.p1) * t;
    }
};
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
