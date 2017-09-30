
#define EPS (1e-10)
#define equals(a,b) (fabs((a) - (b)) < EPS)
static const int COUNTER_CLOCKWISE = 1; //反时钟方向
static const int CLOCKWISE = -1; //时钟方向
static const int ONLINE_BACK = 2; //顺直线方向
static const int ONLINE_FRONT = -2; //反直线方向
static const int ON_SEGMENT = 0; //线段中间
namespace Vector2d
{
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
    inline double norm(Point p);

    inline double abs(Point p);
    //点乘
    inline double dot(Vector a,Vector b);
    //X乘
    inline double cross(Vector a,Vector b);
    //是否垂直
    inline bool isOrthogonal(Vector a,Vector b);
    inline bool isOrthogonal(Point a1,Point a2,Point b1,Point b2);
    inline bool isOrthogonal(Segment s1,Segment s2);
    //是否平行
    inline bool isParallel(Vector a,Vector b);
    inline bool isParallel(Point a1,Point a2,Point b1,Point b2);
    inline bool isParallel(Segment s1,Segment s2);
    //投影点
    inline Point project(Segment s,Point P);
    inline Point reflect(Segment s,Point p);
    inline double getDistance(Point a,Point b);
    //点到线段垂直距离
    inline double getDistanceLP(Line l, Point p);
    //点到线段最短距离
    inline double getDistanceSP(Segment s,Point p);
    //线段之间最短距离
    inline double getDistance(Segment s1,Segment s2);
    //判断点跟直线位置关系
    inline int ccw(Point p0,Point p1,Point p2);
    inline bool intersect(Point p1,Point p2,Point p3,Point p4);
    inline bool intersect(Segment s1,Segment s2);
    Point getCrossPoint(Segment s1,Segment s2);


}