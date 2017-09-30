#include "Vector.hpp"
namespace Vector2d
{
    double norm(Point p){return p.x*p.x +p.y*p.y;}

    double abs(Point p){ return sqrt(norm(p));}
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
        return equals(dot(s1.p2-s1.p1,s2.p2-s2.p1),0.0);;
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
    Point project(Segment s,Point p){
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
    //线段之间最短距离
    double getDistance(Segment s1,Segment s2){
        if (intersect(s1,s2)) return 0.0;
        return min(min(getDistanceSP(s1,s2.p1),getDistanceSP(s1,s2.p2),
                   min(getDistanceSP(s2,s1.p1),getDistanceSP(s2,s1.p2))));
    }
    //相交点
    Point getCrossPoint(Segment s1,Segment s2){
        Vector base = s2.p2 - s2.p1;
        double d1 = abs(cross(base, s1.p1 - s2.p1));
        double d2 = abs(cross(base, s1.p2 - s2.p1));
        double t = d1/(d1 + d2);
        return s1.p1 + (s1.p2 - s1.p1) * t;
    }
    //圆和直线交点
    pair<Point,Point> getCrossPoints(Circle c,Line l){
        assert(intersect(c,1));
        Vector pr = project(1,c.c);
        Vector e = (l.p2 - l.p1) / abs(l.p2 - l.p1);
        double base = sqrt(c.r*c.r - norm(pr - c.c);
        return make_pair(pr + e*base,pr-e*base);
    }
    double arg(Vector p) { return atan2(p.y,p.x);}
    Vector polar(double a, double r) {return Point(cos(r)*a,sin(r)*a);}
    //圆和圆交点
    pair<Point,Point> getCrossPoints(Circle c1,Circle c2){
        assert(intersect(c1,c2));
        double d = abs(c1.c - c2.c);
        double a = acos((c1.r*c1.r + d*d - c2.r*c2.r)/(2*c1.r*d));
        double t = arg(c2.c - c1.c);
        return make_pair(c1.c + polar(c1.r,t+a), c1.c + polar(c1.r,t-a));
    }
    //点的内包
    /*
     IN 2
     ON 1
     OUT 0
     */
    int contains(Polygon g,Point p){
        int n= g.size();
        bool x = false;
        for(int i=0; i < n; i++){
            Point a = g[i] - p,b= g[(i+1)%n] - p;
            if(abs(cross(a,b)) < EPS && dot(a,b) < EPS) return 1;
            if( a.y < EPS && EPS < b.y && cross(a,b)> EPS) x = !x;
        }
        return (x ？ 2 : 0);
    }

}
