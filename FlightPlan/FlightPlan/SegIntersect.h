//calculate the intersection point of two lines
#pragma once

#include <iostream>
#include <math.h>
using namespace std;


class Point2{
public:
    float x, y ;
    void set(float x, float y){this->x=x ;this->y=y;}
    void set(Point2&p){x=p.x;y=p.y;}
};

typedef  Point2 Vector2;

//这是计算机图形学中的，利用两条线段的参数化形式求交点的问题，书籍参见<<计算机图形学(OpenGL版)(第3版>> P155 练习4.6.1 设计求交算法
//自己花了点时间实现了，因为这很基础.
//书籍作者是  Francis S Hill,Jr.       Stephen M Kelley  清华大学出版社
//用四个点表示两条线段，如果线段没有焦点（直线有相交），返回0 ，如果有，则返回1.
//如果相交，那么焦点的值会存入InterPt.如果父直线不相交，则返回-1
int segIntersect(Point2 A,Point2  B,Point2 C,Point2 D,Point2 &InterPtr) ;
int SystemOfLinearEquationOfTwo(Vector2 b,Vector2 c,Vector2 d,float& t,float& u);
int SystemOfLinearEquationOfTwo(float bx,float by, float cx,float cy,float dx,float dy,float& t,float& u);
