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

//���Ǽ����ͼ��ѧ�еģ����������߶εĲ�������ʽ�󽻵�����⣬�鼮�μ�<<�����ͼ��ѧ(OpenGL��)(��3��>> P155 ��ϰ4.6.1 ������㷨
//�Լ����˵�ʱ��ʵ���ˣ���Ϊ��ܻ���.
//�鼮������  Francis S Hill,Jr.       Stephen M Kelley  �廪��ѧ������
//���ĸ����ʾ�����߶Σ�����߶�û�н��㣨ֱ�����ཻ��������0 ������У��򷵻�1.
//����ཻ����ô�����ֵ�����InterPt.�����ֱ�߲��ཻ���򷵻�-1
int segIntersect(Point2 A,Point2  B,Point2 C,Point2 D,Point2 &InterPtr) ;
int SystemOfLinearEquationOfTwo(Vector2 b,Vector2 c,Vector2 d,float& t,float& u);
int SystemOfLinearEquationOfTwo(float bx,float by, float cx,float cy,float dx,float dy,float& t,float& u);
