#include <iostream>
#include <math.h>
#include "SegIntersect.h"

using namespace std;

int segIntersect(Point2 A,Point2  B,Point2 C,Point2 D,Point2 &InterPtr) 
{
    float t,u ;
    int result ;
    Vector2 b;
    b.set(B.x-A.x,B.y-A.y) ;
    Vector2 c  ;
    c.set(C.x-A.x,C.y-A.y) ;
    Vector2 d;
    d.set(D.x-C.x,D.y-C.y) ;		
    result=SystemOfLinearEquationOfTwo(b,c,d,t,u);
    switch(result)
    {
    case 1:
        //both t or u are ok,here we use t
        InterPtr.x=A.x+t*b.x;
        InterPtr.y=A.y+t*b.y;
#ifdef _DEBUG
        cout<<"/t"<<"t="<<t<<" ,u="<<u<<"/n";
#endif // _DEBUG

        break;
    case -1:
        break;
    case 0:
        //we also calculate the intersect point
        InterPtr.x=A.x+t*b.x;
        InterPtr.y=A.y+t*b.y;
#ifdef _DEBUG
        cout<<"/t"<<"t="<<t<<" ,u="<<u<<"/n";
#endif // _DEBUG
        break;
    default:
        break;
    }

    return result;
}

//������Ҫ��дһ����Ԫһ�η��������⺯����
//��ʽ�磺6t=5/2+7/2u
//              -7/2t=-3/2+5/2u

//�����ڼ�����в��ؿ���������Ӻ���˵ķ���,��Ϊ��������ҪѰ����С������
//��������   ��uʱ���� ������ʽ��ߵ�t��ϵ������Ϊ1����ô��ʽ1�ұߺ͵�ʽ2�ұ���ȣ����������u
//��tʱ��ͬ���

//parameter ax,ay represents VectorA(ax,ay), and so on.
int SystemOfLinearEquationOfTwo(float bx,float by, float cx,float cy,float dx,float dy,float& t,float& u)
{
    Vector2 b  ;
    b.set(bx,by) ;
    Vector2 c  ;
    c.set(cx,cy) ;
    Vector2 d;
    d.set(dx,dy) ;
    return SystemOfLinearEquationOfTwo(b,c,d,t,u);
}

//����Ԫһ�η�����
int SystemOfLinearEquationOfTwo(Vector2 b,Vector2 c,Vector2 d,float& t,float& u)
{
    //�ж�����a��d�Ƿ�ƽ��,�ȹ�һ������
    Vector2 unia ;
    Vector2 unid;
    unia.x=b.x/(sqrt(b.x*b.x+b.y*b.y)); //��ĸ�Ȳ�Ϊ0����Ϊ���Ϊ0����0�����ˣ�������ʵ�����û������
    unia.y=b.y/(sqrt(b.x*b.x+b.y*b.y));
    unid.x=d.x/(sqrt(d.x*d.x+d.y*d.y));
    unid.y=d.y/(sqrt(d.x*d.x+d.y*d.y));
    if (unia.x == unid.x && unia.y==unid.y)
    {
        return -1 ; //û�н��㣬�߶�ƽ��,�������Ľ�
    }

    //����u
    //a.x==0��a.y==0����ͬʱ����
    if(b.x==0)// �����ж�if(a.x==0 && d.x!=0)����Ϊ���a.x,d.x��Ϊ0����Ȼ��ƽ�е���� 
    {
        u=-c.x/d.x ;
    }
    else if (b.y==0)// && d.y!=0)
    {
        u=-c.y/d.y;
    }
    else if (b.x !=0 && b.y != 0)
    {
        //������1����ͬʱ����a.x
        //������2����ͬʱ����a.y
        //��������ȼ�:c.x/a.x+d.x/a.x *u = c.y/a.y+d.y*u /a.y;
        u=(c.y/b.y-c.x/b.x)/(d.x/b.x-d.y/b.y) ;
    }

    //����t
    //d.x==0 && d.y==0�ǲ���ͬʱ���ֵģ���Ϊ(0,0)����������û������
    if (d.x==0)// && a.x!=0)
    {
        t=c.x/b.x;
    }
    else if (d.y ==0 )//&& a.y!=0)
    {
        t=c.y/b.y;
    }
    else if(d.x !=0 && d.y != 0)
    {
        //����t
        //������1����ͬʱ����d.x
        //����������ͬʱ����d.y
        // ͬ�ϵ�:a.x/d.x*t-c.x/d.x=a.y/d.y*t-c.y/d.y
        t=(c.x/d.x-c.y/d.y)/(b.x/d.x-b.y/d.y) ;
    }

    //����߶��н��㣬��ôt��u��Ȼ����0��1֮��
    if(0<=t&&t<=1 && 0<=u&&u<=1) 
    {
        return 1;//�߶��н���
    }
    else if (t<0 || t>1 || u<0 || u>1)
    {
        return 0;//�߶�û�н��㣬�����߶����ڵ�ֱ���н���
    }
    return -9999 ;
}

/*
void test(Point2 A,Point2 B,Point2 C,Point2 D,Point2& InterPtr)
{
    int result ;
    cout<<"A("<<A.x<<","<<A.y<<")/t";
    cout<<"B("<<B.x<<","<<B.y<<")/n";
    cout<<"C("<<C.x<<","<<C.y<<")/t";
    cout<<"D("<<D.x<<","<<D.y<<")/n";
    result=segIntersect(A,B,C,D,InterPtr);
    switch(result)
    {
    case -1:
        cout <<"ֱ��ƽ�У�û�н���/n";
        break;
    case 0:
        cout <<"�߶��޽��㣬�߶����ڵ�ֱ���н���,����Ϊ:";
        cout <<"("<<InterPtr.x<<","<<InterPtr.y<<")/n";
        break;
    case 1:
        cout <<"�߶��н��㣬����Ϊ";
        cout <<"("<<InterPtr.x<<","<<InterPtr.y<<")/n";
        break;
    default:
        cout <<"�쳣����������/n";
        break;
    }

}

int main()
{
    Point2 A,B,C,D;
    Point2 InterPtr;

    //test1,�߶��н��� pass
    A.set(1,4);
    B.set(7,1/2.0);
    C.set(7/2.0,5/2.0);
    D.set(7,5);
    test(A,B,C,D,InterPtr); // t=1/109,u=46/109

    //test2,�߶��н��� pass
    cout<<"/n";
    A.set(1.0,4.0);
    B.set(7.0,1.0/2);
    C.set(5.0,0);
    D.set(0,7.0);
    test(A,B,C,D,InterPtr); //	t=16/49,u=20/49

    //test3,�߶�ƽ�� pass
    cout<<"/n";
    A.set(0,7);
    B.set(7,0);
    C.set(8,-1);
    D.set(10,-3);
    test(A,B,C,D,InterPtr);

    //test4,�߶��޽��㣬������ֱ���н��� pass
    cout<<"/n";
    A.set(0,2);
    B.set(5,0);
    C.set(0,5);
    D.set(3,3);
    test(A,B,C,D,InterPtr);// t=9/4, t=15/4  

    //test5,ƽ�е������ͬʱa.x=0,,d.x=0�����  pass

    cout<<"/n";
    A.set(0,2);
    B.set(0,4);
    C.set(1,2);
    D.set(1,4);
    test(A,B,C,D,InterPtr); 

    return 0 ;
}
*/

/*output:

A(1,4)  B(7,0.5)
C(3.5,2.5)      D(7,5)
t=0.422018 ,u=0.00917431
�߶��н��㣬����Ϊ(3.53211,2.52294)

A(1,4)  B(7,0.5)
C(5,0)  D(0,7)
t=0.326531 ,u=0.408163
�߶��н��㣬����Ϊ(2.95918,2.85714)

A(0,7)  B(7,0)
C(8,-1) D(10,-3)
ֱ��ƽ�У�û�н���

A(0,2)  B(5,0)
C(0,5)  D(3,3)
t=2.25 ,u=3.75
�߶��޽��㣬�߶����ڵ�ֱ���н���,����Ϊ:(11.25,-2.5)

A(0,2)  B(0,4)
C(1,2)  D(1,4)
ֱ��ƽ�У�û�н���
Press any key to continue . . .


*/