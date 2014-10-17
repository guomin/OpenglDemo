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

//看来我要编写一个二元一次方程组的求解函数啊
//形式如：6t=5/2+7/2u
//              -7/2t=-3/2+5/2u

//我想在计算机中不必考虑那种相加和相乘的方法,因为这样还需要寻找最小公倍数
//可以这样   求u时，另 两个等式左边的t的系数都化为1，那么等式1右边和等式2右边相等，很容易求出u
//求t时，同理吧

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

//求解二元一次方程组
int SystemOfLinearEquationOfTwo(Vector2 b,Vector2 c,Vector2 d,float& t,float& u)
{
    //判断向量a和d是否平行,先归一化向量
    Vector2 unia ;
    Vector2 unid;
    unia.x=b.x/(sqrt(b.x*b.x+b.y*b.y)); //分母比不为0，因为如果为0就是0向量了，而这在实际情况没有意义
    unia.y=b.y/(sqrt(b.x*b.x+b.y*b.y));
    unid.x=d.x/(sqrt(d.x*d.x+d.y*d.y));
    unid.y=d.y/(sqrt(d.x*d.x+d.y*d.y));
    if (unia.x == unid.x && unia.y==unid.y)
    {
        return -1 ; //没有交点，线段平行,有无数的解
    }

    //计算u
    //a.x==0和a.y==0不会同时出现
    if(b.x==0)// 无需判断if(a.x==0 && d.x!=0)，因为如果a.x,d.x都为0，必然是平行的情况 
    {
        u=-c.x/d.x ;
    }
    else if (b.y==0)// && d.y!=0)
    {
        u=-c.y/d.y;
    }
    else if (b.x !=0 && b.y != 0)
    {
        //方程组1两边同时除以a.x
        //方程组2两边同时除以a.y
        //另方程组相等即:c.x/a.x+d.x/a.x *u = c.y/a.y+d.y*u /a.y;
        u=(c.y/b.y-c.x/b.x)/(d.x/b.x-d.y/b.y) ;
    }

    //计算t
    //d.x==0 && d.y==0是不会同时出现的，因为(0,0)向量在这里没有意义
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
        //计算t
        //方程组1两边同时除以d.x
        //方程组两边同时除以d.y
        // 同上得:a.x/d.x*t-c.x/d.x=a.y/d.y*t-c.y/d.y
        t=(c.x/d.x-c.y/d.y)/(b.x/d.x-b.y/d.y) ;
    }

    //如果线段有交点，那么t和u必然介于0和1之间
    if(0<=t&&t<=1 && 0<=u&&u<=1) 
    {
        return 1;//线段有交点
    }
    else if (t<0 || t>1 || u<0 || u>1)
    {
        return 0;//线段没有交点，但是线段所在的直线有交点
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
        cout <<"直线平行，没有交点/n";
        break;
    case 0:
        cout <<"线段无交点，线段所在的直线有交点,交点为:";
        cout <<"("<<InterPtr.x<<","<<InterPtr.y<<")/n";
        break;
    case 1:
        cout <<"线段有交点，交点为";
        cout <<"("<<InterPtr.x<<","<<InterPtr.y<<")/n";
        break;
    default:
        cout <<"异常发生，请检查/n";
        break;
    }

}

int main()
{
    Point2 A,B,C,D;
    Point2 InterPtr;

    //test1,线段有交点 pass
    A.set(1,4);
    B.set(7,1/2.0);
    C.set(7/2.0,5/2.0);
    D.set(7,5);
    test(A,B,C,D,InterPtr); // t=1/109,u=46/109

    //test2,线段有交点 pass
    cout<<"/n";
    A.set(1.0,4.0);
    B.set(7.0,1.0/2);
    C.set(5.0,0);
    D.set(0,7.0);
    test(A,B,C,D,InterPtr); //	t=16/49,u=20/49

    //test3,线段平行 pass
    cout<<"/n";
    A.set(0,7);
    B.set(7,0);
    C.set(8,-1);
    D.set(10,-3);
    test(A,B,C,D,InterPtr);

    //test4,线段无交点，但所在直线有交点 pass
    cout<<"/n";
    A.set(0,2);
    B.set(5,0);
    C.set(0,5);
    D.set(3,3);
    test(A,B,C,D,InterPtr);// t=9/4, t=15/4  

    //test5,平行的情况，同时a.x=0,,d.x=0的情况  pass

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
线段有交点，交点为(3.53211,2.52294)

A(1,4)  B(7,0.5)
C(5,0)  D(0,7)
t=0.326531 ,u=0.408163
线段有交点，交点为(2.95918,2.85714)

A(0,7)  B(7,0)
C(8,-1) D(10,-3)
直线平行，没有交点

A(0,2)  B(5,0)
C(0,5)  D(3,3)
t=2.25 ,u=3.75
线段无交点，线段所在的直线有交点,交点为:(11.25,-2.5)

A(0,2)  B(0,4)
C(1,2)  D(1,4)
直线平行，没有交点
Press any key to continue . . .


*/