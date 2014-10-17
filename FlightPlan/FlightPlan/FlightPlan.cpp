// ToDo: 1、不同飞机机型，ab不同； 2、垂直方向绕行； 3、 先调整时间，再调整位置

#include "glut.h" 
#include <stdio.h>
#include <vector>
#include <time.h>
#include <atltime.h>
#include "SegIntersect.h"
#include <float.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ShellAPI.h>

#include "log4cplus/logger.h"
#include "log4cplus/fileappender.h"
#include "log4cplus/layout.h"
#include "log4cplus/ndc.h"
#include "log4cplus/helpers/loglog.h"

#ifdef _DEBUG
#pragma comment(lib, "log4cplusUD.lib")
#else
#pragma comment(lib, "log4cplusU.lib")
#endif

#pragma warning(disable:4996)
#pragma warning(disable:4819)
#pragma warning(disable:4244)

using namespace log4cplus;
using namespace std;

//
#define NUM 100 // 折线的最大折线段

int Flag = 0;	// 标记是否已经开始绘制折线
int RFlag = 0;	// 标记是否已经完成一个矩形
int Function = 1;	// 标记选择的功能是画折线还是矩形
int winWidth = 5000, winHeight = 5000;	// 窗口的宽度和高度
int Mousex, Mousey;	// 用于记录当前鼠标的位置
int n = 0;			// 用于记录折线有几段
int m = 0;			// 用于记录矩形个数
int tx = 0, ty = 0;
int oX = winWidth/2;
int oY = winHeight/2;

//旋转的变量
float angle = 0.0f;
float lx = 0.0f, lz=-1.0f;
float x = 0.0f, z= 5.0f;
float fraction=1.0f;

//////////////////////////////////////////////////////////////////////
const GLfloat PI = 3.14;

/// record the state of mouse
GLboolean mouserdown = GL_FALSE;
GLboolean mouseldown = GL_FALSE;
GLboolean mousemdown = GL_FALSE;

/// when a mouse-key is pressed, record current mouse position 
static GLint mousex = 0, mousey = 0;

static GLfloat center[3] = {0.0f, 0.0f, -1.0f}; /// center position
static GLfloat eye[3] = {0.0f, 0.0f, 0.0f}; /// eye's position

static GLfloat yrotate = PI/4; /// angle between y-axis and look direction
static GLfloat xrotate = PI/4; /// angle between x-axis and look direction
static GLfloat celength = 5.0f;//20.0f;/// lenght between center and eye

static GLfloat mSpeed = 0.4f; /// center move speed
static GLfloat rSpeed = 0.02f; /// rotate speed
static GLfloat lSpeed = 0.4f; /// reserved

/// calculate the eye position according to center position and angle,length
void CalEyePosition()
{
	if(yrotate > PI/2.2) yrotate = PI/2.2;   /// 限制看得方向
	if(yrotate < 0.01)  yrotate = 0.01;
	if(xrotate > 2*PI)   xrotate = 0.01;
	if(xrotate < 0)   xrotate = 2 * PI;
	if(celength > 50)  celength = 50;     ///  缩放距离限制
	if(celength < 5)   celength = 5;
	/// 下面利用球坐标系计算 eye 的位置，
	/*eye[0] = center[0] + celength * sin(yrotate) * cos(xrotate);  
	eye[2] = center[2] + celength * sin(yrotate) * sin(xrotate);
	eye[1] = center[1] + celength * cos(yrotate);*/
	center[0] = celength * sin(yrotate) * cos(xrotate);  
	center[2] = celength * sin(yrotate) * sin(xrotate);
	center[1] = celength * cos(yrotate);
}

/// center moves
void MoveBackward()              /// center 点沿视线方向水平向后移动
{
	center[0] += mSpeed * cos(xrotate);
	center[2] += mSpeed * sin(xrotate);
	CalEyePosition();
}

void MoveForward()
{
	center[0] -= mSpeed * cos(xrotate);
	center[2] -= mSpeed * sin(xrotate);
	CalEyePosition();
}

/// visual angle rotates
void RotateLeft()
{
	xrotate -= rSpeed;
	CalEyePosition();
}

void RotateRight()
{
	xrotate += rSpeed;
	CalEyePosition();
}

void RotateUp()
{
	yrotate += rSpeed;
	CalEyePosition();
}

void RotateDown()
{
	yrotate -= rSpeed;
	CalEyePosition();
}

/// CALLBACK func for keyboard presses
void KeyFunc(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'a': RotateLeft(); break;
	case 'd': RotateRight();break;
	case 'w': MoveForward(); break;
	case 's': MoveBackward(); break;
	case 'q': RotateUp(); break;
	case 'e': RotateDown(); break;
	}
	glutPostRedisplay();
}

/// CALLBACK func for mouse kicks
void MouseFunc(int button, int state, int x, int y)
{
	if(state == GLUT_DOWN)
	{
		if(button == GLUT_RIGHT_BUTTON) mouserdown = GL_TRUE;
		if(button == GLUT_LEFT_BUTTON) mouseldown = GL_TRUE;
		if(button == GLUT_MIDDLE_BUTTON)mousemdown = GL_TRUE;
	}
	else
	{
		if(button == GLUT_RIGHT_BUTTON) mouserdown = GL_FALSE;
		if(button == GLUT_LEFT_BUTTON) mouseldown = GL_FALSE;
		if(button == GLUT_MIDDLE_BUTTON)mousemdown = GL_FALSE;
	}
	mousex = x, mousey = y;
}

/// CALLBACK func for mouse motions
void MouseMotion(int x, int y)
{
	if(mouserdown == GL_TRUE)
	{       /// 所除以的数字是调整旋转速度的，随便设置，达到自己想要速度即可
		xrotate += (x - mousex) / 80.0f;     
		yrotate -= (y - mousey) / 120.0f;
	}

	if(mouseldown == GL_TRUE)
	{
		celength += (y - mousey) / 25.0f;
	}
	mousex = x, mousey = y;
	CalEyePosition();
	glutPostRedisplay();
}
bool is3d = false;
void LookAt()            /// 调用 gluLookAt(), 主要嫌直接调用要每次都写好几个参数。。
{
	if (is3d)
	{
		CalEyePosition();
		gluLookAt(eye[0], eye[1], eye[2],
			center[0], center[1], center[2],
			0, 1, 0);
	}
	else
	{
		gluLookAt(0.0, 0.0, 0.0,
		0.0, 0.0, -1.0f,
		0, 1, 0);
	}
	
}

//////////////////////////////////////////////////////////////////////

char *conFileName = "conFile.txt";

// 线性结构体
struct LineNode {
	int x1;
	int y1;
	int x2;
	int y2;
}Line[NUM];
// 矩形结构体
struct Rectangle {
	int x1;
	int y1;
	int x2;
	int y2;
}Rect[NUM];


#define  POS_FRONT 1	// 前
#define  POS_SIDE  0	// 两侧
#define  POS_BACK -1	// 后

bool isAdjustTimeFirst = false;// 是否先调整速度时间

// 空间离散化步长
const double DX = 10.0;
const double DY = 10.0;
const double DZ = 10.0;

// 时空半径
const double TIME_RADIUS = 60.0;
const double SPACE_RADIUS = 1000.0;

// 冲突调整控制参数
const double DELTA = 0.5;

// 日志
Logger static _logger = Logger::getInstance(LOG4CPLUS_TEXT("test.subtestof_filelog"));



// 时空点坐标，xyz可能需要根据经度、纬度和高度进行换算
// 空间点为离散化数值，x = xi * DX
struct Point
{
    int xi;
    int yi;
    int zi;
    time_t t;
	int ptype;	// 临时定义在这里，表示机型

    Point()
    {
        xi = 0;
        yi = 0;
        zi = 0;
        t = 0;
		ptype = 0;
    }

    Point(int xi, int yi, int zi, int ptype=0)
    {
        this->xi = xi;
        this->yi = yi;
        this->zi = zi;
        t = 0;
		this->ptype = ptype;
    }
};

// 如果时间为0，则表示不区分时间

// 白名单
typedef vector<Point> WhiteList;

// 黑名单
typedef vector<Point> BlackList;

WhiteList g_wl;
BlackList g_bl;
BlackList g_bl2;

// 动态添加航线需要的结构
WhiteList g_wlByConf;
BlackList g_blByConf;


// 路径
typedef vector<Point> Route;

// 激活路径，即正在飞行的路径，若current_point = -1表示该路径变为非激活路径
class ActiveRoute
{
public:
    int current_point;
    double a;   // 椭球范围的长轴
    double b;   // 椭球范围的短轴
    Route route;
	Point dep;
	Point arr;
	float r;
	float g;
	float blue;
};

// 激活路径组
typedef vector<ActiveRoute> ActiveRouteGroup;

struct ControlPara
{
    double vx;	// 速度
    double vy;
    double vz;
    double ax;	// 加速度
    double ay;
    double az;
};

// 错误代码
enum ErrorCode
{
    SUCCESS = 0,
    DDA_ERROR,
    NN_OUT_RADIUS,
    NO_COLLISION_FOUND,
    NO_SOLUTION,
    OUT_MAX_ITER,
    SOME_ERROR
};

// OpenGL
static time_t g_OpenGL_time = 0;
const time_t g_OpenGL_time_step = 60;
const float g_OpenGL_draw_delta = 5e-5;


ActiveRouteGroup route_group;

inline int DiscreteX(double x)
{
    return (int)(x / DX + 0.5);
}

inline int DiscreteY(double y)
{
    return (int)(y / DY + 0.5);
}

inline int DiscreteZ(double z)
{
    return (int)(z / DZ + 0.5);
}

int PrintRoute(FILE *fp, const Route &route, int level)
{
    if (level == 0)
    {
        for (Route::const_iterator it = route.begin(); it != route.end(); ++it)
        {
            fprintf(fp, "%d, %d, %d, %d\n", it->xi, it->yi, it->zi, it->t);
        }
    }
    else if (level == 1)
    {
        for (Route::const_iterator it = route.begin(); it != route.end(); ++it)
        {
            fprintf(fp, "%f, %f, %f, %s", it->xi * DX, it->yi * DY, it->zi * DZ, ctime(&(it->t)));
        }
    }

    return SUCCESS;

}

// 数字化/离散化路径各点
int DDA3D(int x0, int y0, int z0, int x1, int y1, int z1, Route &route)
{
    route.clear();
    int dx = x1 - x0;
    int dy = y1 - y0;
    int dz = z1 - z0;
    double x = x0;
    double y = y0;
    double z = z0;

    if (abs(dx) >= abs(dy) && abs(dx) >= abs(dz))
    {
        if (dx == 0)
        {
            return DDA_ERROR;
        }

        double my = (double)dy / dx;
        double mz = (double)dz / dx;

        if (x0 < x1)
        {
            for (int xi = x0; xi <= x1; xi++)
            {
                int yi = (int)y;
                int zi = (int)z;

                route.push_back(Point(xi, yi, zi));

                y += my;
                z += mz;
            }
        }
        else
        {
            for (int xi = x0; xi >= x1; xi--)
            {
                int yi = (int)y;
                int zi = (int)z;

                route.push_back(Point(xi, yi, zi));

                y -= my;
                z -= mz;
            }
        }
    }
    else if (abs(dy) >= abs(dz) && abs(dy) >= abs(dx))
    {
        if (dy == 0)
        {
            return DDA_ERROR;
        }

        double mz = (double)dz / dy;
        double mx = (double)dx / dy;

        if (y0 < y1)
        {
            for (int yi = y0; yi <= y1; yi++)
            {
                int zi = (int)z;
                int xi = (int)x;

                route.push_back(Point(xi, yi, zi));

                z += mz;
                x += mx;
            }
        }
        else
        {
            for (int yi = y0; yi >= y1; yi--)
            {
                int zi = (int)z;
                int xi = (int)x;

                route.push_back(Point(xi, yi, zi));

                z -= mz;
                x -= mx;
            }
        }
    }
    else // (abs(dz) >= abs(dx) && abs(dz) >= abs(dy))
    {
        if (dz == 0)
        {
            return DDA_ERROR;
        }

        double mx = (double)dx / dz;
        double my = (double)dy / dz;

        if (z0 < z1)
        {
            for (int zi = z0; zi <= z1; zi++)
            {
                int xi = (int)x;
                int yi = (int)y;

                route.push_back(Point(xi, yi, zi));

                x += mx;
                y += my;
            }
        }
        else
        {
            for (int zi = z0; zi >= z1; zi--)
            {
                int xi = (int)x;
                int yi = (int)y;

                route.push_back(Point(xi, yi, zi));

                x -= mx;
                y -= my;
            }
        }
    }

    return SUCCESS;
}

int LineRoute(const Point &dep, const Point &arr, Route &route)
{
    int ret = DDA3D(dep.xi, dep.yi, dep.zi, arr.xi, arr.yi, arr.zi, route);
    if (ret != SUCCESS)
    {
        return ret;
    }

    int n_point = route.size();
    double mt = (double)(arr.t - dep.t) / (n_point - 1);
    double t = (double)dep.t;

    for (int i = 0; i < n_point; i++)
    {
        route[i].t = (time_t)t;
        t += mt;
    }

    return SUCCESS;
}


inline double pow2(double x)
{
    return x * x;
}

int TimeDistance(const Point &p1, const Point &p2)
{
    return (int)abs((int)(p1.t - p2.t));
}

//  空间距离：3维
double SpaceDistance(const Point &p1, const Point &p2)
{
    double d2 = pow2((p1.xi - p2.xi) * DX) + pow2((p1.yi - p2.yi) * DY) + pow2((p1.zi - p2.zi) * DZ);
    return sqrt(d2);
}

//  空间距离：2维
double SpaceDistance2D(const Point &p1, const Point &p2)
{
    double d2 = pow2((p1.xi - p2.xi) * DX) + pow2((p1.yi - p2.yi) * DY);
    return sqrt(d2);
}

// 时空距离
double Distance(const Point &p1, const Point &p2)
{
    return SpaceDistance(p1, p2) / SPACE_RADIUS + TimeDistance(p1, p2) / TIME_RADIUS;
}

// 各向同性判断是否冲突
bool IsCollision(const Point &p1, const Point &p2)
{
    return Distance(p1, p2) < 2;
}

// p2是路径中某一点，p1是其前一点，根据椭圆范围判断p0是否冲突
bool IsCollisionEllipse(const Point &p1, const Point &p2, const Point &p0, double a, double b)
{
    double num = -((p0.xi - p1.xi)*(p1.xi - p2.xi)*DX*DX + (p0.yi - p1.yi)*(p1.yi - p2.yi)*DY*DY + (p0.zi - p1.zi)*(p1.zi - p2.zi)*DZ*DZ);
    double den = pow2(p0.xi - p1.xi)*DX*DX +  pow2(p0.yi - p1.yi)*DY*DY + pow2(p0.zi - p1.zi)*DZ*DZ;

    if (den < 1e-6) // p1和p2重合，退化为圆的情况，为安全起见，取其半径为长轴a
    {
        if (SpaceDistance(p2, p0) >= a)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        double u = num / den;
        Point pf;
        pf.xi = DiscreteX((p1.xi + u * (p2.xi - p1.xi)) * DX);
        pf.yi = DiscreteX((p1.yi + u * (p2.xi - p1.xi)) * DY);
        pf.zi = DiscreteX((p1.zi + u * (p2.xi - p1.xi)) * DZ);
        double disPfP2 = SpaceDistance(pf, p2);
        double disP0Pf = SpaceDistance(p0, pf);
        double c = sqrt(pow2(a) - pow2(b));

        if (u > 0)  // 右侧
        {
            if (disPfP2 >= a - c)
            {
                return false;
            }
            else
            {
                double x = disPfP2 + c;
                double y = sqrt(1 - pow2(x / a)) * b;

                if (disP0Pf >= y)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }
        else  // 左侧
        {
            if (disPfP2 >= a + c)
            {
                return false;
            }
            else
            {
                double x = disPfP2 - c;
                double y = sqrt(1 - pow2(x / a)) * b;

                if (disP0Pf >= y)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }
    }
}

// p1和p2分别是两条不同路径中某一点，pre_p1和pre_p2分别是它们的前一点，根据椭圆范围恶法判断p1和p2是否冲突
// 根据p1和pre_p1判p2，再根据p2和pre_p2判p1
bool IsCollisionAni(const Point &p1, const Point &pre_p1, double a1, double b1, const Point &p2, const Point &pre_p2, double a2, double b2)
{
    if (TimeDistance(p1, p2) >= TIME_RADIUS)
    {
        return false;
    }

    if (IsCollisionEllipse(pre_p1, p1, p2, a1, b1) || IsCollisionEllipse(pre_p2, p2, p1, a2, b2))
    {
        return true;
    }
    else
    {
        return false;
    }
}

int AlignTime(const Route &route1, const Route &route2, vector<int> &route2_align_route1)
{

    Route::const_iterator min_it = route2.begin();
    double min_dis = DBL_MAX;
    for (Route::const_iterator it1 = route1.begin(); it1 != route1.end(); ++it1)
    {
        int min_time_dis = INT_MAX;

        for (Route::const_iterator it2 = min_it; it2 != route2.end(); ++it2)
        {
            int time_dis = TimeDistance(*it1, *it2);
            if (time_dis < min_time_dis)
            {
                min_it = it2;
            }
            else
            {
                break;
            }
        }

        double dis = Distance(*it1, *min_it);
        if (dis < min_dis)
        {
            min_dis = dis;
        }
    }

    return SUCCESS;
}

//int NearestNeighbor(const Route &route1, const Route &route2, Route::const_iterator &nn_it1, Route::const_iterator &nn_it2, double &nn_dis)
//{
//	int n1 = route1.size();
//
//	if (route1.front().t - route2.back().t > TIME_RADIUS || route2.front().t - route1.back().t > TIME_RADIUS)
//	{
//		return NN_OUT_RADIUS;
//	}
//
//	Route::const_iterator min_it = route2.begin();
//	double min_dis = DBL_MAX;
//	for (Route::const_iterator it1 = route1.begin(); it1 != route1.end(); ++it1)
//	{
//		int min_time_dis = INT_MAX;
//
//		for (Route::const_iterator it2 = min_it; it2 != route2.end(); ++it2)
//		{
//			int time_dis = TimeDistance(*it1, *it2);
//			if (time_dis <= min_time_dis)
//			{
//				min_time_dis = time_dis;
//				min_it = it2;
//			}
//			else
//			{
//				break;
//			}
//		}
//
//		double dis = Distance(*it1, *min_it);
//		if (dis < min_dis)
//		{
//			min_dis = dis;
//			nn_it1 = it1;
//			nn_it2 = min_it;
//		}
//	}
//
//	nn_dis = min_dis;
//
//	return SUCCESS;
//}

// 找一条路径和一点的最近邻
int NearestNeighbor(const Point &point, const Route &route, Route::const_iterator &nn_it, double &nn_dis)
{
    if (point.t - route.back().t > TIME_RADIUS || route.front().t - point.t > TIME_RADIUS)
    {
        return NN_OUT_RADIUS;
    }

    double min_dis = DBL_MAX;
    for (Route::const_iterator it = route.begin(); it != route.end(); ++it)
    {
        double dis = Distance(point, *it);

        if (dis < min_dis)
        {
            min_dis = dis;
            nn_it = it;
        }
    }

    nn_dis = min_dis;

    return SUCCESS;
}

// 找第一个冲突点
int FindFirstCollision(const Route &route1, double a1, double b1, const Route &route2, double a2, double b2, Route::const_iterator &nn_it1, Route::const_iterator &nn_it2, double &nn_dis)
{
    if (route1.empty() || route2.empty())
    {
        return NN_OUT_RADIUS;
    }
    if (route1.front().t - route2.back().t > TIME_RADIUS || route2.front().t - route1.back().t > TIME_RADIUS)
    {
        return NN_OUT_RADIUS;
    }

    Route::const_iterator min_it = route2.begin();
    for (Route::const_iterator it1 = route1.begin(); it1 != route1.end(); ++it1)
    {
        int min_time_dis = INT_MAX;

        for (Route::const_iterator it2 = min_it; it2 != route2.end(); ++it2)
        {
            int time_dis = TimeDistance(*it1, *it2);
            if (time_dis <= min_time_dis)
            {
                min_time_dis = time_dis;
                min_it = it2;
            }
            else
            {
                break;
            }
        }

        if (it1->t == 1370096700)
        {
            int a = 0;
        }

        Route::const_iterator pre_it1 = it1;
        Route::const_iterator pre_min_it = min_it;
        if (pre_it1 != route1.begin())
        {
            pre_it1--;
        }
        if (pre_min_it != route2.begin())
        {
            pre_min_it--;
        }

		// 此处加上根据机型获取a1、b1、a2、b2的代码

        //if (IsCollision(*it1, *min_it))   // 各向同性判决冲突
        if (IsCollisionAni(*pre_it1, *it1, a1, b1, *pre_min_it, *min_it, a2, b2))   // 各向异性判决冲突
        {

            nn_dis = Distance(*it1, *min_it);
            nn_it1 = it1;
            nn_it2 = min_it;
            return SUCCESS;
        }
    }

    return NO_COLLISION_FOUND;
}

// 规划路径，根据出发地、目的地以及所有占用的路径进行路径规划
int FlightPlanRoute(const Point &dep, const Point &arr, const ActiveRouteGroup &route_group, Route &route, double a, double b)
{
    // dep 到 arr 连点成线生成 route
    LineRoute(dep, arr, route);
    bool no_collision = true;

    // 进行50次循环，如果仍没有规划成功，则把第一个冲突点提前10个点重新递归规划
    for (int i = 0; i < 50; i++)
    {
        no_collision = true;
        for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
        {
            if (it->current_point >= 0)	// 激活路径
            {
                // 逐条路径判断是否与刚生成的 route 存在冲突，如果存在，则修正 route
                Route::const_iterator fc_it1;
                Route::const_iterator fc_it2;
                double nn_dis;

                if (FindFirstCollision(it->route, it->a, it->b, route, a, b, fc_it1, fc_it2, nn_dis) == SUCCESS)
                {
                    double delta_x = 0;
                    double delta_y = 0;
                    double delta_z = 0;
                    double delta_t = 0;
                    for (ActiveRouteGroup::const_iterator it1 = route_group.begin(); it1 != route_group.end(); ++it1)
                    {
                        Route::const_iterator nn_it;
                        double nn_dis;
                        NearestNeighbor(*fc_it2, it1->route, nn_it, nn_dis);

						if (isAdjustTimeFirst)	// 如果想先调整速度时间，此处为true；默认为false
						{
							if (i<10)	// 前10次，只从时间/速度方面调整
							{
								double scale = 5 / (pow2(nn_dis)  + DELTA);
								delta_x = 0;
								delta_y = 0;
								delta_z = 0;
								delta_t += (fc_it2->t - nn_it->t) * scale;
							} 
							else
							{
								double scale = 5 / (pow2(nn_dis)  + DELTA);
								delta_x += (fc_it2->xi - nn_it->xi) * scale;
								delta_y += (fc_it2->yi - nn_it->yi) * scale;
								delta_z += (fc_it2->zi - nn_it->zi) * scale;
								delta_t += (fc_it2->t - nn_it->t) * scale;
							}

						}
						else
						{
							double scale = 5 / (pow2(nn_dis)  + DELTA);
							delta_x += (fc_it2->xi - nn_it->xi) * scale;
							delta_y += (fc_it2->yi - nn_it->yi) * scale;
							delta_z += (fc_it2->zi - nn_it->zi) * scale;
							delta_t += (fc_it2->t - nn_it->t) * scale;
						}
                    }

                    Point mid;
                    mid.xi = fc_it2->xi + (int)(delta_x + 0.5);
                    mid.yi = fc_it2->yi + (int)(delta_y + 0.5);
                    mid.zi = fc_it2->zi + (int)(delta_z + 0.5);
                    mid.t = fc_it2->t + (int)(delta_t + 0.5);

                    Route route1;
                    Route route2;
                    int ret = FlightPlanRoute(dep, mid, route_group, route1, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }
                    ret = FlightPlanRoute(mid, arr, route_group, route2, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }

                    route.resize(route1.size() + route2.size() - 1);
                    copy(route1.begin(), route1.end(), route.begin());
                    copy(route2.begin() + 1, route2.end(), route.begin() + route1.size());

                    no_collision = false;
                }
                else
                {
                    no_collision = true;
                }
            }
        }

        if (no_collision)
        {
            break;
        }
    }

    if (!no_collision)
    {
        for (int i = 0; i < 50; i++)
        {
            bool no_collision = true;
            for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
            {
                if (it->current_point >= 0)	// 激活路径
                {
                    // 逐条路径判断是否与刚生成的 route 存在冲突，如果存在，则修正 route
                    Route::const_iterator fc_it1;
                    Route::const_iterator fc_it2;
                    double nn_dis;

                     if (FindFirstCollision(it->route, it->a, it->b, route, a, b, fc_it1, fc_it2, nn_dis) == SUCCESS)
                    {
                        Route::const_iterator fc_it2_less_10 = fc_it2;
                        
                        while (fc_it2_less_10 != route.begin() + 10 && fc_it2_less_10 != fc_it2 - 10)    // 向前回溯10个点，同时要保证不能太接近起点
                        {
                            fc_it2_less_10--;
                        }
                        fc_it2 = fc_it2_less_10;

                        double delta_x = 0;
                        double delta_y = 0;
                        double delta_z = 0;
                        double delta_t = 0;
                        for (ActiveRouteGroup::const_iterator it1 = route_group.begin(); it1 != route_group.end(); ++it1)
                        {
                            Route::const_iterator nn_it;
                            double nn_dis;
                            NearestNeighbor(*fc_it2, it1->route, nn_it, nn_dis);

                            double scale = 5 / (pow2(nn_dis)  + DELTA);
                            delta_x += (fc_it2->xi - nn_it->xi) * scale;
                            delta_y += (fc_it2->yi - nn_it->yi) * scale;
                            delta_z += (fc_it2->zi - nn_it->zi) * scale;
                            delta_t += (fc_it2->t - nn_it->t) * scale;
                        }

                        Point mid;
                        mid.xi = fc_it2->xi + (int)(delta_x + 0.5);
                        mid.yi = fc_it2->yi + (int)(delta_y + 0.5);
                        mid.zi = fc_it2->zi + (int)(delta_z + 0.5);
                        mid.t = fc_it2->t + (int)(delta_t + 0.5);

                        Route route1;
                        Route route2;
                        int ret = FlightPlanRoute(dep, mid, route_group, route1, a, b);
                        if (ret != SUCCESS)
                        {
                            return ret;
                        }
                        ret = FlightPlanRoute(mid, arr, route_group, route2, a, b);
                        if (ret != SUCCESS)
                        {
                            return ret;
                        }

                        route.resize(route1.size() + route2.size() - 1);
                        copy(route1.begin(), route1.end(), route.begin());
                        copy(route2.begin() + 1, route2.end(), route.begin() + route1.size());

                        no_collision = false;
                    }
                    else
                    {
                        no_collision = true;
                    }
                }
            }
        }
    }

    if (!no_collision)
    {
        return OUT_MAX_ITER;
    }

    return SUCCESS;
}


int WLPlanRoute(const Point &dep, const Point &arr, const ActiveRouteGroup &route_group, const WhiteList &wl, Route &route, double a, double b)
{
    WhiteList wl_c2;  // 第二类白名单
	wl_c2.reserve(wl.size() + 2);	// 加上起点和终点
	wl_c2.push_back(dep);
	wl_c2.insert(wl_c2.begin() + 1, wl.begin(), wl.end());
	wl_c2.push_back(arr);

	Route route_part;

    vector<double> dis_cum;

	route.push_back(dep);

    WhiteList::iterator pre_it = wl_c2.begin();
    for (WhiteList::iterator it = wl_c2.begin() + 1; it != wl_c2.end(); ++it)
    {
        if (it->t != 0)
        {
            int n = dis_cum.size();
            if (n == 0) // 之前没有第一类白名单
            {
                int ret = FlightPlanRoute(*pre_it, *it, route_group, route_part, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }
				route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
            }
            else
            {
                double dis_sum = 0;
                dis_sum = dis_cum[n - 1] + SpaceDistance(*(it - 1), *it);

                for (int i = 0; i < n; i++)
                {
                    (pre_it + i + 1)->t  = pre_it->t + (it->t - pre_it->t) * dis_cum[i] / dis_sum;

                    int ret = FlightPlanRoute(*(pre_it + i), *(pre_it + i + 1), route_group, route_part, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }
					route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
                }

				int ret = FlightPlanRoute(*(it - 1), *(it), route_group, route_part, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }
				route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
            }

			pre_it = it;

            dis_cum.clear();
        }
        else
        {
			double dis = SpaceDistance(*(it - 1), *it);
			int n = dis_cum.size();
			if (n == 0)
			{
				dis_cum.push_back(dis);
			}
			else
			{
				dis_cum.push_back(dis + dis_cum.back());
			}
        }
    }

    return SUCCESS;
}

struct CrossPtr
{
    Point pt;  // 路径穿过禁忌区，路径上的点
    BlackList::const_iterator bl_it;    // 路径穿过禁忌区，禁忌区边的起点
};

int GetCrossPoint(const BlackList &bl, const Point &dep, const Point &arr, vector<CrossPtr> &vcp)
{
    vcp.clear();
    Point2 A;
    A.set(dep.xi, dep.yi);

    Point2 B;
    B.set(arr.xi, arr.yi);

    for (BlackList::const_iterator bl_it = bl.begin(); bl_it != bl.end() - 1; ++bl_it)
    {
        Point2 C;
        C.set(bl_it->xi, bl_it->yi);

        Point2 D;
        D.set((bl_it + 1)->xi, (bl_it + 1)->yi);

        Point2 X;

        if (segIntersect(A, B, C, D, X) == 1)
        {
            if (abs(X.x - D.x) > 0.5 || abs(X.y - D.y) > 0.5)   // 由于BL各段首尾相接，为了避免重复计算，将尾处的相交点去除
            {
                CrossPtr cp;
                cp.pt.xi = X.x;
                cp.pt.yi = X.y;

                cp.bl_it = bl_it;
                vcp.push_back(cp);
            }
        }
    }

    return SUCCESS;
}

/**
 * 针对禁飞区（黑名单）
 */
int BLPlanRoute(const Point &dep, const Point &arr, const ActiveRouteGroup &route_group, const BlackList &bl, Route &route, double a, double b)
{
    bool bl_plan = false;
    if (bl.front().t == 0)
    {
        bl_plan = true;
    }
    else
    {
        if (bl.front().t > arr.t || bl.back().t < dep.t)
        {
            bl_plan = false;
        }
        else
        {
            bl_plan = true;
        }
    }

    // 黑名单规划
    if (bl_plan)
    {
        // 先按正常规划航线
        vector<CrossPtr> vcp; // 路径穿过禁忌区的点
        GetCrossPoint(bl, dep, arr, vcp);

        bool noCross = false;
        bool noSolution = false;
        if (vcp.size() == 0)
        {
            noSolution = false;
            noCross = true;
        }
        else if (vcp.size() == 1)
        {
            if (vcp[0].pt.xi == vcp[0].bl_it->xi && vcp[0].pt.yi == vcp[0].bl_it->yi) // 掠过顶点的情况
            {
                noSolution = false;
                noCross = true;
            }
            else
            {
                noSolution = true;
            }
        }
        else if (vcp.size() % 2 == 1)   // 奇数次穿过代表无解
        {
            noSolution = true;
        }
        else
        {
            noSolution = false;
            noCross = false;
        }

        if (noSolution)
        {
            return NO_SOLUTION;
        }
        else if (noSolution == false && noCross == false)       // 有解但有穿过
        {
            // 路径遍历经过起点-若干顶点（0到n个）-终点的距离，从中找最短的一个
            vector<double> dep_vtx_dis(bl.size());  // 起点到每个顶点的距离
            vector<double> vtx_arr_dis(bl.size());  // 每个顶点到终点的距离
            vector<double> vtx_vtx_dis(bl.size());  // 相邻的顶点间的距离
            vector<int>    dep_vtx_ncross(bl.size());  // 起点到每个顶点的路径与禁忌区的交点数
            vector<int>    vtx_arr_ncross(bl.size());  // 每个顶点到终点的路径与禁忌区的交点数

            int i = 0;
            for (BlackList::const_iterator bl_it = bl.begin(); bl_it != bl.end() - 1; ++bl_it)
            {
                dep_vtx_dis[i] = SpaceDistance2D(dep, *bl_it);
                vtx_arr_dis[i] = SpaceDistance2D(*bl_it, arr);
                vtx_vtx_dis[i] = SpaceDistance2D(*bl_it, *(bl_it + 1));

                vector<CrossPtr> vcp_tmp1;
                GetCrossPoint(bl, dep, bl[i], vcp_tmp1);
                int nCross1 = 0;
                for (vector<CrossPtr>::const_iterator vit = vcp_tmp1.begin(); vit != vcp_tmp1.end(); ++vit)
                {
                    if (abs(vit->pt.xi - bl[i].xi) > 2 || abs(vit->pt.yi - bl[i].yi) > 2)
                    {
                        nCross1++;
                    }
                }
                dep_vtx_ncross[i] = nCross1;


                vector<CrossPtr> vcp_tmp2;
                GetCrossPoint(bl, bl[i], arr, vcp_tmp2);
                int nCross2 = 0;
                for (vector<CrossPtr>::const_iterator vit = vcp_tmp2.begin(); vit != vcp_tmp2.end(); ++vit)
                {
                    if (abs(vit->pt.xi - bl[i].xi) > 2 || abs(vit->pt.yi - bl[i].yi) > 2)
                    {
                        nCross2++;
                    }
                }
                vtx_arr_ncross[i] = nCross2;

                i++;
            }

            double min_dis = DBL_MAX;
            int vtx_in = -1;
            int vtx_out = -1;
            int vtx_dir = -1;   // 0 表示i直接到j，1表示跨过最大值

            int nVtx = bl.size() - 1;
            for (int i = 0; i < nVtx; i++)
            {
                if (dep_vtx_ncross[i] > 0)
                {
                    continue;
                }

                for (int j = 0; j < nVtx; j++)
                {
                    if (vtx_arr_ncross[j] > 0)
                    {
                        continue;
                    }

                    // 起点到入点i + 出点j到终点
                    double dis = dep_vtx_dis[i] + vtx_arr_dis[j];
                    int direction = -1;

                    // 中间i到j或者j到i选较小者
                    if (i > j)
                    {
                        double dis_mid1 = 0;
                        for (int k = j; k < i; k++)
                        {
                            dis_mid1 += vtx_vtx_dis[k];
                        }

                        double dis_mid2 = 0;
                        for (int k = i; k < nVtx; k++)
                        {
                            dis_mid2 += vtx_vtx_dis[k];
                        }
                        for (int k = 0; k < j; k++)
                        {
                            dis_mid2 += vtx_vtx_dis[k];
                        }

                        if (dis_mid1 <= dis_mid2)
                        {
                            direction = 0;
                            dis += dis_mid1;
                        }
                        else
                        {
                            direction = 1;
                            dis += dis_mid2;
                        }
                    }
                    else if (j > i)
                    {
                        double dis_mid1 = 0;
                        for (int k = i; k < j; k++)
                        {
                            dis_mid1 += vtx_vtx_dis[k];
                        }

                        double dis_mid2 = 0;
                        for (int k = j; k < nVtx; k++)
                        {
                            dis_mid2 += vtx_vtx_dis[k];
                        }
                        for (int k = 0; k < i; k++)
                        {
                            dis_mid2 += vtx_vtx_dis[k];
                        }

                        if (dis_mid1 <= dis_mid2)
                        {
                            direction = 0;
                            dis += dis_mid1;
                        }
                        else
                        {
                            direction = 1;
                            dis += dis_mid2;
                        }
                    }

                    if (dis < min_dis)
                    {
                        min_dis = dis;
                        vtx_in = i;
                        vtx_out = j;
                        vtx_dir = direction;
                    }
                } 
            }

            ////////////////////////////////////////////////
            // 择优以后生产路径
            if (min_dis == DBL_MAX)
            {
                return -1;
            }
            else
            {
                Route route_tmp1;
                Route route_tmp2;
                time_t cum_t = 0;
                int g_bl_zi = 0;

                g_bl_zi = (dep.zi + arr.zi) / 2;
                Point pt_vtx_in(bl[vtx_in]);
                pt_vtx_in.zi = g_bl_zi;
                pt_vtx_in.t = (arr.t - dep.t) / min_dis * dep_vtx_dis[vtx_in] + dep.t;
                cum_t = pt_vtx_in.t;
                Point pt_vtx_out(bl[vtx_out]);
                pt_vtx_out.zi = g_bl_zi;
                pt_vtx_out.t = arr.t - (arr.t - dep.t) / min_dis * vtx_arr_dis[vtx_out];
                int ret = FlightPlanRoute(dep, pt_vtx_in, route_group, route_tmp1, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }
                ret = FlightPlanRoute(pt_vtx_out, arr, route_group, route_tmp2, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }

                if (vtx_in == vtx_out)
                {
                    route.resize(route_tmp1.size() + route_tmp2.size() - 1);
                    copy(route_tmp1.begin(), route_tmp1.end(), route.begin());
                    copy(route_tmp2.begin() + 1, route_tmp2.end(), route.begin() + route_tmp1.size());
                }
                else if (vtx_in > vtx_out)
                {
                    vector<Route> vroute;
                    if (vtx_dir == 0)
                    {
                        for (int k = vtx_in; k > vtx_out; --k)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k-1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k-1] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            int ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }
                    }
                    else
                    {
                        for (int k = vtx_in; k < nVtx; k++)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k+1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            int ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }

                        for (int k = 0; k < vtx_out; k++)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k+1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            int ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }
                    }
                    int n_pt = route_tmp1.size() + route_tmp2.size() - 1;
                    for (vector<Route>::const_iterator vrit = vroute.begin(); vrit != vroute.end(); ++vrit)
                    {
                        n_pt += vrit->size() - 1;
                    }

                    route.resize(n_pt);
                    copy(route_tmp1.begin(), route_tmp1.end(), route.begin());
                    int n_pt_cur = route_tmp1.size();
                    for (vector<Route>::const_iterator vrit = vroute.begin(); vrit != vroute.end(); ++vrit)
                    {
                        copy(vrit->begin() + 1, vrit->end(), route.begin() + n_pt_cur);
                        n_pt_cur += vrit->size() - 1;
                    }
                    copy(route_tmp2.begin() + 1, route_tmp2.end(), route.begin() + n_pt_cur);
                }
                else //(vtx_in < vtx_out)
                {
                    vector<Route> vroute;
                    if (vtx_dir == 0)
                    {
                        for (int k = vtx_in; k < vtx_out; k++)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k+1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }
                    }
                    else
                    {
                        for (int k = vtx_in; k > 0; k--)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k-1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k-1] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            int ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }

                        for (int k = nVtx; k > vtx_out; k--)
                        {
                            Point pt1(bl[k]);
                            pt1.zi = g_bl_zi;
                            pt1.t = cum_t;
                            Point pt2(bl[k-1]);
                            pt2.zi = g_bl_zi;
                            pt2.t = (arr.t - dep.t) / min_dis * vtx_vtx_dis[k-1] + cum_t;
                            cum_t = pt2.t;
                            Route route_tmp;
                            int ret = FlightPlanRoute(pt1, pt2, route_group, route_tmp, a, b);
                            if (ret != SUCCESS)
                            {
                                return ret;
                            }
                            vroute.push_back(route_tmp);
                        }
                    }
                    int n_pt = route_tmp1.size() + route_tmp2.size() - 1;
                    for (vector<Route>::const_iterator vrit = vroute.begin(); vrit != vroute.end(); ++vrit)
                    {
                        n_pt += vrit->size() - 1;
                    }

                    route.resize(n_pt);
                    copy(route_tmp1.begin(), route_tmp1.end(), route.begin());
                    int n_pt_cur = route_tmp1.size();
                    for (vector<Route>::const_iterator vrit = vroute.begin(); vrit != vroute.end(); ++vrit)
                    {
                        copy(vrit->begin() + 1, vrit->end(), route.begin() + n_pt_cur);
                        n_pt_cur += vrit->size() - 1;
                    }
                    copy(route_tmp2.begin() + 1, route_tmp2.end(), route.begin() + n_pt_cur);
                }
            }
        }
        else       // 有解但无穿过
        {
            int ret = FlightPlanRoute(dep, arr, route_group, route, a, b);
            if (ret != SUCCESS)
            {
                return ret;
            }
        }
    }
    else       // 不需要黑名单规划
    {
        int ret = FlightPlanRoute(dep, arr, route_group, route, a, b);
        if (ret != SUCCESS)
        {
            return ret;
        }
    }


    bool no_collision = true;

    // 进行50次循环，如果仍没有规划成功，则把第一个冲突点提前10个点重新递归规划
    for (int i = 0; i < 50; i++)
    {
        no_collision = true;
        for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
        {
            if (it->current_point >= 0)	// 激活路径
            {
                // 逐条路径判断是否与刚生成的 route 存在冲突，如果存在，则修正 route
                Route::const_iterator fc_it1;
                Route::const_iterator fc_it2;
                double nn_dis;

                if (FindFirstCollision(it->route, it->a, it->b, route, a, b, fc_it1, fc_it2, nn_dis) == SUCCESS)
                {
                    double delta_x = 0;
                    double delta_y = 0;
                    double delta_z = 0;
                    double delta_t = 0;
                    for (ActiveRouteGroup::const_iterator it1 = route_group.begin(); it1 != route_group.end(); ++it1)
                    {
                        Route::const_iterator nn_it;
                        double nn_dis;
                        NearestNeighbor(*fc_it2, it1->route, nn_it, nn_dis);

                        double scale = 5 / (pow2(nn_dis)  + DELTA);
                        delta_x += (fc_it2->xi - nn_it->xi) * scale;
                        delta_y += (fc_it2->yi - nn_it->yi) * scale;
                        delta_z += (fc_it2->zi - nn_it->zi) * scale;
                        delta_t += (fc_it2->t - nn_it->t) * scale;
                    }

                    Point mid;
                    mid.xi = fc_it2->xi + (int)(delta_x + 0.5);
                    mid.yi = fc_it2->yi + (int)(delta_y + 0.5);
                    mid.zi = fc_it2->zi + (int)(delta_z + 0.5);
                    mid.t = fc_it2->t + (int)(delta_t + 0.5);

                    Route route1;
                    Route route2;
                    int ret = FlightPlanRoute(dep, mid, route_group, route1, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }
                    ret = FlightPlanRoute(mid, arr, route_group, route2, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }

                    route.resize(route1.size() + route2.size() - 1);
                    copy(route1.begin(), route1.end(), route.begin());
                    copy(route2.begin() + 1, route2.end(), route.begin() + route1.size());

                    no_collision = false;
                }
                else
                {
                    no_collision = true;
                }
            }
        }

        if (no_collision)
        {
            break;
        }
    }

    if (!no_collision)
    {
        for (int i = 0; i < 50; i++)
        {
            bool no_collision = true;
            for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
            {
                if (it->current_point >= 0)	// 激活路径
                {
                    // 逐条路径判断是否与刚生成的 route 存在冲突，如果存在，则修正 route
                    Route::const_iterator fc_it1;
                    Route::const_iterator fc_it2;
                    double nn_dis;

                    if (FindFirstCollision(it->route, it->a, it->b, route, a, b, fc_it1, fc_it2, nn_dis) == SUCCESS)
                    {
                        Route::const_iterator fc_it2_less_10 = fc_it2;

                        while (fc_it2_less_10 != route.begin() + 10 && fc_it2_less_10 != fc_it2 - 10)    // 向前回溯10个点，同时要保证不能太接近起点
                        {
                            fc_it2_less_10--;
                        }
                        fc_it2 = fc_it2_less_10;

                        double delta_x = 0;
                        double delta_y = 0;
                        double delta_z = 0;
                        double delta_t = 0;
                        for (ActiveRouteGroup::const_iterator it1 = route_group.begin(); it1 != route_group.end(); ++it1)
                        {
                            Route::const_iterator nn_it;
                            double nn_dis;
                            NearestNeighbor(*fc_it2, it1->route, nn_it, nn_dis);

                            double scale = 5 / (pow2(nn_dis)  + DELTA);
                            delta_x += (fc_it2->xi - nn_it->xi) * scale;
                            delta_y += (fc_it2->yi - nn_it->yi) * scale;
                            delta_z += (fc_it2->zi - nn_it->zi) * scale;
                            delta_t += (fc_it2->t - nn_it->t) * scale;
                        }

                        Point mid;
                        mid.xi = fc_it2->xi + (int)(delta_x + 0.5);
                        mid.yi = fc_it2->yi + (int)(delta_y + 0.5);
                        mid.zi = fc_it2->zi + (int)(delta_z + 0.5);
                        mid.t = fc_it2->t + (int)(delta_t + 0.5);

                        Route route1;
                        Route route2;
                        int ret = FlightPlanRoute(dep, mid, route_group, route1, a, b);
                        if (ret != SUCCESS)
                        {
                            return ret;
                        }
                        ret = FlightPlanRoute(mid, arr, route_group, route2, a, b);
                        if (ret != SUCCESS)
                        {
                            return ret;
                        }

                        route.resize(route1.size() + route2.size() - 1);
                        copy(route1.begin(), route1.end(), route.begin());
                        copy(route2.begin() + 1, route2.end(), route.begin() + route1.size());

                        no_collision = false;
                    }
                    else
                    {
                        no_collision = true;
                    }
                }
            }
        }
    }


    return SUCCESS;
}

// 实现方式：可以认为起始点、白名单上的点，将整个航线分成多段；
// 针对每一段进行规划
int BL_WL_PlanRoute(const Point &dep, const Point &arr, const ActiveRouteGroup &route_group, const BlackList &bl, const WhiteList &wl, Route &route, double a, double b)
{
    WhiteList wl_c2;  // 第二类白名单
    wl_c2.reserve(wl.size() + 2);	// 加上起点和终点
    wl_c2.push_back(dep);
    wl_c2.insert(wl_c2.begin() + 1, wl.begin(), wl.end());
    wl_c2.push_back(arr);

    Route route_part;

    vector<double> dis_cum;

    route.push_back(dep);

    WhiteList::iterator pre_it = wl_c2.begin();
    for (WhiteList::iterator it = wl_c2.begin() + 1; it != wl_c2.end(); ++it)
    {
        if (it->t != 0)
        {
            int n = dis_cum.size();
            if (n == 0) // 之前没有第一类白名单
            {
                int ret = BLPlanRoute(*pre_it, *it, route_group, bl, route_part, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }
                route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
            }
            else
            {
                double dis_sum = 0;
                dis_sum = dis_cum[n - 1] + SpaceDistance(*(it - 1), *it);

                for (int i = 0; i < n; i++)
                {
                    (pre_it + i + 1)->t  = pre_it->t + (it->t - pre_it->t) * dis_cum[i] / dis_sum;

                    int ret = BLPlanRoute(*(pre_it + i), *(pre_it + i + 1), route_group, bl, route_part, a, b);
                    if (ret != SUCCESS)
                    {
                        return ret;
                    }
                    route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
                }

                int ret = BLPlanRoute(*(it - 1), *(it), route_group, bl, route_part, a, b);
                if (ret != SUCCESS)
                {
                    return ret;
                }
                route.insert(route.end() - 1, route_part.begin() + 1, route_part.end());
            }

            pre_it = it;

            dis_cum.clear();
        }
        else
        {
            double dis = SpaceDistance(*(it - 1), *it);
            int n = dis_cum.size();
            if (n == 0)
            {
                dis_cum.push_back(dis);
            }
            else
            {
                dis_cum.push_back(dis + dis_cum.back());
            }
        }
    }

    return SUCCESS;
}


int FlightAddRoute(ActiveRouteGroup &route_group, const Route &route, double a, double b)
{
    ActiveRoute aroute;
    aroute.current_point = 0;
    aroute.a = a;
    aroute.b = b;
    aroute.route = route;

    route_group.push_back(aroute);

    return SUCCESS;
}

int FlightAddRoute(ActiveRouteGroup &route_group, const Route &route, double a, double b, Point dep, Point arr, float r, float g, float blue)
{
	ActiveRoute aroute;
	aroute.current_point = 0;
	aroute.a = a;
	aroute.b = b;
	aroute.route = route;

	aroute.arr.xi = arr.xi;
	aroute.arr.yi = arr.yi;
	aroute.arr.zi = arr.zi;
	aroute.dep.xi = dep.xi;
	aroute.dep.yi = dep.yi;
	aroute.dep.zi = dep.zi;
	aroute.r = r;
	aroute.g = g;
	aroute.blue = blue;

	route_group.push_back(aroute);

	return SUCCESS;
}

// 根据当前点计算飞机控制参数
int GetControlPara(const ActiveRoute &route, ControlPara &para)
{

    return SUCCESS;
}

// 计算一架飞机是否在另一架飞机的位置，计算p0->p，与p->nextp两个矢量的夹角，夹角为锐角，则在后方；否则在前方或左右。
// 余弦定理：cos(c) = (a^2+b^2-c^2)/2ab
int GetThePositionOfOnePlanToAnother(const Point &p, const Point &nextp, const Point &p0)
{
	double a2 = pow2(p0.xi-p.xi) + pow2(p0.yi-p.yi) + pow2(p0.zi-p.zi);
	double b2 = pow2(nextp.xi - p.xi) + pow2(nextp.yi-p.yi) + pow2(nextp.zi-p.zi);
	double c2 = pow2(p0.xi-p.xi) + pow2(p0.yi-p.yi) + pow2(p0.zi - p.zi);

	// 其实只计算a^2+b^2-c^2的>0 =0 <0
	double d2 = a2 + b2 - c2;
	if (d2 > 1e-6)	// >0
	{
		return POS_FRONT;
	}
	else if (d2 < -(1e-6))
	{
		return POS_BACK;
	}
	else
	{
		return POS_SIDE;
	}

}

// 根据机型、相对位置获取a、b值
bool GetABByPositionAndPlaneType(int mType, int nType, int posType, double &a, double &b)
{
	// 此处规则，根据实际a、b的值

	return true;
}


void drawTwoLineWithArray();
// OpenGL
float delta = 0.02;
void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1, 0, 0);
    glPushMatrix();
	LookAt();

	int count = 0;
    for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
    {
		glPointSize(5);
		glBegin(GL_POINTS);
		glColor3f(it->r, it->g, it->blue);
		/*if (count%3==0)
		{
			glColor3f(it->r, it->g, it->blue);
		}
		else if (count%3==1)
		{
			glColor3f(0.0, 1.0, 0.0);
		}
		else
		{
			glColor3f(0.0, 0.0, 1.0);
		}*/
		float tx = it->route.begin()->xi*g_OpenGL_draw_delta;
		float ty = it->route.begin()->yi*g_OpenGL_draw_delta;
		glVertex2f(tx, ty);
		tx = it->route[it->route.size()-2].xi*g_OpenGL_draw_delta;
		ty = it->route[it->route.size()-2].yi*g_OpenGL_draw_delta;
		glVertex2f(tx, ty);
		glEnd();

		count++;
        if (it->current_point >= 0)	// 激活路径
        {
            for (size_t i = it->current_point; i < it->route.size(); i++)
            {
                if (it->route[i].t >= g_OpenGL_time)
                {
                    glTranslatef(it->route[i].xi * g_OpenGL_draw_delta, it->route[i].yi * g_OpenGL_draw_delta, it->route[i].zi * g_OpenGL_draw_delta);
                   
					glutWireSphere(0.002, 10, 8);
                    glTranslatef(-it->route[i].xi * g_OpenGL_draw_delta, -it->route[i].yi * g_OpenGL_draw_delta, -it->route[i].zi * g_OpenGL_draw_delta);
                    break;
                }
            }
        }
    }

	// 禁飞区
	glEnable(GL_LINE_STIPPLE);
	glLineWidth(1.0f);
	/*glBegin(GL_LINES);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_bl.size(); i++)
	{
		glVertex2f(g_bl[i].xi*g_OpenGL_draw_delta, g_bl[i].yi*g_OpenGL_draw_delta);
		if (i > 0 && i<g_bl.size()-1)
		{
			glVertex2f(g_bl[i].xi*g_OpenGL_draw_delta, g_bl[i].yi*g_OpenGL_draw_delta);
		}
	}
	glEnd();*/

	/*glBegin(GL_LINES);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_bl.size(); i++)
	{
		glVertex2f(g_bl2[i].xi*g_OpenGL_draw_delta, g_bl2[i].yi*g_OpenGL_draw_delta);
		if (i > 0 && i<g_bl2.size()-1)
		{
			glVertex2f(g_bl2[i].xi*g_OpenGL_draw_delta, g_bl2[i].yi*g_OpenGL_draw_delta);
		}
	}
	glEnd();*/
	/*glBegin(GL_LINES);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_blByConf.size(); i++)
	{
		glVertex2f(g_blByConf[i].xi*g_OpenGL_draw_delta, g_blByConf[i].yi*g_OpenGL_draw_delta);
		if (i > 0 && i<g_bl2.size()-1)
		{
			glVertex2f(g_blByConf[i].xi*g_OpenGL_draw_delta, g_blByConf[i].yi*g_OpenGL_draw_delta);
		}
	}
	glEnd();*/

	Point preP;

	glBegin(GL_POLYGON);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_bl.size(); i++)
	{
		if (i==0)
		{
			preP.xi = g_bl[i].xi;
			preP.yi = g_bl[i].yi;
			continue;
		}
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_bl[i].xi*g_OpenGL_draw_delta, g_bl[i].yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_bl[i].xi*g_OpenGL_draw_delta, g_bl[i].yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		preP.xi = g_bl[i].xi;
		preP.yi = g_bl[i].yi;
	}
	glEnd();

	glBegin(GL_POLYGON);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_bl2.size(); i++)
	{
		if (i==0)
		{
			preP.xi = g_bl2[i].xi;
			preP.yi = g_bl2[i].yi;
			continue;
		}
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_bl2[i].xi*g_OpenGL_draw_delta, g_bl2[i].yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_bl2[i].xi*g_OpenGL_draw_delta, g_bl2[i].yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		preP.xi = g_bl2[i].xi;
		preP.yi = g_bl2[i].yi;
	}
	glEnd();

	glBegin(GL_POLYGON);
	glColor3f(0.5, 0.5, 0.5);
	for (int i = 0; i < g_blByConf.size(); i++)
	{
		if (i==0)
		{
			preP.xi = g_blByConf[i].xi;
			preP.yi = g_blByConf[i].yi;
			continue;
		}
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_blByConf[i].xi*g_OpenGL_draw_delta, g_blByConf[i].yi*g_OpenGL_draw_delta, DiscreteZ(20000)*g_OpenGL_draw_delta);
		glVertex3f(g_blByConf[i].xi*g_OpenGL_draw_delta, g_blByConf[i].yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		glVertex3f(preP.xi*g_OpenGL_draw_delta, preP.yi*g_OpenGL_draw_delta, DiscreteZ(0)*g_OpenGL_draw_delta);
		preP.xi = g_blByConf[i].xi;
		preP.yi = g_blByConf[i].yi;
	}
	glEnd();
	
	// 
	//int tm, tn;
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// 线性模式画图
	//glColor3f(1.0f, 0.0f, 0.0f);				// 指定当前的绘图颜色
	//if (Function == 1)
	//{
	//	for (tm = 0; tm < n; tm++)
	//	{
	//		glBegin(GL_LINES);	// 绘制直线段
	//		glVertex2i(Line[tm].x1, Line[tm].y1);
	//		glVertex2i(Line[tm].x2, Line[tm].y2);
	//		glEnd();
	//	}
	//	// 动态绘制鼠标动作
	//	if (Flag == 1)
	//	{
	//		glBegin(GL_LINES);
	//		glVertex2i(Line[tm].x1, Line[tm].y1);
	//		glVertex2i(Mousex, Mousey);
	//		glEnd();
	//	}
	//}
	//else
	//{
	//	// 绘制矩形
	//	for (tn = 0; tn < m; tn++)
	//	{
	//		glRecti(Rect[tn].x1, Rect[tn].y1, Rect[tn].x2, Rect[tn].y2);
	//	}
	//	// 动态绘制鼠标动作
	//	if (RFlag == 1)
	//	{
	//		glRecti(Rect[tn].x1, Rect[tn].y1, Mousex, Mousey);
	//	}
	//}
	//
	// end

	// 白名单
	glPointSize(3.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);
	for (int i=0; i< g_wl.size(); i++)
	{
		glVertex2f(g_wl[i].xi*g_OpenGL_draw_delta, g_wl[i].yi*g_OpenGL_draw_delta);
	}
	//glVertex2f(tx*g_OpenGL_draw_delta*27, ty*g_OpenGL_draw_delta*27);
	LOG4CPLUS_INFO(_logger, tx+"\t"+ty);
	glEnd();

	glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);
	for (int i=0; i< g_wlByConf.size(); i++)
	{
		glVertex2f(g_wlByConf[i].xi*g_OpenGL_draw_delta, g_wlByConf[i].yi*g_OpenGL_draw_delta);
	}
	//glVertex2f(tx*g_OpenGL_draw_delta*27, ty*g_OpenGL_draw_delta*27);
	LOG4CPLUS_INFO(_logger, tx+"\t"+ty);
	glEnd();
	// end
	//drawTwoLineWithArray();

    glPopMatrix();
    glutSwapBuffers();

	//
	/*viewing transformation */  
	
	 
	//delta += -0.02;
	//glScalef(1.0,2.0,1.0); /*modeling transformation */  
	//glutWireCube(1.0);  
	//
    glFlush();
}


void init(void)
{
    glClearColor(0,0,0,0);
    glShadeModel(GL_FLAT);
}

void reshape(int w, int h)
{
	// 保存当前窗口的大小
	winWidth = w;
	winHeight = h;	
	glViewport(0, 0, w, h);// 指定窗口显示区域
	//glMatrixMode(GL_PROJECTION);	// 指定设置投影参数
	//glLoadIdentity();		// 调用单位矩阵，去掉以前的投影参数设置
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);
	// 保存当前窗口的大小
	winWidth = w;
	winHeight = h;	
	//glViewport(0, 0, w, h);// 指定窗口显示区域
	
	//gluOrtho2D(0.0, winWidth, 0.0, winHeight);

    glMatrixMode(GL_PROJECTION);
	LookAt();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW); 
    glLoadIdentity();

	glViewport(0, 0, w, h);// 指定窗口显示区域
	//gluOrtho2D(0.0, winWidth, 0.0, winHeight);
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 'f':
        g_OpenGL_time = g_OpenGL_time + g_OpenGL_time_step;
        //glutPostRedisplay();
        break;
    case 'b':
        g_OpenGL_time = g_OpenGL_time - g_OpenGL_time_step;
        //glutPostRedisplay();
        break;
	case 'a': RotateLeft(); break;
	case 'd': RotateRight();break;
	case 'w': MoveForward(); break;
	case 's': MoveBackward(); break;
	case 'q': RotateUp(); break;
	case 'e': RotateDown(); break;
	/*case GLUT_KEY_LEFT :
		angle -= 0.01f;
		lx = sin(angle);
		lz = -cos(angle);
		glutPostRedisplay();
		break;
	case GLUT_KEY_RIGHT :
		angle += 0.01f;
		lx = sin(angle);
		lz = -cos(angle);
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP :
		x += lx * fraction;
		z += lz * fraction;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN :
		x -= lx * fraction;
		z -= lz * fraction;
		glutPostRedisplay();*/
    default:
        break;
    }
	glutPostRedisplay();
}

void GetFlightPlanFromConf();
void StaticFlightPlanByHand();

void ProcessMenu(int value)
{
	Function = value;
	n = 0;
	Flag = 0;
	m = 0;
	RFlag = 0;
	if (Function == 1)
	{
		ShellExecute(NULL, L"open", L"NOTEPAD.EXE", L"conFile.txt", NULL, SW_SHOWNORMAL);
	}
	else if(Function == 2)
	{
		is3d = false;

	// 初始化路径
	route_group.clear();
	g_wl.clear();
	g_wlByConf.clear();
	g_bl.clear();
	g_bl2.clear();
	g_blByConf.clear();
	StaticFlightPlanByHand();
	GetFlightPlanFromConf();
	int time_min = INT_MAX;
	for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
	{
		if (it->current_point >= 0)	// 激活路径
		{
			if (it->route[it->current_point].t < time_min)
			{
				time_min = it->route[it->current_point].t;
			}
		}
	}
	g_OpenGL_time = time_min;
	}
	else if (Function == 3)
	{
		is3d = true;
	}
	else if (Function == 4)
	{
		is3d = false;
	}
	

	glutPostRedisplay();
}

void MousePlot(GLint button, GLint action, GLint xMouse, GLint yMouse)
{
	if (Function == 1)
	{
		if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)
		{
			tx = -winWidth/2 + xMouse;
			ty = winHeight/2 - yMouse;

			if (Flag == 0)
			{
				Flag = 1;
				Line[n].x1 = xMouse;
				Line[n].y1 = winHeight - yMouse;
			}
			else
			{
				Line[n].x2 = xMouse;
				Line[n].y2 = winHeight - yMouse;
				n++;
				// 折线的第二点作为下一段线的第一个的点
				Line[n].x1 = Line[n-1].x2;
				Line[n].y1 = Line[n-1].y2;
			}
		}
	}
	else
	{
		// 矩形处理
		if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)
		{
			tx = -winWidth/2 + xMouse;
			ty = winHeight/2 - yMouse;

			if (RFlag == 0)
			{
				RFlag = 1;
				Rect[m].x1 = xMouse;
				Rect[m].y1 = winHeight - yMouse;
			}
			else
			{
				RFlag = 0;
				Rect[m].x2 = xMouse;
				Rect[m].y2 = winHeight - yMouse;
				m++;
				glutPostRedisplay();
			}
		}
	}
}

void PassiveMouseMove(GLint xMouse, GLint yMouse)
{
	Mousex = xMouse;
	Mousey = winHeight - yMouse;
	glutPostRedisplay();
}


int FlightSim(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);   //缓存模式
    glutInitWindowSize(5000, 5000);    //显示框的大小
    //glutInitWindowPosition(400,400); //确定显示框左上角的位置

    int time_min = INT_MAX;
    for (ActiveRouteGroup::const_iterator it = route_group.begin(); it != route_group.end(); ++it)
    {
        if (it->current_point >= 0)	// 激活路径
        {
            if (it->route[it->current_point].t < time_min)
            {
                time_min = it->route[it->current_point].t;
            }
        }
    }
    g_OpenGL_time = time_min;

    glutCreateWindow("飞行模拟程序（f：向前；b：向后）");
	//
	glutCreateMenu(ProcessMenu);
	glutAddMenuEntry("打开配置文件", 1);
	glutAddMenuEntry("重新加载配置", 2);
	glutAddMenuEntry("三维视图", 3);
	glutAddMenuEntry("截面图", 4);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	glutMouseFunc(MousePlot);
	glutPassiveMotionFunc(PassiveMouseMove);	// 指定鼠标移动响应函数
	//
    init();
    glutDisplayFunc(display);	// display/reshape/keyboard 三个函数需要自己实现，让opengl库调用。
	//glutDisplayFunc(drawTwoLineWithArray);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMainLoop(); //进人GLUT事件处理循环
    return 0;
}

// 手工设定的飞行计划（航线）
void StaticFlightPlanByHand()
{
	// 白名单，一定要经过的点
	WhiteList wl;
	{
		Point p;
		p.xi = DiscreteX(30000.0);
		p.yi = DiscreteY(90000.0);
		p.zi = DiscreteZ(10000.0);
		p.t = 0;

		wl.push_back(p);
		g_wl.push_back(p);

		p.xi = DiscreteX(70000.0);
		p.yi = DiscreteY(70000.0);
		p.zi = DiscreteZ(10000.0);
		p.t = 0;

		wl.push_back(p);
		g_wl.push_back(p);
	}

	// 第一类黑名单，一个正方形区域。 各边首尾顺次相接，可以任意多条边，凸多边形闭合区域。
	BlackList bl;
	{
		Point p1;
		p1.xi = DiscreteX(30000.0);
		p1.yi = DiscreteY(30000.0);
		p1.zi = DiscreteZ(10000.0);
		p1.t = 0;

		Point p2;
		p2.xi = DiscreteX(70000.0);
		p2.yi = DiscreteY(30000.0);
		p2.zi = DiscreteZ(10000.0);
		p2.t = 0;

		Point p3;
		p3.xi = DiscreteX(70000.0);
		p3.yi = DiscreteY(70000.0);
		p3.zi = DiscreteZ(10000.0);
		p3.t = 0;

		Point p4;
		p4.xi = DiscreteX(30000.0);
		p4.yi = DiscreteY(90000.0);
		p4.zi = DiscreteZ(10000.0);
		p4.t = 0;

		//Point p5; // new
		//p5.xi = DiscreteX(50000.0);
		//p5.yi = DiscreteY(00000.0);
		//p5.zi = DiscreteZ(10000.0);
		//p5.t = 0;

		//Point p6; // new
		//p6.xi = DiscreteX(50000.0);
		//p6.yi = DiscreteY(100000.0);
		//p6.zi = DiscreteZ(10000.0);
		//p6.t = 0;

		bl.push_back(p1);
		//bl.push_back(p5); //new
		bl.push_back(p2);
		bl.push_back(p3);
		//bl.push_back(p6); //new
		bl.push_back(p4);
		bl.push_back(p1);

		g_bl.push_back(p1);
		g_bl.push_back(p2);
		g_bl.push_back(p3);
		g_bl.push_back(p4);
		g_bl.push_back(p1);
	}

	BlackList bl2;
	{
		Point p1;
		p1.xi = DiscreteX(10000.0);
		p1.yi = DiscreteY(170000.0);
		p1.zi = DiscreteZ(10000.0);
		p1.t = 0;

		Point p2;
		p2.xi = DiscreteX(30000.0);
		p2.yi = DiscreteY(170000.0);
		p2.zi = DiscreteZ(10000.0);
		p2.t = 0;

		Point p3;
		p3.xi = DiscreteX(30000.0);
		p3.yi = DiscreteY(130000.0);
		p3.zi = DiscreteZ(10000.0);
		p3.t = 0;

		Point p4;
		p4.xi = DiscreteX(10000.0);
		p4.yi = DiscreteY(130000.0);
		p4.zi = DiscreteZ(10000.0);
		p4.t = 0;

		g_bl2.push_back(p1);
		g_bl2.push_back(p2);
		g_bl2.push_back(p3);
		g_bl2.push_back(p4);
		g_bl2.push_back(p1);
	}

	Route route1;
	double route1_a = 2000;
	double route1_b = 1000;
	{
		Point dep;
		dep.xi = DiscreteX(0.0);
		dep.yi = DiscreteY(0.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 0);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(100000.0);
		arr.yi = DiscreteY(100000.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 0);
		arr.t = t2.GetTime();

		int ret = LineRoute(dep, arr, route1);

		FlightAddRoute(route_group, route1, route1_a, route1_b, dep, arr, 1.0, 0.0, 0.0);
		//FILE *fp = fopen("route1.dat", "wt");
		//PrintRoute(fp, route1, 0);
		//fclose(fp);
	}

	Route route2;
	//{
	//    Point dep;
	//    dep.xi = DiscreteX(0.0);
	//    dep.yi = DiscreteY(100000.0);
	//    dep.zi = DiscreteZ(10000.0);
	//    CTime t1(2013, 6, 1, 22, 0, 1);
	//    dep.t = t1.GetTime();

	//    Point arr;
	//    arr.xi = DiscreteX(100000.0);
	//    arr.yi = DiscreteY(0.0);
	//    arr.zi = DiscreteZ(10000.0);
	//    CTime t2(2013, 6, 1, 22, 50, 1);
	//    arr.t = t2.GetTime();

	//    int ret = LineRoute(dep, arr, route2);

	//    FILE *fp = fopen("route2_line.dat", "wt");
	//    PrintRoute(fp, route2, 0);
	//    fclose(fp);
	//}
	//double route1_a = 2000;
	//double route1_b = 1000;

	double route2_a = 1000;
	double route2_b = 500;
	//FlightAddRoute(route_group, route1, route1_a, route1_b);
	{
		Point dep;
		dep.xi = DiscreteX(0.0);
		dep.yi = DiscreteY(100000.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 1);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(100000.0);
		arr.yi = DiscreteY(0.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 1);
		arr.t = t2.GetTime();

		//int ret = FlightPlanRoute(dep, arr, route_group, route2);
		//int ret = WLPlanRoute(dep, arr, route_group, wl, route2);
		//int ret = BLPlanRoute(dep, arr, route_group, bl, route2);
		int ret = BL_WL_PlanRoute(dep, arr, route_group, bl, wl, route2, route2_a, route2_b);

		FlightAddRoute(route_group, route2, route2_a, route2_b, dep, arr, 0.0, 1.0, 0.0);
		//FILE *fp = fopen("route2.dat", "wt");
		//PrintRoute(fp, route2, 0);
		//fclose(fp);
	}

	//FlightAddRoute(route_group, route2, route2_a, route2_b);

	// 添加第3条路径： 
	// a) 声明 Route型变量 route3
	// b) 定义该航线飞机的椭球型防撞区域的长轴route3_a、短轴route3_b
	// c) 定义该航线的起点、终点、及时间
	// d) 在现有航线、黑白名单情况下，规划航线。调用函数BL_WL_PlanRoute()
	// e）将本航线添加到现有航线组中。——FlightAddRoute()
	Route route3;
	double route3_a = 1000;
	double route3_b = 500;
	{
		Point dep;
		dep.xi = DiscreteX(-50000.0);
		dep.yi = DiscreteY(0.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 0);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(0.0);
		arr.yi = DiscreteY(-50000.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 0);
		arr.t = t2.GetTime();

		//int ret = FlightPlanRoute(dep, arr, route_group, route2);
		//int ret = WLPlanRoute(dep, arr, route_group, wl, route2);
		//int ret = BLPlanRoute(dep, arr, route_group, bl, route2);
		int ret = BL_WL_PlanRoute(dep, arr, route_group, bl, wl, route3, route3_a, route3_b);

		FlightAddRoute(route_group, route3, route3_a, route3_b, dep, arr, 0.0, 0.0, 1.0);
	}

	Route route4;
	double route4_a = 1000;
	double route4_b = 500;
	{
		Point dep;
		dep.xi = DiscreteX(-150000.0);
		dep.yi = DiscreteY(120000.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 0);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(120000.0);
		arr.yi = DiscreteY(120000.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 0);
		arr.t = t2.GetTime();

		//int ret = FlightPlanRoute(dep, arr, route_group, route2);
		//int ret = WLPlanRoute(dep, arr, route_group, wl, route2);
		//int ret = BLPlanRoute(dep, arr, route_group, bl, route2);
		int ret = BLPlanRoute(dep, arr, route_group, bl, route4, route4_a, route4_b);

		FlightAddRoute(route_group, route4, route4_a, route4_b, dep, arr, 0.0, 1.0, 1.0);
	}

	Route route5;
	double route5_a = 1000;
	double route5_b = 500;
	{
		Point dep;
		dep.xi = DiscreteX(-150000.0);
		dep.yi = DiscreteY(100000.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 0);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(100100.0);
		arr.yi = DiscreteY(70100.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 0);
		arr.t = t2.GetTime();

		//int ret = FlightPlanRoute(dep, arr, route_group, route2);
		//int ret = WLPlanRoute(dep, arr, route_group, wl, route2);
		//int ret = BLPlanRoute(dep, arr, route_group, bl, route2);
		int ret = BLPlanRoute(dep, arr, route_group, bl, route5, route5_a, route5_b);

		FlightAddRoute(route_group, route5, route5_a, route5_b, dep, arr, 1.0, 0.0, 1.0);
	}

	Route route6;
	double route6_a = 1000;
	double route6_b = 500;
	{
		Point dep;
		dep.xi = DiscreteX(-150000.0);
		dep.yi = DiscreteY(150000.0);
		dep.zi = DiscreteZ(10000.0);
		CTime t1(2013, 6, 1, 22, 0, 0);
		dep.t = t1.GetTime();

		Point arr;
		arr.xi = DiscreteX(150000.0);
		arr.yi = DiscreteY(150000.0);
		arr.zi = DiscreteZ(10000.0);
		CTime t2(2013, 6, 1, 22, 50, 0);
		arr.t = t2.GetTime();

		//int ret = FlightPlanRoute(dep, arr, route_group, route2);
		//int ret = WLPlanRoute(dep, arr, route_group, wl, route2);
		//int ret = BLPlanRoute(dep, arr, route_group, bl, route2);
		int ret = BLPlanRoute(dep, arr, route_group, g_bl2, route6, route6_a, route6_b);

		FlightAddRoute(route_group, route6, route6_a, route6_b, dep, arr, 0.3, 0.7, 0.0);
	}
}

// 根据配置文件生成飞行计划（航线）
void GetFlightPlanFromConf()
{
	ifstream fin(conFileName, std::ios::in);
	char line[1024]={0};
	Point dep;
	Point arr;
	Route route;
	double route_a = 1000;
	double route_b = 500;
	double year, month, day, hour, minite, second;
	double x, y, z;
	int count = 0;
	{
	{
		// 生成白名单
		fin.getline(line, sizeof(line));
		std::stringstream num(line);
		num >> count;
		for (int i=0; i<count; i++)
		{
			Point p;

			fin.getline(line, sizeof(line));
			std::stringstream word(line);
			word >> x; word >> y; word >> z;
			p.xi = DiscreteX(x);
			p.yi = DiscreteY(y);
			p.zi = DiscreteZ(z);
			p.t = 0;
			g_wlByConf.push_back(p);
		}
		fin.getline(line, sizeof(line));	// 读取分隔符行
	}

		// 生成黑名单禁飞区
		{
			fin.getline(line, sizeof(line));
			std::stringstream num(line);
			num >> count;
			Point firstP;
			if (count > 2)	// 最小是3条边
			{
				for (int j=0; j<count; j++)
				{
					Point p0;
					fin.getline(line, sizeof(line));
					std::stringstream word(line);
					word >> x; word >> y; word >> z;
					p0.xi = DiscreteX(x);
					p0.yi = DiscreteY(y);
					p0.zi = DiscreteZ(z);
					p0.t = 0;
					if (j==0)
					{
						firstP.xi = p0.xi;
						firstP.yi = p0.yi;
						firstP.zi = p0.zi;
						firstP.t = p0.t;
					}
					g_blByConf.push_back(p0);
				}
				g_blByConf.push_back(firstP);
			}
			fin.getline(line, sizeof(line));	// 读取分隔符行
		}


		// 读取飞机个数
		{
			fin.getline(line, sizeof(line));
			std::stringstream num(line);
			num >> count;
			for (int k=0; k<count; k++)
			{
				Route tmpRoute;
				// 读取dep
				{
					fin.getline(line, sizeof(line));
					std::stringstream word(line);
					word >> x; word >> y; word >> z;
					dep.xi = DiscreteX(x);
					dep.yi = DiscreteY(y);
					dep.zi = DiscreteZ(z);
					word >> year; word>>month; word >> day; word >> hour; word >> minite; word >> second; 
					CTime t0(year, month, day, hour, minite, second);
					dep.t = t0.GetTime();
				}
				// 读取arr
				{
					fin.getline(line, sizeof(line));
					std::stringstream word(line);
					word >> x; word >> y; word >> z;
					arr.xi = DiscreteX(x);
					arr.yi = DiscreteY(y);
					arr.zi = DiscreteZ(z);
					word >> year; word>>month; word >> day; word >> hour; word >> minite; word >> second; 
					CTime t0(year, month, day, hour, minite, second);
					arr.t = t0.GetTime();
				}

				// 生成路径
				int ret = BL_WL_PlanRoute(dep, arr, route_group, g_blByConf, g_wlByConf, tmpRoute, route_a, route_b);
				//int ret = LineRoute(dep, arr, route);
				if (ret != SUCCESS)
				{
					LOG4CPLUS_INFO(_logger, "guihua chenggong");
				}
				else
				{
					LOG4CPLUS_INFO(_logger, "guihua chenggong");
				}

				FlightAddRoute(route_group, tmpRoute, route_a, route_b, dep, arr, 0.7, 0.4, 0.3);
			}
			//fin.getline(line);
		}
	}
	fin.clear();
	fin.close();
}

/**
 * 入口函数
 */
int main(int argc, char *argv[])
{

	// 初始化日志
	SharedAppenderPtr _append(new FileAppender(LOG4CPLUS_TEXT("fg.log")));
	_append->setName(LOG4CPLUS_TEXT("file log test"));
	std::auto_ptr<Layout> _layout(new PatternLayout(LOG4CPLUS_TEXT("%D [%l] %m%n")));
	_append->setLayout(_layout);
	_logger.addAppender(_append);
	LOG4CPLUS_INFO(_logger, "FlightPlan start.");

	
	StaticFlightPlanByHand();
	GetFlightPlanFromConf();

	// 记录航线轨迹
	if (true)
	{
		if (route_group.size() > 0)
		{
			ofstream frecord("record.txt");
			if (frecord)
			{
				for (int i = 0; i < route_group.size(); i++ )
				{
					Route tmpRoute = route_group[i].route;
					frecord << "第"<<i<<"条航线， 共："<< tmpRoute.size()<< "个点" << endl;
					for (int j = 0; j < tmpRoute.size(); j++)
					{
						if ((j%30==0) || (j==tmpRoute.size()-1))
						{
							frecord << "(" << tmpRoute[j].xi << "," << tmpRoute[j].yi << "," << tmpRoute[j].zi<< ")" << "\t";
						}
					}
					frecord << endl;
				}
				frecord.close();
			}
		}
	}

	// 计算好了各个航线，下边开始模

	//
    FlightSim(argc, argv);

	LOG4CPLUS_INFO(_logger, "FlightPlan finished.");
    return SUCCESS;
}