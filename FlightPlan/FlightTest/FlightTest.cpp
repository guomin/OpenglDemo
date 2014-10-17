// FlightTest.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "glut.h"
#include <cmath>

using namespace std;
/*
#define NUM 100 // 折线的最大折线段

int Flag = 0;	// 标记是否已经开始绘制折线
int RFlag = 0;	// 标记是否已经完成一个矩形
int Function = 1;	// 标记选择的功能是画折线还是矩形
int winWidth = 800, winHeight = 600;	// 窗口的宽度和高度
int Mousex, Mousey;	// 用于记录当前鼠标的位置
int n = 0;			// 用于记录折线有几段
int m = 0;			// 用于记录矩形个数

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

void Initial(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);	// 设置窗口背景颜色
}

void ChangeSize(int w, int h)
{
	// 保存当前窗口的大小
	winWidth = w;
	winHeight = h;	
	glViewport(0, 0, w, h);// 指定窗口显示区域
	glMatrixMode(GL_PROJECTION);	// 指定设置投影参数
	glLoadIdentity();		// 调用单位矩阵，去掉以前的投影参数设置
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);
}

void ProcessMenu(int value)
{
	Function = value;
	n = 0;
	Flag = 0;
	m = 0;
	RFlag = 0;
	glutPostRedisplay();
}

void Display()
{
	int i, j;
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// 线性模式画图
	glClear(GL_COLOR_BUFFER_BIT);				// 用当前背景色填充窗口
	glColor3f(1.0f, 0.0f, 0.0f);				// 指定当前的绘图颜色
	if (Function == 1)
	{
		for (i = 0; i < n; i++)
		{
			glBegin(GL_LINES);	// 绘制直线段
			glVertex2i(Line[i].x1, Line[i].y1);
			glVertex2i(Line[i].x2, Line[i].y2);
			glEnd();
		}
		// 动态绘制鼠标动作
		if (Flag == 1)
		{
			glBegin(GL_LINES);
			glVertex2i(Line[i].x1, Line[i].y1);
			glVertex2i(Mousex, Mousey);
			glEnd();
		}
	}
	else
	{
		// 绘制矩形
		for (j = 0; j < m; j++)
		{
			glRecti(Rect[j].x1, Rect[j].y1, Rect[j].x2, Rect[j].y2);
		}
		// 动态绘制鼠标动作
		if (RFlag == 1)
		{
			glRecti(Rect[j].x1, Rect[j].y1, Mousex, Mousey);
		}
	}
	glutSwapBuffers();	// 交换缓冲区
}

void MousePlot(GLint button, GLint action, GLint xMouse, GLint yMouse)
{
	if (Function == 1)
	{
		if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)
		{
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

int _tmain(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);	// 使用双缓冲及RGB模型
	glutInitWindowSize(400, 300);	// 指定窗口的尺寸
	glutInitWindowPosition(100, 100);	// 指定窗口在屏幕上的位置
	glutCreateWindow("橡皮筋技术");
	glutCreateMenu(ProcessMenu);
	glutAddMenuEntry("画折线", 1);
	glutAddMenuEntry("画矩形", 2);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	glutDisplayFunc(Display);
	glutReshapeFunc(ChangeSize);
	glutMouseFunc(MousePlot);
	glutPassiveMotionFunc(PassiveMouseMove);	// 指定鼠标移动响应函数
	Initial();
	glutMainLoop();
	return 0;
}
*/
const GLfloat PI = 3.14;

/// record the state of mouse
GLboolean mouserdown = GL_FALSE;
GLboolean mouseldown = GL_FALSE;
GLboolean mousemdown = GL_FALSE;

/// when a mouse-key is pressed, record current mouse position 
static GLint mousex = 0, mousey = 0;

static GLfloat center[3] = {0.0f, 0.0f, 0.0f}; /// center position
static GLfloat eye[3] = {0.0f, 0.0f, 5.0f}; /// eye's position

static GLfloat yrotate = PI/4; /// angle between y-axis and look direction
static GLfloat xrotate = PI/4; /// angle between x-axis and look direction
static GLfloat celength = 20.0f;/// lenght between center and eye

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
	eye[0] = center[0] + celength * sin(yrotate) * cos(xrotate);  
	eye[2] = center[2] + celength * sin(yrotate) * sin(xrotate);
	eye[1] = center[1] + celength * cos(yrotate);
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

void LookAt()            /// 调用 gluLookAt(), 主要嫌直接调用要每次都写好几个参数。。
{
	CalEyePosition();
	gluLookAt(eye[0], eye[1], eye[2],
		center[0], center[1], center[2],
		0, 1, 0);
}

void init(void)  
{  
	glClearColor(0.0,0.0,0.0,0.0);  
	glShadeModel(GL_FLAT);  
}  
void display(void)  
{  
	glClear(GL_COLOR_BUFFER_BIT);  
	glColor3f(1.0,1.0,1.0);  
	glLoadIdentity(); /*clear the matrix */  
	/*viewing transformation */  
	//gluLookAt(0.0,0.0,5.0,0.0,0.0,0.0,0.0,1.0,0.0);
	LookAt();
	//glScalef(1.0,2.0,1.0); /*modeling transformation */  
	glutWireCube(1.0);  
	glFlush();  
}  
void reshape(int w,int h)  
{  
	glViewport(0,0,(GLsizei)w,(GLsizei)h);  
	glMatrixMode(GL_PROJECTION);  
	glLoadIdentity();  
	glFrustum(-1.0,1.0,-1.0,1.0,1.5,20.0);  
	glMatrixMode(GL_MODELVIEW);  
	LookAt();
}  
int main(int argc,char**argv)  
{  
	glutInit(&argc,argv);  
	glutInitDisplayMode(GLUT_SINGLE |GLUT_RGB);  
	glutInitWindowSize(500,500);  
	glutInitWindowPosition(100,100);  
	glutCreateWindow(argv [0]);  
	init();  
	glutKeyboardFunc(KeyFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MouseMotion);
	glutDisplayFunc(display);  
	glutReshapeFunc(reshape);  
	glutMainLoop();  
	return 0;  
} 
