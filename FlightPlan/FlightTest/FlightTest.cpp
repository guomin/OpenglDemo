// FlightTest.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "glut.h"
#include <cmath>

using namespace std;
/*
#define NUM 100 // ���ߵ�������߶�

int Flag = 0;	// ����Ƿ��Ѿ���ʼ��������
int RFlag = 0;	// ����Ƿ��Ѿ����һ������
int Function = 1;	// ���ѡ��Ĺ����ǻ����߻��Ǿ���
int winWidth = 800, winHeight = 600;	// ���ڵĿ�Ⱥ͸߶�
int Mousex, Mousey;	// ���ڼ�¼��ǰ����λ��
int n = 0;			// ���ڼ�¼�����м���
int m = 0;			// ���ڼ�¼���θ���

// ���Խṹ��
struct LineNode {
	int x1;
	int y1;
	int x2;
	int y2;
}Line[NUM];
// ���νṹ��
struct Rectangle {
	int x1;
	int y1;
	int x2;
	int y2;
}Rect[NUM];

void Initial(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);	// ���ô��ڱ�����ɫ
}

void ChangeSize(int w, int h)
{
	// ���浱ǰ���ڵĴ�С
	winWidth = w;
	winHeight = h;	
	glViewport(0, 0, w, h);// ָ��������ʾ����
	glMatrixMode(GL_PROJECTION);	// ָ������ͶӰ����
	glLoadIdentity();		// ���õ�λ����ȥ����ǰ��ͶӰ��������
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
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// ����ģʽ��ͼ
	glClear(GL_COLOR_BUFFER_BIT);				// �õ�ǰ����ɫ��䴰��
	glColor3f(1.0f, 0.0f, 0.0f);				// ָ����ǰ�Ļ�ͼ��ɫ
	if (Function == 1)
	{
		for (i = 0; i < n; i++)
		{
			glBegin(GL_LINES);	// ����ֱ�߶�
			glVertex2i(Line[i].x1, Line[i].y1);
			glVertex2i(Line[i].x2, Line[i].y2);
			glEnd();
		}
		// ��̬������궯��
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
		// ���ƾ���
		for (j = 0; j < m; j++)
		{
			glRecti(Rect[j].x1, Rect[j].y1, Rect[j].x2, Rect[j].y2);
		}
		// ��̬������궯��
		if (RFlag == 1)
		{
			glRecti(Rect[j].x1, Rect[j].y1, Mousex, Mousey);
		}
	}
	glutSwapBuffers();	// ����������
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
				// ���ߵĵڶ�����Ϊ��һ���ߵĵ�һ���ĵ�
				Line[n].x1 = Line[n-1].x2;
				Line[n].y1 = Line[n-1].y2;
			}
		}
	}
	else
	{
		// ���δ���
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
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);	// ʹ��˫���弰RGBģ��
	glutInitWindowSize(400, 300);	// ָ�����ڵĳߴ�
	glutInitWindowPosition(100, 100);	// ָ����������Ļ�ϵ�λ��
	glutCreateWindow("��Ƥ���");
	glutCreateMenu(ProcessMenu);
	glutAddMenuEntry("������", 1);
	glutAddMenuEntry("������", 2);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	glutDisplayFunc(Display);
	glutReshapeFunc(ChangeSize);
	glutMouseFunc(MousePlot);
	glutPassiveMotionFunc(PassiveMouseMove);	// ָ������ƶ���Ӧ����
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
	if(yrotate > PI/2.2) yrotate = PI/2.2;   /// ���ƿ��÷���
	if(yrotate < 0.01)  yrotate = 0.01;
	if(xrotate > 2*PI)   xrotate = 0.01;
	if(xrotate < 0)   xrotate = 2 * PI;
	if(celength > 50)  celength = 50;     ///  ���ž�������
	if(celength < 5)   celength = 5;
	/// ��������������ϵ���� eye ��λ�ã�
	eye[0] = center[0] + celength * sin(yrotate) * cos(xrotate);  
	eye[2] = center[2] + celength * sin(yrotate) * sin(xrotate);
	eye[1] = center[1] + celength * cos(yrotate);
}

/// center moves
void MoveBackward()              /// center �������߷���ˮƽ����ƶ�
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
	{       /// �����Ե������ǵ�����ת�ٶȵģ�������ã��ﵽ�Լ���Ҫ�ٶȼ���
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

void LookAt()            /// ���� gluLookAt(), ��Ҫ��ֱ�ӵ���Ҫÿ�ζ�д�ü�����������
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
