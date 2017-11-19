/**************************************************
 *						  *
 * Program that implements Hermit, Interpolation  *
 * and Bezier curves with OpenGL.		  *
 *						  *	
 *       by					  *
 *						  *
 * Vlasis Gogousis [vgogousis@gmail.com]    	  *
 * Kiki Hatzistavrou [kikihatzistavrou@gmail.com] *
 *						  *
 * Aristotle University of Thessaloniki, 2014	  *
 *						  *
 **************************************************/

#include <math.h>
#include <stdio.h>
#include <GL/glut.h>
#include <iostream>

using namespace std;

#define T 100 // number of u subintervals
#define HERMIT 1
#define CUBIC_INTERPOLATION 2
#define BEZIER_P0_EQ_P6 3
#define BEZIER_DOUBLE 4

static GLfloat **ctrlPoints; //control points
int typeMode = 1; //type of curve
static int N = 0; //number of control points

int i;
int j;
int uInt;
// window size
int ww = 1000;
int wh = 800;
// remember the moving control name
int MOVENAME = -1;
// set up pick radius for detecting movement of a control point
int pickRadius = 50;

void createCtrlPoints()
{
	if ((typeMode == HERMIT) || (typeMode == CUBIC_INTERPOLATION)) N = 4;
	else if ((typeMode == BEZIER_P0_EQ_P6) || (typeMode == BEZIER_DOUBLE)) N = 7;

	ctrlPoints = new GLfloat*[N];
	for (i = 0; i<N; i++) ctrlPoints[i] = new GLfloat[3];
	for (i = 0; i<N; i++)
	{
		ctrlPoints[i][0] = 100 * (i + 1);
		ctrlPoints[i][1] = 50 * (i + 1) + 300;
		ctrlPoints[i][2] = 0;
	}
}

void display2DControlPoints()
{
	glPointSize(5.0);
	glColor3f(0.0f, 0.0f, 0.0f);

	glBegin(GL_POINTS);
	for (i = 0; i < N; i++)
	{
		glVertex2i(ctrlPoints[i][0], ctrlPoints[i][1]);
	}
	glEnd();

	glFlush();
}
//factorial of n
int fact(int n)
{
	if (n == 0) return 1;
	if (n == 1) return 1;
	if (n>1) return n*fact(n - 1);
	else return 0;
}
//binomial coeffecient n over k
int binomial(int n, int k)
{
	return fact(n) / (fact(k)*(fact(n - k)));
}

void interpolate()
{
	glLineWidth(1.0);
	glColor3f(1.0f, 0.0f, 0.0f);

	GLfloat *bu = new GLfloat[N];

	glBegin(GL_LINE_STRIP);

	for (uInt = 0; uInt <= T; uInt++)
	{

		GLfloat u = uInt / (GLfloat)T;
		bu[0] = (GLfloat)((-9.0 / 2.0)*(u - (1.0 / 3.0))*(u - (2.0 / 3.0))*(u - 1.0));
		bu[1] = (27.0 / 2.0)*u*(u - (2.0 / 3.0))*(u - 1.0);
		bu[2] = (-27.0 / 2.0)*u*(u - (1.0 / 3.0))*(u - 1.0);
		bu[3] = (9.0 / 2.0)*u*(u - (1.0 / 3.0))*(u - (2.0 / 3.0));

		GLfloat x = 0.0;
		GLfloat y = 0.0;

		for (i = 0; i < N; i++)
		{
			x += bu[i] * ctrlPoints[i][0];
			y += bu[i] * ctrlPoints[i][1];
		}
		glVertex2i(x, y);
	}

	glEnd();
	glFlush();
	display2DControlPoints();
}

void hermite()
{
	glLineWidth(1.0);
	glColor3f(0.0f, 1.0f, 0.0f);

	glBegin(GL_LINE_STRIP);
	glVertex2i(ctrlPoints[0][0], ctrlPoints[0][1]);
	glVertex2i(ctrlPoints[1][0], ctrlPoints[1][1]);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex2i(ctrlPoints[2][0], ctrlPoints[2][1]);
	glVertex2i(ctrlPoints[3][0], ctrlPoints[3][1]);
	glEnd();

	glLineWidth(1.0);
	glColor3f(1.0f, 0.0f, 0.0f);

	GLfloat *fu = new GLfloat[N];

	glBegin(GL_LINE_STRIP);

	for (uInt = 0; uInt <= T; uInt++)
	{
		GLfloat u = uInt / (GLfloat)T;
		GLfloat u2 = u*u;
		GLfloat u3 = u2*u;
		fu[0] = 2.0*u3 - 3.0*u2 + 1.0;
		fu[1] = -2.0*u3 + 3.0*u2;
		fu[2] = u3 - 2.0*u2 + u;
		fu[3] = u3 - u2;

		GLfloat x = 0.0;
		GLfloat y = 0.0;

		// p0 = ctrlPoints[0]
		x += fu[0] * ctrlPoints[0][0];
		y += fu[0] * ctrlPoints[0][1];

		// p1 = ctrlPoints[3]
		x += fu[1] * ctrlPoints[3][0];
		y += fu[1] * ctrlPoints[3][1];

		//  tangent at p0 = ctrlPoints[1]-ctrlPoints[0]
		x += fu[2] * (1.0*(ctrlPoints[1][0] - ctrlPoints[0][0]));
		y += fu[2] * (1.0*(ctrlPoints[1][1] - ctrlPoints[0][1]));

		//  tangent at p1 = ctrlPoints[3]-ctrlPoints[2]
		x += fu[3] * (1.0*(ctrlPoints[2][0] - ctrlPoints[3][0]));
		y += fu[3] * (1.0*(ctrlPoints[2][1] - ctrlPoints[3][1]));

		glVertex2i(x, y);
	}
	glEnd();
	glFlush();
	display2DControlPoints();
}
//b(u) coeffecient
GLfloat bernstein(GLfloat u, int i, int n)
{
	GLfloat coeff = binomial(n, i);
	GLfloat ui = pow(u, i);
	GLfloat oneMinusu = pow(1.0 - u, n - i);

	return (coeff*ui*oneMinusu);
}
//bezier curve from point P0 to Pk (rank = Pk - P0)
void bezier(int P0, int Pk)
{
	glLineWidth(1.0);
	glColor3f(1.0f, 0.0f, 0.0f);

	glBegin(GL_LINE_STRIP);

	for (uInt = 0; uInt <= T; uInt++)
	{
		GLfloat u = uInt / (GLfloat)T;

		GLfloat x = 0.0;
		GLfloat y = 0.0;

		for (i = P0; i <= Pk; i++)
		{
			GLfloat b = bernstein(u, i - P0, Pk - P0);
			x += b*ctrlPoints[i][0];
			y += b*ctrlPoints[i][1];
		}
		glVertex2i(x, y);
	}

	glEnd();
	glFlush();
	display2DControlPoints();
}

void doubleBezier()
{
	//draw curve from point P0 to P3
	bezier(0, 3);
	//draw curve from point P3 to P6
	bezier(3, 6);
}
//make control points co-linear, based on point Pk being moved
void linearControlPoints(int Pk)
{
	if (Pk == 2)
	{
		float a = (ctrlPoints[3][1] - ctrlPoints[2][1]) / (ctrlPoints[3][0] - ctrlPoints[2][0]);
		int b = -a*ctrlPoints[3][0] + ctrlPoints[3][1];
		ctrlPoints[4][0] = 2 * ctrlPoints[3][0] - ctrlPoints[2][0];
		ctrlPoints[4][1] = a*ctrlPoints[4][0] + b;
	}
	if (Pk == 4)
	{
		float a = (ctrlPoints[4][1] - ctrlPoints[3][1]) / (ctrlPoints[4][0] - ctrlPoints[3][0]);
		int b = -a*ctrlPoints[3][0] + ctrlPoints[3][1];
		ctrlPoints[2][0] = 2 * ctrlPoints[3][0] - ctrlPoints[4][0];
		ctrlPoints[2][1] = a*ctrlPoints[2][0] + b;
	}
	if (Pk == 3)
	{
		ctrlPoints[4][0] = 2 * ctrlPoints[3][0] - ctrlPoints[2][0];

	}
}

void mainMenu(int id)
{
	typeMode = id;
	if (typeMode > 5) exit(0);
	//in any case create N control points, depending on application (typeMode)
	createCtrlPoints();
	//for typeMode BEZIER_P0_EQ_P6 make start and end points overlap (P0=P6)
	if (typeMode == BEZIER_P0_EQ_P6)
	{
		ctrlPoints[0][0] = ctrlPoints[N - 1][0];
		ctrlPoints[0][1] = ctrlPoints[N - 1][1];
	}
	//for typeMode BEZIER_DOUBLE make points near intersection co-linear (P2,P3,P4)
	if (typeMode == BEZIER_DOUBLE)
	{
		linearControlPoints(2);
		if ((ctrlPoints[0][0] == ctrlPoints[6][0]) && (ctrlPoints[0][1] == ctrlPoints[6][1]))
		{
			ctrlPoints[6][0] += 100;
			ctrlPoints[6][1] += 100;
		}
	}
	glutPostRedisplay();
}

void createMenu()
{
	glutCreateMenu(mainMenu);
	glutAddMenuEntry("+ Hermite  ", HERMIT);
	glutAddMenuEntry("+ Cubic Interpolating Form  ", CUBIC_INTERPOLATION);
	glutAddMenuEntry("+ 6th Degree Bezier Form  ", BEZIER_P0_EQ_P6);
	glutAddMenuEntry("+ Two Interval Bezier Form  ", BEZIER_DOUBLE);
	glutAddMenuEntry("+ Exit  ", 5);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void init()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ww, 0.0, wh);

	createCtrlPoints();
}
// mouse function
void myPick(int button, int state, int xPosition, int yPosition)
{
	// left mouse button down
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		GLuint newX = xPosition;
		GLuint newY = wh - yPosition;
		// determine which control point is picked
		int choiceFound = 0;
		for (i = 0; i < N && !choiceFound; i++)
		{
			// Use globally defined pickRadius
			if ((abs(ctrlPoints[i][0] - newX) <= pickRadius) &&
				(abs(ctrlPoints[i][1] - newY) <= pickRadius))
			{
				MOVENAME = i;
				choiceFound = 1;
			}
		}
	}
	// left mouse button up
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
		MOVENAME = -1;
	}
	glutPostRedisplay();
}
// mouse motion function
void myMouseMove(int xPosition, int yPosition)
{
	if (MOVENAME > -1)
	{
		GLuint newX = xPosition;
		GLuint newY = wh - yPosition;
		//calculate difference with new mouse position
		float dX = newX - ctrlPoints[MOVENAME][0];
		float dY = newY - ctrlPoints[MOVENAME][1];

		ctrlPoints[MOVENAME][0] = newX;
		ctrlPoints[MOVENAME][1] = newY;
		//make sure start and end points intersect/overlap
		if (typeMode == BEZIER_P0_EQ_P6)
		{
			if (MOVENAME == 0)
			{
				ctrlPoints[6][0] = newX;
				ctrlPoints[6][1] = newY;
			}
			else if (MOVENAME == 6)
			{
				ctrlPoints[0][0] = newX;
				ctrlPoints[0][1] = newY;
			}
		}
		//make sure union points are co-linear
		if ((typeMode == BEZIER_DOUBLE))
		{
			if (MOVENAME == 2)
			{
				linearControlPoints(2);
			}
			else if (MOVENAME == 4)
			{
				linearControlPoints(4);
			}
			else if (MOVENAME == 3)
			{
				ctrlPoints[2][0] += dX; ctrlPoints[4][0] += dX;
				ctrlPoints[2][1] += dY; ctrlPoints[4][1] += dY;
			}
		}

		glutPostRedisplay();
	}
}

void drawCurve()
{
	switch (typeMode) {
	case HERMIT:
		hermite();
		break;
	case CUBIC_INTERPOLATION:
		interpolate();
		break;
	case BEZIER_P0_EQ_P6:
		bezier(0, 6);
		break;
	case BEZIER_DOUBLE:
		doubleBezier();
		break;
	default:
		exit(0);
	}
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);
	drawCurve();
	glFlush();
	glutSwapBuffers();
}

void reshape(GLsizei w, GLsizei h)
{
	/* adjust clipping box */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, (GLdouble)w, 0.0, (GLdouble)h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/* adjust viewport and clear */
	glViewport(0, 0, w, h);
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
	/* set global size for use by drawing routine */
	ww = w;
	wh = h;
}

void credits()
{
	printf("*********************************\n");
	printf("*     Graphics Curve Solver     *\n");
	printf("*********************************\n\n");

	printf("Program implements:\n\n   + Cubic Hermit Curves\n   + Cubic Interpolation Curves\n   + 6th Degree Bezier Curves\n   + Two-segment Bezier Curves with C0, C1 continuity  \n\n");
	printf("Instructions:\n\n     Right click anywhere on the pop-up window to select any of the above\n\n");

	printf("------------------------------------------------------------------------------\n\n");
	printf(" by Vlasis Gogousis & Kiki Hatzistavrou   \n");
	printf(" Aristotle University of Thessaloniki, 2014\n\n");
	printf("------------------------------------------------------------------------------\n");
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(ww, wh);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Grapics Curve Solver");
	glutDisplayFunc(myDisplay);
	glutMouseFunc(myPick);
	glutMotionFunc(myMouseMove);
	glutReshapeFunc(reshape);
	init();

	credits();

	createMenu();
	glutMainLoop();

	system("pause");
	return 0;
}
