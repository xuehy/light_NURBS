// NURBS tangent and point inversion
#include "../source/nurbs.h"
#include <GL/glut.h>
#include <GL/gl.h>


cv::Point2f P[20];
cv::Point3f Pw[20];
double U[] = {0,0,0,0,0.1,0.5,1,1,1,1,1};//{0,0.5/9.0,1.5/9.0,2.5/9.0,3.5/9.0,4.5/9.0,5.5/9.0,6.5/9.0,7.5/9.0,1};
cv::Point2f Q[20];
double u[10];
int n = 6;
int p = 3;
int N = 12;
void myInit(void)
{
	glClearColor(0.0,0.0,0.0,0.0);//设置背景色
	
	glEnable(GL_POINT_SMOOTH);  
 glEnable(GL_LINE_SMOOTH);  
 glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points  
 glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines  

	Q[0] = cv::Point2f(0,0);
	Q[1] = cv::Point2f(1,1);
	Q[2] = cv::Point2f(2,0.5);
	Q[3] = cv::Point2f(3,-1);

	Q[4] = cv::Point2f(4,-0.8);
	Q[5] = cv::Point2f(5.5,0);
	Q[6] = cv::Point2f(6,0.5);
	Q[7] = cv::Point2f(7,0.9);
	Q[8] = cv::Point2f(8,1);
	Q[9] = cv::Point2f(9,0.8);
	Q[10] = cv::Point2f(10,-0.8);
	Q[11] = cv::Point2f(11,-1.9);
	Q[12] = cv::Point2f(14,-9);
	nurbs::chord_length_param(N,Q,u);	
	
	
	u[0] = 0.0;u[N-1] = 1;
	for(int i = 0; i < N-1; i++) cout<<u[i]<<endl;
	for(int i = 0; i < n + p + 2; i++) cout<<"U"<<U[i]<<endl;
	//int span = nurbs::FindSpan(5,3,u[2],U);
	//cout<<span<<endl;
	nurbs::NurbsCurveFit(n,p,N-1 ,U,u,Q,P);
	for(int i= 0; i < n+1; i++) 
	{
		Pw[i].x = P[i].x; Pw[i].y = P[i].y; Pw[i].z = 1;
		cout<<P[i]<<endl;
	}
}
void mykey(unsigned char key,int x,int y)
{
	
	 /*if(key == 'a')
	 {
		glClear(GL_COLOR_BUFFER_BIT);
		if( glutGetModifiers() == GLUT_ACTIVE_ALT)
		{
			X -= 0.1;
		}
		else 
			X -= 0.01;
		glutPostRedisplay();
	 }
	 else if(key == 'd')
	 {
		glClear(GL_COLOR_BUFFER_BIT);
		if( glutGetModifiers() == GLUT_ACTIVE_ALT)
			X += 0.1;
		else
			X += 0.01;
		glutPostRedisplay();
	 }
	 else if(key == 'w')
	 {
		glClear(GL_COLOR_BUFFER_BIT);
		if( glutGetModifiers() == GLUT_ACTIVE_ALT)
			Y += 0.1;
		else
			Y += 0.01;
		glutPostRedisplay();
	 }
	 else if(key == 's')
	 {
		glClear(GL_COLOR_BUFFER_BIT);
		if( glutGetModifiers() == GLUT_ACTIVE_ALT)
			Y -= 0.1;
		else 
			Y -= 0.01;
		glutPostRedisplay();
	 }*/
}
void myDisplay(void)
{
	int i;
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glColor3f(1.0,1.0,1.0);
	glLineWidth(2.0);
	
	cv::Point2f C;
	glBegin(GL_LINE_STRIP);
		for(double i = 0.0; i < 1.0; i+=0.01)
		{
			nurbs::CurvePoint(n,p,U,Pw,i,&C);
			glVertex2f(C.x,C.y);
		}
	glEnd();
	glColor3f(1.0,1.0,0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	for(i = 0;i < N;i++)
		glVertex2f(Q[i].x,Q[i].y);
	glEnd();
	glColor3f(1.0,0.0,0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	for(i = 0;i < n + 1;i++)
		glVertex2f(P[i].x,P[i].y);
	glEnd();
	glutSwapBuffers();
}
void myReshape(GLsizei w,GLsizei h)
{
	glViewport(0,0,w,h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if(w <=h)
		glOrtho(-5.0,10.0,-5.0*(GLfloat)h/(GLfloat)w,10.0*(GLfloat)h/

(GLfloat)w,-5.0,10.0);
	else
		glOrtho(-5.0*(GLfloat)w/(GLfloat)h,10.0*(GLfloat)w/(GLfloat)h,-

5.0,10.0,-5.0,10.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0,0.0,-4.0);
}
int main(int argc,char ** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowSize(600,400);
	glutInitWindowPosition(200,200);
	glutCreateWindow("NURBS curve");
	/*绘制与显示*/
	myInit();
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	 glutKeyboardFunc(mykey);
	glutMainLoop();
	return(0);
}