//
//  main.cpp
//  Minggu2
//
//  Created by Student on 9/8/16.
//  Copyright (c) 2016 Student. All rights reserved.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <gl/freeglut.h>
#include <math.h>

#include "Object3D.h"


using namespace std;
using namespace cv;


static int tick = 0;
int off;


static int DRAW_TYPE_REALTIME = 0;
static int WINDOW_HEIGHT = 640;
static int WINDOW_WIDTH = 640;
static bool is_capturing = false;
static bool is_sobel_enabled = false;
static bool is_canny_enabled = false;
static bool is_lbp_enabled = false;
static bool is_export_enabled = false;
static bool is_write_enabled = false;

static int global_n = 0;
static int global_lbp_radius = 0;
static int flobal_lbp_neighbors = 0;

static int sudut = 0;
static int sudut2 = 0;
float sudutVector = 0;
bool vectorReturn = false;
bool writeOnce = false;
bool running = false;
string BASE_PATH = "C:/Users/Rizal Fahmi/Google Drive/Tugas/Tugas PENS/Tugas/Semester 5/Grafika Komputer/Projek/Files";
string BASE_OUTPUT_PATH = "C:/Users/Rizal Fahmi/Google Drive/Tugas/Tugas PENS/Tugas/Semester 5/Grafika Komputer/Projek/Output";


typedef struct {
	float v[3];
} vector2D_t;

typedef struct {
	float x;
	float y;
} point2d_t;

typedef struct {
	float x, y, z;
}point3d_t;

typedef struct {
	float r;
	float g;
	float b;
} color_t;

typedef struct {
	float v[3];
}Vector3D_t;

typedef struct {
	float m[3][3];
}matrix3D_t;

typedef struct {
	int numberOfVertices;
	short int pnt[32];
}face_t;

typedef struct {
	int numberOfVertices;
	point3d_t pnt[100];
	int numberofFaces;
	face_t fc[32];
}object3D_t;

typedef struct {
	float x, y, z, r, g, b;
}Point3D_color_t;

typedef struct {
	int numberOfVertices;
	Point3D_color_t* pnt = new Point3D_color_t[80000];
	int numberofFaces;
	face_t* fc = new face_t[80000];
}object3D_color_t;

typedef struct {
	float m[3][3];
} matrix2D_t;



void setColor(color_t col) {
	glColor3f(col.r, col.g, col.b);
}

void drawPoint(float x, float y) {
	glColor3f(0, 0, 1.0);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
}

void drawPolygon(point2d_t pnt[], int n) {
	int i;
	glBegin(GL_LINE_LOOP);
	glColor3f(1.0, 1.0, 1.0);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

void drawPolygon(Point3D_color_t pnt[], int n) {
	int i;
	glBegin(GL_LINE_LOOP);
	for (i = 0; i<n; i++) {
		glColor3f(pnt[i].r, pnt[i].g, pnt[i].b);
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

void createInvisible_color(Point3D_color_t P[]) {
	//glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i<3; i++) {
		P[i].z = 0.0;
		glVertex3f(P[i].x, P[i].y, P[i].z);
	}
	glEnd();
}

void createVisible_color(Point3D_color_t P[]) {
	//interpolasi warna dari tiga titik
	glBegin(GL_POLYGON);
	float r = ((P[0].r + P[1].r + P[2].r) / 3)/255;
	float g = ((P[0].g + P[1].g + P[2].g) / 3)/255;
	float b = ((P[0].b + P[1].b + P[2].b) / 3)/255;
	//printf("r : %f\n", r);
	//printf("g : %f\n", g);
	//printf("b : %f\n\n", b);
	color_t col = { r,g,b };
	setColor(col);
	for (int i = 0; i<3; i++) {
		P[i].z = 0.0;
		glVertex3f(P[i].x, P[i].y, P[i].z);
	}
	glEnd();
}

Mat createVisible_color2DImage(Point3D_color_t P[], Mat image) {
	float r = ((P[0].r + P[1].r + P[2].r) / 3) / 255;
	float g = ((P[0].g + P[1].g + P[2].g) / 3) / 255;
	float b = ((P[0].b + P[1].b + P[2].b) / 3) / 255;

	//Vec3b color(r, g, b);
	Point points[1][3];
	//vector<Point> data;
	for (int i = 0; i<3; i++) {
		//data.push_back(Point(P[i].x, P[i].y));
		points[0][i] = Point(P[i].x, P[i].y);
	}

	const Point* ppt[1] = { points[0] };
	int npt[] = { 3 };
	fillPoly(image, ppt, npt, 1, Scalar(r, g, b));

	return image;

}

void draw3DPolygon(Vector3D_t vec[], int n) {
	int i;
	glBegin(GL_LINE_LOOP);
	glColor3f(1.0, 1.0, 1.0);
	for (i = 0; i<n; i++) {
		glVertex3f(vec[i].v[0], vec[i].v[1], vec[i].v[2]);
	}
	glEnd();
}

void drawSinglePolygon(float x, float y) {
	glBegin(GL_LINE_LOOP);
	glColor3f(1.0, 1.0, 1.0);
	glVertex2f(x, y);
	glEnd();
}

void drawPolyline(point2d_t pnt[], int n)
{
	int i;
	glBegin(GL_LINE_STRIP);
	glColor3f(1.0, 1.0, 1.0);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}
void fillPolygon(point2d_t pnt[], int n,
	color_t color)
{
	int i;
	setColor(color);
	glBegin(GL_POLYGON);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}

	glEnd();
}

void fillPolygonVector(vector2D_t pnt[], int n,
	color_t color)
{
	int i;
	setColor(color);
	glBegin(GL_POLYGON);
	for (i = 0; i<n; i++) {
		glVertex3f(pnt[i].v[0], pnt[i].v[1], pnt[i].v[2]);
	}

	glEnd();
}

void GradatePolygon(point2d_t pnt[], int
	n, color_t color)
{
	int i;
	glBegin(GL_POLYGON);
	for (i = 0; i<n; i++) {
		setColor(color);
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}


void drawUsingMathEquation() {

	point2d_t shape[360];
	double srad, r;
	for (int s = 0; s<360; s++) {
		srad = (s)*3.14 / 180;
		r = sin(6 * srad);
		shape[s].x = (float)(r*cos(srad));
		shape[s].y = (float)(r*sin(srad));
	}
	color_t col = { 1.,1.,0.0 };
	fillPolygon(shape, 360, col);

}

void drawSinFucntion(void) {
	point2d_t p[360];
	for (int i = 0; i<360; i++) {
		p[i].x = (float)i;
		p[i].y = (float)sin((float)i / 57.3);
	}
	float x = sudut;
	float y = (float)sin((float)sudut / 57.3);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
	drawPolyline(p, 360);
}

void drawOsiloscope(void) {
	point2d_t p[360];

	for (int i = 0; i<360; i++) {
		p[i].x = (float)(i);
		p[i].y = (float)sin((float)2 * 180 * 0.1*0.2)*(1 + sudut*sin(2 * 180 * 0.1*0.3));
	}
	drawPolyline(p, 360);
}

void drawMmovingSin(void) {
	point2d_t p[360];

	for (int i = 0; i<360; i++) {
		p[i].x = (float)(i);
		p[i].y = (float)sin((float)(i + sudut) / 15);
	}

	float x = sudut;
	float y = (float)sin((float)(sudut + sudut) / 15);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
	drawPolyline(p, 360);
}

void createClock(point2d_t p, int n, float x2, float y2, float x1, float y1) {
	glBegin(GL_LINE_STRIP);
	glColor3f(1.0, 1.0, 1.0);

	glVertex2f(x2, y2);
	glVertex2f(x1, y1);

	glEnd();
}

void drawCircle(float r, int n) {
	point2d_t p[360];
	float a = 6.28 / n;
	for (int i = 0; i<n; i++) {
		p[i].x = r*(float)cos((float)i*a);
		p[i].y = r*(float)sin((float)i*a);
	}
	color_t col = { 0.0,1.0,0.0 };
	fillPolygon(p, n, col);
	//drawPolygon(p, n);
}



void drawClock(float r, int n) {
	point2d_t p[360];
	float a = 6.28 / n;
	for (int i = 0; i<n; i++) {
		p[i].x = r*(float)cos((float)i*a);
		p[i].y = r*(float)sin((float)i*a);
	}
	float x2;
	float y2;
	drawPolygon(p, n);
	float jari = 0.75;
	for (int i = 0; i<360; i++) {
		double teta = (float)(i*tick)*3.14 / 100;
		p[i].x = (float)(jari*cos(teta));
		p[i].y = (float)(jari*sin(teta));

	}
	//float xPoint = (float) (r*cos(0));
	//float yPoint = (float) (r*sin(0));

	int tick = 0;
	double teta = (float)(sudut / 57.3);
	x2 = (float)(jari*cos(teta));
	y2 = (float)(jari*sin(teta));

	double teta1 = (float)(sudut / 57.3);
	float x1 = (float)(jari*cos(teta1));
	float y1 = (float)(jari*sin(teta1));
	sudut--; if (sudut <= -360.0) sudut = 0.0;
	drawCircle(0.005, 40);
	//createClock(p, 360, x2, y2, x1, y1);
	glBegin(GL_LINES);
	glColor3f(1.0, 1.0, 1.0);

	glVertex2f(x2, y2);
	//glEnd();
	glVertex2f(x1, y1);

	glEnd();

}
void drawMovingCircle2(point2d_t p0, float r, int n) {
	point2d_t p[360];

	float a = 6.28 / n;
	for (int i = 0; i<360; i++) {
		p[i].x = (p0.x + r*(float)cos((float)i*a));
		p[i].y = (p0.y + r*(float)sin((float)i*a));
	}
	color_t col = { 1.0,0.0,1.0 };
	fillPolygon(p, n, col);

}
void drawMovingCircle(point2d_t p0, float r1, float r2, int n) {
	point2d_t p[360];
	point2d_t p2[360];
	float a = 6.28 / n;
	double teta1 = (float)(sudut / 0.01);
	float x1 = (float)(r1*cos(teta1));
	float y1 = (float)(r1*sin(teta1));
	for (int i = 0; i<360; i++) {
		p[i].x = (p0.x + r1*(float)cos((float)i*a));
		p[i].y = (p0.y + r1*(float)sin((float)i*a));
	}
	//sudut--; if(sudut<=-360.0) sudut=0.0;
	color_t col = { 1.0,0.0,1.0 };
	fillPolygon(p, n, col);

	for (int i = 0; i<360; i++) {
		p2[i].x = (p0.x + r2*(float)cos((float)i*a));
		p2[i].y = (p0.y + r2*(float)sin((float)i*a));
	}

	drawPolygon(p2, 40);

	drawMovingCircle2(p2[sudut2], 10, 40);
	sudut2++;
	if (sudut2 == 360) {
		sudut2 = 0;
	}

}



void drawPlanetRevolution(float r, int n) {
	point2d_t p1[360];
	point2d_t p2[360];
	point2d_t p3[360];
	float a = 6.28 / n;
	float r2 = 50;
	for (int i = 0; i<360; i++) {
		p1[i].x = r*(float)cos(((float)i*a));
		p1[i].y = r*(float)sin((float)i*a);


	}
	float x2;
	float y2;
	drawPolygon(p1, n);
	float jari = 40;
	for (int i = 0; i<n; i++) {
		p2[i].x = jari*(float)cos((float)i*a);
		p2[i].y = jari*(float)sin((float)i*a);
	}


	color_t col = { 1.0,1.0,1.0 };
	fillPolygon(p2, n, col);

	drawMovingCircle(p1[sudut], 30, 60, 40);

}

vector2D_t point2Vector(point2d_t pnt) {
	vector2D_t vec;
	vec.v[0] = pnt.x;
	vec.v[1] = pnt.y;
	vec.v[2] = 1.;
	return vec;
}

Vector3D_t point2Vector3D(Point3D_color_t pnt) {
	Vector3D_t vec;
	vec.v[0] = pnt.x;
	vec.v[1] = pnt.y;
	vec.v[2] = pnt.z;

	return vec;
}

Point3D_color_t Vector3D2Point(Vector3D_t vec, color_t col) {
	Point3D_color_t p;
	p.x = vec.v[0];
	p.y = vec.v[1];
	p.z = vec.v[2];
	p.r = col.r;
	p.g = col.g;
	p.b = col.b;
	return p;
}

matrix2D_t createIdentity(void) {
	matrix2D_t u;
	int i, j;
	for (i = 0; i<3; i++) {
		for (j = 0; j<3; j++)
			u.m[i][j] = 0.;
		u.m[i][i] = 1.;
	}
	return u;
}

matrix2D_t translationMTX(float dx, float dy) {
	matrix2D_t trans = createIdentity();
	trans.m[0][2] = dx;
	trans.m[1][2] = dy;
	return trans;
}

matrix2D_t rotationMTX(float theta)
{
	matrix2D_t rotate = createIdentity();
	float cs = cos(theta);
	float sn = sin(theta);
	rotate.m[0][0] = cs; rotate.m[0][1] = -sn;
	rotate.m[1][0] = sn; rotate.m[1][1] = cs;
	return rotate;
}




matrix3D_t createIdentity3D(void) {
	matrix3D_t rotate;
	rotate.m[0][0] = 0.0;
	rotate.m[0][1] = 0.0;
	rotate.m[0][2] = 0.0;
	rotate.m[1][0] = 0.0;
	rotate.m[1][1] = 0.0;
	rotate.m[1][2] = 0.0;
	rotate.m[2][0] = 0.0;
	rotate.m[2][1] = 0.0;
	rotate.m[2][1] = 0.0;
	return rotate;
}

matrix3D_t rotationX(float teta) {
	matrix3D_t rotate = createIdentity3D();
	rotate.m[0][0] = 1.0;
	rotate.m[0][1] = 0.0;
	rotate.m[0][2] = 0.0;
	rotate.m[1][0] = 0.0;
	rotate.m[1][1] = cos(teta / 57.3);
	rotate.m[1][2] = -sin(teta / 57.3);
	rotate.m[2][0] = 0.0;
	rotate.m[2][1] = sin(teta / 57.3);
	rotate.m[2][2] = cos(teta / 57.3);
	return rotate;
}

matrix3D_t rotationY(float teta) {
	matrix3D_t rotate = createIdentity3D();
	rotate.m[0][0] = cos(teta / 57.3);
	rotate.m[0][1] = 0.0;
	rotate.m[0][2] = sin(teta / 57.3);

	rotate.m[1][0] = 0.0;
	rotate.m[1][1] = 1.0;
	rotate.m[1][2] = 0.0;

	rotate.m[2][0] = -sin(teta / 57.3);
	rotate.m[2][1] = 0.0;
	rotate.m[2][2] = cos(teta / 57.3);
	return rotate;
}

matrix3D_t rotationZ(float teta) {
	matrix3D_t rotate = createIdentity3D();
	rotate.m[0][0] = cos(teta / 57.3);
	rotate.m[0][1] = -sin(teta / 57.3);
	rotate.m[0][2] = 0.0;

	rotate.m[1][0] = sin(teta / 57.3);
	rotate.m[1][1] = cos(teta / 57.3);
	rotate.m[1][2] = 0.0;

	rotate.m[2][0] = 0.0;
	rotate.m[2][1] = 0.0;
	rotate.m[2][2] = 1.0;
	return rotate;
}

vector2D_t operator * (matrix2D_t a, vector2D_t b)
{
	vector2D_t c;
	//c=a*b
	int i, j;
	for (i = 0; i<3; i++) {
		c.v[i] = 0;
		for (j = 0; j<3; j++)
			c.v[i] += a.m[i][j] * b.v[j];
	}
	return c;
}

matrix2D_t operator * (matrix2D_t a, matrix2D_t b)
{
	matrix2D_t c;
	//c=a*b
	int i, j, k;
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) {
		c.m[i][j] = 0;
		for (k = 0; k<3; k++)
			c.m[i][j] += a.m[i][k] * b.m[k][j];
	}
	return c;
}

matrix3D_t operator * (matrix3D_t a, matrix3D_t b) {
	matrix3D_t c;
	int i, j, k;
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) {
		c.m[i][j] = 0;
		for (k = 0; k<3; k++)
			c.m[i][j] += a.m[i][k] * b.m[k][j];
	}
	return c;
}


Vector3D_t operator + (Vector3D_t a, Vector3D_t b) {
	Vector3D_t c;
	for (int i = 0; i<3; i++) {
		c.v[i] = a.v[i] + b.v[i];
	}
	return c;
}

Vector3D_t operator - (Vector3D_t a, Vector3D_t b) {
	Vector3D_t c;
	for (int i = 0; i<3; i++) {
		c.v[i] = a.v[i] - b.v[i];
	}

	return c;
}

Vector3D_t operator * (matrix3D_t a, Vector3D_t b) {
	Vector3D_t c;
	for (int i = 0; i<3; i++) {
		c.v[i] = 0;
		for (int j = 0; j<3; j++) {
			c.v[i] += a.m[i][j] * b.v[j];
		}
	}
	return c;
}

Vector3D_t operator ^ (Vector3D_t a, Vector3D_t b) {
	Vector3D_t v;
	v.v[0] = a.v[1] * b.v[2] - a.v[2] * b.v[1];
	v.v[1] = a.v[2] * b.v[0] - a.v[0] * b.v[2];
	v.v[2] = a.v[0] * b.v[1] - a.v[1] * b.v[0];
	return v;
}

vector2D_t setToNeutral(vector2D_t vec, float dx, float dy) {
	vector2D_t v;
	v.v[0] = vec.v[0] - dx;
	v.v[1] = vec.v[1] - dy;
	v.v[2] = vec.v[2];
	return v;
}

vector2D_t setToNormal(vector2D_t vec, float dx, float dy) {
	vector2D_t v;
	v.v[0] = vec.v[0] + dx;
	v.v[1] = vec.v[1] + dy;
	v.v[2] = vec.v[2];
	return v;
}

vector2D_t setRotateToNeutral(vector2D_t vec, float theta) {
	vector2D_t v;
	v.v[0] = vec.v[0] * cos(theta) - vec.v[1] * sin(theta);
	v.v[1] = vec.v[0] * sin(theta) + vec.v[1] * cos(theta);
	v.v[2] = vec.v[2];
	return v;
}

void drawLine(Point3D_color_t pnt1, Point3D_color_t pnt2, color_t col) {
	setColor(col);
	glLineWidth(2);
	glBegin(GL_LINES);
	glVertex3f(pnt1.x, pnt1.y, pnt1.z);
	glVertex3f(pnt2.x, pnt2.y, pnt2.z);
	glEnd();
}

void setBackground(void) {
	color_t warna = { 1.0, 1.0, 1.0 };

	Point3D_color_t kiri = { -640.0, 0, 0 };
	Point3D_color_t kanan = { 640.0, 0, 0 };
	drawLine(kiri, kanan, warna);

	Point3D_color_t atas = { 0, -640, 0 };
	Point3D_color_t bawah = { 0, 640, 0 };
	drawLine(atas, bawah, warna);

	Point3D_color_t pojok1 = { -640.0, -640.0, 0 };
	Point3D_color_t pojok2 = { 640.0, 640.0, 0 };
	drawLine(pojok1, pojok2, warna);
}


void drawBattery() {
	point2d_t rectangle1[4] = { { -75.0,125.0 },{ 75.0,125.0 },{ 75.0,-125.0 },{ -75.0,-125.0 } };
	point2d_t rectangle2[4] = { { -75.0,125.0 },{ 75.0,125.0 },{ 75.0,150.0 },{ -75.0,150.0 } };
	point2d_t rectangle3[4] = { { -75.0,-125.0 },{ 75.0,-125.0 },{ 75.0,-150.0 },{ -75.0,-150.0 } };
	color_t col1 = { 0.0,1.0,1.0 };
	color_t col2 = { 1.0,1.0,1.0 };
	//    fillPolygon(rectangle1, 4, col1);
	//    fillPolygon(rectangle2, 4, col2);
	//drawSinglePolygon(-50, 50);

	vector2D_t rectangle_vec1[4];
	vector2D_t rectangle_vec2[4];
	vector2D_t rectangle_vec3[4];
	matrix2D_t rectangle2TransMTX[4];
	matrix2D_t rectangle2RotMTX[4];
	matrix2D_t rectangle2MTX[4];

	for (int i = 0; i<4; i++) {
		rectangle_vec1[i] = point2Vector(rectangle1[i]);
		rectangle_vec2[i] = setToNeutral(point2Vector(rectangle2[i]), 75., 125.);
		rectangle_vec3[i] = point2Vector(rectangle3[i]);
		rectangle2TransMTX[i] = translationMTX(75., 125.);
		rectangle2RotMTX[i] = rotationMTX(sudutVector);
		rectangle2MTX[i] = operator*(rectangle2TransMTX[i], rectangle2RotMTX[i]);
		rectangle_vec2[i] = setRotateToNeutral(rectangle_vec2[i], sudutVector);
		rectangle_vec2[i] = operator*(rectangle2MTX[i], rectangle_vec2[i]);
		//rectangle_vec2[i] = setToNormal(rectangle_vec2[i], 75., 125.);

	}

	fillPolygonVector(rectangle_vec1, 4, col1);
	fillPolygonVector(rectangle_vec2, 4, col2);
	fillPolygonVector(rectangle_vec3, 4, col2);

}

void timer(int value) {
	sudut+=10;
	if (sudut == 360) {
		sudut = 0;
	}
	glutPostRedisplay();
	glutTimerFunc(300, timer, 0);

}

void timerVector(int value) {
	//sudut2++;
	if (!vectorReturn && sudutVector>-0.9) {
		sudutVector -= 0.02;
	}
	else {
		vectorReturn = true;
	}
	if (vectorReturn && sudutVector<0) {
		sudutVector += 0.02;
	}
	else {
		vectorReturn = false;
	}
	glutPostRedisplay();
	glutTimerFunc(75, timerVector, 0);
}



void drawOff(object3D_color_t ob, float tetha, int type) {

	int i, j;

	Vector3D_t vec[3];
	Point3D_color_t p[6];

	if (type == 0)
		tetha = sudut;


	matrix3D_t matrix_X = rotationX(tetha);
	matrix3D_t matrix_Y = rotationY(tetha);
	matrix3D_t matrix_Z = rotationZ(tetha);

	matrix3D_t matrix_Trans;

	Vector3D_t normalVector;
	Vector3D_t vecbuff[10];


	matrix3D_t tilting = rotationX(tetha)*rotationY(tetha);
	//printf("vertices : %d\n", ob.numberOfVertices);
	//printf("faces : %d\n",ob.numberofFaces);

	for (i = 0; i<ob.numberofFaces; i++) {
		//if (!writeOnce)
		//	printf("face ke-%d, vert:%d   ==>\n", i, ob.fc[i].numberOfVertices);
		for (j = 0; j<ob.fc[i].numberOfVertices; j++) {
		//	if (!writeOnce)
		//		printf(" %d : %f,%f,%f, color : %f,%f,%f\n", ob.fc[i].pnt[j], ob.pnt[ob.fc[i].pnt[j]].x, ob.pnt[ob.fc[i].pnt[j]].y, ob.pnt[ob.fc[i].pnt[j]].z,
		//			ob.pnt[ob.fc[i].pnt[j]].r, ob.pnt[ob.fc[i].pnt[j]].g, ob.pnt[ob.fc[i].pnt[j]].b);
			p[j].x = ob.pnt[ob.fc[i].pnt[j]].x;
			p[j].y = ob.pnt[ob.fc[i].pnt[j]].y;
			p[j].z = ob.pnt[ob.fc[i].pnt[j]].z;
			p[j].r = ob.pnt[ob.fc[i].pnt[j]].r;
			p[j].g = ob.pnt[ob.fc[i].pnt[j]].g;
			p[j].b = ob.pnt[ob.fc[i].pnt[j]].b;
			vec[j] = point2Vector3D(p[j]);

			vec[j] = operator*(tilting, vec[j]);
			color_t c = { p[j].r,p[j].g,p[j].b };
			p[j] = Vector3D2Point(vec[j], c);
		}
		normalVector = (vec[1] - vec[0]) ^ (vec[0] - vec[2]);

		Point3D_color_t pNormal[3];
		if (normalVector.v[2]<0) {
			int p0 = ob.fc[i].pnt[0];
			int p1 = ob.fc[i].pnt[1];
			int p2 = ob.fc[i].pnt[2];
			pNormal[0] = p[0];
			pNormal[1] = p[1];
			pNormal[2] = p[2];
			//createInvisible_color(pNormal);
		}
		if (normalVector.v[2]>0) {
			int p0 = ob.fc[i].pnt[0];
			int p1 = ob.fc[i].pnt[1];
			int p2 = ob.fc[i].pnt[2];
			pNormal[0] = p[0];
			pNormal[1] = p[1];
			pNormal[2] = p[2];
			createVisible_color(pNormal);
		}
		printf("\n");
	}


	//writeOnce = true;

}






object3D_color_t readOffFile() {


	ifstream ifs1(BASE_PATH + "/off14/"+to_string(off)+".off");

	if (!ifs1) {
		std::cout << "cannot open file_readoff";
		exit(1);
	}

	char str[255];
	char *tok;
	char delims[] = " ";
	int x, y = 0;
	int numVertices = 0;
	int numFaces = 0;
	int num[10];
	float numTof[10];
	int pntIndices = 0;
	int fcIndices = 0;
	object3D_color_t ob;
	char *next_token1 = NULL;

	while (!ifs1.eof() && (y < 2 + numVertices + numFaces)) {
		ifs1.getline(str, 255);
		if (y == 1) {//mendapatkan nilai numVertices dan numFaces
			tok = strtok_s(str, delims, &next_token1);
			x = 0;
			while (tok && x < 2) {
				num[x] = atof(tok);
				tok = strtok_s(NULL, delims, &next_token1);
				x++;
			}
			numVertices = (int)num[0];
			numFaces = (int)num[1];
			ob.numberOfVertices = (int)num[0];
			ob.numberofFaces = (int)num[1];
			//printf("number of vertices = %d\n", numVertices);
			//printf("number of faces = %d\n", numFaces);

		}
		else if (y > 1 && y<2 + numVertices) {
			tok = strtok_s(str, delims, &next_token1);

			x = 0;
			while (tok && x < 6) {
				numTof[x] = atof(tok);
				tok = strtok_s(NULL, delims, &next_token1);
				x++;
			}
			ob.pnt[pntIndices] = { (float)numTof[0],(float)numTof[1],(float)numTof[2],(float)numTof[3],(float)numTof[4],(float)numTof[5] };
			//printf("index : %d\n", pntIndices);
			//printf("vertices x = %f\n", ob.pnt[pntIndices].x);
			pntIndices++;

		}
		else if (y > (numVertices + 1 + 2)) {
			tok = strtok_s(str, delims, &next_token1);
			x = 0;
			while (tok && x < 4) {
				num[x] = atof(tok);
				tok = strtok_s(NULL, delims, &next_token1);
				x++;
			}
			ob.fc[fcIndices].numberOfVertices = (int)num[0];
			for (int i = 0; i < ob.fc[fcIndices].numberOfVertices; i++)
				ob.fc[fcIndices].pnt[i] = (int)num[i + 1];
			//printf("num ver : %f\n", ob.pnt[113].x);
			fcIndices++;
		}
		y++;
	}

	//drawOff(ob);

	return ob;
	
}


void writePNGFromExport(Mat image, int n, string fileName) {
	int whichFile = off;
	for (int i = 0; i < n; i++) {
		cout << "Writing png into file " + fileName + "_" + to_string(i + 1) + ".png" "\n";
		imwrite(BASE_OUTPUT_PATH + fileName + "_" + to_string(i + 1) + ".png", image);
	}
}

Mat getSobel(Mat originalMat) {


	Mat grey;

	cvtColor(originalMat, grey, CV_BGR2GRAY);

	Mat sobelX;
	Sobel(grey, sobelX, CV_32F, 1, 0);

	double minVal, maxVal;
	minMaxLoc(sobelX, &minVal, &maxVal);

	Mat drawSobelMat;
	sobelX.convertTo(drawSobelMat, CV_8U, 255.0 / (maxVal - minVal), -minVal*255.0 / (maxVal - minVal));


	return drawSobelMat;

}

Mat getCanny(Mat originalMat) {


	int edgeThresh = 1;
	int lowThreshold = 0;
	int const max_lowThreshold = 100;
	int ratio = 3;
	int kernel_size = 3;
	char* windows_name = "3D Object By Rizal Fahmi";

	Mat original_gray;
	Mat dst, detected_edges;

	dst.create(originalMat.size(), originalMat.type());
	cvtColor(originalMat, original_gray, CV_BGR2GRAY);

	blur(original_gray, detected_edges, Size(3, 3));

	Canny(detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size);

	dst = Scalar::all(0);

	originalMat.copyTo(dst, detected_edges);

	/*for (int i = 0; i < n; i++)
	imwrite(BASE_PATH + "/Canny/" + "canny_" + to_string((i + 1)) + ".png", dst);*/

	//Canny(detected_edges,detected_edges)

	return dst;

}

Mat captureWindow() {
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);

	Mat temp = Mat::zeros(height, width, CV_8UC3);
	Mat tempImage;


	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, temp.data);

	cvtColor(temp, tempImage, CV_BGR2RGB);

	flip(tempImage, temp, 0);

	/*for (int i = 0; i < global_n; i++) {
		string fileName = "/2D PNG/" + to_string(off) + "_off_" + to_string(i + 1) + ".png";
		cout << "Writing png into file " + fileName + "\n";
		imwrite(BASE_OUTPUT_PATH + fileName, temp);
	}*/

	return temp;
}


void writeFeature(Mat dst, string path) {
	//ofstream fout(BASE_OUTPUT_PATH + "/LBP/feature/features.txt");
	ofstream fout(BASE_OUTPUT_PATH +  path);
	int unsigned *ptr;
	ptr = (unsigned int *)dst.data;
	cout << "Writing features into file " + path << endl;
	for (int i = 0; i < (dst.rows); i++) {
		for (int j = 0; j < (dst.cols); j++) {
			fout << *ptr << " ";
			ptr++;
		}
		fout << endl;
	}

	fout.close();
}

template <typename _Tp>
Mat OLBP_(Mat originalMat, int radius, int neighbors, int n) {

	Mat dst;

	neighbors = max(min(neighbors, 31), 1);
	dst = Mat::zeros(originalMat.rows - 2 * radius, originalMat.cols - 2 * radius, CV_32SC1);

	for (int n = 0; n < neighbors; n++) {
		// sample points
		float x = static_cast<float>(radius) * cos(2.0*M_PI*n / static_cast<float>(neighbors));
		float y = static_cast<float>(radius) * -sin(2.0*M_PI*n / static_cast<float>(neighbors));

		// relative indices
		int fx = static_cast<int>(floor(x));
		int fy = static_cast<int>(floor(y));
		int cx = static_cast<int>(ceil(x));
		int cy = static_cast<int>(ceil(y));

		// fractional part
		float ty = y - fy;
		float tx = x - fx;

		// set interpolation weights
		float w1 = (1 - tx) * (1 - ty);
		float w2 = tx *(1 - ty);
		float w3 = (1 - tx) * ty;
		float w4 = tx * ty;

		// iterate obtained data
		for (int i = radius; i < originalMat.rows - radius; i++) {
			for (int j = radius; j < originalMat.cols - radius; j++) {
				float t = w1*originalMat.at<_Tp>(i + fy, j + fx) + w2*originalMat.at<_Tp>(i + fy, j + cx) + w3*originalMat.at<_Tp>(i + cy, j + fx) + w4*originalMat.at<_Tp>(i + cy, j + cx);
				// add tolerance to destinated mat, floating point need precise value
				dst.at<unsigned int>(i - radius, j - radius) += ((t>originalMat.at<_Tp>(i, j)) && (abs(t - originalMat.at<_Tp>(i, j)) > numeric_limits<float>::epsilon())) << n;

			}
		}


		

	}

	return dst;

}

void userdraw(void) {
	//static int tick=0;
	//point2d_t bintang[10] = {{1.1284, 175.72847},{61.04209, 57.358},{208.57805, 47.92467},{91.88762, -37.27787},{154.86341, -165.08167},{-0.72383, -102.10588},{-148.90215, -165.08167},{-87.77859, -35.42564},{-221.13908, 44.22021},{-72.96076, 59.03804}};

	//drawPolygon(bintang, 10);
	color_t col = { 0.0,0.1,0.0 };
	//fillPolygon(bintang, 10, col);
	//GradatePolygon(bintang, 10, col);

	//drawSinFucntion();
	//drawCircle(1., 40);
	//drawUsingMathEquation(tick);
	//drawClock(1., 40);
	//tick=tick+1;
	//drawMmovingSin();
	//drawOsiloscope();
	//drawPlanetRevolution(200,40);
	//drawBattery();
	//drawPrism();

	Mat capturedMat;


	drawOff(readOffFile(), 0.0, DRAW_TYPE_REALTIME);

	if (is_capturing) {
		string fileName = "";
		string fileFeatureName = "";
		capturedMat = captureWindow();
		if (is_export_enabled) {
			fileName = "/2D PNG/" + to_string(off) + "_off_";
		}
		if (is_sobel_enabled) {
			capturedMat = getSobel(capturedMat);
			fileName = "/Sobel/" + to_string(off) + "_off_sobel";
		}
		if (is_canny_enabled) {
			capturedMat = getCanny(capturedMat);
			fileName = "/Canny/" + to_string(off) + "_off_canny";
		}
		if (is_lbp_enabled) {
			capturedMat = OLBP_<float>(capturedMat, global_lbp_radius, flobal_lbp_neighbors, global_n);
			fileName = "/LBP/" + to_string(off) + "_off_lbp";
			fileFeatureName = "/LBP/feature/features.txt";
			writeFeature(capturedMat, fileFeatureName);
		}

		writePNGFromExport(capturedMat, global_n, fileName);

		is_capturing = false;
		is_export_enabled = false;
		is_sobel_enabled = false;
		is_canny_enabled = false;
		is_lbp_enabled = false;
	}
		

}

void keyboardInput(unsigned char key, int x, int y) {
	if (key == 13) {
		//printf("return is pressed");
		glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
		is_capturing = false;
		is_export_enabled = false;
		is_sobel_enabled = false;
		is_canny_enabled = false;
		is_lbp_enabled = false;
		glutLeaveMainLoop();
	}
}

void display(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	userdraw();
	glutSwapBuffers();
	glutKeyboardFunc(keyboardInput);
}

void initOpenGLWindow() {

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow("3D Object By Rizal Fahmi");
	//glClearColor(1.0, 1.0, 1.0, 1.0);
	//gluOrtho2D(-320., 320., -240., 240.);
	//gluOrtho2D(0., 360., -1.1, 1.1);
	//gluOrtho2D(-300.0, 300.0, -300.0, 300.0);
	glClearColor(0, 0, 0, 0);
	gluOrtho2D(-2.0, 2.0, -2.0, 2.0);
	//gluOrtho2D(-180., 180., -135., 135.);
	glutTimerFunc(1, timer, 0);

	//glutTimerFunc(1, timerVector, 0);
	glutIdleFunc(display);
	glutDisplayFunc(display);

	glutMainLoop();

}



void openOpenGLAndcaptureWindow() {

	is_capturing = true;
	is_export_enabled = true;

	initOpenGLWindow();

}


/*
	Export file onto png format
*/

Mat exportFile(float tetha) {
	int whichFile = off;
	object3D_color_t ob;
	ob = readOffFile();

	vector<Point2f> points;

	Vector3D_t vec[3];
	Point3D_color_t p[6];

	matrix3D_t matrix_X = rotationX(tetha);
	matrix3D_t matrix_Y = rotationY(tetha);
	matrix3D_t matrix_Z = rotationZ(tetha);

	matrix3D_t matrix_Trans;

	matrix3D_t tilting = rotationX(tetha)*rotationY(tetha);

	Vector3D_t normalVector;
	Vector3D_t vecbuff[10];
	Mat image = Mat::zeros(WINDOW_WIDTH, WINDOW_HEIGHT, CV_8UC3);;

	for(int i = 0; i < ob.numberofFaces; i++) {
		for (int j = 0; j<ob.fc[i].numberOfVertices; j++) {

			p[j].x = ob.pnt[ob.fc[i].pnt[j]].x;
			p[j].y = ob.pnt[ob.fc[i].pnt[j]].y;
			p[j].z = ob.pnt[ob.fc[i].pnt[j]].z;
			p[j].r = ob.pnt[ob.fc[i].pnt[j]].r;
			p[j].g = ob.pnt[ob.fc[i].pnt[j]].g;
			p[j].b = ob.pnt[ob.fc[i].pnt[j]].b;
			vec[j] = point2Vector3D(p[j]);

			vec[j] = operator*(tilting, vec[j]);
			color_t c = { p[j].r,p[j].g,p[j].b };
			p[j] = Vector3D2Point(vec[j], c);
		}
		normalVector = (vec[1] - vec[0]) ^ (vec[0] - vec[2]);

		Point3D_color_t pNormal[3];
		if (normalVector.v[2]<0) {
			int p0 = ob.fc[i].pnt[0];
			int p1 = ob.fc[i].pnt[1];
			int p2 = ob.fc[i].pnt[2];
			pNormal[0] = p[0];
			pNormal[1] = p[1];
			pNormal[2] = p[2];
			//createInvisible_color(pNormal);
		}
		if (normalVector.v[2]>0) {
			int p0 = ob.fc[i].pnt[0];
			int p1 = ob.fc[i].pnt[1];
			int p2 = ob.fc[i].pnt[2];
			pNormal[0] = p[0];
			pNormal[1] = p[1];
			pNormal[2] = p[2];
			//image = createVisible_color2DImage(pNormal, image);
			float r = ((pNormal[0].r + pNormal[1].r + pNormal[2].r) / 3);
			float g = ((pNormal[0].g + pNormal[1].g + pNormal[2].g) / 3);
			float b = ((pNormal[0].b + pNormal[1].b + pNormal[2].b) / 3);

			Vec3b color(r, g, b);
			Point points[1][3];
			//vector<Point> data;
			for (int i = 0; i<3; i++) {
				//data.push_back(Point(P[i].x, P[i].y));
				points[0][i] = Point(pNormal[i].x*(-280)+300, pNormal[i].y*280+300);
				//float x = (float)(pNormal[i].x*(-280)+300);
				//float y = (float)(pNormal[i].y*280+300);
				//image.at<Vec3b>(Point(x, y)) = color;
				//printf("point %d : %f, %f\n", i, x, y);
			}

			const Point* ppt[1] = { points[0] };
			int npt[] = { 3 };
			fillPoly(image, ppt, npt, 1, Scalar(r, g, b));
			
		}
	}

	/*Vec3b color(3, 54, 120);

	for (int i = 0; i < 640; i++) {
		for (int j = 0; j < 640; j++)
		{
			image.at<Vec3b>(Point(i, j)) = color;
		}
	}*/

	return image;
}

void processAllOff(float tetha, int radius, int neighbors, int n) {
	ifstream ifs2(BASE_PATH + "/read_off.txt");
	
	if (!ifs2) {
		cout << "Cannot open file " + BASE_PATH + "/read_off.txt";
		exit(1);
	}

	Mat image, imageSobel, imageCanny, imageLBP;
	char str[255];
	char *next_token = NULL;
	char delims[] = {'/','.'};
	char *tok;
	char str_off[255];
	int x = 0;
	string fileNameExport, fileNameSobel, fileNameCanny, fileNameLBP, fileFeatureName;

	while (!ifs2.eof()) {
		ifs2.getline(str, 255);
		tok = strtok_s(str, delims, &next_token);
		x = 0;
		while (tok && x < 2) {
			if (x == 1) {
				strcpy(str_off, tok);
				off = atoi(str_off);
				cout << "Processing file " + to_string(off) + ".off" << endl;
				image = exportFile(tetha);
				fileNameExport = "/All Off/2D PNG/" + to_string(off) + "_off";
				writePNGFromExport(image, n, fileNameExport);
				imageSobel = getSobel(image);
				fileNameSobel = "/All Off/Sobel/" + to_string(off) + "_off_sobel";
				writePNGFromExport(imageSobel, n, fileNameSobel);
				imageCanny = getCanny(image);
				fileNameCanny = "/All Off/Canny/" + to_string(off) + "_off_canny";
				writePNGFromExport(imageCanny, n, fileNameCanny);
				imageLBP = OLBP_<float>(image, radius, neighbors, n);
				fileNameLBP = "/All Off/LBP/" + to_string(off) + "_off_lbp";
				writePNGFromExport(imageLBP, n, fileNameLBP);
				fileFeatureName = "/All Off/LBP/features/features_" + to_string(off) + "_off.txt";
				writeFeature(imageLBP, fileFeatureName);
				
			}
			tok = strtok_s(NULL, delims, &next_token);
			x++;
		}
	}
}



void readOff(){

	int numOfSobel = 0;
	int numOfCanny = 0;
	int numOfLBP = 0;
	string strtetha;

	int pil=0;

	float tetha = 0;

	while (pil != 7) {

		cout << "Pilih salah satu menu untuk memulai\n";
		cout << "=======================================================\n";

		cout << "1. Tampilkan object 3D dari file off\n";
		cout << "2. Eksport object 3D ke gambar PNG\n";
		cout << "3. Sobel\n";
		cout << "4. Canny\n";
		cout << "5. LBP\n";
		cout << "6. Proses semua file off\n";
		cout << "7. Keluar\n";
		cout << "Masukkan pilihan : ";

		cin >> pil;

		if (pil == 1) {
			cout << "Masukkan file off ke berapa yang ingin ditampilkan : ";
			cin >> off;
			initOpenGLWindow();
		}

		if (pil == 2) {

			float tetha = 0.0;
			int n = 1;
			int metode = 0;

			cout << "Masukkan file off ke berapa yang ingin dieksport : ";
			cin >> off;

			cout << "Masukkan metode : ";
			cin >> metode;

			if (metode == 1) {
				cout << "Masukkan tetha : ";
				cin >> tetha;

				cout << "Masukkan jumlah copy gambar : ";
				cin >> n;
				Mat image = exportFile(tetha);
				string fileName = "/2D PNG/" + to_string(off) + "_off_";
				writePNGFromExport(image, n, fileName);
			}
			else if (metode == 2) {
				cout << "Masukkan jumlah copy gambar : ";
				cin >> global_n;
				openOpenGLAndcaptureWindow();
			}
			
		}

		if (pil == 3) {
			float tetha = 0.0;
			int n = 1;
			int metode = 0;

			cout << "Masukkan file off ke berapa yang ingin dieksport dengan Sobel : ";
			cin >> off;

			cout << "Masukkan metode konversi : ";
			cin >> metode;

			if (metode == 1) {
				cout << "Masukkan tetha : ";
				cin >> tetha;

				cout << "Masukkan jumlah copy gambar : ";
				cin >> n;
				Mat image = exportFile(tetha);
				image = getSobel(image);
				string fileName = "/Sobel/" + to_string(off) + "_off_sobel";
				writePNGFromExport(image, n, fileName);

			}
			else if (metode == 2) {
				cout << "Masukkan jumlah copy gambar : ";
				cin >> global_n;

				is_sobel_enabled = true;

				openOpenGLAndcaptureWindow();
			}

		}

		if (pil == 4) {
			float tetha = 0.0;
			int n = 1;
			int metode = 0;

			cout << "Masukkan file off ke berapa yang ingin dieksport dengan Canny : ";
			cin >> off;

			cout << "Masukkan metode konversi : ";
			cin >> metode;

			if (metode == 1) {
				cout << "Masukkan tetha : ";
				cin >> tetha;

				cout << "Masukkan jumlah copy gambar : ";
				cin >> n;
				Mat image = exportFile(tetha);
				image = getCanny(image);
				string fileName = "/Canny/" + to_string(off) + "_off_canny";
				writePNGFromExport(image, n, fileName);

			}
			else if (metode == 2) {
				cout << "Masukkan jumlah copy gambar : ";
				cin >> global_n;

				is_canny_enabled = true;

				openOpenGLAndcaptureWindow();
			}

		}

		if (pil == 5) {
			float tetha = 0.0;
			int n = 1;
			int metode = 0;
			int radius = 0;
			int neighbors = 0;

			cout << "Masukkan file off ke berapa yang ingin dieksport dengan LBP : ";
			cin >> off;


			cout << "Masukkan metode konversi : ";
			cin >> metode;


			if (metode == 1) {
				cout << "Masukkan tetha : ";
				cin >> tetha;

				cout << "Masukkan jumlah copy gambar : ";
				cin >> n;

				cout << "Masukkan radius : ";
				cin >> radius;

				cout << "Masukkan neighbors : ";
				cin >> neighbors;


				Mat image = exportFile(tetha);
				image = OLBP_<float>(image, radius, neighbors, n);
				string fileName = "/LBP/" + to_string(off) + "_off_lbp";
				writePNGFromExport(image, n, fileName);
				string fileFeatureName = "/LBP/feature/features.txt";
				writeFeature(image, fileFeatureName);

			}
			else if (metode == 2) {
				cout << "Masukkan jumlah copy gambar : ";
				cin >> global_n;


				cout << "Masukkan radius : ";
				cin >> global_lbp_radius;

				cout << "Masukkan neighbors : ";
				cin >> flobal_lbp_neighbors;

				is_lbp_enabled = true;

				openOpenGLAndcaptureWindow();
			}

		}

		if (pil == 6) {
			float tetha = 0.0;
			int n = 1;
			int metode = 0;
			int radius = 0;
			int neighbors = 0;

			/*cout << "Masukkan metode konversi : ";
			cin >> metode;


			if (metode == 1) {*/
				cout << "Masukkan tetha : ";
				cin >> tetha;

				cout << "Masukkan jumlah copy gambar : ";
				cin >> n;

				cout << "Masukkan radius untuk LBP : ";
				cin >> radius;

				cout << "Masukkan neighbors LBP : ";
				cin >> neighbors;


				processAllOff(tetha,radius,neighbors,n);

			/*}
			else if (metode == 2) {
				cout << "Masukkan jumlah copy gambar : ";
				cin >> global_n;


				cout << "Masukkan radius untuk LBP : ";
				cin >> global_lbp_radius;

				cout << "Masukkan neighbors untuk LBP : ";
				cin >> flobal_lbp_neighbors;

				is_lbp_enabled = true;

				openOpenGLAndcaptureWindow();
			}*/

		}

		if (pil == 7) {
			exit(0);
		}



		//setBackground();


		/*getSobel(numOfSobel);
		getCanny(numOfCanny);
		OLBP_<float>(1, 8, numOfLBP);*/
	}
}



//void drawPrism(void) {
//	int i, j;
//	object3D_color_t prism = { 5,{ { 0,100,0,0.5,0.6,0.7 },{ 100,0,0,0.2,0.4,0.6 },{ 0,0,100,0.8,1.0,0.3 },{ -100,0,0,0.6,0.2,0.9 },{ 0,0,-100,0.1,1.0,0.7 } },
//		6,
//		{
//			{ 3,{ 0,2,1 } },{ 3,{ 0,3,2 } },
//			{ 3,{ 0,4,3 } },{ 3,{ 0,1,4 } },
//			{ 3,{ 1,2,4 } },{ 3,{ 3,4,2 } },
//		}
//	};
//	Vector3D_t vec[3];
//	Point3D_color_t p[6];
//
//	float teta = 40.0;
//
//	matrix3D_t matrix_X = rotationX(teta);
//	matrix3D_t matrix_Y = rotationY(teta);
//	matrix3D_t matrix_Z = rotationZ(teta);
//
//	matrix3D_t matrix_Trans;
//
//	Vector3D_t normalVector;
//	Vector3D_t vecbuff[3];
//
//
//	matrix3D_t tilting = rotationX(teta)*rotationY(teta);
//
//	for (i = 0; i<prism.numberofFaces; i++) {
//	//	if (!writeOnce)
//	//		printf("face ke-%d, vert:%d   ==>\n", i, prism.fc[i].numberOfVertices);
//		for (j = 0; j<prism.fc[i].numberOfVertices; j++) {
//	//		if (!writeOnce)
//	//			printf(" %d : %0.f,%0.f,%0.f, color : %f,%f,%f\n", prism.fc[i].pnt[j], prism.pnt[prism.fc[i].pnt[j]].x, prism.pnt[prism.fc[i].pnt[j]].y, prism.pnt[prism.fc[i].pnt[j]].z,
//	//				prism.pnt[prism.fc[i].pnt[j]].r, prism.pnt[prism.fc[i].pnt[j]].g, prism.pnt[prism.fc[i].pnt[j]].b);
//			//vec[i] = point2Vector3D(prism.pnt[prism.fc[i].pnt[j]]);
//			p[j].x = prism.pnt[prism.fc[i].pnt[j]].x;
//			p[j].y = prism.pnt[prism.fc[i].pnt[j]].y;
//			p[j].z = prism.pnt[prism.fc[i].pnt[j]].z;
//			p[j].r = prism.pnt[prism.fc[i].pnt[j]].r;
//			p[j].g = prism.pnt[prism.fc[i].pnt[j]].g;
//			p[j].b = prism.pnt[prism.fc[i].pnt[j]].b;
//			vec[j] = point2Vector3D(p[j]);
//
//			vec[j] = operator*(tilting, vec[j]);
//			color_t c = { p[j].r,p[j].g,p[j].b };
//			p[j] = Vector3D2Point(vec[j], c);
//		}
//		normalVector = (vec[1] - vec[0]) ^ (vec[0] - vec[2]);
//
//		Point3D_color_t pNormal[3];
//		if (normalVector.v[2]<0) {
//			int p0 = prism.fc[i].pnt[0];
//			int p1 = prism.fc[i].pnt[1];
//			int p2 = prism.fc[i].pnt[2];
//			pNormal[0] = p[0];
//			pNormal[1] = p[1];
//			pNormal[2] = p[2];
//			createInvisible_color(pNormal);
//		}
//		if (normalVector.v[2]>0) {
//			int p0 = prism.fc[i].pnt[0];
//			int p1 = prism.fc[i].pnt[1];
//			int p2 = prism.fc[i].pnt[2];
//			pNormal[0] = p[0];
//			pNormal[1] = p[1];
//			pNormal[2] = p[2];
//			createVisible_color(pNormal);
//		}
//		//printf("Normal Vector ke %d = %f\n",i);
//
//		// drawPolygon(p, 3);
//
//		printf("\n");
//	}
//
//	writeOnce = true;
//
//}









int main(int argc, char ** argv)
{
	
	glutInit(&argc, argv);

	readOff();

	waitKey(0);
	//glutMainLoop();
	return 0;
}

