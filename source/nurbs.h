#ifndef __NURBS_H
#define __NURBS_H
#include <opencv2/opencv.hpp>
#include <iostream>
#include <memory>
using namespace std;

namespace cv
{
class Point4f
{
public:
    double w,x,y,z;
    Point4f()
    {
        w=x=y=z=0;
    }
    Point4f(double a, double b, double c, double d);
    Point4f(const Point4f &P);
    friend Point4f operator +(const Point4f &P1, const Point4f &P2);
    friend Point4f operator -(const Point4f &P1, const Point4f &P2);
    friend Point4f operator *(const Point4f &P,  double a);
    friend Point4f operator *(double a, const Point4f &P);
    friend Point4f operator /(const Point4f &P,  double a);
    Point4f &operator = (const Point4f &);
    friend ostream& operator<<(ostream& os,Point4f& P)
    {
        os<<"["<<P.w<<", "<<P.x<<", "<<P.y<<", "<<P.z<<"]";
        return os;
    }
};
}

namespace nurbs
{
/*Functions for manipulating Bezier curve*/
void Horner1(double *a, int n, double u0, double *C);
void Bernstein(int i, int n, double u, double *B); 
void AllBernstein(int n, double u, double *B);
void PointOnBezierCurve(cv::Point2f *P, int n, double u, cv::Point2f *C);
void deCasteljau1(cv::Point2f *P, int n, double u, cv::Point2f *C);
void deCasteljau1(cv::Point3f *P, int n, double u, cv::Point3f *C);

/*Functions for NURBS curve and surface*/
//input argument u, find its knot span index
//U is the array of knots, n is number of control points - 1, p is the degree of B-spline curve
int FindSpan(int n, int p, double u, double *U);
//Generate p + 1 basis B-spline funcs in N[] given knot span index i, argument u, degree p
void BasisFuns(int i, double u, int p, double *U, double *N);
// calculate derivatives of B-Spline basis functions
void DersBasisFuns(int i, double u, int p, int n, double *U, vector <vector <double> > &ders);
//find Point on B-spline curve given argument u. P[] is array of control points, C is the result
void CurvePoint(int n, int p, double *U, cv::Point2f *P, double u, cv::Point2f *C);
//find Point on B-spline surface
void SurfacePoint(int n, int p, double *U, int m, int q, double *V, vector <vector <cv::Point3f> > &P,double u, double v, cv::Point3f *S);
//find point on NURBS curve, Pw is 3D control points.
void CurvePoint(int n, int p, double *U, cv::Point3f *Pw, double u, cv::Point2f *C);
//find point on NURBS surface
//n,p,U one direction
//m,q,V the other direction
void SurfacePoint(int n, int p, double *U, int m, int q, double *V, const vector <vector <cv::Point4f> > &Pw, double u, double v, cv::Point3f *S);

// Calculate derivatives on B-Spline Curve
// calculate derivatives until order d
// result is in CK. CK[k] is kth derivative vector
// n + 1 is number of control points. 
// p is degree of the curve. 
// U is knots. P is control points
void CurveDerivsAlg(int n, int p, double *U, cv::Point2f *P, double u, int d, cv::Point2f *CK);

// Calculate derivative vectors for B-Spline surface
void SurfaceDerivsAlg(int n, int p, double *U, int m, int q, double *V, 
					  vector <vector <cv::Point3f> > &P,
					  double u, double v, int d, vector <vector <cv::Point3f> > &SKL);

// Calculate derivative vector of NURBS Curve given Aders and wders
void RatCurveDerivs(cv::Point2f *Aders, double *wders, int d, cv::Point2f *CK);

// Calculate derivative vector of NURBS Surface given Aders and wders
void RatSurfaceDerivs( vector <vector <cv::Point3f> > &Aders, vector <vector <double> > &wders, 
					int d, vector <vector <cv::Point3f> > &SKL);

//Calculate derivative vector for NURBS Curve
void NurbsCurveDerivs(int n, int p, double *U, const cv::Point3f *Pw, double u,
					  int d, cv::Point2f *CK);

//Calculate derivative vector for NURBS Surface
void NurbsSurfaceDerivs(int n, int p, double *U,
						int m, int q, double *V,
						const vector <vector <cv::Point4f> > &P,
						double u, double v, int d, vector <vector <cv::Point3f> > &SKL);

// This function uses 弦长参数化方法 to get parameterization of given data point set Q.
void chord_length_param(int n, cv::Point2f *Q, double *u);

// find parameters of a point on the curve
void CurvePointInv(int n, int p, double *U, cv::Point3f *Pw, cv::Point2f *C, double *u);

// find parameters of a point on the surface
void SurfacePointInv(int n, int p, double *U, 
					 int m, int q, double *V, 
					 const vector <vector <cv::Point4f> > &Pw, cv::Point3f *C, 
					 double *u, double *v);

// Input: n = control points num - 1, knots U, parameters u, points Q, m is number of Q - 1
// Output: control points P
void NurbsCurveFit(int n, int p, int m, double *U, double *u, cv::Point2f *Q, cv::Point2f *P);
}
#endif
