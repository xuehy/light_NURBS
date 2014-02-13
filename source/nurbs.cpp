#include "nurbs.h"
const double eps = 1.0e-13;
namespace cv
{
Point4f :: Point4f(double a, double b, double c, double d)
{
    w = a;
    x = b;
    y = c;
    z = d;
}
Point4f :: Point4f(const Point4f &P)
{
    w = P.w;
    x = P.x;
    y = P.y;
    z = P.z;
}

Point4f operator +(const Point4f &P1, const Point4f &P2)
{
    Point4f R(P1.w + P2.w,P1.x+P2.x,P1.y+P2.y,P1.z+P2.z);
    return R;
}
Point4f operator -(const Point4f &P1, const Point4f &P2)
{
    Point4f R(P1.w - P2.w, P1.x - P2.x, P1.y - P2.y ,P1.z - P2.z);
    return R;
}
Point4f operator *(const Point4f &P, double a)
{
    Point4f R(a*P.w, a*P.x, a*P.y, a*P.z);
    return R;
}
Point4f operator *(double a, const Point4f &P)
{
    Point4f R(a*P.w, a*P.x, a*P.y, a*P.z);
    return R;
}
Point4f operator /(const Point4f &P, double a)
{
    Point4f R(P.w / a, P.x / a, P.y / a, P.z / a);
    return R;
}
Point4f & Point4f :: operator =(const Point4f &P)
{
    this->w = P.w;
    this->x = P.x;
    this->y = P.y;
    this->z = P.z;
}

}

namespace nurbs
{
double factorial(int n)
{
	double r = 1.0;
	while( n > 1 )
	{
		r = r * 1.0 * n;
		n--;
	}
	return r;
}
	
double Bin(int k, int i)
{
	return factorial(k) / (factorial(i) * factorial(k - i)) ;
}

void Horner1(double *a, int n, double u0, double *C)
{
    *C = a[n];
    for(int i = n - 1; i >= 0;  i--) *C = *C * u0 + a[i];
}

void Bernstein(int i, int n, double u, double *B)
{
    double *temp = new double[n+1];
    for(int j = 0; j <= n; j++)
    {
        temp[j] = 0;
    }
    temp[n - i] = 1.0;
    double u1 = 1.0 - u;
    for(int k = 1; k <= n; k++)
        for(int j = n; j >= k; j--)
            temp[j] = u1 * temp[j] + u * temp[j - 1];
    *B = temp[n];
    delete []temp;
}

void AllBernstein(int n, double u, double *B)
{
    B[0] = 1.0;
    double u1 = 1.0 - u;
    for(int j = 1; j <= n; j++)
    {
        double saved = 0.0;
        for(int k = 0; k < j; k++)
        {
            double temp = B[k];
            B[k] = saved + u1 * temp;
            saved = u * temp;
        }
        B[j] = saved;
    }
}

void PointOnBezierCurve(cv::Point2f *P, int n, double u, cv::Point2f *C)
{
    double *B = new double [n + 1];
    AllBernstein(n,u,B);
    for(int k = 0; k <= n; k++)
    {
        *C = *C + B[k] * P[k];
    }
    delete []B;
}

void deCasteljau1(cv::Point2f *P, int n, double u, cv::Point2f *C)
{
    cv::Point2f * Q = new cv::Point2f [n + 1];
    for(int i = 0; i <= n; i++)
    {
        Q[i] = P[i];
    }
    for(int k = 1; k <= n; k++)
        for(int i = 0; i <= n - k; i++)
            Q[i] = (1.0-u)*Q[i] + u*Q[i+1];
    *C = Q[0];
    delete []Q;
}

void deCasteljau1(cv::Point3f *P, int n, double u, cv::Point3f *C)
{
    cv::Point3f * Q = new cv::Point3f [n + 1];
    for(int i = 0; i <= n; i++)
    {
        Q[i] = P[i];
    }
    for(int k = 1; k <= n; k++)
        for(int i = 0; i <= n - k; i++)
            Q[i] = (1.0-u)*Q[i] + u*Q[i+1];
    *C = Q[0];
    delete []Q;
}

int FindSpan(int n, int p, double u, double *U)
{
    if(abs(u - U[n+1]) <= eps) return n;
    int low = p;
    int high = n + 1;
    int mid = (low + high) / 2;
    while(u < U[mid] || u >= U[mid + 1])
    {
        if(u < U[mid]) high = mid;
        else	low = mid;
        mid = (low + high) / 2;
    }

    return mid;
}

void BasisFuns(int i, double u, int p, double *U, double *N)
{
    N[0] = 1.0;
    double * left = new double [p+1];
    double * right = new double [p+1];
    for(int j = 1; j <= p; j++)
    {
        left[j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for(int r = 0; r < j; r ++)
        {
            double temp = N[r] / (right[r + 1] + left[j - r]);
			
            N[r] = saved + right[ r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    delete []left;
    delete []right;
}

void DersBasisFuns(int i, double u, int p, int n, double *U, vector <vector <double> > &ders)
{
	vector <vector <double> > ndu;
	ndu.resize(p+1);
	for(int i = 0; i < p + 1; i++) ndu[i].resize(p+1);
	ndu[0][0] = 1.0;
    double * left = new double [p+1];
    double * right = new double [p+1];
	
	for(int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		double saved = 0.0;
		
		for(int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			double temp = ndu[r][j - 1]/ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
		
	for(int j = 0; j <= p; j++)
	{
		ders[0][j] = ndu[j][p];
	}
	vector <vector <double> >a;
	a.resize(2);
	a[0].resize(p+1);a[1].resize(p+1);
	
	int s1,s2;
	for(int r = 0; r <= p; r++)
	{
		s1 = 0, s2 = 1;
		a[0][0] = 1.0;
		int j1,j2;
		int j;
		for(int k = 1; k <= n; k++)
		{
			double d = 0.0;
			int rk = r - k;
			int pk = p - k;
			if( r >= k )
			{
				a[s2][0] = a[s1][0]/ndu[pk+1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
				
			if( rk >= -1 ) j1 = 1;
			else j1 = -rk;
			if( r-1 <= pk) j2 = k - 1;
			else j2 = p - r;
			for(j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
				d += a[s2][j] * ndu[rk+j][pk];
			}
			if(r <= pk)
			{
				a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			j = s1; s1 = s2; s2 = j;
		}
	}
	
	int r = p;
	for(int k = 1; k <= n; k++)
	{
		for(int j = 0; j <= p; j++) ders[k][j] = ders[k][j] * r;
		r *= (p-k);
	}
    delete []left;
    delete []right;
}

void CurvePoint(int n, int p, double *U, cv::Point2f *P, double u, cv::Point2f *C)
{
    int span = FindSpan(n,p,u,U);
    double *N = new double[p + 1];
    BasisFuns(span,u,p,U,N);
    *C = cv::Point2f(0,0);
    for(int i = 0; i <= p; i++)
        *C = *C + N[i] * P[span - p + i];
    delete []N;
}

void SurfacePoint(int n, int p, double *U, int m, int q, double *V, vector <vector <cv::Point3f> > &P,double u, double v, cv::Point3f *S)
{
    int uspan = FindSpan(n,p,u,U);
    int vspan = FindSpan(m,q,v,V);
    double *Nu = new double[p + 1];
    double *Nv = new double[q + 1];
    BasisFuns(uspan,u,p,U,Nu);
    BasisFuns(vspan,v,q,V,Nv);
    int uind = uspan - p;
    int vind;
    cv::Point3f * temp = new cv::Point3f [q + 1];
    for(int l = 0; l <= q; l++)
    {
        temp[l] = cv::Point3f(0,0,0);
        vind = vspan - q + l;
        for(int k = 0; k <= p; k++)
            temp[l] = temp[l] + Nu[k] * P[uind+k][vind];
    }
    *S = cv::Point3f(0,0,0);
    for(int l = 0; l <= q; l++)
    {
        *S = *S + Nv[l] * temp[l];
    }
    delete []Nu;
    delete []Nv;
	delete []temp;
}

void CurvePoint(int n, int p, double *U, cv::Point3f *Pw, double u, cv::Point2f *C)
{
    double *N = new double [p + 1];
    int span = FindSpan(n,p,u,U);
    BasisFuns(span,u,p,U,N);
    cv::Point3f Cw(0,0,0);
    for(int j = 0; j <= p; j++)
        Cw = Cw + N[j] * Pw[span - p + j];
    *C = cv::Point2f(Cw.x/Cw.z,Cw.y/Cw.z);
    delete []N;
}

void SurfacePoint(int n, int p, double *U, int m, int q, double *V, const vector <vector <cv::Point4f> > &Pw, double u, double v, cv::Point3f *S)
{
    int uspan = FindSpan(n,p,u,U);
    int vspan = FindSpan(m,q,v,V);
    double *Nu = new double[p + 1];
    double *Nv = new double[q + 1];
    BasisFuns(uspan,u,p,U,Nu);
    BasisFuns(vspan,v,q,V,Nv);
    cv::Point4f * temp = new cv::Point4f [q + 1];
    for(int l = 0; l <= q; l++)
    {
        temp[l] = cv::Point4f(0,0,0,0);
        for(int k = 0; k <= p; k++)
            temp[l] = temp[l] + Nu[k] * Pw[uspan-p+k][vspan-q+l];
    }
    cv::Point4f Sw(0,0,0,0);
    for(int l = 0; l <= q; l++)
        Sw = Sw + Nv[l] * temp[l];
    *S = cv::Point3f(Sw.w/Sw.z,Sw.x/Sw.z,Sw.y/Sw.z);
    delete []Nu;
    delete []Nv;
    delete []temp;
}

void chord_length_param(int n, cv::Point2f *Q, double *u)
{
	double d = 0;
	for(int i = 1; i < n; i++)
	{
		d += norm(Q[i] - Q[i-1]);
	}
	u[0] = 0;
	for(int i = 1; i < n - 1; i++)
	{
		u[i] = u[i-1] + norm(Q[i]-Q[i-1]) / d;
	}
	u[n - 1] = 1;
}

void CurveDerivsAlg(int n, int p, double *U, cv::Point2f *P, double u, int d, cv::Point2f *CK)
{
	int du = min(p,d);
	for(int k = p + 1; k <= d; k++) CK[k] = cv::Point2f(0,0);
	int span = FindSpan(n,p,u,U);
	vector <vector <double> > nders;
	nders.resize(n+1);
	for(int i = 0; i < n + 1; i++) nders[i].resize(p + 1);
	DersBasisFuns(span, u, p, du, U, nders);
	for(int k = 0; k <= du; k++) 
	{
		CK[k] = cv::Point2f(0,0);
		for(int j = 0; j <= p; j++)
			CK[k] = CK[k] + nders[k][j] * P[span - p + j];
	}
}

void SurfaceDerivsAlg(int n, int p, double *U, int m, int q, double *V, 
					  vector <vector <cv::Point3f> > &P,
					  double u, double v, int d, vector <vector <cv::Point3f> > &SKL)
{
	int du = min(d,p);
	for(int k = p + 1; k <= d; k++)
		for(int l = 0; l <= d - k; l++) SKL[k][l] = cv::Point3f(0,0,0);
	int dv = min(d,q);
	for(int l = q + 1; l <= d; l++)
		for(int k = 0; k <= d - 1; k++) SKL[k][l] = cv::Point3f(0,0,0);
	int uspan = FindSpan(n,p,u,U);
	int vspan = FindSpan(n,q,v,V);
	vector <vector <double> > Nu;
	Nu.resize(n+1);
	for(int i = 0; i < n + 1; i++) Nu[i].resize(p + 1);
	DersBasisFuns(uspan, u, p, du, U, Nu);
	vector <vector <double> > Nv;
	Nv.resize(m+1);
	for(int i = 0; i < m + 1; i++) Nv[i].resize(q + 1);
	DersBasisFuns(vspan, v, q, dv, V, Nv);
	
    cv::Point3f * temp = new cv::Point3f [q + 1];
	for(int k = 0; k <= du; k++)
	{
		for(int s = 0; s <= q; s++)
		{
			temp[s] = cv::Point3f(0,0,0);
			for(int r = 0; r <= p; r++)
				temp[s] = temp[s] + Nu[k][r] * P[uspan - p + r][vspan - q + s];
		}
		
		int dd = min(d - k, dv);
		for(int l = 0; l <= dd; l++)
		{
			SKL[k][l] = cv::Point3f(0,0,0);
			for(int s = 0; s <= q; s++)
				SKL[k][l] = SKL[k][l] + Nv[l][s] * temp[s];
		}
	}
    delete []temp;
}

void RatCurveDerivs(cv::Point2f *Aders, double *wders, int d, cv::Point2f *CK)
{
	for(int k = 0; k <= d; k++)
	{
		cv::Point2f v = Aders[k];
		for(int i = 1; i <= k; i++)
			v = v - Bin(k,i) * wders[i] * CK[k - i];
		CK[k] = 1.0 / wders[0] * v;
	}
}

void RatSurfaceDerivs( vector <vector <cv::Point3f> > &Aders, vector <vector <double> > &wders, 
					int d, vector <vector <cv::Point3f> > &SKL)
{
	for(int k = 0; k <= d; k++)
		for(int l = 0; l <= d-k; l++)
		{
			cv::Point3f v = Aders[k][l];
			for(int j = 1; j <= l; j++)
				v = v - Bin(l,j) * wders[0][j] * SKL[k][l-j];
			for(int i = 1; i <= k; i++)
			{
				v = v - Bin(k,i) * wders[i][0] * SKL[k-i][l];
				cv::Point3f v2 = cv::Point3f(0,0,0);
				for(int j = 1; j <= l; j++)
					v2 = v2 + Bin(l,j) * wders[i][j] * SKL[k-i][l-j];
				v = v - Bin(k,i) * v2;
			}
			SKL[k][l] = 1.0 / wders[0][0] * v;
		}
}

void NurbsCurveDerivs(int n, int p, double *U, const cv::Point3f *Pw, double u,
					  int d, cv::Point2f *CK)
{
	cv::Point2f *Aders = new cv::Point2f [d + 1];
	double *wders = new double [d + 1];
	cv::Point2f *P = new cv::Point2f [n + 1];
	for(int i = 0; i < n + 1; i++)
	{
		P[i].x = Pw[i].x; P[i].y = Pw[i].y;
	}
	
	CurveDerivsAlg(n,p,U,P,u,d,Aders);
	
	cv::Point2f *A = new cv::Point2f [d + 1];
	for(int i = 0; i < n + 1; i++)
	{
		P[i].x = Pw[i].z; P[i].y = 0;
	}
	CurveDerivsAlg(n,p,U,P,u,d,A);
	for(int k = 0; k < d + 1; k++)
	{
		wders[k] = A[k].x;
	}
	/*vector <vector <double> > ders;
	ders.resize(n + 1);
	for(int i = 0; i < n + 1; i++) ders[i].resize(p + 1);
	
	
	int span = FindSpan(n,p,u,U);
	DersBasisFuns(span,u,p,n,U,ders);
	
	for(int k = 0; k < d + 1; k++)
	{
		wders[k] = 0;
		for(int j = p - span; j <= p - span + n; j++)
		{
			wders[k] = wders[k] + ders[k][j] * Pw[span - p + j].z;
		}
	}*/
	
	RatCurveDerivs(Aders,wders,d,CK);
	delete []Aders;
	delete []wders;
	delete []P;
	delete []A;
}

void NurbsSurfaceDerivs(int n, int p, double *U,
						int m, int q, double *V,
						const vector <vector <cv::Point4f> > &P,
						double u, double v, int d, vector <vector <cv::Point3f> > &SKL)
{
	vector <vector <cv::Point3f> > Aders;
	Aders.resize(d + 1);
	for(int i = 0 ; i < d + 1; i++) Aders[i].resize(d + 1);
	
	vector <vector <cv::Point3f> > Pw;
	Pw.resize(n + 1);
	for(int i = 0; i < n + 1; i++) Pw[i].resize(m + 1);
	
	vector <vector <double> > wders;
	wders.resize(d + 1);
	for(int i = 0 ; i < d + 1; i++) wders[i].resize(d + 1);
		
	for(int i = 0; i < n + 1; i++)
		for(int j = 0; j < m + 1; j++)
		{
			Pw[i][j].x = P[i][j].w;
			Pw[i][j].y = P[i][j].x;
			Pw[i][j].z = P[i][j].y;
		}

	SurfaceDerivsAlg(n, p, U, m, q, V, 
					 Pw,
					 u, v, d, Aders);

	for(int i = 0; i < n + 1; i++)
		for(int j = 0; j < m + 1; j++)
		{
			Pw[i][j].x = P[i][j].z;
			Pw[i][j].y = 0;
			Pw[i][j].z = 0;
		}
	vector <vector <cv::Point3f> > tAders;
	tAders.resize(d + 1);
	for(int i = 0 ; i < d + 1; i++) tAders[i].resize(d + 1);
	SurfaceDerivsAlg(n, p, U, m, q, V, 
					 Pw,
					 u, v, d, tAders);
	for(int i = 0; i < d + 1; i++)
		for(int j = 0; j < d + 1; j++)
		{
			wders[i][j] = tAders[i][j].x;
		}
	
	RatSurfaceDerivs(Aders,wders,d,SKL);
}

void CurvePointInv(int n, int p, double *U, cv::Point3f *Pw, cv::Point2f *C, double *u)
{
	int knots_num = n + p + 2;
	int N = 20;
	vector <double> U0;
	
	double eps1 = 1e-8;
	double eps2 = 1e-8;
	for(int i = 0; i < knots_num - 1; i++)
	{
		if( abs(U[i + 1] - U[i]) < eps ) continue;
		
		for(int j = 0; j < N; j++)
		{
			U0.push_back(U[i] + (double)j * (U[i + 1]- U[i]) / (double)N);
		}
	}
	
	double dist = 1000;
	double u0;
	double u1;
	for (int tt = 0; tt < U0.size(); tt++)
	{
		cv::Point2f p1;
		CurvePoint(n,p,U,Pw,U0[tt],&p1);
		double x = norm(*C - p1);
		if( x < dist )
		{
			dist = x;
			u0 = U0[tt];
		}
	}
	
	double a, b;

	for(int i = 0; i < n + p + 1; i++)
	{
		if( U[i] <= u0 && u0 < U[i+1] )
		{
			a = U[i];
			b = U[i+1];
			break;
		}
	}
	// now u0 is the best initilization parameter
	cv::Point2f CK[3];
	cv::Point2f CP;
	int cnt = 0;
	while(1)
	{
		cnt ++;
		if(cnt > 2000) 
		{ 
			break;			
		}
		NurbsCurveDerivs(n,p,U,Pw,u0,2,CK); //CK[0] = C(ui), CK[1] = C'(ui), CK[2] = C''(ui)
		CP = CK[0] - *C;
		u1 = u0 - (CK[1].x * CP.x + CK[1].y * CP.y)/(CK[2] .x + CP.x + CK[2].y * CP.y + norm(CK[1])*norm(CK[1]));
	
		if (norm(CP) <= eps1 || abs( (CK[1].x * CP.x + CK[1].y * CP.y) ) / norm(CK[1]) / norm(CP) <= eps2)
		{
			break; 
		}
		else
		{
			if ( u1 < a ) u1 = a;
			if ( u1 > b ) u1 = b;
			if ( abs(u1-u0) * norm(CK[1]) <= eps1 ) 
			{
				u0 = u1;
				break;
			}
		}
		u0 = u1;
	}
	*u = u0;
}

void SurfacePointInv(int n, int p, double *U, int m, int q, double *V, const vector <vector <cv::Point4f> > &Pw, cv::Point3f *C, double *u, double *v)
{
    double eps1 = 1e-14;
    double eps2 = 1e-14;
	
	double u0,u1,v0,v1;
	double a,b,c,d;
	//find initial value of u,v
	int knots_num_u = n + p + 2;
	int knots_num_v = m + q + 2;
	vector <double> U0;
	vector <double> V0;
    int N = 50;
	
	for(int i = 0; i < knots_num_u - 1; i++)
	{
		if( abs(U[i + 1] - U[i]) < eps ) continue;
		
		for(int j = 0; j < N; j++)
		{
			U0.push_back(U[i] + (double)j * (U[i + 1]- U[i]) / (double)N);
		}
	}
	
	for(int i = 0; i < knots_num_v - 1; i++)
	{
		if( abs(V[i + 1] - V[i]) < eps ) continue;
		
		for(int j = 0; j < N; j++)
		{
			V0.push_back(V[i] + (double)j * (V[i + 1]- V[i]) / (double)N);
		}
	}
	
    double dist = +1000000;

	for (int ut = 0; ut < U0.size(); ut++)
	{
		for (int vt = 0; vt < V0.size(); vt++)
		{
			cv::Point3f P;
			SurfacePoint(n,p,U,m,q,V,Pw,U0[ut],V0[vt],&P);
			double x = norm(*C - P);
			if( x < dist )
			{
				dist = x;
				u0 = U0[ut];
				v0 = V0[vt];
			}
		}
	}
	// u0,v0 are initlal values of u,v
	for(int i = 0; i < n + p + 1; i++)
	{
		if( U[i] <= u0 && u0 < U[i+1] )
		{
			a = U[i];
			b = U[i+1];
			break;
		}
	}
	
	for(int i = 0; i < m + q + 1; i++)
	{
		if( V[i] <= v0 && v0 < V[i+1] )
		{
			c = V[i];
			d = V[i+1];
			break;
		}
	}
	
	vector <vector <cv::Point3f> > SKL;
	SKL.resize(3);
	for(int i = 0; i < 3; i++) SKL[i].resize(3);
	double fu,fv,gu,gv,du,dv,fuv,guv;
	cv::Point3f r;
	
	int cnt = 0;

	while(1)
	{
		cnt++;
		if(cnt > 1000) break;
		NurbsSurfaceDerivs(n,p,U,m,q,V,Pw,u0,v0,2,SKL);
		r = SKL[0][0] - *C;
		fu = norm(SKL[1][0]) * norm(SKL[1][0]) + r.x * SKL[2][0].x + 
			r.y * SKL[2][0].y + r.z * SKL[2][0].z;
		fv = SKL[1][0].x * SKL[0][1].x + SKL[1][0].y * SKL[0][1].y + SKL[1][0].z * SKL[0][1].z +
			r.x * SKL[1][1].x + r.y * SKL[1][1].y + r.z * SKL[1][1].z;
		gu = fv;
		gv = norm(SKL[0][1]) * norm(SKL[0][1]) + r.x * SKL[0][2].x + 
			r.y * SKL[0][2].y + r.z * SKL[0][2].z;
		
		fuv = r.x * SKL[1][0].x + r.y * SKL[1][0].y + r.z * SKL[1][0].z;
		guv = r.x * SKL[0][1].x + r.y * SKL[0][1].y + r.z * SKL[0][1].z;
		
		du = -(gv * fuv - fv * guv) / (fu * gv - fv * gu);
		dv = -(-gu * fuv + fu * guv) / (fu * gv - fv * gu);
        if(du != du) du = 0;

        if(dv != dv) dv = 0;
		u1 = u0 + du;
		v1 = v0 + dv;

		double er1 = abs(fuv) / norm(SKL[1][0]) / norm(r);
		double er2 = abs(guv) / norm(SKL[0][1]) / norm(r);
		if( norm(SKL[0][0] - *C) <= eps1 || (er1 <= eps2 && er2 <= eps2) )
		{
			break;
		}
		else
		{
			if(u1 < a) u1 = a;
			if(u1 > b) u1 = b;
			if(v1 < c) v1 = c;
			if(v1 > d) v1 = d;
			if( norm( du*SKL[1][0] + dv*SKL[0][1] ) < eps1 )
			{
				u0 = u1;
				v0 = v1;
				break;
			}
		}
		u0 = u1;
		v0 = v1;
	}
	
	*u = u0;
	*v = v0;
}
/*void NurbsCurveFit(int n, int p, int k, double *U, double *u, cv::Point2f *Q, cv::Point2f *P)
{
	int m = n + p + 1;
	cv::Mat A(k, n + 1, CV_64F);
	int span;
	double *N = new double [p + 1];
	
	double d = (double)(m + 1) / (double)(n - p + 1);
	for(int j = 1; j <= n - p; j++)
	{
		int i = d * 1.0 * j;
		double alpha = 1.0 * j * d - 1.0 * i;
		U[p + j] = (1-alpha) * u[i - 1] + alpha * u[i];
		cout<<"U["<<p+j<<"]= "<<U[p+j]<<endl;
	}
	// for point Q[i] with param u[i]
	for(int i = 0; i < k; i++)
	{
		span = FindSpan(n,p,u[i],U);
		for(int i = 0; i < p + 1; i++) N[i] = 0;
		BasisFuns(span, u[k], p, U, &N[span - p]);   //N_i,p(u[k])
		
		double temp = 0.0;
		for(int j = 0; j < n + 1; j++) temp += N[j];
		for(int j = 0; j < n + 1; j++)
		{
			A.at<double>(i,j) = N[j] / temp;
		}
		
	}
	cout<<A<<endl;
	cv::Mat Qx(k, 1, CV_64F);
	cv::Mat Qy(k, 1, CV_64F);
	cv::Mat Px(n + 1 , 1, CV_64F);
	cv::Mat Py(n + 1 , 1, CV_64F);
	for(int i = 0; i < k; i++)
	{
		Qx.at<double>(i,0) = Q[i].x;
		Qy.at<double>(i,0) = Q[i].y;
	}

	cv::Mat B = A.t() * A;	
	// solve the problem of sigular matrix A^TA
	if(abs(cv::determinant(B)) < eps)
	{
		cv::Mat penalty = cv::Mat::eye(n + 1, n + 1, CV_64F);
		B = B + 0.000001 * penalty;
	}
	
	B = B.inv();

	B = B * A.t();

	Px  = B * Qx;
	Py  = B * Qy;
	
	for(int i = 0; i < n + 1; i++)
	{
		P[i].x = Px.at<double>(i,0);
		P[i].y = Py.at<double>(i,0);
	}
	delete []N;
}*/
void NurbsCurveFit(int n, int p, int m, double *U, double *u, cv::Point2f *Q, cv::Point2f *P)
{
	cv::Mat A(m - 1, n - 1, CV_64F);
	cv::Mat R_mat = cv::Mat::zeros(n - 1, 2, CV_64F);
	
	double d = (double)(m + 1) / (double)(n - p + 1);
	for(int i = 0; i < p + 1; i++) U[i] = 0;
	for(int i = n + 1; i < n + p + 2; i++) U[i] = 1;
	for(int j = 1; j <= n - p; j++)
	{
		int i = d * 1.0 * j;
		double alpha = 1.0 * j * d - 1.0 * i;
		U[p + j] = (1-alpha) * u[i - 1] + alpha * u[i];
	}
	
	cv::Point2f *R = new cv::Point2f [m + 1];
	int span;
	double *N = new double [n + 1];
	
	for(int k = 1; k < m; k++)
	{
		span = FindSpan(n,p,u[k],U);
		for(int i = 0; i < p + 1; i++) N[i] = 0;
		BasisFuns(span, u[k], p, U, &N[span - p]);   //N_i,p(u[k])
		
		R[k] = Q[k] - N[0] * Q[0] - N[n] * Q[m];
		for(int j = 0; j < n - 1; j++)
		{
			A.at<double>(k - 1, j) = N[j + 1];
			R_mat.at<double>(j, 0) = R_mat.at<double>(j, 0) + N[j + 1] * R[k].x;
			R_mat.at<double>(j, 1) = R_mat.at<double>(j, 1) + N[j + 1] * R[k].y;
		}
	}

	cv::Mat Rx(n - 1, 1, CV_64F);
	cv::Mat Ry(n - 1, 1, CV_64F);
	cv::Mat Px(n - 1 , 1, CV_64F);
	cv::Mat Py(n - 1 , 1, CV_64F);
	Rx = R_mat.col(0);
	Ry = R_mat.col(1);
	cv::Mat B = A.t() * A;
	if(abs(cv::determinant(B)) < eps)
	{
		cv::Mat penalty = cv::Mat::eye(n + 1, n + 1, CV_64F);
		B = B + 0.000001 * penalty;
	}
	Px = B.inv() * Rx;
	Py = B.inv() * Ry;
	
	for(int i = 1; i < n; i++)
	{
		P[i].x = Px.at<double>(i-1,0);
		P[i].y = Py.at<double>(i-1,0);
	}
	P[0] = Q[0];P[n] = Q[m];
	delete []N;
	delete []R;
}
}
