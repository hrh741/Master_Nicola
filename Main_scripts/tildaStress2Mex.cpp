#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <mex.h>
#include <omp.h>

/*    
 *   tildaStress2Mex.cpp
 *   Calculate tilda stress field
 *   To calculate vdis using method of Chakravarthy and Curtin, 2010 MSMSE
 *   Sana Waheed, Imperial College London
 *   11 Nov 2018
*/
// Compile Command for Windows 
// mex tildaStress2Mex.cpp COMPFLAGS="/openmp $COMPFLAGS"

// Compile Command for Mac 
// 
void tildaStress(double *, double *, int *, double *, int, int, int, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
void matrix22_multi(const double A[][2],const double B[][2], double C[][2]);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *xdis, *ydis, *alpha;
	int *type;
	int ndis, flag, npts;
	double *xSource, *ySource;
	double edgeType;
    double *b;
	int i;
	double *s11, *s22, *s12, *s21, *s11x, *s11y, *s12x, *s12y, *s22x, *s22y;

	xdis = (double *) mxGetPr(prhs[0]);
	ydis = (double *) mxGetPr(prhs[1]);
	type = (int *) mxGetPr(prhs[2]);
	alpha = (double *) mxGetPr(prhs[3]);
	ndis = (int) mxGetScalar(prhs[4]);
	flag = (int) mxGetScalar(prhs[5]);
	xSource = (double *) mxGetPr(prhs[6]);
	ySource = (double *) mxGetPr(prhs[7]);
	npts = (int) mxGetScalar(prhs[8]);
	edgeType = mxGetScalar(prhs[9]);
	b = (double *) mxGetPr(prhs[10]);

	plhs[0] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(npts, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[7] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[8] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[9] = mxCreateDoubleMatrix(npts, 1, mxREAL);   
    
	s11 = (double *) mxGetPr(plhs[0]);
	s22 = (double *) mxGetPr(plhs[1]);
	s12 = (double *) mxGetPr(plhs[2]);
	s21 = (double *) mxGetPr(plhs[3]);
    
    s11x = (double *) mxGetPr(plhs[4]);
	s11y = (double *) mxGetPr(plhs[5]);
	s12x = (double *) mxGetPr(plhs[6]);
	s12y = (double *) mxGetPr(plhs[7]);
    s22x = (double *) mxGetPr(plhs[8]);
    s22y = (double *) mxGetPr(plhs[9]);
    
    
	#pragma omp parallel for
	for (i=0; i<npts; i++)
	{
		tildaStress(xdis, ydis, type, alpha, ndis, i, flag, xSource[i], ySource[i], edgeType, b, &s11[i], &s22[i], &s12[i], &s21[i], &s11x[i], &s11y[i], &s12x[i], &s12y[i], &s22x[i], &s22y[i]);
	}
}

void tildaStress(double *xdis, double *ydis, int *type, double *alpha, int ndis, int i, int flag, double xNode, double yNode, double edgeType, double *b, double *s11, double *s22, double *s12, double *s21, double *s11x, double *s11y, double *s12x, double *s12y, double *s22x, double *s22y)
{
	int j, k;
	int ng, isys, sum_isys;
	double *dx, *dy, *r2, *r4;
	dx = new double [ndis];
	dy = new double [ndis];
	r2 = new double [ndis];
	r4 = new double [ndis];
	double dx_temp, dy_temp, r4inv, dx_temp2, dy_temp2;
    double r6inv, Dsigma11x_temp, Dsigma11y_temp, Dsigma12x_temp, Dsigma12y_temp, Dsigma22x_temp, Dsigma22y_temp;
    double Dsigma11x, Dsigma11y, Dsigma12x, Dsigma12y, Dsigma22x, Dsigma22y;
	double delta, delta2, rinv;
	double R[2][2]={0}, RT[2][2]={0};
	double multi1[2][2]={0}, multi2[2][2]={0};
	double sigma[2][2]={0};
	int *typeL;
	typeL = new int [ndis];
	
	if (flag==1 || flag==2)
	/* field point is a (1)-source or (2)-dislocation */
	{
		delta = 0.0005; // *2.0*b
	}
	else if (flag==3)
	/* field point is a (3)-node */
	{
		delta = 0.025; //* 100.0*b
	}
	delta2 = delta*delta;

	s11[0] = 0;
	s22[0] = 0;
	s12[0] = 0;
	s21[0] = 0;
    s11x[0] = 0; s11y[0] = 0;
    s12x[0] = 0; s12y[0] = 0;
    s22x[0] = 0; s11y[0] = 0;

	if (ndis>0)
	{
		if (flag==1 || flag==3)
		/* field point is a (1)-source or (3)-node */
		{
			for (j=0; j<ndis; j++)
			{
				dx[j] = xNode-xdis[j];
				dy[j] = yNode-ydis[j];
				typeL[j] = type[j];
			}
		}
		else if (flag==2)
		/* field point is a (2)-dislocation */
		{
			for (j=0; j<ndis; j++)
			{
				dx[j] = xdis[i]-xdis[j];
				dy[j] = ydis[i]-ydis[j];
				typeL[j] = type[j];
			}
			typeL[i] = 0;
		}
		for (j=0; j<ndis; j++)
		{
			r2[j] = dx[j]*dx[j]+dy[j]*dy[j];
		/* ======================================================== */
		/* Cut-off Distance (same as minDistCheck.m in MATLAB code) */
			if (r2[j]<delta2)
			{
				if (r2[j]==0)
				{
					dx[j] = delta;
					dy[j] = 0;
					r2[j] = delta2;
				}
				else
				{
					rinv = delta/sqrt(r2[j]);
					dx[j] = dx[j]*rinv;
					dy[j] = dy[j]*rinv;
					r2[j] = dx[j]*dx[j]+dy[j]*dy[j];
				}
			}
		/* ======================================================== */
			r4[j] = r2[j]*r2[j];
		}

		/* sum_isys = 0;
		for (ng=0; ng<ngr; ng++)
		{
			for (isys=0; isys<nSystems[ng]; isys++)
			{
				R[0][0] = RC[sum_isys*4];
				R[1][0] = RC[sum_isys*4+1];
				R[0][1] = RC[sum_isys*4+2];
				R[1][1] = RC[sum_isys*4+3];
				//* Transpose Rotation Matrix
				RT[0][0] = R[0][0];
				RT[1][0] = R[0][1];
				RT[0][1] = R[1][0];
				RT[1][1] = R[1][1];
				for (k=0; k<ndis; k++)
				{
					if (alpha[k]==lambdaC[sum_isys] && ngsource[k]==(ng+1))
					{
						dx_temp = dx[k]*R[0][0] + dy[k]*R[0][1];
						dy_temp = dx[k]*R[1][0] + dy[k]*R[1][1];
						r4inv = 1/r4[k];
						dx_temp2 = dx_temp*dx_temp;
						dy_temp2 = dy_temp*dy_temp;

						sigma[0][0] = -(dy_temp*(3*dx_temp2+dy_temp2)*r4inv)*typeL[k];
						sigma[1][1] = (dy_temp*(dx_temp2-dy_temp2)*r4inv)*typeL[k];
						sigma[0][1] = (dx_temp*(dx_temp2-dy_temp2)*r4inv)*typeL[k];
						sigma[1][0] = sigma[0][1];

						matrix22_multi(RT,sigma,multi1);
						matrix22_multi(multi1,R,multi2);

						s11[0] += multi2[0][0]; //* sigma(1,1)
						s22[0] += multi2[1][1]; //* sigma(2,2)
						s12[0] += multi2[0][1]; //* sigma(1,2)
						s21[0] += multi2[1][0]; //* sigma(2,1)
                        
                        
                        // Find grad.sigma for Chakravarth and Curtin 2011 correction
                        r6inv = 1/(r2[k]*r4[k]);
                        Dsigma11x_temp = (-6*dx_temp*dy_temp*(dx_temp2+dy_temp2) + 4*dx_temp*dy_temp*(3*dx_temp2+dy_temp2))*r6inv*typeL[k];
                        Dsigma11y_temp = (-3*(dx_temp2+dy_temp2)*(dx_temp2+dy_temp2) + 4*dy_temp2*(3*dx_temp2+dy_temp2))*r6inv*typeL[k];
                        Dsigma12x_temp = ((3*dx_temp2-dy_temp2)*(dx_temp2+dy_temp2) - 4*dx_temp2*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma12y_temp = (-2*dx_temp*dy_temp*(dx_temp2+dy_temp2) - 4*dx_temp*dy_temp*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma22x_temp = (2*dx_temp*dy_temp*(dx_temp2+dy_temp2) - 4*dx_temp*dy_temp*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma22y_temp = ((dx_temp2-3*dy_temp2)*(dx_temp2+dy_temp2) - 4*dy_temp2*(dx_temp2-dy_temp2))*r6inv*typeL[k];
						
						Dsigma11x = Dsigma11x_temp*RT[0][0] + Dsigma11y_temp*RT[1][0];
						Dsigma11y = Dsigma11x_temp*RT[0][1] + Dsigma11y_temp*RT[1][1];
						Dsigma12x = Dsigma12x_temp*RT[0][0] + Dsigma12y_temp*RT[1][0];
						Dsigma12y = Dsigma12x_temp*RT[0][1] + Dsigma12y_temp*RT[1][1];
						Dsigma22x = Dsigma22x_temp*RT[0][0] + Dsigma22y_temp*RT[1][0];
						Dsigma22y = Dsigma22x_temp*RT[0][1] + Dsigma22y_temp*RT[1][1];
                        
                        s11x[0] += Dsigma11x; s11y[0] += Dsigma11y;
                        s12x[0] += Dsigma12x; s12y[0] += Dsigma12y;
                        s22x[0] += Dsigma22x; s11y[0] += Dsigma22y;   
						
					}
					//* next dislocation
				}
				//* next slip system
				sum_isys += 1;
			}
			//* next grain
		}*/
        for (k=0; k<ndis; k++)
        {
            R[0][0] = cos(alpha[k]);
            R[1][0] = -sin(alpha[k]);
            R[0][1] = sin(alpha[k]);
            R[1][1] = cos(alpha[k]);
            //* Transpose Rotation Matrix
            RT[0][0] = R[0][0];
            RT[1][0] = R[0][1];
            RT[0][1] = R[1][0];
            RT[1][1] = R[1][1];
            
            dx_temp = dx[k]*R[0][0] + dy[k]*R[0][1];
            dy_temp = dx[k]*R[1][0] + dy[k]*R[1][1];
            r4inv = 1/r4[k];
            dx_temp2 = dx_temp*dx_temp;
            dy_temp2 = dy_temp*dy_temp;
            
            sigma[0][0] = -(dy_temp*(3*dx_temp2+dy_temp2)*r4inv)*typeL[k]*b[k];
            sigma[1][1] = (dy_temp*(dx_temp2-dy_temp2)*r4inv)*typeL[k]*b[k];
            sigma[0][1] = (dx_temp*(dx_temp2-dy_temp2)*r4inv)*typeL[k]*b[k];
            sigma[1][0] = sigma[0][1];

            matrix22_multi(RT,sigma,multi1);
            matrix22_multi(multi1,R,multi2);

            s11[0] += multi2[0][0]; //* sigma(1,1)
            s22[0] += multi2[1][1]; //* sigma(2,2)
            s12[0] += multi2[0][1]; //* sigma(1,2)
            s21[0] += multi2[1][0]; //* sigma(2,1)
            
             // Find grad.sigma for Chakravarth and Curtin 2011 correction
                        r6inv = 1/(r2[k]*r4[k]);
                        Dsigma11x_temp = (-6*dx_temp*dy_temp*(dx_temp2+dy_temp2) + 4*dx_temp*dy_temp*(3*dx_temp2+dy_temp2))*r6inv*typeL[k];
                        Dsigma11y_temp = (-3*(dx_temp2+dy_temp2)*(dx_temp2+dy_temp2) + 4*dy_temp2*(3*dx_temp2+dy_temp2))*r6inv*typeL[k];
                        Dsigma12x_temp = ((3*dx_temp2-dy_temp2)*(dx_temp2+dy_temp2) - 4*dx_temp2*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma12y_temp = (-2*dx_temp*dy_temp*(dx_temp2+dy_temp2) - 4*dx_temp*dy_temp*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma22x_temp = (2*dx_temp*dy_temp*(dx_temp2+dy_temp2) - 4*dx_temp*dy_temp*(dx_temp2-dy_temp2))*r6inv*typeL[k];
                        Dsigma22y_temp = ((dx_temp2-3*dy_temp2)*(dx_temp2+dy_temp2) - 4*dy_temp2*(dx_temp2-dy_temp2))*r6inv*typeL[k];
						
						Dsigma11x = Dsigma11x_temp*RT[0][0] + Dsigma11y_temp*RT[1][0];
						Dsigma11y = Dsigma11x_temp*RT[0][1] + Dsigma11y_temp*RT[1][1];
						Dsigma12x = Dsigma12x_temp*RT[0][0] + Dsigma12y_temp*RT[1][0];
						Dsigma12y = Dsigma12x_temp*RT[0][1] + Dsigma12y_temp*RT[1][1];
						Dsigma22x = Dsigma22x_temp*RT[0][0] + Dsigma22y_temp*RT[1][0];
						Dsigma22y = Dsigma22x_temp*RT[0][1] + Dsigma22y_temp*RT[1][1];
                        
                        s11x[0] += Dsigma11x; s11y[0] += Dsigma11y;
                        s12x[0] += Dsigma12x; s12y[0] += Dsigma12y;
                        s22x[0] += Dsigma22x; s11y[0] += Dsigma22y;  
        }
        
		s11[0] *=edgeType;
		s22[0] *=edgeType;
		s12[0] *=edgeType;
		s21[0] *=edgeType;
        
        s11x[0] *=edgeType; s11y[0] *=edgeType;
        s12x[0] *=edgeType; s12y[0] *=edgeType;
        s22x[0] *=edgeType; s11y[0] *=edgeType;
	}
	//* Deallocate dynamic memory
	delete[] dx;
	delete[] dy;
	delete[] r2;
	delete[] r4;
	delete[] typeL;
}

void matrix22_multi(const double A[][2],const double B[][2], double C[][2])
{
    double a,b,c,d,e,f,g,h;
    
    a=A[0][0];
    b=A[0][1];
    c=A[1][0];
    d=A[1][1];
    
    e=B[0][0];
    f=B[0][1];
    g=B[1][0];
    h=B[1][1];
    
    C[0][0]=a*e+b*g;
    C[0][1]=a*f+b*h;
    C[1][0]=c*e+d*g;
    C[1][1]=c*f+d*h;
}