#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <mex.h>
#include <omp.h>

/*    
 *   displacementMex.cpp
 *   Calculate tilda displacement field
 *   Burgers vector is using vector instead of scalar
 *   Zebang Zheng, Imperial College London
 *   16 March 2017
*/

// Compile Command for Windows 
// mex displacementMex.cpp COMPFLAGS="/openmp $COMPFLAGS"

void displacement(int *, double, double, double *, double, double, double *, int, int, double *, double *, double *, int *, double *, double *, double *, int *, double *, double *);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* 
prhs	Array of right-side input arguments.
plhs	Array of left-side output arguments.
nrhs	Number of right-side arguments, or the size of the prhs array.
nlhs	Number of left-side arguments, or the size of the plhs array.(i.e 1 i.e. utilda)
*/

{
	int *node, *nSystems;
	double *xnode, *ynode, *lambdaC;
	double uConst, nu;
	int mno, npts, ndis, nout;
	double *xdis, *ydis, *alpha, *b, *bOut;
	int *type;
	double *xdisOut, *ydisOut, *alphaOut;
	int *typeOut;
	double *utilda;
	int i;
	
	node = (int *) mxGetPr(prhs[0]);
	nSystems = (int *) mxGetPr(prhs[1]);
	xnode = (double *) mxGetPr(prhs[2]);
	ynode = (double *) mxGetPr(prhs[3]);
	lambdaC = (double *) mxGetPr(prhs[4]);
	uConst = mxGetScalar(prhs[5]);
	nu = mxGetScalar(prhs[6]);
	b = (double *) mxGetPr(prhs[7]);
	mno = (int) mxGetScalar(prhs[8]);       // total number of nodes  
	npts = (int) mxGetScalar(prhs[9]);      // no of points of mesh on which u_node is defined 
	ndis = (int) mxGetScalar(prhs[10]);
	nout = (int) mxGetScalar(prhs[11]);
	xdis = (double *) mxGetPr(prhs[12]);
	ydis = (double *) mxGetPr(prhs[13]);
	alpha = (double *) mxGetPr(prhs[14]);
	type = (int *) mxGetPr(prhs[15]);
	xdisOut = (double *) mxGetPr(prhs[16]);
	ydisOut = (double *) mxGetPr(prhs[17]);
	alphaOut = (double *) mxGetPr(prhs[18]);
	typeOut = (int *) mxGetPr(prhs[19]);
    bOut = (double *) mxGetPr(prhs[20]);
    
	plhs[0] = mxCreateDoubleMatrix((2*mno), 1, mxREAL); // Create 2-D, double-precision, floating-point mxArray initialized to 0 
	// mxCreateDoubleMatrix(rows,columns, real/imaginary)
	utilda = (double *) mxGetPr(plhs[0]);
	
	for (i=0; i<mno; i++)
	{
		utilda[2*i] = 0.0;
		utilda[2*i+1] = 0.0;
	}
	#pragma omp parallel for
	for (i=0; i<npts; i++)
	{
		displacement(nSystems, xnode[node[i]-1], ynode[node[i]-1], lambdaC, uConst, nu, b, ndis, nout, xdis, ydis, alpha, type, xdisOut, ydisOut, alphaOut, typeOut, bOut, &utilda[2*(node[i]-1)]);
	}
}

void displacement(int *nSystems, double xNode, double yNode, double *lambdaC, double uConst, double nu, double *b, int ndis, int nout, double *xdis, double *ydis, double *alpha, int *type, double *xdisOut, double *ydisOut, double *alphaOut, int *typeOut, double *bOut, double *u)
{
	int j, k;
	int isys, sum_isys;
	double delta, delta2, rinv;
	double u_temp[2]={0};
	double dx_temp, dy_temp, r2inv, dx_temp2, dy_temp2;
	// double R[2][2]={0}, RT[2][2]={0};
	double *dx, *dy, *r2, *bL;
	dx = new double [ndis+nout];
	dy = new double [ndis+nout];
	r2 = new double [ndis+nout];
	int *typeL;
	typeL = new int [ndis+nout];
    bL = new double [ndis+nout];
    
	u[0] = 0;
	u[1] = 0;
	if ((ndis+nout)>0)
	{
		//delta = 100.0*b; /* field point is a node */
		delta = 0.025; // radius of the core /* field point is a node */
        delta2 = delta*delta; // just for comparing with r2[j]
		for (j=0; j<ndis; j++)
		{
			dx[j] = xNode-xdis[j]; // distance from each node(s) to the each dislocation(s) inside
			dy[j] = yNode-ydis[j];
			typeL[j] = type[j];
            bL[j] = b[j];
		}
// 		for (j=0; j<nout; j++)
// 		{
// 			dx[ndis+j] = xNode-xdisOut[j]; // distance from each node(s) to the each dislocation(s) out
// 			dy[ndis+j] = yNode-ydisOut[j];
// 			typeL[ndis+j] = typeOut[j];
//          bL[ndis+j] = bOut[j];
// 		}
		for (j=0; j<(ndis+nout); j++)
		{
			r2[j] = dx[j]*dx[j]+dy[j]*dy[j];
		/* ======================================================== */
		/* Cut-off Distance (same as minDistCheck.m in MATLAB code) */
			if (r2[j]<delta2) // if dislocation is in the core
			{
				if (r2[j]==0) // if dislocation is on the source
				{
					dx[j] = delta;
					dy[j] = 0;
					r2[j] = delta2;
				}
				else
				{
					rinv  = delta/sqrt(r2[j]);
					dx[j] = dx[j]*rinv;
					dy[j] = dy[j]*rinv;
					r2[j] = dx[j]*dx[j]+dy[j]*dy[j];
				}
			}
		/* ======================================================== */
		}

		for (k=0; k<ndis; k++)
        {
            dx_temp = +dx[k]*cos(alpha[k]) + dy[k]*sin(alpha[k]);
            dy_temp = -dx[k]*sin(alpha[k]) + dy[k]*cos(alpha[k]);

            r2inv = 1/r2[k];
            dx_temp2 = dx_temp*dx_temp;
            dy_temp2 = dy_temp*dy_temp;
            // x and y component of velocity (Needleman)
            u_temp[0] = (0.5*dx_temp*dy_temp*r2inv-(1-nu)*atan(dx_temp/dy_temp))*typeL[k]*bL[k];
            u_temp[1] = (0.5*dy_temp2*r2inv-0.25*(1-2*nu)*log(r2[k]/bL[k]/bL[k]))*typeL[k]*bL[k];
            // update new position of dislocation
            u[0] += cos(alpha[k])*u_temp[0] - sin(alpha[k])*u_temp[1];
            u[1] += sin(alpha[k])*u_temp[0] + cos(alpha[k])*u_temp[1];
            //* next in dislocation
        }
        
//         for (k=0; k<nout; k++)
//         {
//             dx_temp = +dx[ndis+k]*cos(alphaOut[k]) + dy[ndis+k]*sin(alphaOut[k]);
//             dy_temp = -dx[ndis+k]*sin(alphaOut[k]) + dy[ndis+k]*cos(alphaOut[k]);
// 
//             r2inv = 1/r2[ndis+k];
//             dx_temp2 = dx_temp*dx_temp;
//             dy_temp2 = dy_temp*dy_temp;
//  // if dx_temp = infinity?
//             u_temp[0] = (0.5*dx_temp*dy_temp*r2inv-(1-nu)*atan(dx_temp/dy_temp))*typeL[ndis+k]*bL[ndis+k];
//             u_temp[1] = 0.0;
//             u[0] += cos(alphaOut[k])*u_temp[0] - sin(alphaOut[k])*u_temp[1];
//             u[1] += sin(alphaOut[k])*u_temp[0] + cos(alphaOut[k])*u_temp[1];
//             //* next out surface dislocation
//         }
        
		u[0] *= uConst;  // uConst = 1./(2*pi*(1-nu)) defined in InComputational.m
		u[1] *= uConst;
	}
	//* Deallocate dynamic memory
	delete[] dx;
	delete[] dy;
	delete[] r2;
	delete[] typeL;
    delete[] bL;
}