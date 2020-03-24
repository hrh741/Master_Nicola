#include <stdio.h>
#include <math.h>
#include <mex.h>
#include <omp.h>

/*    
 *   hatStressMex.cpp
 *   Calculate hat stress field for 3 noded triangular element
 *   15 June 2017
*/

// Compile Command for Windows 
// mex hatStressMex.cpp COMPFLAGS="/openmp $COMPFLAGS"

// Compile Command for Mac 
// 

void hatStress(double *, int *, int , double *, double *, double *, int, int , double *, double *, double *, double *);
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *uhat, *DC;
	int *ncC;
	double *xnode, *ynode;
    int *elegrainindex, *belong;
	double *xdis, *ydis;
	int npts, mel;
	int i;
	double *s11, *s22, *s12, *s21;

	uhat = (double *) mxGetPr(prhs[0]);
	ncC = (int *) mxGetPr(prhs[1]);
	xnode = (double *) mxGetPr(prhs[2]);
	ynode = (double *) mxGetPr(prhs[3]);
	DC = (double *) mxGetPr(prhs[4]);
    elegrainindex = (int *) mxGetPr(prhs[5]);
    belong = (int *) mxGetPr(prhs[6]);
	npts = (int) mxGetScalar(prhs[7]);
    mel = (int) mxGetScalar(prhs[8]);

	plhs[0] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	s11 = (double *) mxGetPr(plhs[0]);
	s22 = (double *) mxGetPr(plhs[1]);
	s12 = (double *) mxGetPr(plhs[2]);
	s21 = (double *) mxGetPr(plhs[3]);
	#pragma omp parallel for
	for (i=0; i<npts; i++)
	{
		hatStress(uhat, ncC, mel, xnode, ynode, DC, elegrainindex[belong[i]-1], belong[i]-1, &s11[i], &s22[i], &s12[i], &s21[i]);
	}
}

void hatStress(double *uhat, int *ncC, int mel, double *xnode, double *ynode, double *DC, int elegrainindexL, int belongL, double *s11, double *s22, double *s12, double *s21)
{
    double Be[3][6]={0};    
    double DCL[9]={0};
	double U[6]={0};
	double BU[3]={0};
	double BU_temp, detJ;
	double x1, y1, x2, y2, x3, y3;
	int i, j, k;
	int node[3] = {0};
       
	
    // nodal displacements
    for (k=0; k<3; k++)
    {
        node[k] = ncC[belongL+mel*k]-1;
        U[2*k] = uhat[2*node[k]];
        U[2*k+1] = uhat[2*node[k]+1];
    }
	
    x1 = xnode[node[0]]; y1 = ynode[node[0]];
    x2 = xnode[node[1]]; y2 = ynode[node[1]];
    x3 = xnode[node[2]]; y3 = ynode[node[2]];
    
	// detJ = 2*Area
	detJ = fabs((x1-x3)*(y2-y1)-(x1-x2)*(y3-y1));
	  
    Be[0][0] = y2-y3;
    Be[0][2] = y3-y1;
    Be[0][4] = y1-y2;

    Be[1][1] = x3-x2;
    Be[1][3] = x1-x3;
    Be[1][5] = x2-x1;

    Be[2][0] = Be[1][1];
    Be[2][1] = Be[0][0];
    Be[2][2] = Be[1][3];
    Be[2][3] = Be[0][2];
    Be[2][4] = Be[1][5];
    Be[2][5] = Be[0][4];
    
    // local material matrix
    for (i=0; i<9; i++)
    {
        DCL[i] = DC[i+(elegrainindexL-1)*9];
    }
	  
    //* sigmaHat = D*(B*U);
	for (i=0; i<3; i++)
	{
		BU_temp = 0.0;
		for (j=0; j<6; j++)
		{
			BU_temp += Be[i][j]*U[j]/detJ;
		}
		BU[i] = BU_temp;
	}
	
	s11[0] = DCL[0]*BU[0]+DCL[3]*BU[1]+DCL[6]*BU[2];
	s22[0] = DCL[1]*BU[0]+DCL[4]*BU[1]+DCL[7]*BU[2];
	s12[0] = DCL[2]*BU[0]+DCL[5]*BU[1]+DCL[8]*BU[2];
	s21[0] = s12[0];
    
    
//     delete[] Be;
//     delete[] U;
//     delete[] BU;
//     delete[] node;
}
