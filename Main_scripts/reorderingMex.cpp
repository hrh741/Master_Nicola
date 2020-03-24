#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <mex.h>
#include <omp.h>
#include <algorithm>

/*    
 *   reorderingMex.cpp
 *   Reorder the dislocations on the same slip plane
 *   Zebang Zheng, Imperial College London
 *   18 Jan 2016
*/
// Compile Command for Windows 
// mex reorderingMex.cpp COMPFLAGS="/openmp $COMPFLAGS"

// Compile Command for Mac 
// 
void reordering(int *, int *, int, int, double *, double *, double, double, double *, int *, int *, double);
typedef std::pair<double,int> mypair;
bool comparator ( const mypair& l, const mypair& r)
    { return l.first < r.first; }

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int *activeplanes, *plane, *pinned, *irmbound;
    int ndis;
    double *rdis, *vdis, *eta_obs;
    double reps, dt, damping;
    int i, npts;
    double *pinnedout, *irmboundout, *vdisout, *eta_obsout;
    activeplanes = (int *) mxGetPr(prhs[0]);
    plane = (int *) mxGetPr(prhs[1]);
    pinned = (int *) mxGetPr(prhs[2]);
    irmbound = (int *) mxGetPr(prhs[3]);
    ndis = (int) mxGetScalar(prhs[4]);
    rdis = (double *) mxGetPr(prhs[5]);
    vdis = (double *) mxGetPr(prhs[6]);
    eta_obs = (double *) mxGetPr(prhs[7]);
    reps = mxGetScalar(prhs[8]);
    dt = mxGetScalar(prhs[9]);
    damping = mxGetScalar(prhs[10]);
    npts = (int) mxGetScalar(prhs[11]);
    
    plhs[0] = mxCreateDoubleMatrix(ndis, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(ndis, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(ndis, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(ndis, 1, mxREAL);
    vdisout = (double *) mxGetPr(plhs[0]);
	irmboundout = (double *) mxGetPr(plhs[1]);
	pinnedout = (double *) mxGetPr(plhs[2]);
	eta_obsout = (double *) mxGetPr(plhs[3]);
    
    #pragma omp parallel for
	for (i=0; i<npts; i++)
    {
        reordering(activeplanes, plane, ndis, activeplanes[i], rdis, vdis, dt, reps, eta_obs, irmbound, pinned, damping);
    }
    for (i=0; i<ndis; i++)
    {
        vdisout[i] = vdis[i];
        irmboundout[i] = irmbound[i];
        pinnedout[i] = pinned[i];
        eta_obsout[i] = eta_obs[i];
    }
}

void reordering(int *activeplanes, int *plane, int ndis, int actplane, double *rdis, double *vdis, double dt, double reps, double *eta_obs, int *irmbound, int *pinned, double damping)
{
    int i, j, *m, ndisp, r_trigger, n_reorder, x_reorder;
    double *r, *v, *vtrial;
    int dis_indexi, dis_indexj;
    
    m = new int [1000]; // assume the maximum dislocations on one slip plane is 1000
    r = new double [1000]; // assume the maximum dislocations on one slip plane is 1000
    v = new double [1000]; // assume the maximum dislocations on one slip plane is 1000
    vtrial = new double [1000]; // assume the maximum dislocations on one slip plane is 1000
    
    r_trigger = 0;
    n_reorder = 0;
    ndisp = 0;
    for (j=0; j<ndis; j++)
    {
        if ( plane[j] == actplane )
        {
            m[ndisp] = j;
            ndisp++;
        }
    }
    if( ndisp>=2 )
    {
        mypair rr[1000]; // assume the maximum dislocations on one slip plane is 1000
        for (j=0; j<ndisp; j++)
        {
            r[j] = rdis[m[j]];
            vtrial[j] = vdis[m[j]];
        }
        
        for (j=0; j<ndisp; j++)
        {
            rr[j].first = r[j];
            rr[j].second = j;
        }
        std::sort(rr,rr+ndisp);
        std::sort(r,r+ndisp);
        for (j=0; j<ndisp; j++)
            v[j] = vtrial[rr[j].second];
        while (r_trigger == 0 && n_reorder<=10)
        {
            n_reorder++;
            x_reorder = 0;
            for (j=0; j<ndisp-1; j++)
            {
                if ((r[j]+v[j]*dt)-(r[j+1]+v[j+1]*dt)>=0)
                {
                    x_reorder = 1;
                    i = j+1;
                    dis_indexi = m[rr[i].second];
                    dis_indexj = m[rr[j].second];
                    if (v[j]==0 && v[i]>=0)
                    {
                        if (irmbound[dis_indexj]==1 || irmbound[dis_indexj]==1) // modify v(i)
                        {
                            v[i] = (r[i]+reps)/dt;
                            irmbound[dis_indexi] = 0;
                            pinned[dis_indexi] = 0;
                            eta_obs[dis_indexi] = 0;
                        }
                        else if (irmbound[dis_indexj]==2 || irmbound[dis_indexj]==2) // modify v(j)
                        {
                            v[j] = (r[j]-reps)/dt;
                            irmbound[dis_indexj] = 0;
                            pinned[dis_indexj] = 0;
                            eta_obs[dis_indexj] = 0;
                        }
                        else // modify velocity of both
                        {
                            v[i] = (r[i]+reps/2)/dt;
                            irmbound[dis_indexi] = 0;
                            pinned[dis_indexi] = 0;
                            eta_obs[dis_indexi] = 0;
                            v[j] = (r[j]-reps/2)/dt;
                            irmbound[dis_indexj] = 0;
                            pinned[dis_indexj] = 0;
                            eta_obs[dis_indexj] = 0;
                        }
                    }
                    else if (v[j]>0 && v[i]>=0) // modify v(j)
                    {
                        v[j] = v[j]*damping;
                        irmbound[dis_indexj] = 0;
                        pinned[dis_indexj] = 0;
                        eta_obs[dis_indexj] = 0;
                    }
                    else if (v[j]<=0 && v[i]<0) // modify v(i)
                    {
                        v[i] = v[i]*damping;
                        irmbound[dis_indexi] = 0;
                        pinned[dis_indexi] = 0;
                        eta_obs[dis_indexi] = 0;
                    }
                    else if (v[j]>0 && v[i]<0) // modify velocity of both
                    {
                        v[i] = v[i]*damping;
                        irmbound[dis_indexi] = 0;
                        pinned[dis_indexi] = 0;
                        eta_obs[dis_indexi] = 0;
                        v[j] = v[j]*damping;
                        irmbound[dis_indexj] = 0;
                        pinned[dis_indexj] = 0;
                        eta_obs[dis_indexj] = 0;
                    }
                }
            }
            if (x_reorder == 0)
            {
                r_trigger = 1;
            }
        }
        for (j=0; j<ndisp; j++)
        {
            vdis[m[rr[j].second]] = v[j];
        }
    } 
    delete[] m;
    delete[] r;
    delete[] v;
    delete[] vtrial;
}