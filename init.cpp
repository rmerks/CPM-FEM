// file: init.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
VOX* init_voxels(void)
{
    //VOX* pv;
    int v, vx, vy;
    int i;

    VOX *pv = new VOX[NV];

    // set voxel information
    for(vy=0; vy<par.NVY; vy++)
    for(vx=0; vx<(par.NVX); vx++)
    {
        v = vx + vy*(par.NVX);

        //pv[v].x = vx * (par.VOXSIZE); pv[v].y = vy * (par.VOXSIZE);
        pv[v].ctag = 0;
    }
    return pv;
}

////////////////////////////////////////////////////////////////////////////////
NOD* init_nodes(void)
{
    //NOD* pn;
    int n, nx, ny;

    NOD *pn = new NOD[NN];

    // set node information
    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        //pn[n].x = nx * (par.VOXSIZE); pn[n].y = ny * (par.VOXSIZE);
        pn[n].fx = .0; pn[n].fy = .0;
        pn[n].ux = .0; pn[n].uy = .0;
        pn[n].restrictx = FALSE; pn[n].restricty = FALSE;
    }
    return pn;
}

////////////////////////////////////////////////////////////////////////////////
int init_cells(VOX* pv)
{
    int v, vx, vy;
    int NRc;
    double r01;
    double d; int dx, dy; // distance to center

    NRc = 0;
    for(vy=0; vy<par.NVY; vy++)
    for(vx=0; vx<(par.NVX); vx++)
    {
        v = vx + vy*(par.NVX);

        if((vx>0+par.BOUNDARYDIS)&&(vx<(par.NVX)-1-par.BOUNDARYDIS)&&(vy>0+par.BOUNDARYDIS)&&(vy<par.NVY-1-par.BOUNDARYDIS)) // exclude outer rim
        {
            r01 = rand()/(double)RAND_MAX;
            if(r01<par.CELLDENSITY/par.TARGETVOLUME) //if(r01<.1/par.TARGETVOLUME)
            //if((vx==(par.NVX)/2)&&(vy==par.NVY/2))
            //if(((vx==(par.NVX)/2-7)||(vx==(par.NVX)/2+7))&&(vy==par.NVY/2))
            //dx=vx-(par.NVX)/2; dy=vy-par.NVY/2; d=sqrt(dx*dx+dy*dy); if((d<(par.NVX)/8.0) && (r01<1.5/par.TARGETVOLUME))
            {
                NRc++;
                pv[v].ctag = NRc;

            }
        }

    }

    return NRc;
}





////////////////////////////////////////////////////////////////////////////////
void set_forces(NOD* pn)
{
    int n, nx, ny;
    double a;

    //a = (0.0/6.0) * 3.1416;
	a=par.LOADANGLE/180*3.1416;

    for(n=0; n<NN; n++)
    {
        pn[n].fx = .0;
        pn[n].fy = .0;
    }

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        // lower plate (iy==0) loading
        if(ny==0)
        {
            pn[n].fx +=  sin(a)*cos(a)*FORCE;
            pn[n].fy += -cos(a)*cos(a)*FORCE;
        }
        // upper plate (iy==NNY-1) loading
        if(ny==NNY-1)
        {
            pn[n].fx += -sin(a)*cos(a)*FORCE;
            pn[n].fy +=  cos(a)*cos(a)*FORCE;
        }
        // left plate (ix==0) loading
        if(nx==0)
        {
            pn[n].fx += -sin(a)*sin(a)*FORCE;
            pn[n].fy +=  sin(a)*cos(a)*FORCE;
        }
        // right plate (ix==NNX-1) loading
        if(nx==NNX-1)
        {
            pn[n].fx +=  sin(a)*sin(a)*FORCE;
            pn[n].fy += -sin(a)*cos(a)*FORCE;
        }
    }

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        // for loading on the side of a plate, forces are lower
        if(((nx==0)||(nx==NNX-1)) && ((ny==0)||(ny==NNY-1)))
        {
		//WHY???????!!!! volume behoud???
            pn[n].fx *= .5;
            pn[n].fy *= .5;
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
void set_restrictions(NOD* pn)
{

    /*
    int fixn1, fixn2, fixn3;

    // impose boundary conditions
    fixn1=0;
    fixn2=NNX-1;
    //fixn3=NNX*(NNY-1);


    pn[fixn1].restrictx=TRUE;
    pn[fixn1].restricty=TRUE;
    pn[fixn2].restricty=TRUE;
    //pn[fixn3].restrictx=TRUE;
    */

    int n, nx, ny;

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        if((nx==0)||(nx==NNX-1)||(ny==0)||(ny==NNY-1))
        {
            pn[n].restrictx=TRUE;
            pn[n].restricty=TRUE;
        }
    }

}

