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

//Initialize cells that are only present on a cirkel (comment when working on square)
int init_cells(VOX* pv)
{
    int v, vx, vy;
    int NRc;
    double r01, R, Rloc;
    //double d; int dx, dy; // distance to center
    R=min(par.NVX,par.NVY)/2+0.5;
    NRc = 0;
    for(vy=0; vy<par.NVY; vy++)
    for(vx=0; vx<(par.NVX); vx++)
    {
        v = vx + vy*(par.NVX);
        Rloc = sqrt((vx+0.5-R)*(vx+0.5-R)+(vy+0.5-R)*(vy+0.5-R));
        if((Rloc<R-par.BOUNDARYDIS) && (Rloc >30)) // exclude outer rim
        {
            r01 = rand()/(double)RAND_MAX;
            //if(r01<par.CELLDENSITY/par.TARGETVOLUME)
            if(r01<par.CELLDENSITY) 
            {
                NRc++;
                pv[v].ctag = NRc;

            }
        }

    }

    return NRc;
}

/*

//Initialize cells on square (Comment when working on cirkel):
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
*/




////////////////////////////////////////////////////////////////////////////////
/*
//Initialize forces pulling outwards on a cirkel: NEW VERSION!
void set_forces(NOD* pn)
{
    int n, nx, ny;
    double R, mp, Rloc, sina, cosa;
    R = min(NNX,NNY)/2-0.6;
    mp = min(NNX,NNY)/2;
    //cout << endl << "set forces happens " << endl;

    for(n=0; n<NN; n++)
    {
        pn[n].fx = .0;
        pn[n].fy = .0;
    }

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;
        Rloc = sqrt((nx+0.5-mp)*(nx+0.5-mp)+(ny+0.5-mp)*(ny+0.5-mp));

        if((Rloc > R))
        {
            sina = (nx+0.5-mp)/Rloc;
            cosa = (ny+0.5-mp)/Rloc;
            pn[n].fx += FORCE*sina;
            pn[n].fy += FORCE*cosa;
        }
    }
}


*/
/*
//Initialize forces pulling outwards on a cirkel (comment when working on square)
void set_forces(NOD* pn)
{
    int n, m, nx, R, ny, Tr, Tc;//Tr en Tc zijn toggles voor rij en column
    double Rloc, sina, cosa; 
    //Ik ga er nu even vanuit dat NNX en NNY even zijn!
    R=min(NNX,NNY)/2; //straal van cirkel, middelpunt van de cirkel is (R,R) dus de grootste cirkel die past in de linker bovenhoek
    Tr = NN;
    Tc = NN; 
    //a = (0.0/6.0) * 3.1416;
	//a=par.LOADANGLE/180*3.1416;

    for(n=0; n<NN; n++)
    {
        pn[n].fx = .0;
        pn[n].fy = .0;
    }
    //Left lower quadrant
    for(ny=0; ny<R; ny++)
    for(nx=0; nx<R; nx++)
    {
        n = nx + ny*NNX;
        Rloc = sqrt((nx+0.5-R)*(nx+0.5-R)+(ny+0.5-R)*(ny+0.5-R));
        sina = (nx+0.5-R)/Rloc;
        cosa = (ny+0.5-R)/Rloc;

        if((Rloc<R) && (Rloc>(R-1)))
        {
            if((Tr!=ny) && (Tc!=nx))//hoekjes
            {

                pn[n].fx += 0.5*FORCE*(-sina*cosa-sina*sina);
                pn[n].fy += 0.5*FORCE*(-sina*cosa-cosa*cosa);
            }
            else if((ny==0) && (Tr==0))//rechte stuk onderkant
            {
                pn[n].fx += -sina*cosa*FORCE;
                pn[n].fy += -cosa*cosa*FORCE;
            }    
            else if((nx==0) && (Tc==0))//rechte stuk linkerkant
            {
                pn[n].fx += -sina*sina*FORCE;
                pn[n].fy += -sina*cosa*FORCE;
            }   
        Tr=ny;
        Tc=nx;
        }
    }
    //Left upper quadrant 
    for(ny=R; ny<2*R-1; ny++)
    for(nx=0; nx<R; nx++)
    {
        n = nx + ny*NNX;
        m = nx + (2*R-1-ny)*NNX;
        pn[n].fx=pn[m].fx;
        pn[n].fy=-pn[m].fy;
    }
    //Right lower quadrant 
    for(ny=0; ny<R; ny++)
    for(nx=R; nx<2*R-1; nx++)
    {
        n = nx + ny*NNX;
        m = 2*R-1-nx + ny*NNX;
        pn[n].fx=-pn[m].fx;
        pn[n].fy=pn[m].fy;
    }
    //Right upper quadrant 
    for(ny=R; ny<2*R-1; ny++)
    for(nx=R; nx<2*R-1; nx++)
    {
        n = nx + ny*NNX;
        m = 2*R-1-nx + (2*R-1-ny)*NNX;
        pn[n].fx=-pn[m].fx;
        pn[n].fy=-pn[m].fy;
    }
}
*/


void set_forces(NOD* pn)
{ /*
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
            pn[n].fx += sin(a)*cos(a)*FORCE;
            pn[n].fy += -cos(a)*cos(a)*FORCE;
        }
        // upper plate (iy==NNY-1) loading
        if(ny==NNY-1)
        {
            pn[n].fx += -sin(a)*cos(a)*FORCE;
            pn[n].fy += cos(a)*cos(a)*FORCE;
        }
        // left plate (ix==0) loading
        
        if(nx==9)
        {
            pn[n].fx += -FORCE; //-sin(a)*sin(a)*FORCE;
            pn[n].fy += 0; //sin(a)*cos(a)*FORCE;
        }
        
        if(nx==0)
        {
            pn[n].fx += FORCE; //-sin(a)*sin(a)*FORCE;
            pn[n].fy += 0; //sin(a)*cos(a)*FORCE;
        }
        // right plate (ix==NNX-1) loading
      
        if(nx==NNX-10)
        {
            pn[n].fx += FORCE; //sin(a)*sin(a)*FORCE;
            pn[n].fy += 0; //-sin(a)*cos(a)*FORCE;
        } 
        if(nx==NNX-1)
        {
            pn[n].fx += -FORCE; //sin(a)*sin(a)*FORCE;
            pn[n].fy += 0; //-sin(a)*cos(a)*FORCE;
        }
    }

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        // for loading on the side of a plate, forces are lower
        if(((nx==9)||(nx==NNX-10)) && ((ny==0)||(ny==NNY-1)))
        {
		//WHY???????!!!! volume behoud???
            pn[n].fx *= .5;
            pn[n].fy *= .5;
        }
    }
    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        // for loading on the side of a plate, forces are lower
        if(((nx==9)||(nx==NNX-10)) && ((ny==0)||(ny==NNY-1)))
        {
		//WHY???????!!!! volume behoud???
            pn[n].fx = 0;
            pn[n].fy = 0;
        }
    }*/
}

////////////////////////////////////////////////////////////////////////////////
//Set restricitions op bijzondere plekken: 
/*
void set_restrictions(NOD* pn)
{

    
    //int fixn1, fixn2, fixn3;
    // impose boundary conditions
    //fixn1=0;
    //fixn2=NNX-1;
    //fixn3=NNX*(NNY-1);
    //pn[fixn1].restrictx=TRUE;
    //pn[fixn1].restricty=TRUE;
    //pn[fixn2].restricty=TRUE;
    //pn[fixn3].restrictx=TRUE;
    

    int n, nx, ny;

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        if((nx == 10)||(nx == 11)||(nx == 12))
        {
            pn[n].restrictx=TRUE;
            pn[n].restricty=TRUE;
        }
    }
}
*/
/*
//Set restrictions square:
void set_restrictions(NOD* pn)
{

    
    //int fixn1, fixn2, fixn3;
    // impose boundary conditions
    //fixn1=0;
    //fixn2=NNX-1;
    //fixn3=NNX*(NNY-1);
    //pn[fixn1].restrictx=TRUE;
    //pn[fixn1].restricty=TRUE;
    //pn[fixn2].restricty=TRUE;
    //pn[fixn3].restrictx=TRUE;
    

    int n, nx, ny;

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;

        if((nx==0)||(nx==NNX-1)||(ny==0)||(ny==NNY-1))
        {
         //   pn[n].restrictx=TRUE;
           // pn[n].restricty=TRUE;
        }
    }
}*/

//Set restrictions cirkel 2.0-nu ik weet wat nodes ookalweer zijn
void set_restrictions(NOD* pn)
{

    
    //int fixn1, fixn2, fixn3;
    // impose boundary conditions
    //fixn1=0;
    //fixn2=NNX-1;
    //fixn3=NNX*(NNY-1);
    //pn[fixn1].restrictx=TRUE;
    //pn[fixn1].restricty=TRUE;
    //pn[fixn2].restricty=TRUE;
    //pn[fixn3].restrictx=TRUE;
    

    int n, nx, ny;
    double R, mp, Rloc;
    R = min(NNX,NNY)/2-0.1;
    mp = min(NNX,NNY)/2;

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;
        Rloc = sqrt((nx-mp)*(nx-mp)+(ny-mp)*(ny-mp));
        if(Rloc > R)
        {
            pn[n].restrictx=TRUE;
            pn[n].restricty=TRUE;
        }
    
    }
}
/*

//Set restrictions cirkel, we gaan deze functie pas in als global strain false is..:
void set_restrictions(NOD* pn)
{

    
    //int fixn1, fixn2, fixn3;
    // impose boundary conditions
    //fixn1=0;
    //fixn2=NNX-1;
    //fixn3=NNX*(NNY-1);
    //pn[fixn1].restrictx=TRUE;
    //pn[fixn1].restricty=TRUE;
    //pn[fixn2].restricty=TRUE;
    //pn[fixn3].restrictx=TRUE;
    

    int n, nx, ny;
    double R, mp, Rloc;
    R = min(NNX,NNY)/2-0.6;
    mp = min(NNX,NNY)/2;

    for(ny=0; ny<NNY; ny++)
    for(nx=0; nx<NNX; nx++)
    {
        n = nx + ny*NNX;
        Rloc = sqrt((nx+0.5-mp)*(nx+0.5-mp)+(ny+0.5-mp)*(ny+0.5-mp));
        if(Rloc > R)
        {
            pn[n].restrictx=TRUE;
            pn[n].restricty=TRUE;
        }
    
    }
}*/
