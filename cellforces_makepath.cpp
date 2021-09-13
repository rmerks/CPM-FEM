// file: par.CELLFORCEs.c

#include <functions.h>
#include <math.h>


////////////////////////////////////////////////////////////////////////////////
void cell_forces(VOX* pv, NOD* pn, int* csize, int NRc)
{
	int c;
	int n,nx,ny;
	int v,vx,vy, cnttag;
	int NRcelln,cellnodes[NN];
	int i,j, n2;
	double dnx,dny,forcex,forcey;

	for(ny=1; ny<NNY-1; ny++)
   	for(nx=1; nx<NNX-1; nx++)
   	{
   		n = nx + ny*NNX;
		pn[n].fx = 0;
		pn[n].fy = 0;
	}

	if(par.CELLFORCES){

		for(c=0;c<NRc;c++)
		{
			// determine which nodes belong to cell c
			NRcelln = 0;
			for(ny=1; ny<NNY-1; ny++)
	   		for(nx=1; nx<NNX-1; nx++)
	   		{
	   			n = nx + ny*NNX;
				cnttag = 0;
				for(vy=ny-1; vy<ny+1; vy++)
				for(vx=nx-1; vx<nx+1; vx++)
				{
					v = vx + vy*(par.NVX);
					if(pv[v].ctag == c+1)
						cnttag++;
				}
				if(cnttag>0) // all cell nodes
				{
					cellnodes[NRcelln] = n;
					NRcelln++;
				}
			}
			// forces between cellnodes
			for(i=0;i<NRcelln;i++)
			{
				n = cellnodes[i];

				for(j=0;j<NRcelln;j++)
				{

					n2 = cellnodes[j];
					ny=n/NNX; nx=n%NNX;
					dny=(n2/NNX-ny)*(par.VOXSIZE); // y distance between n and n2
					dnx=(n2%NNX-nx)*(par.VOXSIZE); // x distance between n and n2
					else{makepath=0;}
					//check if cell nodes n and n2 are connected by a straight line that stays within the cell c
					if(par.NODECONNEXTION}{
						if(CheckCellNodeConnection(pv,pn,n,n2,c,pathx,pathy,makepath))
						{

							forcex = par.CELLFORCE*dnx;
							forcey = par.CELLFORCE*dny;

							pn[n].fx += forcex;
							pn[n].fy += forcey;
						}
					}
					else
					{
						forcex = par.CELLFORCE*dnx;
						forcey = par.CELLFORCE*dny;

						pn[n].fx += forcex;
						pn[n].fy += forcey;
					}
				
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////



BOOL CheckCellNodeConnection(VOX* pv, NOD* pn, int n1, int n2, int c){

	int ny1 = n1/NNX; int nx1=n1%NNX;
	int ny2 = n2/NNX; int nx2 = n2%NNX;

	int diffx = nx2-nx1; int diffy = ny2-ny1;
	int vx1, vy1, vx, vy, v, n11, n22, nx, ny, n;
	double yideal, xideal;
	int bound;
	BOOL nodesconnected=TRUE;
	int nocellvox = 0;
	double slope =0;
	int counter = 0;

	if(!diffx==0){slope = diffy/diffx;}

	if((abs(diffx)>=abs(diffy))&(!diffy==0))
	{
			if(diffx<0){n22=n1;n11=n2;} //swap the nodes
			else{n11=n1;n22=n2;}

			ny1 = n11/NNX; nx1=n11%NNX;
			ny2 = n22/NNX; nx2 = n22%NNX;
			diffx = nx2-nx1; diffy = ny2-ny1;
			if(diffy>0){vx1=n11%NNX; vy1 = n11/NNX;}
			else{vx1=n11%NNX; vy1 = n11/NNX-1;}
			vx=vx1;
			vy=vy1;
			ny = ny1;
			nx=nx1;
			if(makepath==1){
			pathx[counter]= vx;
			pathy[counter]= vy;}
			if(diffx==0){

				ny++;
				while(ny<ny2){
					n = nx + ny*(NNX);
					//this node should not be on a corner, so we check all four voxels
					//first lower right corner
					vy = n/NNX;
					vx = n%NNX;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX-1;
					vx = n%NNX;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX-1;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					if(nocellvox>2){nodesconnected=FALSE;}
					ny++;
				}}
				else{
		                v = vx + vy*(par.NVX);
				while(vx<nx2-1){
					if(!pv[v].ctag==c+1){nodesconnected=FALSE;}
					vx++;
					yideal= ny1+slope*(vx+0.5-nx1);
					vy = ceil(yideal)-1;
					if(vy<0){vy=0;}
		                        v = vx + vy*(par.NVX);
					counter++;
					if(makepath==1){
					pathx[counter]=vx;
					pathy[counter]=vy;}
				}}

	}
	else
	{
			if(diffy<0){n22=n1;n11=n2;} //swap the nodes
			else{n11=n1;n22=n2;}
			ny1 = n11/NNX; nx1=n11%NNX;
			ny2 = n22/NNX; nx2 = n22%NNX;
			diffx = nx2-nx1; diffy = ny2-ny1;
			if(diffx>0){vx1=n11%NNX; vy1 = n11/NNX;}
			else{	vx1=n11%NNX-1; vy1 = n11/NNX;}
			vx=vx1;
			vy=vy1;
			nx=nx1;
			ny=ny1;
			if(makepath==1){
			pathx[counter]= vx;
			pathy[counter]= vy;}
			if(diffy==0){
				nx++;
				while(nx<nx2){
					n = nx + ny*(NNX);
					//this node should not be on a corner, so we check all four voxels
					//first right lower corner
					vy = n/NNX;
					vx = n%NNX;

		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX-1;
					vx = n%NNX;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX-1;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					vy = n/NNX;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
					if(!pv[v].ctag==c+1){nocellvox++;}
					if(nocellvox>2){nodesconnected=FALSE;}
					nx++;
				}}
				else{
		                v = vx + vy*(par.NVX);
				while(vy<ny2-1){	
					if(!pv[v].ctag==c+1){ nodesconnected=FALSE;}
					vy++;
					xideal = nx1 +(1/slope)*(vy+0.5-ny1);
					vx = ceil(xideal)-1;
					if(vx<0){vx=0;}
		                        v = vx + vy*(par.NVX);
					counter++;
					if(makepath==1){
					pathx[counter]=vx;
					pathy[counter]=vy;}
				}}

	}
	return nodesconnected;
}
