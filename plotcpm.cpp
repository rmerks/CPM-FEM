#include "plotcpm.h"
//#include <iostream>
//using namespace std;

#include "functions.h"
#include <string>
#include <sstream>


PlotCPM::PlotCPM(int xfield, int yfield, char *moviefile) :
    QtGraphics(xfield, yfield, moviefile)
{

    int vecsize=3;
    // allocate field for strain magnitude visualization
    strain = new double**[par.NVX];
    strain[0] = new double*[par.NVX*par.NVY];
    strain[0][0] = new double[vecsize*par.NVX*par.NVY];

    for (int i=1;i<par.NVX;i++) {
        strain[i]=strain[i-1]+par.NVY;
    }
    for (int i=1;i<par.NVX*par.NVY;i++) {
        strain[0][i]=strain[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*par.NVX*par.NVY;i++) {
        strain[0][0][i]=0.;
    }

    force = new double**[NNX];
    force[0] = new double*[NNX*NNY];
    force[0][0] = new double[vecsize*NNX*NNY];

    for (int i=1;i<NNX;i++) {
        force[i]=force[i-1]+NNY;
    }
    for (int i=1;i<NNX*NNY;i++) {
        force[0][i]=force[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*NNX*NNY;i++) {
        force[0][0][i]=0.;
    }


}

void PlotCPM::Plot(VOX *pv) {
    BeginScene();
    for (int vx = 0; vx < par.NVX-1; vx++ ) {
        for (int vy = 0; vy < par.NVY-1; vy++ ) {
            int v=vx + vy*par.NVX;

	int colour;

            if (pv[v].ctag<=0) {
                   colour=-1; //-1
            } else {
	      if(par.CELLCOLOUR){colour =Qt::green;} // make cells grey (yes...Qt::green is grey?)
		else{colour=-1;}
            }
            
	    for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
            		Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }
		

            if ( pv[v].ctag != pv[v+1].ctag )  // if cellborder  etc. etc.
            {
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1-l, (par.PIXPERVOX)*vy +pixy -1*l);
		}}

            }
            else{
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
                Point( colour, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1, (par.PIXPERVOX)*vy +pixy );
		}
	   }


            if ( pv[v].ctag != pv[v+par.NVX].ctag ) {
		for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+pixx-1*l, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1-l );
		}}
            } else
		{
		for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
                Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1 );
		}
		}

            // Cells that touch eachother's corners are NO neighbours //

            if (pv[v].ctag!=pv[v+par.NVX+1].ctag
                    || pv[v+1].ctag!=pv[v+par.NVX].ctag ) {
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1-l, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1-l );
            }}
            else
                Point( colour, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1 );

        }
    }
    EndScene();
}

/*
void write_pstrain(VOX* pv, NOD* pn, int increment)
{
    int v;
    double estrains[3],L1,L2,v1[2],v2[2];
    char filename[20],filename2[20];
    char astring[20];
    FILE *ofp;

    myitostr(increment, astring);
    strcpy(filename, "pstrain");
    strcat(filename, astring);
    strcat(filename, ".out");

    ofp = fopen(filename,"w");
    for(v=0;v<NV;v++)
    {
        get_estrains(pn,v,estrains);
        L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);
        if(L1>L2)
        {
            fprintf(ofp,"%d ", (int)(1000000*L1));
            fprintf(ofp,"%d ", (int)(1000*v1[0]));
            fprintf(ofp,"%d ", (int)(1000*v1[1]));
            fprintf(ofp,"%d\n",(int)(1000000*L2));
        }
        else
        {
            fprintf(ofp,"%d ", (int)(1000000*L2));
            fprintf(ofp,"%d ", (int)(1000*v2[0]));
            fprintf(ofp,"%d ", (int)(1000*v2[1]));
            fprintf(ofp,"%d\n",(int)(1000000*L1));
        }
    }
    fflush(ofp); fclose(ofp);
}
*/


void PlotCPM::PlotHueStrainField(bool STRAINFIELD) {


    BeginScene();
QColor c;
	    double max_strain;
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR = MaxStrainMagnitude();}
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR=0.01;}
	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBAR;
		if(par.MAXCOLORBAR==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

c=Qt::white;

    // Plot vector field

    // Draw strain magnitude

    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {

	if(par.STRAINFIELD)
	{	
            int hue=240-240*strain[vx][vy][2]/max_strain;
	    if(strain[vx][vy][2]>max_strain){hue=0;}
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                cerr << "Strain at ( " << vx << ", " << vy << ") is " << strain[vx][vy][2] << ", max strain: " << max_strain << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	}
        picture->setPen( c );

	
	    for(int pixx = 0; pixx < (par.PIXPERVOX); pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX); pixy++){ 	
            		picture->drawPoint((par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }

	    
        }
    }
EndScene();
}

void PlotCPM::PlotPrincipleStrainField(VOX *pv) {

BeginScene();


    // calculate strains
    //CalculateStrainField(pn);

    // draw strain lines
    for (int vx=0;vx<par.NVX;vx++) { 
        for (int vy=0;vy<par.NVY;vy++) {


	    double max_strain;
	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBAR;
		if(par.MAXCOLORBAR==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

            double linelength = sqrt(2)*(par.PIXPERVOX)/max_strain; 
            double x1= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2-linelength*strain[vx][vy][0]/2);
            double y1= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2-linelength*strain[vx][vy][1]/2);
            double x2= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2+linelength*strain[vx][vy][0]/2);
            double y2= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2+linelength*strain[vx][vy][1]/2);


            if (x1<0) x1=0;
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<0) y1=0;
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<0) x2=0;
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<0) y2=0;
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            // And draw it :-)

            //strain_magnitude[vx][vy] = sqrt(strx*strx + stry*stry);
            int v=vx + vy*par.NVX;

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,1);
			Line(x1-l,y1,x2-l,y2,1);
		}


        }
    }


    EndScene();
}

double PlotCPM::MaxStrainMagnitude(void) const {
    double max=0.;
    for (int i=0;i<par.NVX*par.NVY;i++) {
        //for (int vx=0;vx<par.NVX-1;vx++) {
        //for (int vy=0;vy<par.NVY-1;vy++) {


        double magn=strain[0][i][2];
        if (magn>max) {
            max=magn;
        }
    }
    // }
    return max;
}


void PlotCPM::CalculateForceField(NOD *pn) const{
    for (int nx=1;nx<NNX-1;nx++) {
        for (int ny=1;ny<NNX-1;ny++) {
            int n=nx + ny*NNX;

                // calculate force vector

                force[nx][ny][0] = pn[n].fx;
                force[nx][ny][1] = pn[n].fy;

            // magnitude

            double fx=force[nx][ny][0], fy=force[nx][ny][1];
            force[nx][ny][2]=sqrt(fx*fx + fy*fy);
	}
   }
}

double PlotCPM::MaxForceMagnitude(void) const {
    double max=0.;
    for (int i=0;i<NNX*NNY;i++) {
        double magn=force[0][i][2];
        if (magn>max) {
            max=magn;
        }
    }
    return max;
}

void PlotCPM::CalculateStrainField(NOD *pn) const{
    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {
            int v=vx + vy*par.NVX;

            double estrains[3],L1,L2,v1[2],v2[2];

            get_estrains(pn,v,estrains);

            L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);

            if(L1>L2)
            {

                // calculate strain vector
			
                strain[vx][vy][0] = L1 * v1[0];
                strain[vx][vy][1] = L1 * v1[1];
            }
            else
            {

                strain[vx][vy][0] = L2 * v2[0];
                strain[vx][vy][1] = L2 * v2[1];

            }
            // magnitude

            double strx=strain[vx][vy][0], stry=strain[vx][vy][1];
            strain[vx][vy][2]=sqrt(strx*strx + stry*stry);
	}
   }
}

void PlotCPM::PlotNodalForces(NOD *pn){


BeginScene();

    // draw force lines
    for (int nx=1;nx<NNX-1;nx++) { 
        for (int ny=1;ny<NNY-1;ny++) {

            double linelength = sqrt(2)*(par.PIXPERVOX)/MaxForceMagnitude(); 
            double x1= (par.PIXPERVOX)*nx;
            double y1= (par.PIXPERVOX)*ny;
            double x2= ((par.PIXPERVOX)*nx+linelength*force[nx][ny][0]);
            double y2= ((par.PIXPERVOX)*ny+linelength*force[nx][ny][1]);

            if (x1<(par.PIXPERVOX)) x1=(par.PIXPERVOX);
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<(par.PIXPERVOX)) y1=(par.PIXPERVOX);
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<(par.PIXPERVOX)) x2=(par.PIXPERVOX);
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<(par.PIXPERVOX)) y2=(par.PIXPERVOX);
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            int n=nx + ny*NNX;

if(!(x1==x2) || !(y1==y2)){

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,2);
			Line(x1-l,y1,x2-l,y2,2);
		}
		
	int diffx = x2-x1;
	int diffy = y2-y1;
	double length = sqrt(diffx*diffx+diffy*diffy);

	if(length>0)
	{
		double alpha = acos(diffx/length);
		double gamma = (15.0/180)*3.1416;
		double beta = 3.1416-gamma;
		double lengtha = length/5;
		double lengthv = (length-lengtha*cos(gamma))/cos(beta);
		if(lengthv<0){lengthv=-lengthv;}
		if(y2>y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 - lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 - lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,2);
		Line(x2-l,y2,x4-l,y4,2);}}
		if(y2<y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 + lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 + lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,2);
		Line(x2-l,y2,x4-l,y4,2);}}
	}
	
        }
    }
}


    EndScene();


}

void PlotCPM::PlotNodeConnection(NOD *pn, int *pathx, int* pathy){

    BeginScene();

		//for loop plot the path from node to node calculated in checknodeconnection
		int colour = 5;
		int length = par.NVX;
		for(int i = 0; i < length ; i++)
		{
			int vx = pathx[i]; int vy = pathy[i];
			if(vx>0&&vy>0)
			{
				    for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
					for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
				    		Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
					}
				    }
			}
		}
EndScene();

}

void PlotCPM::StrainColorBar(void){
BeginScene();
QColor c;

double strainhue=0;
double max_strain;

	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBAR;
		if(par.MAXCOLORBAR==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

double cstep = max_strain/(par.NVY*par.PIXPERVOX);

cout << "max_strain" << max_strain << endl;

for(int i = 0; i < par.NVY*par.PIXPERVOX; i++)
{
	    strainhue = i*cstep;
            int hue=240-240*strainhue/max_strain;
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	
        picture->setPen( c );
	
	double x = par.NVX*par.PIXPERVOX + par.WIDTHCOLORBAR;
	double y = par.NVY*par.PIXPERVOX-i;
	for(int pixx = 0; pixx < par.WIDTHCOLORBAR; pixx++)
	{
        	picture->drawPoint(x+pixx,y);
	}
    if(!(i%par.NVY)&&(i>0))
    {
        double val = strainhue;
        if(val>max_strain){val=par.MAXCOLORBAR;}
        DrawColorBarLabel(i,val);
    }

}
EndScene();

}

void PlotCPM::DrawColorBarLabel(int yval, double val){
string valstr;
ostringstream convert;
convert.precision(4);
convert << fixed;
convert << val;
valstr = convert.str();
const char * c = valstr.c_str();
const QString & s = QString(c);
picture->setPen(QColor::fromRgb(0,0,0));
QFont* font = new QFont("Arial");
font->setPixelSize(100*par.NVX/300);
font->setBold(true);
picture->setFont(*font);
this->picture->drawText(QPointF(par.NVX*par.PIXPERVOX+2*par.WIDTHCOLORBAR,par.NVY*par.PIXPERVOX-yval), s);

}







 
