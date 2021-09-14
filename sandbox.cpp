#include "myparameters.h"
extern "C" {
#include <stdlib.h>
}
#include "functions.h"
#include "plotcpm.h"
#include <QApplication>
#include <QMainWindow>
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <string>
#include <fstream>
#include <sstream>


Parameter par;

int main(int argc, char *argv[])
{
    
    try {
        if (argc <= 1)
        {
            cout << "Usage: " << argv[0] << " <ParametersDirectory>" << endl;
            exit(1);
        }
        int nrrdof, d;
        int *dofpos;
        double *f, *u;
        double **klocal;
        int *kcol;
        double *kval;
        NOD *pn;
        VOX *pv;
        int NRc,c,v;
        int *csize,*csumx,*csumy;
        int *cper;
        double *clength, *ecc, *cangle;
        int incr, startincr;
        int *cmxi, *cmyi;
        double* sqdis;
        int contact = 0;
        
        //FOR SAVING FILES
        stringstream currentdir;
        string thisdir;
        thisdir=argv[1];
        currentdir << thisdir;
        string scurrentdir = currentdir.str();
        const char * ccurrentdir = scurrentdir.c_str();
        
        stringstream lstr;
        string lstr1 = scurrentdir;
        string lstr2;
        lstr2 = "/length.txt";
        lstr<<lstr1<<lstr2;
        string ls = lstr.str();
        const char * lfcstr = ls.c_str();
        FILE *ofpl;
        ofpl = fopen(lfcstr,"w");
        
        stringstream pastr;
        string pastr1 = scurrentdir;
        string pastr2;
        pastr2 = "/ratiopa.txt";
        pastr<<pastr1<<pastr2;
        string pas = pastr.str();
        const char * pafcstr = pas.c_str();
        FILE *ofppa;
        ofppa = fopen(pafcstr,"w");
        
        stringstream astr;
        string astr1 = scurrentdir;
        string astr2;
        astr2 = "/area.txt";
        astr<<astr1<<astr2;
        string as = astr.str();
        const char * afcstr = as.c_str();
        FILE *ofpa;
        ofpa = fopen(afcstr,"w");
        
        stringstream sdstr;
        string sdstr1 = scurrentdir;
        string sdstr2;
        sdstr2 = "/sqdis.txt";
        sdstr<<sdstr1<<sdstr2;
        string sds = sdstr.str();
        const char * sdfcstr = sds.c_str();
        FILE *ofpsd;
        ofpsd = fopen(sdfcstr,"w");
        
        stringstream estr;
        string estr1 = scurrentdir;
        string estr2;
        estr2 = "/ecc.txt";
        estr<<estr1<<estr2;
        string es = estr.str();
        const char * efcstr = es.c_str();
        FILE *ofpe;
        ofpe = fopen(efcstr,"w");
        
        stringstream anstr;
        string anstr1 = scurrentdir;
        string anstr2;
        anstr2 = "/angle.txt";
        anstr<<anstr1<<anstr2;
        string ans = anstr.str();
        const char * anfcstr = ans.c_str();
        FILE *ofpan;
        ofpan = fopen(anfcstr,"w");
        
        
        stringstream tccstr;
        string tccstr1 = scurrentdir;
        string tccstr2;
        tccstr2 = "/twocellcontact.txt";
        tccstr<<tccstr1<<tccstr2;
        string tccs = tccstr.str();
        const char * tccfcstr = tccs.c_str();
        FILE *ofptcc;
        ofptcc = fopen(tccfcstr,"w");
        
        stringstream ssreadcells;
        string sreadcells1 = scurrentdir;
        string sreadcells2;
        sreadcells2 = "/ctags0.out";
        ssreadcells<<sreadcells1<<sreadcells2;
        string sreadcells = ssreadcells.str();
        
        stringstream sssigmas;
        string ssigmas1 = scurrentdir;
        string ssigmas2;
        ssigmas2 = "/sigmas/";
        sssigmas<<ssigmas1<<ssigmas2;
        string ssigmas = sssigmas.str();
        const char * sigmasdir = ssigmas.c_str();
        mkdir(sigmasdir,0777);
        
        
        
        
        
        stringstream str;
        string str1, str2;
        str1=argv[1];
        str2="/parameters.txt";
        str<<str1<<str2;
        string s = str.str();
        const char * cstr = s.c_str();
        //First argument is directory
        par.Read(cstr);
        
        /// GRAPHICS //
        
        
        
        
        QApplication a(argc, argv);
        
        
        
        
        setlocale(LC_NUMERIC,"C");
        
        
        
        if(!par.WIDTHCOLORBAR){par.WIDTHCOLORBAR=par.NVX*par.PIXPERVOX/10;}
        if(!par.COLORBAR){par.WIDTHCOLORBAR=0;}
        
        /// INITIALIZE ///
        if ((par.SEED)) { srand((par.SEED));
        }
        // weggelaten omdat deze niet werkt in Unix systemen
        // else { sranddev(); }
        
        mt_init();
        pv = init_voxels();
        pn = init_nodes();
        
        startincr = 0;
        if(startincr==0) {
            
            if(par.CELLCOL){
                NRc = init_cells(pv);
            }
            
            
            if(par.TWOCELL){
                // two cells in the middle
                int vx1=(par.NVX)/2-par.DISTWOCELLS, vx2=(par.NVX)/2+par.DISTWOCELLS;
                int vy=par.NVY/2;
                
                v = vx1 + vy*(par.NVX);
                pv[v].ctag=1;
                v = vx2 + vy*(par.NVX);
                pv[v].ctag=2;
                NRc=2;
            }
            
            
            if(par.ONECELL){
                // one cell in the middle
                int vx=(par.NVX)/2;
                int vy=par.NVY/2;
                
                v = vx + vy*(par.NVX);
                pv[v].ctag=1;
                NRc=1;
            }
            
            
            //USE /and* ....  *and/  to make a big section into commentary
            
            //write_cells(pv,0);
        }
        if(par.READCELLS) {NRc = read_cells(pv,startincr,sreadcells);}
        
        
        
        if(NRc==2){par.TWOCELL=true;}
        
        //write NRC
        stringstream nrcstr;
        string nrcstr1 = scurrentdir;
        string nrcstr2;
        nrcstr2 = "/nrc.txt";
        nrcstr<<nrcstr1<<nrcstr2;
        string nrcs = nrcstr.str();
        const char * nrcfcstr = nrcs.c_str();
        FILE *ofpnrc;
        ofpnrc = fopen(nrcfcstr,"w");
        fprintf(ofpnrc ,"%d ",NRc);
        fclose(ofpnrc);
        
        
        
        
        
        if(par.WRATIOPA){cper = new int[NRc];
            for(c=0;c<NRc;c++) {cper[c]=0;}}
        if(par.WLENGTH | par.WECC | par.WANGLE)
        {
            clength = new double[NRc];
            for(c=0;c<NRc;c++) {clength[c]=0;}
            ecc = new double[NRc];
            for(c=0;c<NRc;c++) {ecc[c]=0;}
            if(!par.WTWOCELLANGLECM){cangle = new double[NRc]; for(c=0;c<NRc;c++) {cangle[c]=0;}}
            if(par.TWOCELL & par.WTWOCELLANGLECM){cangle = new double[3];for(c=0;c<3;c++) {cangle[c]=0;}}
            
        }
        
        cout << "NRc" << NRc << endl;
        
        if(par.WSQDIS){sqdis = new double[NRc];
            for(c=0;c<NRc;c++) {sqdis[c]=0;}}
        csize = new int[NRc];
        for(c=0;c<NRc;c++) {csize[c]=0;}
        csumx = new int[NRc];
        for(c=0;c<NRc;c++) {csumx[c]=0;}
        csumy = new int[NRc];
        for(c=0;c<NRc;c++) {csumy[c]=0;}
        cmxi = new int[NRc];
        cmyi = new int[NRc];
        
        
        
        // set initial volumes and center of mass tags
        for (v=0;v<NV;v++) {
            int y = v/(par.NVX); int x = v%(par.NVX);
            
            if (pv[v].ctag) {
                
                csize[pv[v].ctag-1]++;
                csumx[pv[v].ctag-1]+=x;
                csumy[pv[v].ctag-1]+=y;
                
                
                
            }
        }
        
        
        
        //initial center of mass
        for(c=0;c<NRc;c++)
        {
            
            
            cmxi[c]=csumx[c]/csize[c];
            cmyi[c]=csumy[c]/csize[c];
            
            
        }
        
        
        // set initial volumes
        //for(v=0;v<NV;v++) {if(pv[v].ctag) {csize[pv[v].ctag-1]++;}}
        
        
        
        set_forces(pn);
        //dont set restrictions if static stress rene
        if(!par.GLOBALSTRAIN){set_restrictions(pn);}
        
        // local K matrix
        klocal = set_klocal();
        
        // global K matrix:
        kcol = new int[10*NDOF];
        kval = new double[10*NDOF];
        assembly(kcol,kval,klocal,pv);
        dofpos = new int[NDOF];
        nrrdof = arrange_dofpos(dofpos,pn);
        reduce_K(kcol,kval,dofpos);
        
        FILE *of=fopen("com.dat","w");
        
        /// START SIMULATION ///
        for(incr=startincr; incr<(par.NRINC); incr++)
        {
            printf("\nSTART INCREMENT %d",incr);
            //write_cells(pv,incr);
            
            
            
            cell_forces(pv,pn,csize,NRc);
            
            
            
            // FEA part // parts of this can go out the loop depending on what changes
            f = new double[nrrdof];
            place_node_forces_in_f(pn,f);
            u = new double[nrrdof];
            
            set_disp_of_prev_incr(pn,u);
            
            solvePCG(kcol,kval,u,f,nrrdof);
            
            disp_to_nodes(pn,u);
            
            
            
            //for(int i = 0; i < NN; i++){pn[i].ux = 1;pn[i].uy=3;}
            //free(u); free(f);
            delete [] u; delete [] f;
            
            stringstream fstr;
            string fstr2;
            fstr2="cpmfem%05d.png";
            fstr<<str1<<fstr2;
            string fs = fstr.str();
            const char * fcstr = fs.c_str();
            
            
            if (!(incr%(par.STRIDE))) {
                char *threrror;
                
                PlotCPM plot((par.NVX)*(par.PIXPERVOX)+3*par.WIDTHCOLORBAR,par.NVY*(par.PIXPERVOX));
                
                
                char fname[200];
                sprintf(fname,fcstr,incr);
                
                if(par.PRINCFIELD | par.STRAINFIELD){
                    plot.CalculateStrainField(pn);}
                
                plot.PlotHueStrainField(par.STRAINFIELD);
                
                if(par.PRINCFIELD){
                    plot.PlotPrincipleStrainField(pv);}
                
                plot.Plot(pv); //for art figure
                
                if(par.FORCEFIELD){
                    plot.CalculateForceField(pn);
                    plot.PlotNodalForces(pn);}
                
                if(par.COLORBAR){
                    plot.StrainColorBar();}
                
                
                
                cout << endl <<"Max Strain "<<plot.MaxStrainMagnitude() << endl;
                //for (int n=0;n<NN;n++){cout << "n " << n << " fx "<< pn[n].fx << endl;}
                //for (int n=0;n<NN;n++){cout << "n " << n << " fy "<< pn[n].fy << endl;}
                
                //for (int n=0;n<NN;n++){cout << "n " << n << " ux "<< pn[n].ux << endl;}
                //for (int n=0;n<NN;n++){cout << "n " << n << " uy "<< pn[n].uy << endl;}
                cout << fname << endl;
                plot.Write(fname);
                
                //write_forces(pn,incr);
                //write_disps(pn,incr);
                //write_pstrain(pv,pn,incr);
                
                
            }
            for (int c=0;c<NRc;c++) {
                fprintf(of, "%d %lf %lf\n", c, (double)csumx[c]/(double)csize[c], (double)csumy[c]/(double)csize[c]);
            }
            
            
            if (!(incr%(par.WSTRIDE)))
            {
                
                
                
                if(par.WRATIOPA)
                {
                    CalcPerimeters(pv,cper,NRc);
                    write_ratiopa(cper,csize,NRc,incr,ofppa);
                    
                }
                
                if(par.WLENGTH)
                {
                    CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
                    write_length(clength,NRc,incr,ofpl);
                }
                if(par.WECC)
                {
                    CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
                    write_eccentricity(ecc,NRc,incr,ofpe);
                }
                
                if(par.WANGLE)
                {
                    CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
                    write_cangle(cangle,NRc,incr,ofpan);
                }
                if(par.WAREA)
                {
                    write_area(csize,NRc,incr,ofpa);
                }
                if(par.WSIGMA)
                {
                    write_cells(pv,incr,scurrentdir);
                }
                if(par.WSQDIS)
                {
                    for(int c=0;c<NRc;c++)
                    {
                        sqdis[c] = (cmxi[c]-csumx[c]/csize[c])*(cmxi[c]-csumx[c]/csize[c])+(cmyi[c]-csumy[c]/csize[c])*(cmyi[c]-csumy[c]/csize[c]);
                    }
                    write_sqdis(sqdis,NRc,incr,ofpsd);
                }
                if(par.TWOCELL & par.WTWOCELLCONTACT)
                {
                    contact=check_contact(pv);
                    write_twocellcontact(contact,incr,ofptcc);
                }
            }
            
            CPM_moves(pv,pn,csize, csumx, csumy, incr);
            
            
            
        }
        fclose(of);
        /// END ///
        printf("\nSIMULATION FINISHED!");
        //free(pv); free(pn); free(klocal); free(kcol); free(kval); free(dofpos);
        delete [] pv; delete [] pn; delete [] klocal; delete [] kcol; delete [] kval; delete [] dofpos;
        delete [] csize; delete [] csumx; delete [] csumy;
        if(par.WRATIOPA){delete [] cper;} if(par.WLENGTH){delete [] clength;}
        fclose(ofpl); fclose(ofppa); fclose(ofpe); fclose(ofpan); fclose(ofpa); fclose(ofpsd); fclose(ofptcc);
    } catch (const char *threrror) {
        fprintf(stderr,"%s\n", threrror);
        exit(1);
    }
    
    return 0;
}


