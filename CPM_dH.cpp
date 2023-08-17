// file: cpm_dh.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
double calcdH(VOX* pv, NOD* pn, int* csize, int xt, int xs, int pick, int ttag, int stag, int incr)
{
	double dH, dHcontact, dHvol, dHstr;


	dHcontact = 0;
	dHcontact = calcdHcontact(pv,xt,ttag,stag);
    //cerr << "dHcontact = " << dHcontact << endl;


	dHvol = 0;
	dHvol = calcdHvol(csize,ttag,stag);
    //cerr << "dHvol = " << dHvol << endl;

	dHstr = 0;// DAPHNE: kunnen deze 3 regels niet beter in de if statement hieronder staan?
	dHstr = calcdHstrain(pn,xt,xs,pick,ttag,stag);
    //cerr << "dHstr = " << dHstr << endl;

	dH = dHcontact+dHvol;
	if(par.DUROTAXIS)
	{
		dH+=dHstr;
	}
    

/*
	double q,vmx[8],vmy[8],vm[2];
	q = SQ05;
	vmx[0]=-q; vmx[1]= 0; vmx[2]= q;
	vmx[7]=-1;            vmx[3]= 1;
	vmx[6]=-q; vmx[5]= 0; vmx[4]= q;

	vmy[0]= q; vmy[1]= 1; vmy[2]= q;
	vmy[7]= 0;            vmy[3]= 0;
	vmy[6]=-q; vmy[5]=-1; vmy[4]=-q;

	vm[0] = vmx[pick];
	vm[1] = vmy[pick];
	if(incr>16 && (vm[0]==0 | vm[1]==0) && stag){
	cout << "dHcontact " << dHcontact << endl;
	cout << "dHvol " << dHvol << endl;
	cout << "dHstr " << dHstr << endl;}
*/

	//cout << "dHcontact " << dHcontact << endl;

	return dH;


}




////////////////////////////////////////////////////////////////////////////////
double calcdHcontact(VOX* pv, int xt, int ttag, int stag)
{
	double dHcontact, Hcontact, Hcontactn;
	int nbs[8],n,nbtag;

	Hcontact=0; Hcontactn=0;
	nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
	nbs[7]=xt-1;                    nbs[3]=xt+1;
	nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);
	for(n=0;n<8;n++)
	{
		nbtag = pv[nbs[n]].ctag;
		Hcontact += contactenergy(ttag,nbtag);
		Hcontactn += contactenergy(stag,nbtag);
	}
	dHcontact = Hcontactn-Hcontact;


	return dHcontact;
}

////////////////////////////////////////////////////////////////////////////////
double contactenergy(int tag1, int tag2)
{
	double J;

	J = 0;
	if(tag1!=tag2)
	{
    	if((tag1==0)||(tag2==0))
        	J = JCM;
    	else
        	J = JCC;
	}
	return J;
}





////////////////////////////////////////////////////////////////////////////////
double calcdHvol(int* csize, int ttag, int stag)
{
	double dHvol,dHvolA,dHvolB,V0,eV,eVn;
	int V;

	// assume 2 cells, A (ttag) and B (stag) are involved
	dHvolA=0; dHvolB=0; V0=par.TARGETVOLUME;
	if(ttag) // cell ttag retracts
	{
		V=csize[ttag-1]; eV=(V-V0)/V0; eVn=(V-1-V0)/V0;
		dHvolA = (par.INELASTICITY)*(eVn*eVn-eV*eV);
	}
	if(stag) // cell stag expands
	{
		V=csize[stag-1]; eV=(V-V0)/V0; eVn=(V+1-V0)/V0;
		dHvolB = (par.INELASTICITY)*(eVn*eVn-eV*eV);

	}
	dHvol = dHvolA+dHvolB;

	return dHvol;
}



////////////////////////////////////////////////////////////////////////////////
double calcdHstrain(NOD* pn, int xt, int xs, int pick, int ttag, int stag)
{
	double dHstrain;
	double q,vmx[8],vmy[8],vm[2];
	double estrains[3],L1,L2,v1[2],v2[2];
	double vmv1,vmv2;
	double E1,E2;

	dHstrain = 0;

	// unitvectors for move: vm
	q = SQ05;
	vmx[0]=-q; vmx[1]= 0; vmx[2]= q;
	vmx[7]=-1;            vmx[3]= 1;
	vmx[6]=-q; vmx[5]= 0; vmx[4]= q;

	vmy[0]= q; vmy[1]= 1; vmy[2]= q;
	vmy[7]= 0;            vmy[3]= 0;
	vmy[6]=-q; vmy[5]=-1; vmy[4]=-q;

	vm[0] = vmx[pick];
	vm[1] = vmy[pick];


	if(stag) // expansion
	{
		get_estrains(pn,xt,estrains);

//cout << "xs " << xs%par.NVX << " " << xs/par.NVY << endl;
//cout << "xt " << xt%par.NVX << " " << xt/par.NVY << endl;
//cout << "estrainsx  " << estrains[0] << endl;
//cout << "estrainsy  " << estrains[1] << endl;
//cout << "estrainsxy " << estrains[2] << endl << endl;
		
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);
		//if(L1<0 && L2<0){
		//cout << endl <<"both strains negative" << endl;
		//cout << "L1" << L1 << " L2 " << L2 << endl;}

		// inproducts of move vector with pr. strain vectors
		vmv1 = vm[0]*v1[0]+vm[1]*v1[1];
		vmv2 = vm[0]*v2[0]+vm[1]*v2[1];

		//dHstrain -= sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
		//dHstrain -= sige((par.YOUNGS));

		if(!par.COMPRESSINGSTIFF){
			E1 = (par.YOUNGS); if(L1>0) {E1*=(1+L1/(par.STIFFENINGSTIFF));}
			E2 = (par.YOUNGS); if(L2>0) {E2*=(1+L2/(par.STIFFENINGSTIFF));}
			dHstrain -= sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;


//cout << "E1 " << E1 << endl;
//cout << "E2 " << E2 << endl;
		}
		if(par.COMPRESSINGSTIFF){
			E1 = (par.YOUNGS); if(L1>0) {E1*=(1+L1/(par.STIFFENINGSTIFF));}
			E2 = (par.YOUNGS); if(L2>0) {E2*=(1+L2/(par.STIFFENINGSTIFF));}
			if(L1<0) {E1*=(1-L1/(par.STIFFENINGSTIFF));}
			if(L2<0) {E2*=(1-L2/(par.STIFFENINGSTIFF));}
			dHstrain -= sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;
		}


	}



	if(ttag) // retraction
	{
		get_estrains(pn,xt,estrains);
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);

		// inproducts of move vector with pr. strain vectors
		vmv1 = vm[0]*v1[0]+vm[1]*v1[1];
		vmv2 = vm[0]*v2[0]+vm[1]*v2[1];

		//dHstrain += sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
		//dHstrain += sige((par.YOUNGS));
		
		if(!par.COMPRESSINGSTIFF){
			E1 = (par.YOUNGS); if(L1>0) {E1*=(1+L1/(par.STIFFENINGSTIFF));}
			E2 = (par.YOUNGS); if(L2>0) {E2*=(1+L2/(par.STIFFENINGSTIFF));}
			dHstrain += sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;

//cout << "E1 " << E1 << endl;
//cout << "E2 " << E2 << endl;
		}
		if(par.COMPRESSINGSTIFF){
			E1 = (par.YOUNGS); if(L1>0) {E1*=(1+L1/(par.STIFFENINGSTIFF));}
			E2 = (par.YOUNGS); if(L2>0) {E2*=(1+L2/(par.STIFFENINGSTIFF));}
			if(L1<0) {E1*=(1-L1/(par.STIFFENINGSTIFF));}
			if(L2<0) {E2*=(1-L2/(par.STIFFENINGSTIFF));}
			dHstrain += sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;
		}




	}


	return dHstrain;
}



////////////////////////////////////////////////////////////////////////////////
double sige(double L)
{
	double sigL;
	if(par.WHICHSIGMF==1){
		double x;
		x = (par.STIFFSENSITIVITY)*(L - (par.THRESHOLDSTIFF));
		sigL = (par.LAMBDADUR)/(1+exp(-x)); //sigmoid function
	}

	if(par.WHICHSIGMF==2){
		double x;
		x = (par.STIFFSENSITIVITY)*(L - (par.THRESHOLDSTIFF));
		sigL = (par.LAMBDADUR)/2*(1+x/sqrt(1+x*x)); //new sigmoid function
	}
	
	if(par.WHICHSIGMF==3){
		double alpha = 1/par.THRESHOLDSTIFF;
		sigL = (par.LAMBDADUR)*alpha*L/(1+alpha*L);
	}
	return sigL;
}



/*
////////////////////////////////////////////////////////////////////////////////
double calcdHstrain(NOD* pn, int xt, int xs, int pick, int ttag, int stag)
{
	double dHstrain;
	double q,vmx8[8],vmy8[8],vmx,vmy;
	double estrains[3],exx,eyy,exy,em;

	dHstrain = 0;
	//if ((old==0)||(new==0))
	{
		// unitvectors for move: vm
		q = SQ05;
		vmx8[0]=-q; vmx8[1]= 0; vmx8[2]= q;
		vmx8[7]=-1;             vmx8[3]= 1;
		vmx8[6]=-q; vmx8[5]= 0; vmx8[4]= q;

		vmy8[0]= q; vmy8[1]= 1; vmy8[2]= q;
		vmy8[7]= 0;             vmy8[3]= 0;
		vmy8[6]=-q; vmy8[5]=-1; vmy8[4]=-q;

		vmx = vmx8[pick];
		vmy = vmy8[pick];

		if(stag) // expansion
		{
			get_estrains(pn,xt,estrains);
			exx=estrains[0]; eyy=estrains[1]; exy=.5*estrains[2];
			em = vmx*(exx*vmx+exy*vmy)+vmy*(exy*vmx+eyy*vmy);
			dHstrain -= sige(em);
		}
		if(ttag) // retraction
		{
			get_estrains(pn,xs,estrains);
			exx=estrains[0]; eyy=estrains[1]; exy=.5*estrains[2];
			em = vmx*(exx*vmx+exy*vmy)+vmy*(exy*vmx+eyy*vmy);
			dHstrain += sige(em);
		}
	}
	return dHstrain;
}
*/


