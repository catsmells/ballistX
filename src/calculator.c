#include <stdio.h>
#include "../include/balstats.h"
#include <stdlib.h>
#include <math.h>
#define TOO_LOW 1
#define TOO_HIGH 2
#define PBR_ERROR 3
double degmoa(double deg){
	return(deg*60);}
double degrad (double deg){
	return(deg*M_PI/180);}
double moadeg(double moa){
	return(moa/60);}
double moarad(double moa){
	return(moa/60*M_PI/180);}
double raddeg(double rad){
	return(rad*180/M_PI);}
double radmoa(double rad){
	return(rad*60*180/M_PI);}
double calcfr(double temp,double psi,double hg){
	double vpw=4e-6*pow(temp,3)-0.0004*pow(temp,2)+0.0234*temp-0.2517;
	double frh=0.995*(psi/(psi-(0.3783)*(hg)*vpw));
	return(frh);}
double calcfp(double psi){
	double pp=29.53;
	double fp=0;
	fp=(psi-pp)/(pp);
	return(fp);}
double calcft(double temp,double altx){
	double tssx=-0.0036*altx+59;
	double ft=(temp-tssx)/(459.6+tssx);
	return(ft);}
double calcfa(double altx){
	double fa=0;
	fa=-4e-15*pow(altx,3)+4e-10*pow(altx,2)-3e-5*altx+1;
	return(1/fa);}
double acxxxx(double Γ,double altx,double hgmarker,double temp,double hg){
	double fa=calcfa(altx);
	double ft=calcft(temp,altx);
	double fr=calcfr(temp,hgmarker,hg);
	double fp=calcfp(hgmarker);
	double cd=(fa*(1+ft-fp)*fr);
	return(Γ*cd);}
double rttxa (int dfcalc,double Γ,double vy){
	double vp=vy;
	double val=-1;
	double A=-1;
	double M=-1;
	switch(dfcalc){
		case G0:
			if(vp>2600){
				A=1.5366e-03;
				M=1.67;}
			else if(vp>2000){
				A=5.8497e-03;
				M=1.5;}
			else if(vp>1460){
				A=5.9814e-04;
				M=1.8;}
			else if(vp>1190){
				A=9.5408e-08;
				M=3.0;}
			else if(vp>1040){
				A=2.3385e-18;
				M=6.45;}
			else if(vp>840){
				A=5.9939e-08;
				M=3.0;}
			else if(vp>0){
				A=7.4422e-04;
				M=1.6;}break;
		case G1:
			if(vp>4230){
				A=1.477404177730177e-04;
				M=1.9565;}
			else if(vp>3680){
				A=1.920339268755614e-04;
				M=1.925;}
			else if(vp>3450){
				A=2.894751026819746e-04;
				M=1.875;}
			else if(vp>3295){
				A=4.349905111115636e-04;
				M=1.825;}
			else if(vp>3130){
				A=6.520421871892662e-04;
				M=1.775;}
			else if(vp>2960){
				A=9.748073694078696e-04;
				M=1.725;}
			else if(vp>2830){
				A=1.453721560187286e-03;
				M=1.675;}
			else if(vp>2680){
				A=2.162887202930376e-03;
				M=1.625;}
			else if(vp>2460){
				A=3.209559783129881e-03;
				M=1.575;}
			else if(vp>2225){
				A=3.904368218691249e-03;
				M=1.55;}
			else if(vp>2015){
				A=3.222942271262336e-03;
				M=1.575;}
			else if(vp>1890){
				A=2.203329542297809e-03;
				M=1.625;}
			else if(vp>1810){
				A=1.511001028891904e-03;
				M=1.675;}
			else if(vp>1730){
				A=8.609957592468259e-04;
				M=1.75;}
			else if(vp>1595){
				A=4.086146797305117e-04;
				M=1.85;}
			else if(vp>1520){
				A=1.954473210037398e-04;
				M=1.95;}
			else if(vp>1420){
				A=5.431896266462351e-05;
				M=2.125;}
			else if(vp>1360){
				A=8.847742581674416e-06;
				M=2.375;}
			else if(vp>1315){
				A=1.456922328720298e-06;
				M=2.625;}
			else if(vp>1280){
				A=2.419485191895565e-07;
				M=2.875;}
			else if(vp>1220){
				A=1.657956321067612e-08;
				M=3.25;}
			else if(vp>1185){
				A=4.745469537157371e-10;
				M=3.75;}
			else if(vp>1150){
				A=1.379746590025088e-11;
				M=4.25;}
			else if(vp>1100){
				A=4.070157961147882e-13;
				M=4.75;}
			else if(vp>1060){
				A=2.938236954847331e-14;
				M=5.125;}
			else if(vp>1025){
				A=1.228597370774746e-14;
				M=5.25;}
			else if(vp>980){
				A=2.916938264100495e-14;
				M=5.125;}
			else if(vp>945){
				A=3.855099424807451e-13;
				M=4.75;}
			else if(vp>905){
				A=1.185097045689854e-11;
				M=4.25;}
			else if(vp>860){
				A=3.566129470974951e-10;
				M=3.75;}
			else if (vp>810){
				A=1.045513263966272e-08;
				M=3.25;}
			else if(vp>780){
				A=1.291159200846216e-07;
				M=2.875;}
			else if(vp>750){
				A=6.824429329105383e-07;
				M=2.625;}
			else if(vp>700){
				A=3.569169672385163e-06;
				M=2.375;}
			else if(vp>640){
				A=1.839015095899579e-05;
				M=2.125;}
			else if(vp>600){
				A=5.71117468873424e-05;
				M=1.950;}
			else if(vp>550){
				A=9.226557091973427e-05;
				M=1.875;}
			else if(vp>250){
				A=9.337991957131389e-05;
				M=1.875;}
			else if(vp>100){
				A=7.225247327590413e-05;
				M=1.925;}
			else if(vp>65){
				A=5.792684957074546e-05;
				M=1.975;}
			else if(vp>0){
				A=5.206214107320588e-05;
				M=2.000;}break;
		case G2:
			if(vp>1674){
				A=.0079470052136733;
				M=1.36999902851493;}
			else if(vp>1172){
				A=1.00419763721974e-03;
				M=1.65392237010294;}
			else if(vp>1060){
				A=7.15571228255369e-23;
				M=7.91913562392361;}
			else if(vp>949){
				A=1.39589807205091e-10;
				M=3.81439537623717;}
			else if(vp>670){
				A=2.34364342818625e-04;
				M=1.71869536324748;}
			else if(vp>335){
				A=1.77962438921838e-04;
				M=1.76877550388679;}
			else if(vp>0){
				A=5.18033561289704e-05;
				M=1.98160270524632;}break;
		case G5:
			if(vp>1730){
				A=7.24854775171929e-03;
				M=1.41538574492812;}
			else if(vp>1228){
				A=3.50563361516117e-05;
				M=2.13077307854948;}
			else if(vp>1116){
				A=1.84029481181151e-13;
				M=4.81927320350395;}
			else if(vp>1004){
				A=1.34713064017409e-22;
				M=7.8100555281422;}
			else if(vp>837){
				A=1.03965974081168e-07;
				M=2.84204791809926;}
			else if(vp>335){
				A=1.09301593869823e-04;
				M=1.81096361579504;}
			else if(vp>0){
				A=3.51963178524273e-05;
				M=2.00477856801111;}break;
		case G6:
			if(vp>3236){
				A=0.0455384883480781;
				M=1.15997674041274;}
			else if(vp>2065){
				A=7.167261849653769e-02;
				M=1.10704436538885;}
			else if(vp>1311){
				A=1.66676386084348e-03;
				M=1.60085100195952;}
			else if(vp>1144){
				A=1.01482730119215e-07;
				M=2.9569674731838;}
			else if(vp>1004){
				A=4.31542773103552e-18;
				M=6.34106317069757;}
			else if(vp>670){
				A=2.04835650496866e-05;
				M=2.11688446325998;}
			else if(vp>0){
				A=7.50912466084823e-05;
				M=1.92031057847052;}break;
		case G7:
			if(vp>4200){
				A=1.29081656775919e-09;
				M=3.24121295355962;}
			else if(vp>3000){
				A=0.0171422231434847;
				M=1.27907168025204;}
			else if(vp>1470){
				A=2.33355948302505e-03;
				M=1.52693913274526;}
			else if(vp>1265){
				A=7.97592111627665e-04;
				M=1.67688974440324;}
			else if(vp>1110){
				A=5.71086414289273e-12;
				M=4.3212826264889;}
			else if(vp>960){
				A=3.02865108244904e-17;
				M=5.99074203776707;}
			else if(vp>670){
				A=7.52285155782535e-06;
				M=2.1738019851075;}
			else if(vp>540){
				A=1.31766281225189e-05;
				M=2.08774690257991;}
			else if(vp>0){
				A=1.34504843776525e-05;
				M=2.08702306738884;}break;
		case G8:
			if(vp>3571){
				A=.0112263766252305;
				M=1.33207346655961;}
			else if(vp>1841){
				A=.0167252613732636;
				M=1.28662041261785;}
			else if(vp>1120){
				A=2.20172456619625e-03;
				M=1.55636358091189;}
			else if(vp>1088){
				A=2.0538037167098e-16;
				M=5.80410776994789;}
			else if(vp>976){
				A=5.92182174254121e-12;
				M=4.29275576134191;}
			else if(vp>0){
				A=4.3917343795117e-05;
				M=1.99978116283334;}break;
		default:
			break;}
	if(A!=-1&&M!=-1&&vp>0&&vp<10000){
		val=A*pow(vp,M)/Γ;
		return(val);}
	else
		return(-1);}
double grxxx(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx]);}
	else
		return(0);}
double gptxx(double *sln, int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+1]);}
	else
		return(0);}
double gmnxx(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+2]);}
	else
		return(0);}
double gtdxx(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+3]);}
	else
		return(0);}
double gwdxx(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+4]);}
	else
		return(0);}
double gwdxxMOA(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+5]);}
	else
		return(0);}
double Getvy(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+6]);}
	else
		return(0);}
double GetVx(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+7]);}
	else
		return(0);}
double GetVy(double *sln,int ydxxx){
	double size=sln[__BCMXR__*10+1];
	if(ydxxx<size){
		return(sln[10*ydxxx+8]);}
	else
		return(0);}
double wndgtw(double wndspdtw,double Vi,double xx,double t){
	double Vw=wndspdtw*17.60;
	return(Vw*(t-xx/Vi));}
double hwrstxxc(double wndspdtw,double wndngltw){
	double xwng=dgtrdxx(wndngltw);
	return(cos(xwng)*wndspdtw);}
double cwwx(double wndspdtw,double wndngltw){
	double xwng=dgtrdxx(wndngltw);
	return(sin(xwng)*wndspdtw);}
double zrnxgl(int dfcalc,double Γ,double Vi,double wwshlm,double zrnxglb,double ymmrtq){
	double t =0;
	double dt=1/Vi;
	double y=-wwshlm/12;
	double x=0;
	double da;
	double v=0,vx=0,vy=0;
	double vx1=0,vy1=0;
	double dv=0,dvx=0,dvy=0;
	double Gx=0,Gy=0;
	double bqbqbq=0;
	int quit=0;
	da=dgtrdxx(14);
	for(bqbqbq=0;quit==0;bqbqbq=bqbqbq+da){
		vy=Vi*sin(bqbqbq);
		vx=Vi*cos(bqbqbq);
		Gx=κ*sin(bqbqbq);
		Gy=κ*cos(bqbqbq);
		for(t=0,x=0,y=-wwshlm/12;x<=zrnxglb*3;t=t+dt){
			vy1=vy;
			vx1=vx;
			v=pow((pow(vx,2)+pow(vy,2)),0.5);
			dt=1/v;
			dv=rttxa(dfcalc,Γ,v);
			dvy=-dv*vy/v*dt;
			dvx=-dv*vx/v*dt;
			vx=vx+dvx;
			vy=vy+dvy;
			vy=vy+dt*Gy;
			vx=vx+dt*Gx;
			x=x+dt*(vx + vx1)/2;
			y=y+dt*(vy + vy1)/2;
			if(vy<0&&y<ymmrtq){
				break;}
			if(vy>3*vx){
				break;}}
		if(y>ymmrtq&&da>0){
			da=-da/2;}
		if(y<ymmrtq&&da<0){
			da=-da/2;}
		if(fabs(da)<matrdxx(0.01)){
			quit=1;}
		if(bqbqbq>dgtrdxx(45)){
			quit=1;}}
	return(sdfwfkkk(bqbqbq));}
int sxpxxxxm(int dfcalc,double Γ,double Vi,double wwshlm,double sangpxxc,double sxngxc,double wndspdtw,double wndngltw,double **sxxxxxm){
	double *ptr;
	double t=0;
	double dt=0.5/Vi;
	double v=0;
	double vx=0,vx1=0,vy=0,vy1=0;
	double dv=0,dvx=0,dvy=0;
	double x=0,y=0;
	double sxngxchw=hwrstxxc(wndspdtw,wndngltw);
	double bqbqbdq=cwwx(wndspdtw,wndngltw);
	double Gy=κ*cos(dgtrdxx((sangpxxc+sxngxc)));
	double Gx=κ*sin(dgtrdxx((sangpxxc+sxngxc)));
	int n=0;
	ptr=(double *)malloc(10*__BCMXR__*sizeof(double)+2048);
	vx=Vi*cos(dgtrdxx(sxngxc));
	vy=Vi*sin(dgtrdxx(sxngxc));
	y=-wwshlm/12;
	for(t=0;;t=t+dt){
		vx1=vx,vy1=vy;
		v=pow(pow(vx,2)+pow(vy,2),0.5);
		dt=0.5/v;
		dv=rttxa(dfcalc,Γ,v+sxngxchw);
		dvx=-(vx/v)*dv;
		dvy=-(vy/v)*dv;
		vx=vx+dt*dvx+dt*Gx;
		vy=vy+dt*dvy+dt*Gy;
		if(x/3>=n){
			ptr[10*n+0]=x/3;
			ptr[10*n+1]=y*12;
			ptr[10*n+2]=-drdmwqqc(atan(y/x));
			ptr[10*n+3]=t+dt;
			ptr[10*n+4]=wndgtw(bqbqbdq,Vi,x,t+dt);
			ptr[10*n+5]=drdmwqqc(atan((ptr[10*n+4]/12)/(ptr[10*n+0]*3)));
			ptr[10*n+6]=v;
			ptr[10*n+7]=vx;
			ptr[10*n+8]=vy;
			ptr[10*n+9]=0;
			n++;}
		x=x+dt*(vx+vx1)/2;
		y=y+dt*(vy+vy1)/2;
		if(fabs(vy)>fabs(3*vx))break;
		if(n>=__BCMXR__+1)break;}
	ptr[10*__BCMXR__+1]=(double)n;
	*sxxxxxm=ptr;
	return(n);}
int pbr(int dfcalc,double Γ,double Vi,double wwshlm,double ddvsrn,int ddsinyen,int *oresult){
	double t=0;
	double dt=0.5/Vi;
	double v=0;
	double vx=0,vx1=0,vy=0,vy1=0;
	double dv=0,dvx=0,dvy=0;
	double x=0,y=0;
	double sangpxxc=0;
	double sxngxc=0;
	double fyrstp=10;
	int quit=0;
	double fyrzrop=-1;
	double frzx=0;
	int fvzx=0;
	double fyrxxx=0;
	double mnpbrx=0;
	int mnpkp=0;
	double mxpbrx=0;
	int mxpkp=0;
	int tnxxr=0;
	double Gy=κ*cos(dgtrdxx((sangpxxc+sxngxc)));
	double Gx=κ*sin(dgtrdxx((sangpxxc+sxngxc)));
	while(quit==0){
		int keep=0;
		int keep2=0;
		int tinkeep=0;
		int n=0;
		Gy=κ*cos(dgtrdxx((sangpxxc+sxngxc)));
		Gx=κ*sin(dgtrdxx((sangpxxc+sxngxc)));
		vx=Vi*cos(dgtrdxx(sxngxc));
		vy=Vi*sin(dgtrdxx(sxngxc));
		y=-wwshlm/12;
		x=0;
		y=-wwshlm/12;
		mnpkp=0;
		mxpkp=0;
		fvzx=0;
		tnxxr=0;
		tinkeep=0;
		for(t=0;;t=t+dt){
			vx1=vx,vy1=vy;
			v=pow(pow(vx,2)+pow(vy,2),0.5);
			dt=0.5/v;
			dv=rttxa(dfcalc,Γ,v);
			dvx=-(vx/v)*dv;
			dvy=-(vy/v)*dv;
			vx=vx+dt*dvx+dt*Gx;
			vy=vy+dt*dvy+dt*Gy;
			x=x+dt*(vx+vx1)/2;
			y=y+dt*(vy+vy1)/2;
			if(y>0&&keep==0&&vy>=0){
				fyrzrop=x;
				keep=1;}
			if(y<0&&keep2==0&&vy<=0){
				frzx=x;
				keep2=1;}
			if((12*y)>-(ddvsrn/2)&&mnpkp==0){
				mnpbrx=x;
				mnpkp=1;}
			if((12*y)<-(ddvsrn/2)&&mnpkp==1&&mxpkp==0){
				mxpbrx=x;
				mxpkp=1;}
			if(x>=(ddsinyen*3)&&tinkeep==0){
				tnxxr=(int)((float)100*(float)y*(float)12);
				tinkeep=1;}
			if(fabs(vy)>fabs(3*vx)){break;}
			if(n>=__BCMXR__+1){break;}
			if(vy<0&&fvzx==0){
				fyrxxx=y;
				fvzx=1;}
			if(keep==1&&keep2==1&&mnpkp==1&&mxpkp==1&&fvzx==1&&tinkeep==1){break;}}
		if((fyrxxx*12)>(ddvsrn/2)){
			if(fyrstp>0)
				fyrstp=-fyrstp/2;}
		else if((fyrxxx*12)<=(ddvsrn/2)){
			if(fyrstp<0)
				fyrstp=-fyrstp/2;}
		sxngxc+=fyrstp;if(fabs(fyrstp)<(0.01/60))
			quit=1;}
	oresult[0]=(int)(fyrzrop/3);
	oresult[1]=(int)(frzx/3);
	oresult[2]=(int)(mnpbrx/3);
	oresult[3]=(int)(mxpbrx/3);
	oresult[4]=(int)tnxxr;
	return(0);}
