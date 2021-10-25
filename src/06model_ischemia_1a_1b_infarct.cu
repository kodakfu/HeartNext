//MyEmail:xsjshr1108@163.com
#include"base.cuh"
#include<stdio.h>
const int ENDO = 1;//normal 1,2,3;infarct 4,5,6;ischemia1a 7,8,9;ischemia1b 10,11,12;
const int MCELL = 2;
const int EPI = 3;
const int ENDO_S = 4;
const int MCELL_S = 5;
const int EPI_S = 6;
const int ENDO_S_1A = 7;
const int MCELL_S_1A = 8;
const int EPI_S_1A = 9;
const int ENDO_S_1B = 10;
const int MCELL_S_1B = 11;
const int EPI_S_1B = 12;
__global__ void kernel06(double* Cai_buffer, double* CaSR_buffer, double* CaSS_buffer, double* Nai_buffer,
						 double* Ki_buffer,  double* sm_buffer, double* sh_buffer, double* sj_buffer,
						 double* sxr1_buffer, double* sxr2_buffer,double* sxs_buffer,double* sr_buffer,
						 double* ss_buffer, double* sd_buffer, double* sf_buffer, double* sf2_buffer,
						 double* sfcass_buffer, double* sOO_buffer, double* sRR_buffer, double* sml_buffer,
						 double* shl_buffer, double* volt_buffer, double* istim_buffer, double* itot_buffer,
						 int* type_buffer, int cellNumber)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index >= cellNumber) return;
	
	double &V = volt_buffer[index];
	if(V == 15.0) V = 15.001;
	double &Istim = istim_buffer[index];
	double &out_itot = itot_buffer[index];
	const int type = type_buffer[index];
	// cell param begin
	double &Cai =Cai_buffer[index] ; 
	double &CaSR=CaSR_buffer[index]; 
	double &CaSS=CaSS_buffer[index]; 
	double &Nai =Nai_buffer[index] ; 
	double &Ki  =Ki_buffer[index]  ; 

	double &sm   =sm_buffer[index]   ; 
	double &sh   =sh_buffer[index]   ; 
	double &sj   =sj_buffer[index]   ; 

	double &sxr1 =sxr1_buffer[index] ; 
	double &sxr2 =sxr2_buffer[index] ; 
							    
	double &sxs  =sxs_buffer[index]  ; 
	double &sr   =sr_buffer[index]   ; 
	double &ss   =ss_buffer[index]   ; 
							    
	double &sd   =sd_buffer[index]   ; 
	double &sf   =sf_buffer[index]   ; 
	double &sf2  =sf2_buffer[index]  ; 
	double &sfcass =sfcass_buffer[index] ;
	double &sOO  =sOO_buffer[index] ;
	double &sRR  =sRR_buffer[index] ;

	double &sml=sml_buffer[index];
	double &shl=shl_buffer[index];
	// cell param end

	double INa_l=0;  ///////////////////////////////////////////////
	double INaL=0;
	double Gna_l = 0.0065;//////////////////////////////////////////
	double alpham_nal =0;///////////////////////////////////////////
	double betam_nal=0;/////////////////////////////////////////////
	double minf_nal = 0;////////////////////////////////////////////
	double hinf_nal=0;//////////////////////////////////////////////
	double taum_nal = 0;////////////////////////////////////////////
	double tauh_nal=600.0;//////////////////////////////////////////
                
	const double Cao=2.0; ////luo
	const double Nao=140.0; ////luo
	const double Vc= 0.016404; ////luo-change
	const double Vsr=0.001094; ////luo
	const double Vss=0.00005468;////luo
		
	const double Bufc=0.2; ////luo
	const double Kbufc=0.001; ////luo
	const double Bufsr=10; ////luo
	const double Kbufsr=0.3; ////luo
	const double Bufss=0.4;////luo
	const double Kbufss=0.00025;////luo
	const double Vmaxup=0.006375;////luo
	const double Kup=0.00025;////luo
	const double Vrel=0.102;////luo
	const double k1_=0.15;////luo
	const double k2_=0.045;////luo
	const double k3=0.060;////luo
	const double k4=0.005;////luo
	const double EC=1.5;////luo
	const double maxsr=2.5;////luo
	const double minsr=1.;////luo
	const double Vleak=0.00036; ////luo
	const double Vxfer=0.0038;////luo
		
	const double R=8314.472;////luo
	const double F=96486.3415;////luo-change
	const double T=310.0;////luo
	const double RTONF=(R*T)/F;////luo
		    
	const double CAPACITANCE=0.185; ////luo
        double Gkr=0.153;
	const double pKNa=0.03; ////luo
	const double GK1=5.405; ////luo
        double GNa=14.838; 
	const double GbNa=0.00029; ////luo
	const double KmK=1.0; ////luo
	const double KmNa=40.0; ////luo
	const double knak=2.724; ////luo
	double GCaL=0.00003980;
	const double GbCa=0.000592;////luo  
	const double knaca=1000; ////luo
	const double KmNai=87.5;////luo
	const double KmCa=1.38; ////luo
	const double ksat=0.1;////luo
	const double n=0.35;////luo
	const double GpCa=0.1238; ////luo
	const double KpCa=0.0005;////luo
	const double GpK=0.0146;////luo

	const double  dt=0.02;

	const double inverseVcF2=1/(2*Vc*F);////luo
	const double inverseVcF=1./(Vc*F);////luo
	const double inversevssF2=1/(2*Vss*F);////luo

	double Ko=5.4;////luo
	if(type== ENDO_S_1A ||type== EPI_S_1A||type==MCELL_S_1A||type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B) Ko=8;//liang_change2
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) Ko=8;//liang_change2
	double Gks, Gto;////luo

	double IKr=0,IKs=0,IK1=0,Ito=0,INa=0,IbNa=0,ICaL=0,IbCa=0,INaCa=0,IpCa=0,IpK=0,INaK=0,Irel=0,Ileak=0,Iup=0,Ixfer=0,IKATP=0;////luo
	double k1=0,k2=0,kCaSR=0,dNai=0,dKi=0,dCai=0,dCaSR=0,dCaSS=0,dRR=0;////luo
	double Ek=0,Ena=0,Eks=0,Eca=0;////luo
	double CaCSQN=0,bjsr=0,cjsr=0,CaSSBuf=0,bcss=0,ccss=0,CaBuf=0,bc=0,cc=0,Ak1=0,Bk1=0,rec_iK1=0,rec_ipK=0,rec_iNaK=0,AM=0,BM=0,////luo
		AH_1=0,BH_1=0,AH_2=0,BH_2=0,AJ_1=0,BJ_1=0,AJ_2=0,BJ_2=0,M_INF=0,H_INF=0,J_INF=0,TAU_M=0,TAU_H=0,TAU_J=0,axr1=0,bxr1=0,axr2=0,bxr2=0,////luo
		Xr1_INF=0,Xr2_INF=0,TAU_Xr1=0,TAU_Xr2=0,Axs=0,Bxs=0,Xs_INF=0,TAU_Xs=0,R_INF=0,TAU_R=0,S_INF=0,TAU_S=0,Ad=0,Bd=0,Cd=0,Af=0,Bf=0,Cf=0,////luo
		Af2=0,Bf2=0,Cf2=0,TAU_D=0,D_INF=0,TAU_F=0,F_INF=0,TAU_F2=0,F2_INF=0,TAU_FCaSS=0,FCaSS_INF=0;////luo
	// calculate start
	Ek=RTONF*(log((Ko/Ki)));////luo
	Ena=RTONF*(log((Nao/Nai)));////luo
	Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));////luo
	Eca=0.5*RTONF*(log((Cao/Cai)));////luo

	Ak1=0.1/(1.+exp(0.06*(V-Ek-200)));////luo
	Bk1=(3.*exp(0.0002*(V-Ek+100))+exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));////luo

	rec_iK1=Ak1/(Ak1+Bk1);////luo
	rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));////luo
	rec_ipK=1./(1.+exp((25-V)/5.98));////luo

	if(type == EPI||type == ENDO||type == EPI_S||type == ENDO_S)   //////////////////////////////////////////////
		INa_l = Gna_l*sml*sml*sml*shl*(V-Ena);/////////////////////////////////////////////
	if(type == MCELL||type == MCELL_S){////////////////////////////////////////////////////////////////
		INa_l = 1.46*Gna_l*sml*sml*sml*shl*(V-Ena);////////////////////////////////////////
		//GCaL = GCaL * 6.0;
		//Gkr = Gkr * 0.2;
	}
	INa_l =0;

	double mL=0;
	double hL=1;
	double mLss=1.0/(1.0+exp((-(V+42.85))/5.264));
	double tmL=1.0/(6.765*exp((V+11.64)/34.77)+8.552*exp(-(V+77.42)/5.955));
	mL=mLss-(mLss-mL)*exp(-dt/tmL);
	double hLss=1.0/(1.0+exp((V+87.61)/7.488));
	double thL=200.0;
	hL=hLss-(hLss-hL)*exp(-dt/thL);
	double GNaL=0.0075;
	INaL=GNaL*(V-Ena)*mL*hL;
	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) INaL=INaL*1.5;
	
	double ATPi=4.6;              //liang_change
    	double GKATP=3.9;                                               //liang_change
   	double H=2.0;                                                   //liang_change
    	double n1=0.24;                                                  //liang_change
	double Khalf=0.25;  
        IKATP=GKATP*(1/(1+pow(ATPi/Khalf,H)))*pow(Ko/5.4,n1)*(V-Ek);//liang_change
	if(type== ENDO||type== EPI||type==MCELL) IKATP=0;//liang_change2




	INa=GNa*sm*sm*sm*sh*sj*(V-Ena);////luo
        	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) INa=INa*0.887;
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) INa=INa*0.38;
         
	ICaL=GCaL*sd*sf*sf2*sfcass*4*(V-15)*(F*F/(R*T))*(0.25*exp(2*(V-15)*F/(R*T))*CaSS-Cao)/(exp(2*(V-15)*F/(R*T))-1.);////luo
	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) ICaL=ICaL*0.8;
	if(type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B) ICaL=ICaL*0.5;
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) ICaL=ICaL*0.62;

	if(type== ENDO||type== ENDO_S||type== ENDO_S_1A||type== ENDO_S_1B) Gto=0.073;////luo
	if(type== MCELL||type== EPI||type== MCELL_S||type== EPI_S||type== MCELL_S_1A||type== EPI_S_1A||type== MCELL_S_1B||type== EPI_S_1B)  Gto=0.294;////luo
	Ito=Gto*sr*ss*(V-Ek);////luo
    	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) Ito=Ito*0.5;
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) Ito=Ito*0.37;

	IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);////luo
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) IKr=IKr*0.3;
	
        if(type== ENDO ||type== EPI||type== ENDO_S ||type== EPI_S||type== ENDO_S_1A ||type== EPI_S_1A||type== ENDO_S_1B ||type== EPI_S_1B) Gks=0.392;   //liang_change  //liang_change 
	if(type== MCELL||type== MCELL_S||type== MCELL_S_1A||type== MCELL_S_1B) Gks=0.098;   //liang_change  //liang_change 
	IKs=Gks*sxs*sxs*(V-Eks);////luo
	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) IKs=IKs*0.781;
	if(type== ENDO_S||type== EPI_S||type==MCELL_S) IKs=IKs*0.2;

	IK1=GK1*sqrt(Ko/5.4)*rec_iK1*(V-Ek);////luo
		
	INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*(1./(1+ksat*exp((n-1)*V*F/(R*T))))*(exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);////luo
	if(type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B) INaCa=INaCa*0.6;

	INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;////luo
	if(type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B) INaK=INaK*0.46;

	IpCa=GpCa*Cai/(KpCa+Cai);////luo
	IpK=GpK*rec_ipK*(V-Ek);////luo
	IbNa=GbNa*(V-Ena);////luo
	IbCa=GbCa*(V-Eca);////luo
	if(type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B) IbCa=IbCa*1.3;

	
	out_itot=IKr+IKs+IK1+Ito+INa+IbNa+ICaL+IbCa+INaK+INaCa+IpCa+IpK+Istim+INaL+IKATP+INa_l;	////luo//liang_change2
	
	kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR)));////luo
	k1=k1_/kCaSR;////luo
	k2=k2_*kCaSR;////luo
	dRR=k4*(1-sRR)-k2*CaSS*sRR;////luo
	sRR+=dt*dRR;////luo
	sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);////luo

	Irel=Vrel*sOO*(CaSR-CaSS);////luo
	Ileak=Vleak*(CaSR-Cai);////luo
	Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));////luo
	Ixfer=Vxfer*(CaSS-Cai);////luo
	if(type== ENDO_S_1B||type== EPI_S_1B||type==MCELL_S_1B)
	{
	Irel=Irel*0.65;
	Iup=Iup*0.71;	
	}

	CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);////luo
	dCaSR=dt*(Iup-Irel-Ileak);////luo
	bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;////luo
	cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);////luo
	CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2; ////luo

	CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);////luo
	dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));////luo
	bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;////luo
	ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);////luo
	CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;////luo

	CaBuf=Bufc*Cai/(Cai+Kbufc);////luo
	dCai=dt*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);////luo
	bc=Bufc-CaBuf-dCai-Cai+Kbufc;////luo
	cc=Kbufc*(CaBuf+dCai+Cai);////luo
	Cai=(sqrt(bc*bc+4*cc)-bc)/2;////luo
    
	dNai=-(INa+INaL+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;////luo
	Nai+=dt*dNai;////luo
    
	dKi=-(Istim+IK1+Ito+IKr+IKs+IKATP-2*INaK+IpK)*inverseVcF*CAPACITANCE;////luo
	Ki+=dt*dKi;////luo
 
	AM=1./(1.+exp((-60.-V)/5.));////luo
	BM=0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));////luo
	TAU_M=AM*BM;////luo
	M_INF=1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));////luo
	if(type== ENDO_S_1A||type== EPI_S_1A||type==MCELL_S_1A) 
	M_INF=1./((1.+exp((-55.5-V)/9.03))*(1.+exp((-55.5-V)/9.03)));

	if(V>=-40.)
		{
			AH_1=0.;////luo
			BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));////luo
			TAU_H=1.0/(AH_1+BH_1);////luo
		}
	else
		{
			AH_2=(0.057*exp(-(V+80.)/6.8));////luo
			BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));////luo
			TAU_H=1.0/(AH_2+BH_2);////luo
		}
	H_INF=1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));////luo
	if(V>=-40.)
		{
			AJ_1=0.;////luo
			BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));////luo
			TAU_J= 1.0/(AJ_1+BJ_1);////luo
		}
	else
		{
			AJ_2=(((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*exp(-0.04391*V))*(V+37.78)/(1.+exp(0.311*(V+79.23))));////luo
			BJ_2=(0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));////luo
			TAU_J= 1.0/(AJ_2+BJ_2);////luo
		}
	J_INF=H_INF;////luo


	alpham_nal = 0.32*(V+47.13)/(1.0 - exp(-0.1*(V+47.13)));//////////////////////////////////
	betam_nal = 0.08*exp(V/11.0);/////////////////////////////////////////////////////////////
	minf_nal = alpham_nal/(alpham_nal+betam_nal);/////////////////////////////////////////////
	taum_nal = 1.0/(alpham_nal+betam_nal);////////////////////////////////////////////////////
	hinf_nal = 1.0/(1.0 + exp((V+91.0)/6.1));/////////////////////////////////////////////////
	tauh_nal= 600.0;//////////////////////////////////////////////////////////////////////////

	Xr1_INF=1./(1.+exp((-26.-V)/7.));////luo
	axr1=450./(1.+exp((-45.-V)/10.));////luo
	bxr1=6./(1.+exp((V-(-30.))/11.5));////luo
	TAU_Xr1=axr1*bxr1;////luo
	Xr2_INF=1./(1.+exp((V-(-88.))/24.));////luo
	axr2=3./(1.+exp((-60.-V)/20.));////luo
	bxr2=1.12/(1.+exp((V-60.)/20.));////luo
	TAU_Xr2=axr2*bxr2;////luo

	Xs_INF=1./(1.+exp((-5.-V)/14.));////luo
	//Axs = 1400. * rsqrt(1. + exp((5.-V)/6) );
	Axs=(1400./(sqrt(1.+exp((5.-V)/6))));////luo
	Bxs=(1./(1.+exp((V-35.)/15.)));////luo
	TAU_Xs=Axs*Bxs+80;////luo
	
	if(type== EPI||type== EPI_S||type== EPI_S_1B||type== MCELL||type== MCELL_S||type== MCELL_S_1B)
		{
			R_INF=1./(1.+exp((20-V)/6.));////luo
			S_INF=1./(1.+exp((V+20)/5.));////luo
			TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;////luo
			TAU_S=85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;////luo
		}
	if(type== ENDO||type== ENDO_S||type== ENDO_S_1B)
		{
			R_INF=1./(1.+exp((20-V)/6.));////luo
			S_INF=1./(1.+exp((V+28)/5.));////luo
			TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;////luo
			TAU_S=1000.*exp(-(V+67)*(V+67)/1000.)+8.;////luo
		}
	if(type== EPI_S_1A||type== MCELL_S_1A)
		{
			R_INF=1./(1.+exp((27.2-V)/6.));////luo
			S_INF=1./(1.+exp((V+6.3)/5.));////luo
			TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;////luo
			TAU_S=85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;////luo
		}
	if(type== ENDO_S_1A)
		{
			R_INF=1./(1.+exp((27.2-V)/6.));////luo
			S_INF=1./(1.+exp((V+14.3)/5.));////luo
			TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;////luo
			TAU_S=1000.*exp(-(V+67)*(V+67)/1000.)+8.;////luo
		}

	D_INF=1./(1.+exp((-8-V)/7.5));////luo
	Ad=1.4/(1.+exp((-35-V)/13))+0.25;////luo
	Bd=1.4/(1.+exp((V+5)/5));////luo
	Cd=1./(1.+exp((50-V)/20));////luo
	TAU_D=Ad*Bd+Cd;////luo
	F_INF=1./(1.+exp((V+20)/7));////luo
	Af=1102.5*exp(-(V+27)*(V+27)/225);////luo
	Bf=200./(1+exp((13-V)/10.));////luo
	Cf=(180./(1+exp((V+30)/10)))+20;////luo
	TAU_F=(Af+Bf+Cf);////luo                    
	F2_INF=0.67/(1.+exp((V+35)/7))+0.33;////luo
	Af2=600*exp(-(V+25)*(V+25)/170);////luo
	Bf2=31/(1.+exp((25-V)/10));////luo
	Cf2=16/(1.+exp((V+30)/10));////luo
	TAU_F2=Af2+Bf2+Cf2;////luo
	FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;////luo
	TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;////luo



	sml = minf_nal - (minf_nal-sml)*exp(-dt/taum_nal);////////////////////////////////////////////////////
	shl = hinf_nal - (hinf_nal-shl)*exp(-dt/tauh_nal);////////////////////////////////////////////////////
		
	sm=M_INF-(M_INF-sm)*exp(-dt/TAU_M);////luo
	sh=H_INF-(H_INF-sh)*exp(-dt/TAU_H);////luo
	sj=J_INF-(J_INF-sj)*exp(-dt/TAU_J);////luo
	sxr1=Xr1_INF-(Xr1_INF-sxr1)*exp(-dt/TAU_Xr1);////luo
	sxr2=Xr2_INF-(Xr2_INF-sxr2)*exp(-dt/TAU_Xr2);////luo
	sxs=Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);////luo
	ss=S_INF-(S_INF-ss)*exp(-dt/TAU_S);////luo
	sr=R_INF-(R_INF-sr)*exp(-dt/TAU_R);////luo
	sd =D_INF-(D_INF-sd)*exp(-dt/TAU_D);////luo
	sf=F_INF-(F_INF-sf)*exp(-dt/TAU_F);////luo
	sf2=F2_INF-(F2_INF-sf2)*exp(-dt/TAU_F2);////luo
	sfcass=FCaSS_INF-(FCaSS_INF-sfcass)*exp(-dt/TAU_FCaSS);////luo
}

void call_kernel06(double* Cai_buffer, double* CaSR_buffer, double* CaSS_buffer, double* Nai_buffer,
				   double* Ki_buffer,  double* sm_buffer, double* sh_buffer, double* sj_buffer,
				   double* sxr1_buffer, double* sxr2_buffer,double* sxs_buffer,double* sr_buffer,
				   double* ss_buffer, double* sd_buffer, double* sf_buffer, double* sf2_buffer,
				   double* sfcass_buffer, double* sOO_buffer, double* sRR_buffer, double* sml_buffer,
				   double* shl_buffer, double* volt_buffer, double* istim_buffer, double* itot_buffer,
				   int* type_buffer, int cellNumber, int threadCount)
{
	int blockDim = threadCount;
	int gridDim = cellNumber/blockDim;
	gridDim = gridDim * blockDim < cellNumber ? gridDim + 1 : gridDim;
	/*
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	*/
	kernel06<<< gridDim, blockDim >>>
		( Cai_buffer,  CaSR_buffer,  CaSS_buffer,  Nai_buffer,
		  Ki_buffer,  sm_buffer,  sh_buffer,  sj_buffer,
		  sxr1_buffer,  sxr2_buffer,  sxs_buffer,  sr_buffer,
		  ss_buffer,  sd_buffer,  sf_buffer,  sf2_buffer,
		  sfcass_buffer,  sOO_buffer,  sRR_buffer,  sml_buffer,
		  shl_buffer,  volt_buffer,  istim_buffer,  itot_buffer,
		  type_buffer, cellNumber);
	//cudaDeviceSynchronize();
	/*
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &time, start, stop );
	printf("%fms\n", time);
	cudaEventDestroy( start );
	cudaEventDestroy( stop );
	*/
}

