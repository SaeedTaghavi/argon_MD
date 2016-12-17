#include <iostream>
#include <fstream>
#include"Ar.h"
#include<cmath>
#include<algorithm>
using namespace std;

int main()
{
      double X[1000];
      double Y[1000];
      double Z[1000];
      double V_x[1000];
      double V_y[1000];
      double V_z[1000];
      double XT[1000];
      double YT[1000];
      double ZT[1000];
      double V_xT[1000];
      double V_yT[1000];
      double V_zT[1000];
      double MSD[1000];
      double VCF[1000];
      double CNORM[1000];
	int new_config,equal_loop,step_eq,step_main,step_save;
      int x,y;
        int i;

	double del=0.1e-10/Ar::H;
	int maxbin=round(1.0/del);
	int N=int(Ar::N);
	double G[maxbin];
	int nbin;
	cout<<del<<endl;
	cout<<maxbin<<endl;

	ifstream setup("config_file");
    	setup >>new_config>>equal_loop>>step_eq>>step_main>>step_save;
	
	ifstream in("./data/radial", ios::in | ios::binary);
	in.seekg(0,std::ios::beg);
	double n_step;
	double x_ij, y_ij, z_ij;
	double R_ij;	
	int kmax=step_main/step_save;
	double ktmax;
	double ncor;
	int kmin=0;
	int kct;

	    for(int k0=0; k0<kmax; k0++) 
	    {
		
		in.seekg(k0*48004);
		in.read((char *) &x,sizeof x);
		in.seekg(k0*48004+4);
		in.read((char *) &X,sizeof X);
		in.seekg(k0*48004+8004);
		in.read((char *) &Y,sizeof Y);
		in.seekg(k0*48004+16004);
		in.read((char *) &Z,sizeof Z);
		in.seekg(k0*48004+24004);
		in.read((char *) &V_x,sizeof V_x);
		in.seekg(k0*48004+32004);
		in.read((char *) &V_y,sizeof V_y);
		in.seekg(k0*48004+40004);
		in.read((char *) &V_z,sizeof V_z);
		
		ncor=(kmax-kmin)/2;
		ktmax=std::min(kmax,int(k0+ncor));
		
		for(int kt=k0;kt<ktmax;kt++)
		{
			in.seekg(kt*48004);
			in.read((char *) &x,sizeof x);
			in.seekg(kt*48004+4);
			in.read((char *) &XT,sizeof XT);
			in.seekg(kt*48004+8004);
			in.read((char *) &YT,sizeof YT);
			in.seekg(kt*48004+16004);
			in.read((char *) &ZT,sizeof ZT);
			in.seekg(kt*48004+24004);
			in.read((char *) &V_xT,sizeof V_xT);
			in.seekg(kt*48004+32004);
			in.read((char *) &V_yT,sizeof V_yT);
			in.seekg(kt*48004+40004);
			in.read((char *) &V_zT,sizeof V_zT);
			
			kct=kt-k0;
			
		
		
		    for(int k=0;k<N;k++)
		    {
			x_ij=X[k]-XT[k];	
			y_ij=Y[k]-YT[k];	
			z_ij=Z[k]-ZT[k];
			
			MSD[kct]+=x_ij*x_ij+y_ij*y_ij+z_ij*z_ij;
			VCF[kct]+=V_xT[k]*V_x[k]+V_yT[k]*V_y[k]+V_zT[k]*V_z[k];
			//CNORM[kct]+=1.0;
		    }
		
		    CNORM[kct]+=1.0;

		}//kt

		

	    }//k0
	    
	    in.close();
	    double rho=N/8.0;
	    double rlower, rupper,vol,dind;
	    std::ofstream myfile;
	    std::ofstream myfile2;
	    myfile.open("./data/cor",std::ios::in | std::ios::app);
	    myfile2.open("./data/vel",std::ios::in | std::ios::app);
	    //myfile<<0<<"\t"<<0<<endl;
	    for(int a=0;a<kmax/2;a++)
	    {
		MSD[a]/=CNORM[a];
		//VCF[a]/=CNORM[a];
		//cout<<CNORM[a]<<endl; 
		myfile<<a/double(step_save)<<"\t"<<(MSD[a]*Ar::H*Ar::H)/1e-10<<endl;
		myfile2<<a/double(step_save)<<"\t"<<VCF[a]<<endl;
	    }
	    myfile.close();
	    myfile2.close();


		    return 0;
            }
