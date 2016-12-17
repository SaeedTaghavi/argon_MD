#include <iostream>
#include <fstream>
#include"Ar.h"
#include<cmath>
#include<math.h>
using namespace std;

int main()
{
      double X[1000];
      double Y[1000];
      double Z[1000];
      double V_x[1000];
      double V_y[1000];
      double V_z[1000];
      int x,y;
        int i;

	double del=0.1e-10/Ar::H;
	int maxbin=int(1.0/del);
	int N=int(Ar::N);
	double G[200];
	int nbin;
	cout<<del<<endl;
	cout<<maxbin<<endl;
	    ifstream in("./data/radial", ios::in | ios::binary);
	   in.seekg(0,std::ios::beg);
	double x_ij, y_ij, z_ij;
       double R_ij;
	double dt=Ar::dt;
	int step;		

	int a,b;
	int step_eq, step_main;

	ifstream setup("config_file");
	setup>>a>>b>>step_eq>>step_main;

	int nmax=step_main/Ar::step_save;


       for(int s=0;s<maxbin;s++)
       {
	   G[s]=0;
       }

	    for(i=0; i<nmax; i++) // show values read from file
	    {
		
		in.seekg(i*48004);
		in.read((char *) &x,sizeof x);
		in.seekg(i*48004+4);
		in.read((char *) &X,sizeof X);
		in.seekg(i*48004+8004);
		in.read((char *) &Y,sizeof Y);
		in.seekg(i*48004+16004);
		in.read((char *) &Z,sizeof Z);
		in.seekg(i*48004+24004);
		in.read((char *) &V_x,sizeof V_x);
		in.seekg(i*48004+32004);
		in.read((char *) &V_y,sizeof V_y);
		in.seekg(i*48004+40004);
		in.read((char *) &V_z,sizeof V_z);
		
		//cout<<x<<"\t"<<X[1]<<"\t"<<"\t"<<Y[1]<<"\t"<<Z[1]<<"\t"<<V_x[1]<<"\t"<<V_y[1]<<"\t"<<V_z[1]<<endl;

		for(int l=0;l<N;l++)
		{
		    for(int k=0;k<N;k++)
		    {
			if(k!=l){
			x_ij=X[k]-X[l];	
			y_ij=Y[k]-Y[l];	
			z_ij=Z[k]-Z[l];

			x_ij-=2.0*int(x_ij);	
			y_ij-=2.0*int(y_ij);	
			z_ij-=2.0*int(z_ij);


			R_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
			nbin=int(R_ij/del)+1;
			//cout<<nbin<<"\t"<<R_ij<<endl;
			if(nbin<=maxbin){G[nbin]+=1.0;}	
			}
		    }


		}


	    }
	 //   in.seekg(48004);

	   //in.read((char *) &step,sizeof step);
	cout<<Ar::H<<endl;	
	    in.close();
	    double rho=N/8.0;
	    double rlower, rupper,vol,dind;

	    std::ofstream myfile;
	    myfile.open("./data/plot",std::ios::in | std::ios::trunc);
	    myfile<<0<<"\t"<<0<<endl;
	    for(int k=1;k<maxbin+1;k++)
	    {

		rlower=(k-1.)*del;
		rupper=k*del;

		vol=4.0*3.1415926*(rupper*rupper*rupper-rlower*rlower*rlower)/3.0;
		dind=vol*rho;
		//cout<<double(N*nmax)*dind<<endl;
		//cout<<1.0/vol*rho<<endl;
		//cout<<dind<<endl;
		G[k]=G[k]/(double(N));
    
		//cout<<(pow(rupper,3)-pow(rlower,3))<<endl;
		G[k]=G[k]/(double(nmax)*dind);
		//cout<<G[k]<<endl;
		myfile<<k*0.1<<"\t"<<G[k]<<endl;
	    }
	    myfile.close();


		    return 0;
            }
