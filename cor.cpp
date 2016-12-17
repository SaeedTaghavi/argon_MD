#include <iostream>
#include <fstream>
#include"Ar.h"
#include<cmath>
using namespace std;

int main()
{
      double fnum[4] = {9.5, -3.4, 1.0, 2.1};
      double a=1;
      double b=3;
      double tab[2][2]={{5,11},{32,3.14}};
      double tab3[2][2]={{3,24},{9,6.28}};
      double f2[4];
      double tab2[2][2];
      double tab4[2][2];
      double X[1000];
      double Y[1000];
      double Z[1000];
      double V_x[1000];
      double V_y[1000];
      double V_z[1000];
      int x,y;
        int i;

	double del=0.1e-10/Ar::H;
	int maxbin=round(1.0/del);
	int N=int(Ar::N);
	double G[maxbin];
	int nbin;
	cout<<del<<endl;
	cout<<maxbin<<endl;

	ifstream in("./data/radial", ios::in | ios::binary);
	in.seekg(0,std::ios::beg);
	double n_step=Ar::step_save;
	double x_ij, y_ij, z_ij;
	double R_ij;	
	int nmax=500;
	double X_all[nmax][1000];
	double Y_all[nmax][1000];
	double Z_all[nmax][1000];
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
		
		for(int k=0;k<N;k++)
		{
			X_all[i][k]=X[k];
			Y_all[i][k]=Y[k];
			Z_all[i][k]=Z[k];

		}
		//cout<<x<<"\t"<<X[1]<<"\t"<<"\t"<<Y[1]<<"\t"<<Z[1]<<"\t"<<V_x[1]<<"\t"<<V_y[1]<<"\t"<<V_z[1]<<endl;

		/*for(int l=0;l<N-1;l++)
		{
		    for(int k=l+1;k<N;k++)
		    {
			x_ij=X[k]-X[l];	
			y_ij=Y[k]-Y[l];	
			z_ij=Z[k]-Z[l];
			
			R_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
			nbin=round(R_ij/del)+1;
			G[nbin]+=1.0;	

		    }


		}
		*/

	    }
	    
	    in.close();
	    double rho=N/8.0;
	    double rlower, rupper,vol,dind;
		cout<<X_all[10][10]<<endl;

		cout<<Y_all[10][10]<<endl;
		cout<<Z_all[10][10]<<endl;
	    /*std::ofstream myfile;
	    //myfile.open("./data/plot",std::ios::in | std::ios::app);
	    //myfile<<0<<"\t"<<0<<endl;
	    for(int k=1;k<maxbin+1;k++)
	    {

		rlower=(k-1.)*del;
		rupper=k*del;

		//vol=4.0*3.1415926*(rupper*rupper*rupper-rlower*rlower*rlower)/3.0;
		dind=vol*rho;
		//cout<<double(N*nmax)*dind<<endl;
		//cout<<1.0/vol*rho<<endl;
		//cout<<G[k]<<endl;
		G[k]=G[k]/(N);
    
		//cout<<(pow(rupper,3)-pow(rlower,3))<<endl;
		G[k]=G[k]/(60.0);
		//cout<<G[k]<<endl;
		//myfile<<k<<"\t"<<G[k]<<endl;
	    }
	    myfile.close();

*/
		    return 0;
            }
