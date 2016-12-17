#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include"Ar.h"
#include<fstream>

// Symulacja ciekłego argonu prędkościowym algorytmem Verleta oraz potencjałem Lenarda - Jonesa//
// wersja 2.4 (dodanie algorytmu obliczania błędów energii)
// -----------------------------------------------------------//



//------------------------- Obliczanie sił i energii potencjalnej -----------------------------------

void Force(double *X,double *Y, double *Z, double *F_x, double *F_y, double *F_z, double N,double *E_point)
{

	for(int a=0;a<1000;a++)
	{
		F_x[a]=0.;
		F_y[a]=0.;
		F_z[a]=0.;
	}

	*E_point=0.0;

	 
	double x_i, y_i, z_i;
	double x_ij, y_ij, z_ij;
	double r_ij2, rs2, rs6, rs12;
	double D_p;
	

	for(int i=0;i<N-1;i++)
	{
		x_i=X[i];
		y_i=Y[i];
		z_i=Z[i];

		for(int j=i+1;j<N;j++)
		{

			x_ij=x_i - X[j];
			y_ij=y_i - Y[j];
			z_ij=z_i - Z[j];


			x_ij=x_ij - 2.0*int(x_ij);
			y_ij=y_ij - 2.0*int(y_ij);
			z_ij=z_ij - 2.0*int(z_ij);

			r_ij2=x_ij*x_ij+y_ij*y_ij+z_ij*z_ij;

			if( r_ij2 <= 1.0)
			{
				rs2=Ar::sigma_n2/r_ij2;
				rs6=rs2*rs2*rs2;
				rs12=rs6*rs6;

				*E_point=*E_point+rs12-rs6;
				D_p=(2.0*rs12-rs6)/r_ij2;

				F_x[i]=F_x[i]+D_p*x_ij;
				F_y[i]=F_y[i]+D_p*y_ij;
				F_z[i]=F_z[i]+D_p*z_ij;


				F_x[j]=F_x[j]-D_p*x_ij;
				F_y[j]=F_y[j]-D_p*y_ij;
				F_z[j]=F_z[j]-D_p*z_ij;
			}

		}

	}
	
	*E_point*=4;
	for (int l=0;l<N;l++)
	{

	    F_x[l]*=24;
	    F_y[l]*=24;
	    F_z[l]*=24;
	}



}

//-------------------------------------------------------------------------------------------------------------



//--------------------------------- Skalowanie prędkości układu ----------------------------------------------

void Scale_V(double *V_x, double *V_y, double *V_z, double N,double Scale, double *Ek_point)
{
	*Ek_point=0.0;
	for(int p=0;p<N;p++)
	{
		V_x[p]*=Scale;
		V_y[p]*=Scale;
		V_z[p]*=Scale;
		*Ek_point+=V_x[p]*V_x[p]+V_y[p]*V_y[p]+V_z[p]*V_z[p];
	}

}
// -------------------------------------------------------------------------------------------------------------



// --------------------------------Redukowanie pędu układu -----------------------------------------------------

void Reduce_M(double *V_x, double *V_y, double *V_z, double N, double *Ek_point)
{

	double Vcx=0.;
	double Vcy=0.;
	double Vcz=0.;

	*Ek_point=0.0;

	for(int k=0;k<N;k++)
	{
	    Vcx+=V_x[k];
	    Vcy+=V_y[k];
	    Vcz+=V_z[k];
 
	}    

	Vcx/=N;
	Vcy/=N;
	Vcz/=N;

	for(int l=0;l<N;l++)
	{
		V_x[l]-=Vcx;
		V_y[l]-=Vcy;
		V_z[l]-=Vcz;
		*Ek_point+=V_x[l]*V_x[l]+V_y[l]*V_y[l]+V_z[l]*V_z[l];
	}

}

//-----------------------------------------------------------------------------------------------------------



int main()
{
	srand(time(NULL));
	 using namespace std;
	int  M=Ar::M;
	int new_config;
	int equal_loop;
       	double X[1000];
       	double Y[1000];
	double Z[1000];
	double V_x[1000];
	double V_y[1000];
	double V_z[1000];
	double F_x[1000];
	double F_y[1000];
	double F_z[1000];

	double N=Ar::N;
	
       	double Ek_n=0;   // energia kinetyczna 
       	double *Ek_point=&Ek_n;

	 double E_p=0.0; //energia potencjalna 
	 double *Ep_point=&E_p;	
	 double T_r=Ar::epsilon_k*Ek_n/(Ar::dt_n2*3.0*(N-1.0));
	 double Skal_T=sqrt(Ar::T/T_r);	 
	 int ijk=0;

	int step_eq; // Ilość kroków czasowych (dochodzenie do równowagi)	
	int step_main; // Ilość kroków czasowych
	int step_save;

// -------------- Załadowanie parametrów ----------------------------------------
	
    ifstream setup("config_file");
    setup >>new_config>>equal_loop>>step_eq>>step_main>>step_save;
//-------------------------------------------------------------------------

// ------------------Inicjalizacja układu (losowanie położeń i prędkości) ------------------------
if(new_config==1)
{
	for(int i=1;i<=M;i++)
	{
		for(int j=1;j<=M;j++)
		{
			for(int k=1;k<=M;k++)
			{

				X[ijk]=-1+(i-1)*Ar::A;
				Y[ijk]=-1+(k-1)*Ar::A;
				Z[ijk]=-1+(j-1)*Ar::A;

				X[ijk]+=Ar::disp*(double((rand()%1000000))/1000000.-0.5)*2.;
				Y[ijk]+=Ar::disp*(double((rand()%1000000))/1000000.-0.5)*2.;
				Z[ijk]+=Ar::disp*(double((rand()%1000000))/1000000.-0.5)*2.;


				X[ijk]=X[ijk]-2.*int(X[ijk]);
				Y[ijk]=Y[ijk]-2.*int(Y[ijk]);
				Z[ijk]=Z[ijk]-2.*int(Z[ijk]);


				V_x[ijk]=Ar::V_maxn*(double((rand()%1000000))/1000000.-0.5)*2.;
				V_y[ijk]=Ar::V_maxn*(double((rand()%1000000))/1000000.-0.5)*2.;
				V_z[ijk]=Ar::V_maxn*(double((rand()%1000000))/1000000.-0.5)*2.;

				
				ijk+=1;
			}
		}
	}
}
// -------------------------------------------------------------------------------------------


// ----------------------- Wczytywanie określonej konfiguracji ------------------------------
if(new_config==0)
{
   
    ifstream state_load("./data/config/end_state", ios::in | ios::binary);
    for(int k=0;k<5;k++)
    {
	switch(k)
	{
	    case 0:
		state_load.seekg(k*sizeof(X));
		state_load.read((char *) &X, sizeof X);
		break;
	    case 1:
		state_load.seekg(k*sizeof(Y));
		state_load.read((char *) &Y, sizeof Y);
		break;
	    case 2:
		state_load.seekg(k*sizeof(Z));
		state_load.read((char *) &Z, sizeof Z);
		break;
	    case 3:
		state_load.seekg(k*sizeof(V_x));
		state_load.read((char *) &V_x, sizeof V_x);
		break;
	    case 4:
		state_load.seekg(k*sizeof(V_y));
		state_load.read((char *) &V_y, sizeof V_y);
		break;
	    case 5:
		state_load.seekg(k*sizeof(V_z));
		state_load.read((char *) &V_z, sizeof V_z);
		break;
	}
	


    }
    /*ifstream X_load("./data/config/X", ios::in | ios::binary);
    X_load.read((char *) &X, sizeof X);
    X_load.close();

    ifstream Y_load("./data/config/Y", ios::in | ios::binary);
    Y_load.read((char *) &Y, sizeof Y);
    Y_load.close();

    ifstream Z_load("./data/config/Z", ios::in | ios::binary);
    Z_load.read((char *) &Z, sizeof Z);
    Z_load.close();

    ifstream Vx_load("./data/config/Vx", ios::in | ios::binary);
    Vx_load.read((char *) &V_x, sizeof V_x);
    Vx_load.close();

    ifstream Vy_load("./data/config/Vy", ios::in | ios::binary);
    Vy_load.read((char *) &V_y, sizeof V_y);
    Vy_load.close();

    ifstream Vz_load("./data/config/Vz", ios::in | ios::binary);
    Vz_load.read((char *) &V_z, sizeof V_z);
    Vz_load.close();*/
}
//-------------------------------------------------------------------------------------------
	Reduce_M(V_x,V_y,V_z,N,Ek_point);

	T_r=Ar::epsilon_k*Ek_n/(Ar::dt_n2*3.0*(N-1.0));
	Skal_T=sqrt(Ar::T/T_r);
	
	cout.width(10);
	cout<<internal<<"Parametry symulacji"<<"\n";
	cout<<"--------------------------------------------------------------------"<<"\n";
	cout<<"Liczba atomów: "<< N<<endl;
	cout<<"Temperatura żądana: "<< Ar::T<<endl;
	cout<<"Krok czasowy: "<< Ar::dt<<endl;
	cout<<"V_max: "<< Ar::V_maxn<<endl;
	
	Scale_V(V_x,V_y,V_z,N,Skal_T,Ek_point);


	T_r=Ar::epsilon_k*Ek_n/(Ar::dt_n2*3.0*(N-1.0));
	
    
	Force(X,Y,Z,F_x,F_y,F_z,N,Ep_point); // Obliczanie sił w czasie t
    	
	cout<<"Energia całkowita przed dochodzeniem do równowagi: "<<Ek_n+E_p<<endl;
	cout<<"Krok czasowy zredukowany: "<< Ar::dt_n<<endl;
	cout<<"---------------------------------------------------------------------"<<endl;



// --------------- Dochodzenie do równowagi -----------------------------------------------
if(equal_loop==1)
{
	for(int a=0;a<step_eq;a++)
	{

	    for(int p=0;p<N;p++)			    //Aktualizacja połżeń
	    {
		X[p]=X[p]+V_x[p]+0.5*F_x[p]*Ar::dt_n2;
		Y[p]=Y[p]+V_y[p]+0.5*F_y[p]*Ar::dt_n2;
		Z[p]=Z[p]+V_z[p]+0.5*F_z[p]*Ar::dt_n2;
		
		X[p]=X[p]-2.*int(X[p]);
		Y[p]=Y[p]-2.*int(Y[p]);
		Z[p]=Z[p]-2.*int(Z[p]);

		V_x[p]=V_x[p]+0.5*F_x[p]*Ar::dt_n2;		    // Aktualizacja prędkości w czasie t + 1/2 dt
		V_y[p]=V_y[p]+0.5*F_y[p]*Ar::dt_n2;
		V_z[p]=V_z[p]+0.5*F_z[p]*Ar::dt_n2;

	    }

	   
	    Force(X,Y,Z,F_x,F_y,F_z,N,Ep_point);	    // Obliczanie sił w czasie t + dt

	    for(int r=0;r<N;r++)
	    {

		V_x[r]=V_x[r]+0.5*F_x[r]*Ar::dt_n2;		    // Aktualizacja prędkości w czasie t +  dt
		V_y[r]=V_y[r]+0.5*F_y[r]*Ar::dt_n2;
		V_z[r]=V_z[r]+0.5*F_z[r]*Ar::dt_n2;

	    }

	    Reduce_M(V_x,V_y,V_z,N,Ek_point);		    //Redukowanie pędu
	    

	    T_r=Ar::epsilon_k*Ek_n/(Ar::dt_n2*3.0*(Ar::N-1.0));
	    Skal_T=sqrt(Ar::T/T_r);


	    if(abs(Ar::T-T_r)>3.0)
	    {
		
		Scale_V(V_x,V_y,V_z,N,Skal_T,Ek_point);			
	    }

	}
}
// ------------------Główna pętla -------------------------------------------

	double E1=0.;
	double E2=0.;
	double T1=0.;
 	double T2=0.;
	double E_total;
	double errorE,errorT;
	double E10,E20,T10,T20;	

	cout.precision(10);
	ofstream radial("./data/radial",ios::out | ios::binary);
	//cout.width(12);
	cout<<"Krok | "<< " Energia |"<<" sigma E |"<< " Temperatura | "<<" sigma T |"<<endl;	
	for(int a=0;a<step_main;a++)
	{
	    
	    for(int p=0;p<N;p++)			    //Aktualizacja połżeń
	    {
		X[p]=X[p]+V_x[p]+0.5*F_x[p]*Ar::dt_n2;
		Y[p]=Y[p]+V_y[p]+0.5*F_y[p]*Ar::dt_n2;
		Z[p]=Z[p]+V_z[p]+0.5*F_z[p]*Ar::dt_n2;
		
		X[p]=X[p]-2.*int(X[p]);
		Y[p]=Y[p]-2.*int(Y[p]);
		Z[p]=Z[p]-2.*int(Z[p]);

		V_x[p]=V_x[p]+0.5*F_x[p]*Ar::dt_n2;		    // Aktualizacja prędkości w czasie t + 1/2 dt
		V_y[p]=V_y[p]+0.5*F_y[p]*Ar::dt_n2;
		V_z[p]=V_z[p]+0.5*F_z[p]*Ar::dt_n2;

	    }

	   
	    Force(X,Y,Z,F_x,F_y,F_z,N,Ep_point);	    // Obliczanie sił w czasie t + dt

	    for(int r=0;r<N;r++)
	    {

		V_x[r]=V_x[r]+0.5*F_x[r]*Ar::dt_n2;		    // Aktualizacja prędkości w czasie t +  dt
		V_y[r]=V_y[r]+0.5*F_y[r]*Ar::dt_n2;
		V_z[r]=V_z[r]+0.5*F_z[r]*Ar::dt_n2;

	    }

	    Reduce_M(V_x,V_y,V_z,N,Ek_point);		    //Redukowanie pędu
	    
	    E_total=Ek_n+E_p;	
	    T_r=Ar::epsilon_k*Ek_n/(Ar::dt_n2*3.0*(N-1.0));
	    
	    T1=T1+T_r;
	    T2=T2+T_r*T_r;

	    E1=E1+E_total;
	    E2=E2+E_total*E_total;

		
	    E10=E1/double(a+1.);
	    E20=E2/double(a+1.);
	    T10=T1/double(a+1.);
	    T20=T2/double(a+1.);

	    errorT=sqrt(T20-T10*T10);
	    errorE=sqrt(E20-E10*E10);
	
	    if(a%step_save==0)
	    {
		cout.setf(ios::showpoint);
		cout.width(2);
		cout<<internal<<a<< "\t"<< E_total<<"\t"<<Ek_n<<"\t"<<errorE<<"\t"<<T_r<<"\t"<<errorT<<endl;					

		radial.write((char *) &a, sizeof a);
		radial.write((char *) &X, sizeof X);
		radial.write((char *) &Y, sizeof Y);
		radial.write((char *) &Z, sizeof Z);
		radial.write((char *) &V_x, sizeof V_x);
		radial.write((char *) &V_y, sizeof V_y);
		radial.write((char *) &V_z, sizeof V_z);


	    }

	}

    radial.close();
//-----------------------------------------------------------------------

    ofstream state_file("./data/end_state", ios::out | ios::binary);

    state_file.write((char *) &X, sizeof X);
    state_file.write((char *) &Y, sizeof Y);
    state_file.write((char *) &Z, sizeof Z);
    state_file.write((char *) &V_x, sizeof V_x);
    state_file.write((char *) &V_y, sizeof V_y);
    state_file.write((char *) &V_z, sizeof V_z);
    state_file.close();
    /*X_file.write((char *) &X, sizeof X);
    X_file.close();


    ofstream Y_file("./data/Y", ios::out | ios::binary);
    Y_file.write((char *) &Y, sizeof Y);
    Y_file.close();


    ofstream Z_file("./data/Z", ios::out | ios::binary);
    Z_file.write((char *) &Z, sizeof Z);
    Z_file.close();


    ofstream Vx_file("./data/Vx", ios::out | ios::binary);
    Vx_file.write((char *) &V_x, sizeof V_x);
    Vx_file.close();


    ofstream Vy_file("./data/Vy", ios::out | ios::binary);
    Vy_file.write((char *) &V_y, sizeof V_y);
    Vy_file.close();

    ofstream Vz_file("./data/Vz", ios::out | ios::binary);
    Vz_file.write((char *) &V_z, sizeof V_z);
    Vz_file.close();
*/
    return 0;

}
