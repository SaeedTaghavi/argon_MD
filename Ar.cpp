#include<cmath>
namespace Ar
{
	extern const double T(87.3);
	extern const double ro(1395);
	extern const double N_a=6.022045e23;
	extern const double m_g=40.;
	extern const double k_b=1.38066e-23;
	extern const double m_kg=m_g*1e-3/N_a;
	extern const double epsilon_k=119.8;
	extern const double epsilon_j=epsilon_k*k_b;
	extern const double sigma=3.41e-10;

	extern const int  M=10;
	extern const double N=pow(M,3);
	extern const double A=2./double(M);

	extern const double T_n=T/epsilon_k;
	extern const double L=cbrt((N*m_g/N_a)*1e-3/ro);
	extern const double H=L/2.0;
	extern const double sigma_n=sigma/H;
	extern const double sigma_n2=sigma_n*sigma_n;

	extern const double disp=0.02;
	extern const double dt=1e-14;
	extern const double dt_2=dt*dt;

	extern const int step_eq=2000;
	extern const int step_main=2000;
	extern const int step_save=20;

		
	extern const double dt_n=dt/sqrt((H*H*m_g*1e-3)/(N_a*epsilon_j));
	extern const double dt_n2=dt_n*dt_n;
	extern const double V_maxn=dt_n*sqrt(3.0*T_n);
}


