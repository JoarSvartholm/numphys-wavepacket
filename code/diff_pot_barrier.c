/********************************************************************

 Sample file for user defined functions

 To be used with the program "wavepacket"

 Claude Dion, February 2018

********************************************************************/

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "wavediff.h"

#ifdef MPI
#include "mpi.h"
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#define SQUARE(z) ((z) * (z))

#define NORM2(z) (conj(z)*z)

/* Shared variables (if needed)*/

// Potential
double *V;
// For user observables
FILE *fp_user;
double a;


void
initialize_potential (const parameters params, const int argc,
		      char ** const argv)
{
  /* Initialize the potential */

  extern double *V;
  V = (double *) calloc (params.nx_local, sizeof (double));

  /* If the potential is time-independent, it can be precalculated here,
     stored in V, and used by the potential function */
		 extern double a;
		 a = 0.01;
		 double V0 = 2e4;

		 for(size_t i=0;i<params.nx_local;i++){
			 if(fabs(params.x[i])<a) V[i] = V0;
		 }

  return;
}


void
potential (const parameters params, const double t, double * const pot)
{
  /** Potential function **/


  // Use the precalculated potential
  extern double *V;

  for (size_t i = 0; i < params.nx_local; ++i)
    {
      pot[i] = V[i];
    }

  return;
}


void
initialize_wf (const parameters params, const int argc, char ** const argv,
	       double complex *psi)
{
  /** Initialize the  wave function **/

  /* Parameters have to be hard coded, read from a file, or
     read from the command line */
		 double sigma0 = 0.1;
		 double x0 = -1.;
		 double k0 = atof(argv[argc-1]);

		 for(size_t i = 0; i < params.nx_local; i++){

			 double psiReal = pow(1./(M_PI*sigma0*sigma0),0.25)*exp(-(params.x[i]-x0)*(params.x[i]-x0)/(2.*sigma0*sigma0))*cos(k0*params.x[i]);
			 double psiImag = pow(1./(M_PI*sigma0*sigma0),0.25)*exp(-(params.x[i]-x0)*(params.x[i]-x0)/(2.*sigma0*sigma0))*sin(k0*params.x[i]);
			 psi[i] = psiReal + I*psiImag;
		 }



  return;
}


void
initialize_user_observe (const parameters params, const int argc,
			 char ** const argv)
{
  /* Inistialization needed for user_observe */


	char fname[60];
	sprintf(fname,"../data/pot_barrier_RT_k0_%s.dat",argv[argc-1]);
	fp_user = fopen(fname,"w+");

  fprintf(fp_user, "t    R    T\n");

  return;
}


void
user_observe (const parameters params, const double t,
	      const double complex * const psi)
{
  /* User-defined observables, such as transmission and reflection */
	extern double a;
	double R=0;
	double T=0;

	for(int i=0;params.x[i]<=-a;i++){
		if(i==0 || params.x[i]==-a) R+=0.5*NORM2((psi[i]))*params.dx;
		else R+= NORM2(psi[i])*params.dx;
	}
	for(int i=49152/2;i<params.nx;i++){
		if(i==(int)params.nx || params.x[i]==a) T+=0.5*NORM2(psi[i])*params.dx;
		else if(params.x[i]>a) T+= NORM2(psi[i])*params.dx;
	}

	fprintf(fp_user, "%e %e %e\n", t, R, T);



  return;
}

void initialize_hamiltonian(const parameters params, gsl_matrix_complex *H){

	double hbar = params.hbar;
	extern double *V;

	gsl_complex H1;
	gsl_complex H0;
	GSL_SET_COMPLEX(&H1,-hbar*hbar/(2.*params.dx*params.dx),0);


    gsl_matrix_complex_set(H,0,1,H1);
  for(int i=1;i<params.nx-1;i++){
		GSL_SET_COMPLEX(&H0,hbar*hbar/(params.dx*params.dx)+V[i],0);
    gsl_matrix_complex_set(H,i,i,H0);
    gsl_matrix_complex_set(H,i,i+1,H1);
    gsl_matrix_complex_set(H,i,i-1,H1);
  }
    gsl_matrix_complex_set(H,params.nx-1,params.nx-2,H1);

}

void hamiltonian_operator(const parameters params, double complex *psi, double complex *Hpsi){

	double htmx2,th;
	double Hii,Hi1;
	double complex Hp1,Hp2,Hp3;
	extern double *V;
	int N = params.nx_local;

	htmx2 = params.hbar*params.dt/(params.mass*params.dx*params.dx);
	th = params.dt/params.hbar;

	Hi1 = -0.5*htmx2;
	Hii = htmx2 + th*V[0];

	Hp2 =Hii*psi[0];
	Hp3 =Hi1*psi[1];

	//Set first element
	Hpsi[0] = Hp2 + Hp3;

	for(size_t i=1;i<N-1;i++){
		Hii = htmx2 + th*V[i];

		//Compute sum elements
		Hp1 = Hi1*psi[i-1];
		Hp2 = Hii*psi[i];
		Hp3 = Hi1*psi[i+1];

		//Set element i
		Hpsi[i] = Hp1 + Hp2 + Hp3;
	}

	Hii = htmx2 + th*V[N-1];

	Hp1 = Hi1*psi[N-2];
	Hp2 = Hii*psi[N-1];

	//Set last element
	Hpsi[N-1] = Hp1 + Hp2;


		/*for(int i=0;i<params.nx;i++){
			psi0 = gsl_vector_complex_get(Hpsi,i);
			gsl_vector_complex_set(Hpsi,i,gsl_complex_mul(alpha,psi0));
		}*/

}

void take_step(double complex type,double complex *Hpsi, double complex *psi,parameters params){

	for(size_t i=0;i<params.nx_local;i++){
		psi[i] += type*Hpsi[i];
	}

}

void swap_vectors(double complex *psi1,double complex *psi2,parameters params){
	double complex temp;

	for(size_t i=0;i<params.nx_local;i++){
		temp = psi1[i];
		psi1[i] = psi2[i];
		psi2[i] = temp;
	}


}
