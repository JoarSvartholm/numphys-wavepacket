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
	       gsl_vector_complex *psi)
{
  /** Initialize the  wave function **/

  /* Parameters have to be hard coded, read from a file, or
     read from the command line */
		 double sigma0 = 0.1;
		 double x0 = -1.;
		 double k0 = atof(argv[argc-1]);

		 for(size_t i = 0; i < params.nx_local; i++){
			 double complex psiTemp = pow(1./(M_PI*sigma0*sigma0),0.25)*(cos(k0*params.x[i])+ I*sin(k0*params.x[i]))*exp(-(params.x[i]-x0)*(params.x[i]-x0)/(2.*sigma0*sigma0));
			 gsl_complex temp = gsl_complex_rect(creal(psiTemp),cimag(psiTemp));
			 gsl_vector_complex_set(psi,i,temp);
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

void hamiltonian_operator(const parameters params, gsl_vector_complex *psi, gsl_vector_complex *Hpsi, gsl_complex alpha){

	gsl_complex H0,H1,psi0,psi1,psi2,a1,a2,a3;
	double hbar = params.hbar;
	extern double *V;
	GSL_SET_COMPLEX(&H1,-hbar*hbar/(2.*params.dx*params.dx),0);
	GSL_SET_COMPLEX(&H0,hbar*hbar/(params.dx*params.dx)+V[0],0);

	gsl_vector_complex_set(Hpsi,0,gsl_complex_add(gsl_complex_mul(gsl_vector_complex_get(psi,0),H0),gsl_complex_mul(gsl_vector_complex_get(psi,1),H1)));
	for(int i=1;i<params.nx-1;i++){
		GSL_SET_COMPLEX(&H0,hbar*hbar/(params.dx*params.dx)+V[i],0);

		psi0 = gsl_vector_complex_get(psi,i-1);
		psi1 = gsl_vector_complex_get(psi,i);
		psi2 = gsl_vector_complex_get(psi,i+1);

		a1 =gsl_complex_mul(psi0,H1);
		a2 =gsl_complex_mul(psi1,H0);
		a3 =gsl_complex_mul(psi2,H1);

		gsl_vector_complex_set(Hpsi,i,gsl_complex_add(a1,gsl_complex_add(a2,a3)));
	}
		GSL_SET_COMPLEX(&H0,hbar*hbar/(params.dx*params.dx)+V[params.nx-1],0);

		psi0 = gsl_vector_complex_get(psi,params.nx-1);
		psi1 = gsl_vector_complex_get(psi,params.nx-2);
		gsl_vector_complex_set(Hpsi,params.nx-1,gsl_complex_add(gsl_complex_mul(psi0,H0),gsl_complex_mul(psi1,H1)));

		/*for(int i=0;i<params.nx;i++){
			psi0 = gsl_vector_complex_get(Hpsi,i);
			gsl_vector_complex_set(Hpsi,i,gsl_complex_mul(alpha,psi0));
		}*/

}
