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

#include "wavepacket.h"

#ifdef MPI
#include "mpi.h"
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#define SQUARE(z) ((z) * (z))

/* Shared variables (if needed)*/

// Potential
double *V;
// For user observables
FILE *fp_user;

void
initialize_potential (const parameters params, const int argc,
		      char ** const argv)
{
  /* Initialize the potential */

  extern double *V;
  V = (double *) calloc (params.nx_local, sizeof (double));

  /* If the potential is time-independent, it can be precalculated here,
     stored in V, and used by the potential function */
		 double a = 0.02;
		 double V0 = -1e-5;

		 for(size_t i=0;i<params.nx_local;i++){
			 V[i] = 0;
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
		 double x0 = -2.;
		 double k0 = 240;

		 for(size_t i = 0; i < params.nx_local; i++){
			 psi[i] = pow(1./(M_PI*sigma0*sigma0),0.25)*(cos(k0*params.x[i])+ I*sin(k0*params.x[i]))*exp(-(params.x[i]-x0)*(params.x[i]-x0)/(2.*sigma0*sigma0));
		 }



  return;
}


void
initialize_user_observe (const parameters params, const int argc,
			 char ** const argv)
{
  /* Inistialization needed for user_observe */

  return;
}


void
user_observe (const parameters params, const double t,
	      const double complex * const psi)
{
  /* User-defined observables, such as transmission and reflection */

  return;
}
