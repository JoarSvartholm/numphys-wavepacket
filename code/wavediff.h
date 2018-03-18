#ifndef __WAVEPACKET__
#define __WAVEPACKET__

#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_spmatrix.h>

/* Structures */

typedef struct
{
  /* Parameters and grid */
  int size, rank;
  size_t nx, ny, nz, n, nx_local, nx0, n_local;
  double x_min, y_min, z_min, x_max, y_max, z_max, dx, dy, dz;
  double *x, *y, *z, *x2, *y2, *z2;
  double mass, dt, hbar;
  bool is1D;
} parameters;


/* Function declarations */

void abort ();


double expectation1D (const parameters, const int, const double * const,
		      const double complex * const);

void initialize_potential (const parameters, const int, char ** const);

void initialize_user_observe (const parameters, const int, char ** const);

void initialize_wf (const parameters, const int, char ** const,
		    double complex *);

double complex integrate3D (const parameters, const double complex * const,
			    const double complex * const);

double norm (const parameters, const double complex * const);

void potential (const parameters, const double, double * const);

void renormalize (const parameters, double complex *);

void user_observe (const parameters, const double,
		   const double complex * const);

void initialize_hamiltonian(const parameters,gsl_matrix_complex *);

void hamiltonian_operator( const parameters ,double complex *, double complex *);

void take_step(double complex ,double complex *, double complex *,parameters );

void swap_vectors(double complex *,double complex *,parameters );

void
print_wf2 (const parameters params, double complex *psi,
	  const char * const wf_text);

#endif /* __WAVEPACKET__ */
