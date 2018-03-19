#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavediff.h"
#include "wavepacket_private.h"

int main(int argc, char *argv[]) {

  size_t alloc_local;
  unsigned int nt_inner, nt_outer;
  double t;
  parameters params;
  output_flags output;
  char results_file[80], wf_text[80], wf_bin[80];
  FILE *fp_results;

  read_parameters (&params, &output, &nt_inner, &nt_outer, argv[1], results_file, wf_text, wf_bin);

  alloc_local = params.n;
  params.nx_local = params.nx;
  params.nx0 = 0;


  // Total local number of grid points
  params.n_local = params.nx_local * params.ny * params.nz;

  // Make spatial grid
  make_grid (&params);

  //Stepping constant
  double complex HALF_BACK = 0.5*I;
  double complex HALF_FWD = -0.5*I;
  double complex FWD_STEP = -2.*I;

  //Create array for wavefunction
  double complex *psi;
  psi = (double complex *) malloc (params.n_local * sizeof (double complex));
  double complex *Hpsi;
  Hpsi = (double complex *) malloc (params.n_local * sizeof (double complex));
  double complex *psi1;
  psi1 = (double complex *) malloc (params.n_local * sizeof (double complex));

  initialize_user_observe (params, argc, argv);

  initialize_wf(params,argc,argv,psi);
  initialize_wf(params,argc,argv,psi1);


  initialize_potential (params, argc, argv);

  //Generate symmetric IC
  //back euler step
  hamiltonian_operator(params,psi1,Hpsi);
  take_step(HALF_BACK,Hpsi,psi1,params);
  //forward euler step
  hamiltonian_operator(params,psi,Hpsi);
  take_step(HALF_FWD,Hpsi,psi,params);

  //Make initial observation
  t = 0.5*params.dt;
  if(output.yes) user_observe (params, t,psi);
  //Main loop
  for(int i = 0;i<nt_outer;i++){
    for(int j = 0; j<nt_inner;j++){
      //take full step
      hamiltonian_operator(params,psi,Hpsi);
      take_step(FWD_STEP,Hpsi,psi1,params);
      //prepare for next iteration (put psi as current step and psi1 as previous)
      swap_vectors(psi1,psi,params);
      t+=params.dt;
    }
    //Print some results
	  if(output.yes) user_observe (params, t,psi);
  }


  //free variables

// Clean up

  if (params.rank == 0)
    {
      fclose (fp_results);
    }


  free (params.x2);
  free (params.y2);
  free (params.z2);

  // Print out final wave function
  if (output.wf == 1 || output.wf == 3)
    print_wf2 (params, psi, wf_text);


  // Finish clean up
  free (psi);
  free (Hpsi);
  free (psi1);

  free (params.x);
  free (params.y);
  free (params.z);
  return 0;
}


/*** Routines for initialization ***/

void
make_grid (parameters *params)
{
  /** Construct the three-dimensional spatial grid **/

  params->x = (double *) malloc (params->nx_local * sizeof (double));
  params->y = (double *) malloc (params->ny * sizeof (double));
  params->z = (double *) malloc (params->nz * sizeof (double));
  params->x2 = (double *) malloc (params->nx_local * sizeof (double));
  params->y2 = (double *) malloc (params->ny * sizeof (double));
  params->z2 = (double *) malloc (params->nz * sizeof (double));

  if (params->x == NULL || params->y == NULL || params->z == NULL
      || params->x2 == NULL || params->y2 == NULL || params->z2 == NULL)
    {
      fprintf (stderr, "Out of memory while creating spatial grid.\n");
      abort ();
    }

  params->dx = (params->x_max - params->x_min) / (double)(params->nx - 1);

  if (params->ny == 1)
    params->dy = 1.;
  else
    params->dy = (params->y_max - params->y_min) / (double)(params->ny - 1);

  if (params->nz == 1)
    params->dz = 1.;
  else
    params->dz = (params->z_max - params->z_min) / (double)(params->nz - 1);

  for (size_t i = 0; i < params->nx_local; ++i)
    {
      params->x[i] = params->x_min + (double)(i + params->nx0) * params->dx;
      params->x2[i] = SQUARE(params->x[i]);
    }

  for (size_t j = 0; j < params->ny; ++j)
    {
      params->y[j] = params->y_min + (double)j * params->dy;
      params->y2[j] = SQUARE(params->y[j]);
    }

  for (size_t k = 0; k < params->nz; ++k)
    {
      params->z[k] = params->z_min + (double)k * params->dz;
      params->z2[k] = SQUARE(params->z[k]);
    }

  return;
}



/*** Routines for input and output ***/

void
read_parameters (parameters * const params, output_flags * const output,
		 unsigned int * const nt_inner, unsigned int * const nt_outer,
		 const char * const parameter_file,
		 char * const results_file, char * const wf_text,
		 char * const wf_bin)
{
  /** Read parameters from the file specified on the command line **/

  if (params->rank != 0)
    {
      /* (This function should only be called by the rank=0 processor in a
	 multi-processor environment) */
      fprintf (stderr, "Function read_parameters called by a processor other than rank 0 (rank = %d).\n", params->rank);
      abort ();
    }

  bool error = false, have_results_file = false, have_wf_output_text = false,
    have_wf_output_bin = false, have_mass = false,
    have_nx = false, have_ny = false, have_nz = false, have_x_min = false,
    have_x_max = false, have_y_min = false, have_y_max = false,
    have_z_min = false, have_z_max = false, have_dt = false, have_nt = false,
    have_nprint = false;
  int units = 0, nt, nprint;

  // Open file for reading
  FILE *fp_in = fopen (parameter_file, "r");
  if (fp_in == NULL)
    {
      fprintf (stderr, "Could not open parameter file %s\n",
	       parameter_file);
      abort ();
    }

  // Set all output flags to false before reading parameters
  output->yes = false;
  output->norm = false;
  output->energy = false;
  output->x_avg = false;
  output->y_avg = false;
  output->z_avg = false;
  output->sx = false;
  output->sy = false;
  output->sz = false;
  output->autoc = false;
  output->user = false;

  /* Read parameters */
  char key[30], value[80];
  while (fscanf (fp_in, "%s = %s", key, value) != EOF)
    {
      if (!strcmp (key, "results_file"))
	{
	  strcpy (results_file, value);
	  have_results_file = true;
	}
      else if (!strcmp (key, "wf_output_text"))
	{
	  strcpy (wf_text, value);
	  have_wf_output_text = true;
	}
      else if (!strcmp (key, "wf_output_binary"))
	{
	  strcpy (wf_bin, value);
	  have_wf_output_bin = true;
	}
      else if (!strcmp (key, "units"))
	{
	  if (!strcmp (value, "AU"))
	    {
	      units = 0;
	    }
	  else if (!strcmp (value, "SI"))
	    {
	      units = 1;
	    }
	  else
	    {
	      fprintf (stderr, "In parameter file %s, value %s for key units unknown, ignored\n", parameter_file, value);
	    }
	}
      else if (!strcmp (key, "mass"))
	{
	  params->mass = atof (value);
	  have_mass = true;
	}
      else if (!strcmp (key, "nx"))
	{
	  params->nx = atoi (value);
	  have_nx = true;
	}
      else if (!strcmp (key, "ny"))
	{
	  params->ny = atoi (value);
	  have_ny = true;
	}
      else if (!strcmp (key, "nz"))
	{
	  params->nz = atoi (value);
	  have_nz = true;
	}
      else if (!strcmp (key, "x_min"))
	{
	  params->x_min = atof (value);
	  have_x_min = true;
	}
      else if (!strcmp (key, "x_max"))
	{
	  params->x_max = atof (value);
	  have_x_max = true;
	}
      else if (!strcmp (key, "y_min"))
	{
	  params->y_min = atof (value);
	  have_y_min = true;
	}
      else if (!strcmp (key, "y_max"))
	{
	  params->y_max = atof (value);
	  have_y_max = true;
	}
      else if (!strcmp (key, "z_min"))
	{
	  params->z_min = atof (value);
	  have_z_min = true;
	}
      else if (!strcmp (key, "z_max"))
	{
	  params->z_max = atof (value);
	  have_z_max = true;
	}
      else if (!strcmp (key, "dt"))
	{
	  params->dt = atof (value);
	  have_dt = true;
	}
      else if (!strcmp (key, "nt"))
	{
	  nt = atoi (value);
	  have_nt = true;
	}
      else if (!strcmp (key, "nprint"))
	{
	  nprint = atoi (value);
	  have_nprint = true;
	}
      else if (!strcmp (key, "output"))
	{
	  if (!strcmp (value, "norm"))
	    {
	      output->norm = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "energy"))
	    {
	      output->energy = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "x_avg"))
	    {
	      output->x_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "y_avg"))
	    {
	      output->y_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "z_avg"))
	    {
	      output->z_avg = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sx"))
	    {
	      output->sx = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sy"))
	    {
	      output->sy = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "sz"))
	    {
	      output->sz = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "autocorrelation"))
	    {
	      output->autoc = true;
	      output->yes = true;
	    }
	  else if (!strcmp (value, "user_defined"))
	    {
	      output->user = true;
	      output->yes = true;
	    }
	  else
	    {
	      fprintf (stderr, "In parameter file %s, output value '%s' unknown, ignored\n", parameter_file, value);
	    }
	}
      else
	{
	  fprintf (stderr, "In parameter file %s, parameter %s unknown, ignored\n", parameter_file, key);
	}
    }

  fclose (fp_in);

  /* Check for missing parameters */
  if (!have_mass)
    {
      fprintf (stderr, "Parameter mass missing from file %s\n",
	       parameter_file);
      error = true;
    }

  // 1D means x
  if (!have_nx)
    {
      fprintf (stderr, "Parameter nx missing from file %s\n",
	       parameter_file);
      fprintf (stderr, "(1D system must be defined along x, 2D system in xy plane)\n");
      error = true;
    }
  else if (!have_ny && !have_nz)
    {
      fprintf (stdout, "Parameter ny and nz missing from file %s, values set to 1\n",
	       parameter_file);
      params->ny = 1;
      params->nz = 1;
    }
  else if (!have_ny)
    {
      fprintf (stdout, "Parameter ny missing from file %s, value set to 1\n",
	       parameter_file);
      params->ny = 1;
    }
  else if (!have_nz)
    {
      fprintf (stdout, "Parameter nz missing from file %s, value set to 1\n",
	       parameter_file);
      params->nz = 1;
    }

  // 2D means xy
  if (have_nx && params->ny <= 1 && params->nz > 1)
    {
      fprintf (stderr, "2D system must be defined along the x and y axes\n");
      error = true;
    }

  if (!have_x_min)
    {
      fprintf (stderr, "Parameter x_min missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_x_max)
    {
      fprintf (stderr, "Parameter x_max missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_y_min)
    {
      if(params->ny > 1)
	{
	  fprintf (stderr, "Parameter y_min missing from file %s\n", parameter_file);
	  error = true;
	}
      else
	{
	  fprintf (stdout, "Parameter y_min missing from file %s, value set to 0.\n", parameter_file);
	  params->y_min = 0.;
	}
    }

  if (!have_y_max)
    {
      if (params->ny > 1)
	{
	  fprintf (stderr, "Parameter y_max missing from file %s\n",
		   parameter_file);
	  error = true;
	}
      else
	{
	  params->y_max = params->y_min;
	}
    }


  if (!have_z_min)
    {
      if(params->nz > 1)
	{
	  fprintf (stderr, "Parameter z_min missing from file %s\n",
		   parameter_file);
	  error = true;
	}
      else
	{
	  fprintf (stdout, "Parameter z_min missing from file %s, value set to 0.\n", parameter_file);
	  params->z_min = 0.;
	}
    }

  if (!have_z_max)
    {
      if (params->nz > 1)
	{
	  fprintf (stderr, "Parameter z_max missing from file %s\n",
		   parameter_file);
	  error = true;
	}
      else
	{
	  params->z_max = params->z_min;
	}
    }

  if (!have_dt)
    {
      fprintf (stderr, "Parameter dt missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_nt)
    {
      fprintf (stderr, "Parameter nt missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (!have_nprint && output->yes)
    {
      fprintf (stderr, "Parameter nprint missing from file %s\n",
	       parameter_file);
      error = true;
    }

  if (error) abort ();

  if (!have_results_file)
    {
      fprintf (stdout, "No results_file specified, using default value 'results'\n");
      strcpy (results_file, "results");
    }

  /* Check validity of the parameters */
  if (params->mass <= 0.)
    {
      fprintf (stderr, "mass = %g, must be greater than zero.\n", params->mass);
      error = true;
    }

  if (params->nx == 0)
    {
      fprintf (stderr, "nx has to be greater than zero\n");
      fprintf (stderr, "(1D system must be defined along x, 2D system in xy plane))\n");
      error = true;
    }
#ifdef MPI
  else if (params->nx < params->size)
    {
      fprintf (stderr, "fewer grid points in x than number of processors\n");
      error = true;
    }
#endif

  if (params->ny == 0)
    {
      if (params->nz > 1)
	{
	  fprintf (stderr, "ny has to be greater than zero\n");
	  fprintf (stderr, "(2D system must be defined in xy plane))\n");
	  error = true;
	}
      else
	{
	  fprintf (stdout, "ny should be greater than zero\n");
	  fprintf (stdout, "Assuming 1D system and setting ny = 1\n");
	  params->ny = 1;
	}
    }

  if (params->nz == 0)
    {
      fprintf (stdout, "nz should be greater than zero\n");
      fprintf (stdout, "Assuming reduced dimensionality and setting nz = 1\n");
      params->nz = 1;
    }

  if (params->x_max <= params->x_min)
    {
      fprintf (stderr, "x_max = %lg must be greater than x_min = %lg\n",
	       params->x_max, params->x_min);
      error = true;
    }

  if (params->ny > 1 && params->y_max <= params->y_min)
    {
      fprintf (stderr, "y_max = %lg must be greater than y_min = %lg\n",
	       params->y_max, params->y_min);
      error = true;
    }

  if (params->nz > 1 && params->z_max <= params->z_min)
    {
      fprintf (stderr, "z_max = %lg must be greater than z_min = %lg\n",
	       params->z_max, params->z_min);
      error = true;
    }

  if (error) abort ();

  if (output->yes && nprint > nt)
    {
      fprintf (stdout, "* Warning *\n");
      fprintf (stdout, "nprint = %d is greater than nt = %d.\n",
	       nprint, nt);
      output->yes = 0;
    }
  else if (output->yes && nt % nprint != 0)
    {
      fprintf (stdout, "* Warning *\n");
      fprintf (stdout, "nprint = %d is not an integer multiple of nt = %d.\n",
	       nprint, nt);

      nt = (nt / nprint + 1) * nprint;

      fprintf (stdout, "nt is modified to %d.\n", nt);
    }

  // Figure out the final output of the wave function
  output->wf = 0;

  if (have_wf_output_text)
    output->wf = 1;

  if (have_wf_output_bin)
    output->wf += 2;

  if (!output->yes)
    {
      if (output->wf == 0)
	{
	  fprintf (stderr, "Program will produce no output.\nAborting...\n");
	  abort ();
	}
      else
	{
	  fprintf (stdout, "No output during time evolution,\nonly final wave function will be printed out\n");

	  *nt_inner = nt - 1;
	  *nt_outer = 1;
	}
    }
  else
    {
      *nt_inner = nprint - 1;
      *nt_outer = nt/nprint;
    }

  if (output->y_avg && params->ny == 1)
    {
      fprintf (stdout, "No grid along y, so y_avg will not be calculated\n");
      output->y_avg = 0;
    }

  if (output->z_avg && params->nz == 1)
    {
      fprintf (stdout, "No grid along z, so z_avg will not be calculated\n");
      output->z_avg = 0;
    }

  if (output->sy && params->ny == 1)
    {
      fprintf (stdout, "No grid along y, so sy will not be calculated\n");
      output->sy = 0;
    }

  if (output->sz && params->nz == 1)
    {
      fprintf (stdout, "No grid along z, so sz will not be calculated\n");
      output->sz = 0;
    }

  /* Derived parameters */
  params->n = params->nx * params->ny * params->nz;

  if (units == 1)
    {
      params->hbar = 1.054571726e-34; // Planck constant over 2 pi (2010 CODATA)
    }
  else
    {
      params->hbar = 1.;
    }

  /* Print out all parameters */
  fprintf (stdout, "*** Parameters ***\n");
  if (units == 0)
    {
      fprintf (stdout, "units = AU\n");
    }
  else
    {
      fprintf (stdout, "units = SI\n");
    }
  fprintf (stdout, "mass = %.16g\n", params->mass);
  fprintf (stdout, "nx = %zu\n", params->nx);
  fprintf (stdout, "ny = %zu\n", params->ny);
  fprintf (stdout, "nz = %zu\n", params->nz);
  fprintf (stdout, "x_min = %.16g\n", params->x_min);
  fprintf (stdout, "x_max = %.16g\n", params->x_max);
  fprintf (stdout, "y_min = %.16g\n", params->y_min);
  fprintf (stdout, "y_max = %.16g\n", params->y_max);
  fprintf (stdout, "z_min = %.16g\n", params->z_min);
  fprintf (stdout, "z_max = %.16g\n", params->z_max);
  fprintf (stdout, "dt = %.16g\n", params->dt);
  fprintf (stdout, "nt = %d\n", nt);
  fprintf (stdout, "nprint = %d\n\n", nprint);

  if (output->yes)
    fprintf (stdout, "results_file = %s\n", results_file);

  if (output->wf == 1 || output->wf == 3)
    fprintf (stdout, "wf_output_text = %s\n", wf_text);

  if (output->wf == 2 || output->wf == 4)
    fprintf (stdout, "wf_output_bin = %s\n", wf_bin);

  if (output->yes)
    {
      fprintf (stdout, "\n");
      if (output->norm)
	fprintf (stdout, "output = norm\n");
      if (output->energy)
	fprintf (stdout, "output = energy\n");
      if (output->x_avg)
	fprintf (stdout, "output = x_avg\n");
      if (output->y_avg)
	fprintf (stdout, "output = y_avg\n");
      if (output->z_avg)
	fprintf (stdout, "output = z_avg\n");
      if (output->sx)
	fprintf (stdout, "output = sx\n");
      if (output->sy)
	fprintf (stdout, "output = sy\n");
      if (output->sz)
	fprintf (stdout, "output = sz\n");
      if (output->autoc)
	fprintf (stdout, "output = autocorrelation\n");
      if (output->user)
	fprintf (stdout, "output = user_defined\n");
    }

  return;
}

//print functions
void
print_wf2 (const parameters params, double complex *psi,
	  const char * const wf_text)
{
  /** Print out the wave function psi in a text file **/

  size_t index_x, index_xy, index;
  FILE *fp;

  if (params.rank == 0)
    {
      fp = fopen (wf_text, "w");
      if (fp == NULL)
	{
	  fprintf (stderr, "Could not open file %s for output.\n", wf_text);
	  abort ();
	}
    }

  bool print_y = (params.ny > 1);
  bool print_z = (params.nz > 1);

  if (params.rank == 0)
    {
      // Print header
      fprintf (fp, "# x");

      if (print_y)
	fprintf (fp, " y");

      if (print_z)
	fprintf (fp, " z");

      fprintf (fp, " psi_r psi_i psi2 \n");

      // Processor zero prints its wave function
      for (size_t i = 0; i < params.nx_local; ++i)
	{
	  index_x = i * params.nz * params.ny;
	  for (size_t j = 0; j < params.ny; ++j)
	    {
	      index_xy = j * params.nz + index_x;
	      for (size_t k = 0; k < params.nz; ++k)
		{
		  index = k + index_xy;

		  fprintf (fp, "%13.6e ", params.x[i]);

		  if (print_y)
		    fprintf (fp, "%13.6e ", params.y[j]);

		  if (print_z)
		    fprintf (fp, "%13.6e ", params.z[k]);

		  double psiR = creal(psi[index]);
      double psiI = cimag(psi[index]);
      double psi2 = psiR*psiR + psiI*psiI;
		  fprintf (fp, "%13.6e %13.6e %.6e\n",
			   psiR,psiI, psi2);
		}
	    }
	}
    }
  }
