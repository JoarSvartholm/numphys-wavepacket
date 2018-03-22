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

}
