#include "utils.h"
#include "nfft.h"

nfft_plan my_plan;
fftw_complex *f_hat_x, *f_hat_y, *f_hat_z;
fftw_complex *f_x, *f_y, *f_z;

void FElecRecipENUF_INIT_(int *ikspace, int *iLx,
	int *iLy, int *iLz, int *Nchargee){
	int N[3], n[3],kspace, Lx, Ly, Lz, kklimit, precision, mm;

	kspace = *ikspace;
	Lx = *iLx;
	Ly = *iLy;
	Lz = *iLz;

	kklimit = 2*(kspace);
	precision = 2;

	N[0]=kklimit;
	N[1]=kklimit;
	N[2]=kklimit;

	n[0]= kklimit>Lx?kklimit:Lx;
	n[1]= kklimit>Ly?kklimit:Ly;
	n[2]= kklimit>Lz?kklimit:Lz;

	mm = *Nchargee; //n[0]*n[1]*n[2];

	//nfft_init_2d(&my_plan,12,12,19);
	nfft_init_specific(&my_plan, 3, N, mm, n, precision,
		PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT |
		MALLOC_F | FFT_OUT_OF_PLACE | FFTW_INIT, FFTW_MEASURE |
		FFTW_DESTROY_INPUT);

	f_hat_x=(fftw_complex*)fftw_malloc((my_plan.N_L)*sizeof(fftw_complex));
	f_hat_y=(fftw_complex*)fftw_malloc((my_plan.N_L)*sizeof(fftw_complex));
	f_hat_z=(fftw_complex*)fftw_malloc((my_plan.N_L)*sizeof(fftw_complex));

	// force
	f_x =(fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
	f_y =(fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
	f_z =(fftw_complex*)fftw_malloc(my_plan.M*sizeof(fftw_complex));
}



void FElecRecipENUF_(int *ikspace, int *Npmaxx, int *Nchargee, double *alphaa,
	double *dxINV, double *dyINV, double *dzINV,
	int idcharge[*Nchargee], double charge[*Npmaxx],
	double rX[*Npmaxx], double rY[*Npmaxx], double rZ[*Npmaxx],
	double fkx[*Npmaxx],double fky[*Npmaxx],double fkz[*Npmaxx],
	double *nfft_engg){

	int j, jjj, N[3], kspace, Ncharge, kklimit, Npmax;
	double dLxINV,dLyINV,dLzINV,alpha_chg;

	kspace = *ikspace;
	Npmax = *Npmaxx;
	Ncharge = *Nchargee;

	dLxINV = *dxINV;
	dLyINV = *dyINV;
	dLzINV = *dzINV;

	alpha_chg = (*alphaa);
	kklimit = 2*(kspace);
	N[0]=kklimit;
	N[1]=kklimit;
	N[2]=kklimit;

	// init pseudo random nodes
	for(j=0; j<Ncharge; j++){
		jjj = idcharge[j] - 1 ;
		my_plan.x[3*j+0] = rX[jjj]*dLxINV;
		my_plan.x[3*j+1] = rY[jjj]*dLyINV;
		my_plan.x[3*j+2] = rZ[jjj]*dLzINV;
		my_plan.f[j][0]  = charge[jjj];
		my_plan.f[j][1]  = 0.0;
	}

	for(j=Ncharge ; j<my_plan.M; j++){
		my_plan.x[3*j+0]=0.0;
		my_plan.x[3*j+1]=0.0;
		my_plan.x[3*j+2]=0.0;
		my_plan.f[j][0] =0.0;
		my_plan.f[j][1] =0.0;
	}

	// if(my_plan.nfft_flags & PRE_PSI)
	nfft_precompute_psi(&my_plan);

	nfft_transposed(&my_plan);

	double w, gu, gf, uSumF, virSumF;
	int zlim, ylim, xlim;

	w = (M_PI*dLxINV/alpha_chg)*(M_PI*dLxINV/alpha_chg);
	gu = 0.5*dLxINV/M_PI;
	gf = 2.0*dLxINV*dLxINV;

	uSumF = 0.0;
	virSumF = 0.0;

	zlim = N[2]/2;
	ylim = N[1]/2;
	xlim = N[0]/2;

	int nx,ny,nz,nvv;
	int const zhigh = 1;
	for ( nz = -zlim; nz < zhigh; ++nz){
		int const yhigh = (nz == 0 ? 1 : ylim);
		for ( ny = -ylim; ny < yhigh; ++ny){
			int const xhigh = ( nz == 0 && ny == 0 ? 1 : xlim);
			for ( nx = -xlim; nx < xhigh; ++nx){

				nvv = nx*nx + ny*ny + nz*nz;
				if (nvv == 0 || nvv > (kspace*kspace))
					continue;

				double fMult,vnSq;
				vnSq  = 1.0*(nx*nx+ny*ny+nz*nz);
				fMult = exp(-w*vnSq)/vnSq;

				int kk;
				kk = (nx + xlim);
				kk *= 2*ylim;
				kk += (ny + ylim);
				kk *= 2*zlim;
				kk += (nz + zlim);

				// Structure factor components
				double const sumC =  my_plan.f_hat[kk][0];
				double const sumS = -my_plan.f_hat[kk][1];

				// Energy for +-[nx,ny,nz]
				uSumF += fMult*((sumC*sumC) + (sumS*sumS));
			}
		}
	}
	uSumF *= 2.0*gu; // Symmetry gives factor 2.0
	*nfft_engg = uSumF;


	// F-coeff for interactions (symmetry not used)
	for ( nz = -zlim; nz < zlim; ++nz){
		for ( ny = -ylim; ny < ylim; ++ny){
			for ( nx = -xlim; nx < xlim; ++nx){
				nvv = nx*nx + ny*ny + nz*nz;

				int kk;
				kk  = (nx + xlim);
				kk *= 2*ylim;
				kk += (ny + ylim);
				kk *= 2*zlim;
				kk += (nz + zlim);

				double vnSq,fMult;
				vnSq = 1.0*nvv;
				fMult = 0.0;

				if (nvv <= (kspace*kspace) && nvv != 0)
					fMult = exp (- w * vnSq) / vnSq;

				double const sumC = fMult*my_plan.f_hat[kk][0];
				double const sumS = fMult*my_plan.f_hat[kk][1];

				f_hat_x[kk][0] =  nx*sumC;
				f_hat_x[kk][1] =  nx*sumS;
				f_hat_y[kk][0] =  ny*sumC;
				f_hat_y[kk][1] =  ny*sumS;
				f_hat_z[kk][0] =  nz*sumC;
				f_hat_z[kk][1] =  nz*sumS;
			}
		}
	}


	// Components in x, y, and z direction
	SWAPC (my_plan.f_hat, f_hat_x);
	nfft_conjugated(&my_plan);
	SWAPC (my_plan.f, f_x);

	SWAPC (my_plan.f_hat, f_hat_y);
	nfft_conjugated(&my_plan);
	SWAPC (my_plan.f, f_y);

	SWAPC (my_plan.f_hat, f_hat_z);
	nfft_conjugated(&my_plan);
	SWAPC (my_plan.f, f_z);

	double gnf = 0.0;

	for (j=0; j<Ncharge; j++){
		jjj = idcharge[j] - 1;
		gnf = charge[jjj]*gf;
		fkx[jjj] = gnf*f_x[j][1];
		fky[jjj] = gnf*f_y[j][1];
		fkz[jjj] = gnf*f_z[j][1];
	}

}


void FElecRecipENUF_FINALIZE_(){
	nfft_finalize(&my_plan);
	free(f_hat_x);
	free(f_hat_y);
	free(f_hat_z);
	free(f_x);
	free(f_y);
	free(f_z);
}
