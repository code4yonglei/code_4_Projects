#include <cufft.h>
#include "Info.cuh"
#include <math.h>

#ifndef __ENUFFORCE_CUH__
#define __ENUFFORCE_CUH__

typedef struct{
        
    unsigned int number;    
    Real alpha;
    Real Uk; 
    int kc;
     
    int3 N;
    int3 n;
	Real3 minv;
    int N_L;
    int n_L;
    int precision;        
    int spatial_expansion;  

    Real sigma;    
    Real b;    
    Real *c_phi_inv0;
    Real *c_phi_inv1;
    Real *c_phi_inv2;
	Real *c_coeff;
	
    int blocks_M;
    int blocks_N_L;
    Real *uk;
	unsigned int block_size;
} recip_plan;


// parameters used for the calculatio of EE in reciprocal space
typedef struct{
    CUFFTCOMPLEX* gpu_g;
    Real *gpu_uk; 

    Real *gpu_c_phi_inv0; 
    Real *gpu_c_phi_inv1; 
    Real *gpu_c_phi_inv2;
   
    CUFFTCOMPLEX *gpu_g_hatx;
    CUFFTCOMPLEX *gpu_g_haty; 
    CUFFTCOMPLEX *gpu_g_hatz; 

} gpu_malloc;


cudaError_t cuenuf(Real4* d_force,
					Real4* d_pos,
					Real* d_charge,
					const BoxSize& box,
					unsigned int* d_group_members,
					unsigned int group_size,
					cufftHandle plan, 
					recip_plan *recip, 
					gpu_malloc *gpuMalloc);

cudaError_t fix_exclusions2(Real4 *d_force,
                           ForceLog& force_log,
                           const Real4 *d_pos,
                           const Real *d_charge,
                           const BoxSize& box,
                           const unsigned int *d_n_ex,
                           const unsigned int *d_exlist,
                           const Index2D& nex,
                           Real m_kappa,
                           unsigned int *d_group_members,
                           unsigned int group_size,
                           int block_size);

#endif
