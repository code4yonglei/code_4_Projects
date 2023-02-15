#include "ENUFForce.cuh"

// Constant memory for gridpoint weighting
#define CONSTANT_SIZE 2048


#ifndef SINGLE_PRECISION
static __device__ inline double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull,
            assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif


__global__ static void spread(const Real4 *d_pos,
                              const Real *d_charge,
                              unsigned int *d_group_members,
							  int3 n, 
							  int expansion, 
							  int gpu_icon_num, 
							  CUFFTCOMPLEX* gpu_g, 
							  int precision, 
							  Real b,
							  Real binv,
							  Real bpinv,
							  Real3 minv)
{
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;

    int i = bid * blockDim.x + tid;
    if(i < gpu_icon_num){
        unsigned int idx = d_group_members[i];
        Real4 posi = d_pos[idx];

		Real3 c;
		int3 u;
		c.x = posi.x * minv.x;
		u.x = c.x;
		Real dx = c.x-Real(u.x);
		if(c.x < Real(0.0)){
			dx += Real(1.0);
			u.x -= 1;
		}
		u.x = u.x - precision;

		c.y = posi.y * minv.y;
		u.y = c.y;
		Real dy = c.y-Real(u.y);
		if(c.y < Real(0.0)){
			dy += Real(1.0);
			u.y -= 1;
		}
		u.y = u.y - precision;
		
		c.z = posi.z * minv.z;
		u.z = c.z;
		Real dz = c.z-Real(u.z);
		if(c.z < Real(0.0)){
			dz += Real(1.0);	
			u.z -= 1;
		}
		u.z = u.z - precision;

        Real f_real = d_charge[idx];
        Real bax = exp_gala(-dx*(dx+Real(2.0)*precision)*binv)*bpinv;
        Real bay = exp_gala(-dy*(dy+Real(2.0)*precision)*binv)*bpinv;
        Real baz = exp_gala(-dz*(dz+Real(2.0)*precision)*binv)*bpinv;

        Real plx = exp_gala(dx*Real(2.0)*binv);
        Real ply = exp_gala(dy*Real(2.0)*binv);
        Real plz = exp_gala(dz*Real(2.0)*binv);

		Real fax = bax;
        for(int t0=0; t0<expansion; t0++){
			int l0 = u.x + t0;
			if(l0 >= n.x) l0 -= n.x;
			if(l0 < 0)  l0 += n.x;
			Real psi0 = fax*d_coeff[t0];
			int l_plain0 = l0*n.y;
			Real fay = bay;
			for(int t1=0; t1<expansion; t1++){
				int l1 = u.y + t1;
				if(l1 >= n.y) l1 -= n.y;
				if(l1 < 0)  l1 += n.y;
				Real psi1 = fay*d_coeff[t1];
				int l_plain1 = l_plain0+l1;
				Real faz = baz;
				for(int t2=0; t2<expansion; t2++){
					int l2 = u.z + t2;
					if(l2 >= n.z) l2 -= n.z;
					if(l2 < 0)  l2 += n.z;
					Real psi2 = faz*d_coeff[t2];
                    Real phi_product = psi0*psi1*psi2;
                    int l_plain = l_plain1*n.z + l2;                  
                    atomicAdd(&(gpu_g[l_plain].x), phi_product* f_real); 
					faz *= plz;
				}
				fay *= ply;
			}
			fax *= plx;
		}
	}  
}


__global__ static void scale_subdivide(const BoxSize box, 
							Real* uk, 
							Real* c_phi_inv0, 
							Real* c_phi_inv1, 
							Real* c_phi_inv2, 
							int3 N, 
							int3 n, 
							CUFFTCOMPLEX* g_hat, 
                            CUFFTCOMPLEX* g_hatx, 
							CUFFTCOMPLEX* g_haty, 
							CUFFTCOMPLEX* g_hatz,
							Real inner, 
							int kc)
{
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;
    extern __shared__ Real shared[];

    int i = bid * blockDim.x + tid;  
    shared[tid] = 0;
	Real3 L = box.getL();
	Real3 LINV = box.getLINV();
    if(i < N.x*N.y*N.z){
        int t0 = i / (N.y*N.z);
        int t1 = i / N.z - t0 * N.y;
        int t2 = i % N.z;
		int3 k, ks;
        if(t0 < N.x/2){
            k.x = t0;
            ks.x = N.x/2 + t0;
		}else{
    	    k.x = n.x - N.x + t0;
    	    ks.x = t0 - N.x/2;
		}

        if(t1 < N.y/2){
            k.y = t1;
            ks.y = N.y/2 + t1;
		}else{
    	    k.y = n.y - N.y + t1;
    	    ks.y = t1 - N.y/2;
		}

        if(t2 < N.z/2){
            k.z = t2;
            ks.z = N.z/2 + t2;
		}else{
    	    k.z = n.z - N.z + t2;
    	    ks.z = t2 - N.z/2;
		}

        int k0 = ks.x-N.x/2;
        int k1 = ks.y-N.y/2;
        int k2 = ks.z-N.z/2;
        int k0_sq = k0*k0;
        int k1_sq = k1*k1;
        int k2_sq = k2*k2;
        int ks_square = k0_sq + k1_sq + k2_sq;
        
        if((ks_square!=0)&&(ks_square<=kc*kc)){
            Real c_phi_inv_k = c_phi_inv0[ks.x]* c_phi_inv1[ks.y]* c_phi_inv2[ks.z];
//            int ks_plain = (ks.x*N.y + ks.y)*N.z + ks.z;
            int k_plain = (k.x*n.y + k.y)*n.z + k.z;
			
			Real box_yz = L.y*L.z;
			Real box_xz = L.x*L.z;
			Real box_xy = L.x*L.y;
			
            Real omega = k0_sq * box_yz * box_yz + k1_sq * box_xz * box_xz + k2_sq * box_xy * box_xy;
            Real factor = exp_gala(-inner * omega) / omega;
            Real real = g_hat[k_plain].x * c_phi_inv_k;
            Real image = g_hat[k_plain].y * c_phi_inv_k;
            shared[tid] = factor * (real * real + image * image);

            Real f_real = factor * real;
            Real f_image = factor * image;       

            g_hatx[k_plain].x = f_real * k0 * LINV.x * c_phi_inv_k;
            g_hatx[k_plain].y = f_image * k0 * LINV.x * c_phi_inv_k;
            g_haty[k_plain].x = f_real * k1 * LINV.y * c_phi_inv_k;
            g_haty[k_plain].y = f_image * k1 * LINV.y * c_phi_inv_k;
            g_hatz[k_plain].x = f_real * k2 * LINV.z * c_phi_inv_k;
            g_hatz[k_plain].y = f_image * k2 * LINV.z * c_phi_inv_k;
		}
	}  

    int offset = blockDim.x / 2;

    __syncthreads();
    while(offset > 0) {
        if(tid < offset) {
            shared[tid] += shared[tid + offset];
        }
        offset >>= 1;
        __syncthreads();
    } 

    if(tid == 0) 
        uk[bid] = shared[0];       

}


__global__ static void interpolate(const Real4 *d_pos,
								   const Real *d_charge,
                                   Real volume,
                                   unsigned int *d_group_members,
                                   Real4 *d_force,
								   int3 n,							   
								   int expansion, 
								   int icon_num, 
								   CUFFTCOMPLEX* gpu_g_x, 
								   CUFFTCOMPLEX* gpu_g_y, 
								   CUFFTCOMPLEX* gpu_g_z,							  
								   int precision, 
								   Real b,
								   Real binv,
								   Real bpinv,
								   Real3 minv)
{
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;

    int i = bid * blockDim.x + tid;  
    if(i < icon_num){ 
		Real fx_image=0, fy_image=0, fz_image=0;
        unsigned int idx = d_group_members[i];
        Real prefix = Real(2.0)*d_charge[idx]*volume;
        Real4 posi = d_pos[idx];

		Real3 c;
		int3 u;
		c.x = posi.x * minv.x; 
		u.x = c.x;
		Real dx = c.x-Real(u.x);
		if(c.x < Real(0.0)){
			dx += Real(1.0);
			u.x -= 1;
		}
		u.x = u.x - precision;

		c.y = posi.y * minv.y; 
		u.y = c.y;
		Real dy = c.y-Real(u.y);
		if(c.y < Real(0.0)){
			dy += Real(1.0);
			u.y -= 1;
		}
		u.y = u.y - precision;
		
		c.z = posi.z * minv.z;
		u.z = c.z;
		Real dz = c.z-Real(u.z);
		if(c.z < Real(0.0)){
			dz += Real(1.0);	
			u.z -= 1;
		}
		u.z = u.z - precision;	
		
        Real bax = exp_gala(-dx*(dx+Real(2.0)*precision)*binv)*bpinv;
        Real bay = exp_gala(-dy*(dy+Real(2.0)*precision)*binv)*bpinv;
        Real baz = exp_gala(-dz*(dz+Real(2.0)*precision)*binv)*bpinv;

        Real plx = exp_gala(dx*Real(2.0)*binv);
        Real ply = exp_gala(dy*Real(2.0)*binv);
        Real plz = exp_gala(dz*Real(2.0)*binv);

		Real4 force = d_force[idx];
		Real fax = bax;
        for(int t0=0; t0<expansion; t0++){
			int l0 = u.x + t0;
			if(l0 >= n.x) l0 -= n.x;
			if(l0 < 0)  l0 += n.x;
			Real psi0 = fax*d_coeff[t0];
			int l_plain0 = l0*n.y;
			Real fay = bay;
			for(int t1=0; t1<expansion; t1++){
				int l1 = u.y + t1;
				if(l1 >= n.y) l1 -= n.y;
				if(l1 < 0)  l1 += n.y;
				Real psi1 = fay*d_coeff[t1];	
				int l_plain1 = l_plain0 + l1;
				Real faz = baz;
				for(int t2=0; t2<expansion; t2++){
					int l2 = u.z + t2;				
					if(l2 >= n.z) l2 -= n.z;
					if(l2 < 0)  l2 += n.z;
					Real psi2 = faz*d_coeff[t2];					
					Real phi_product = psi0*psi1*psi2;
                    int l_plain = l_plain1*n.z + l2; 		                        
                    fx_image += phi_product* gpu_g_x[l_plain].y;                    
                    fy_image += phi_product* gpu_g_y[l_plain].y;		     
                    fz_image += phi_product* gpu_g_z[l_plain].y;
					faz *= plz;					
				} 
				fay *= ply;
			}
			fax *= plx;
		}
        force.x += prefix * fx_image;       
        force.y += prefix * fy_image;        
        force.z += prefix * fz_image; 
		d_force[idx] = force;
	}
}


cudaError_t cuenuf(Real4* d_force,
					Real4* d_pos,
					Real* d_charge,
					const BoxSize& box,
					unsigned int* d_group_members,
					unsigned int group_size,
					cufftHandle plan, 
					recip_plan *recip, 
					gpu_malloc *gpuMalloc) 
{
    cudaMemcpyToSymbol(d_coeff, &(recip->c_coeff[0]),
		recip->spatial_expansion * sizeof(Real));
	
    cudaMemset(gpuMalloc->gpu_g, 0, recip->n_L * sizeof(CUFFTCOMPLEX));
	
    CUFFTCOMPLEX* gpu_g = (CUFFTCOMPLEX *)gpuMalloc->gpu_g; 

    spread<<<recip->blocks_M, recip->block_size>>>(d_pos, 
		d_charge, 
		d_group_members,
		recip->n, 
		recip->spatial_expansion, 
		recip->number, 
		gpu_g, 
		recip->precision, 
		recip->b,
		1.0/recip->b,
		1.0/sqrt(M_PI*recip->b),
		recip->minv);
	cudaThreadSynchronize();					
	
    CUFFTEXEC(plan, gpu_g, gpu_g, CUFFT_FORWARD);

    Real volume = box.getVolume();     
    Real inner = (Real)pow((M_PI / (volume*recip->alpha)), 2);        

    cudaMemset(gpuMalloc->gpu_g_hatx, 0, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMemset(gpuMalloc->gpu_g_haty, 0, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMemset(gpuMalloc->gpu_g_hatz, 0, recip->n_L * sizeof(CUFFTCOMPLEX)); 

    <<<recip->blocks_N_L, recip->block_size, recip->block_size * sizeof(Real)>>>(box, 
		gpuMalloc->gpu_uk, 
		gpuMalloc->gpu_c_phi_inv0, 
		gpuMalloc->gpu_c_phi_inv1, 
		gpuMalloc->gpu_c_phi_inv2, 
		recip->N, 
		recip->n, 
		gpu_g, 
        gpuMalloc->gpu_g_hatx, 
		gpuMalloc->gpu_g_haty, 
		gpuMalloc->gpu_g_hatz,
		inner, 
		recip->kc); 
	cudaThreadSynchronize();



    CUFFTEXEC(plan, gpuMalloc->gpu_g_hatx, gpuMalloc->gpu_g_hatx, CUFFT_INVERSE);
    CUFFTEXEC(plan, gpuMalloc->gpu_g_haty, gpuMalloc->gpu_g_haty, CUFFT_INVERSE);
    CUFFTEXEC(plan, gpuMalloc->gpu_g_hatz, gpuMalloc->gpu_g_hatz, CUFFT_INVERSE);
				cudaThreadSynchronize();

    interpolate<<<recip->blocks_M, recip->block_size>>>(d_pos,
		d_charge, 
		volume, 
		d_group_members, 
		d_force,
		recip->n,
		recip->spatial_expansion,
		recip->number,
		gpuMalloc->gpu_g_hatx, 
		gpuMalloc->gpu_g_haty, 
		gpuMalloc->gpu_g_hatz,
		recip->precision, 
		recip->b, 
		1.0/recip->b,
		1.0/sqrt(M_PI*recip->b),
		recip->minv);														
    return cudaSuccess;
}



//! The developer has chosen not to document this function
__global__ void gpu_fix_exclusions2_kernel(Real4 *d_force,
                                          ForceLog force_log,
                                          const Real4 *d_pos,
                                          const Real *d_charge,
                                          const BoxSize box,
                                          const unsigned int *d_n_neigh,
                                          const unsigned int *d_nlist,
                                          const Index2D nli,
                                          Real kappa,
                                          unsigned int *d_group_members,
                                          unsigned int group_size)
{
    // start by identifying which particle we are to handle
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx < group_size){
        unsigned int idx = d_group_members[group_idx];
        const Real sqrtpi = sqrt_gala(M_PI);
        unsigned int n_neigh = d_n_neigh[idx];
        Real4 posi = d_pos[idx];
        Real qi = d_charge[idx];
        // initialize the force
        Real4 force = d_force[idx];
		Real virial = Real(0.0);
		Real6 virial_matrix = ToReal6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		if(force_log.virial)
			virial = force_log.d_virial[idx];
		if(force_log.virial_matrix)
			virial_matrix = force_log.d_virial_matrix[idx];
        unsigned int cur_j = 0;

        // prefetch neighbor index
        unsigned int next_j = d_nlist[nli(idx, 0)];
		for (int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++){
			// read the current neighbor index (MEM TRANSFER: 4 bytes)
			// prefetch the next value and set the current one
			cur_j = next_j;
			next_j = d_nlist[nli(idx, neigh_idx+1)];
			
			// get the neighbor's position (MEM TRANSFER: 16 bytes)
			Real4 posj = d_pos[cur_j];
			Real qj = d_charge[cur_j];
			
			Real dx = posi.x - posj.x;
			Real dy = posi.y - posj.y;
			Real dz = posi.z - posj.z;
			
			box.minDisImage(dx, dy, dz);
			Real rsq = dx*dx + dy*dy + dz*dz;
			
			// calculate r squard (FLOPS: 5)
			Real r = sqrt_gala(rsq);
			Real qiqj = qi * qj;
			Real erffac = erf_gala(kappa * r) / r;
			Real force_divr = qiqj * (-Real(2.0) * exp_gala(-rsq * kappa * kappa) * kappa / (sqrtpi * rsq) + erffac / rsq);
			Real pair_eng = qiqj * erffac;

			if(force_log.virial)
				virial -= Real(1.0)/Real(6.0) * rsq * force_divr; 
			if(force_log.virial_matrix){
				Real force_div2r = Real(0.5) * force_divr;
				virial_matrix.x -=  dx * dx * force_div2r;   // press_tensor_xx
				virial_matrix.y -=  dx * dy * force_div2r;   // press_tensor_xy
				virial_matrix.z -=  dx * dz * force_div2r;   // press_tensor_xz
				virial_matrix.w -=  dy * dy * force_div2r;   // press_tensor_yy
				virial_matrix.m -=  dy * dz * force_div2r;   // press_tensor_yz
				virial_matrix.n -=  dz * dz * force_div2r;   // press_tensor_zz			
			}				

			force.x -= dx * force_divr;
			force.y -= dy * force_divr;
			force.z -= dz * force_divr;
			force.w -= pair_eng*Real(0.5);
		}
		d_force[idx] = force;
		if(force_log.virial)	
			force_log.d_virial[idx] = virial;
		if(force_log.virial_matrix)
			force_log.d_virial_matrix[idx] = virial_matrix;
	}
}



//! The developer has chosen not to document this function
cudaError_t fix_exclusions2(Real4 *d_force,
                           ForceLog& force_log,
                           const Real4 *d_pos,
                           const Real *d_charge,
						   const BoxSize& box,
                           const unsigned int *d_n_ex,
                           const unsigned int *d_exlist,
                           const Index2D& nex,
                           Real kappa,
                           unsigned int *d_group_members,
                           unsigned int group_size,
                           int block_size)
{
    dim3 grid( (int)ceil((double)group_size / (double)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);

    gpu_fix_exclusions2_kernel <<< grid, threads >>> (d_force,
		force_log, 
		d_pos,
		d_charge,
		box,
		d_n_ex,
		d_exlist,
		nex,
		kappa,
		d_group_members,
		group_size);
    return cudaSuccess;
}
