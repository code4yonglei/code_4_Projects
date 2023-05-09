#include <cufft.h>
#include "cunfft.h"

__device__ inline REAL phi(REAL x, int n, REAL b)
{
    REAL phi = exp(-pow(x*n,2)/b)/sqrt(PI*b);
    return phi;
}


REAL phi_hut(cunfft_plan* cunfft, int k, int d)
{
    REAL phi_hut = exp(-(pow(PI*(k)/cunfft->n[d],2.0) * cunfft->b));
    return phi_hut;
}


extern "C" void cunfft_init(REAL box_length, int N0, int N1, int N2, REAL sigma, int precision, cunfft_plan* cunfft, gpu_malloc* gpuMalloc)
{
    cunfft->box_length = box_length;

    cunfft->N = (int*)malloc(3 * sizeof(int));
    cunfft->N[0] = N0;
    cunfft->N[1] = N1;
    cunfft->N[2] = N2;

    cunfft->sigma = sigma;
    cunfft->n = (int*)malloc(3 * sizeof(int));
    for(int t = 0; t < 3; t++)
	cunfft->n[t] = cunfft->N[t] * sigma;
    
    cunfft->precision = precision;
    cunfft->N_L = cunfft->N[0] * cunfft->N[1] * cunfft->N[2];
    cunfft->n_L = cunfft->n[0] * cunfft->n[1] * cunfft->n[2];
    
    cunfft->b = (2*sigma*precision) / ((2*sigma-1)*PI);        
    cunfft->spatial_expansion = 2*precision + 2;

    cunfft->c_phi_inv = (REAL**)malloc(3 * sizeof(REAL*));
    for(int t=0; t<3; t++)
	{
		cunfft->c_phi_inv[t] = (REAL*)malloc(cunfft->N[t] * sizeof(REAL));
		for(int ks=0; ks<cunfft->N[t]; ks++)
			cunfft->c_phi_inv[t][ks]= 1.0/(phi_hut(cunfft,ks-cunfft->N[t]/2,t));
    }
    cunfft->position = (REAL*)malloc(3 * cunfft->number * sizeof(REAL));
    cunfft->f = (REAL*)malloc(2 * cunfft->number * sizeof(REAL));
    cunfft->f_hat = (REAL*)malloc(2 * cunfft->N_L * sizeof(REAL));

    cunfft->blocks_M = (cunfft->number + THREAD_NUM - 1) / THREAD_NUM;
    cunfft->blocks_N_L = (cunfft->N_L + THREAD_NUM - 1) / THREAD_NUM;

    cudaMalloc((void**) &gpuMalloc->gpu_position, 3 * cunfft->number * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_spatial_begin, 3 * cunfft->number * sizeof(int));
    cudaMalloc((void**) &gpuMalloc->gpu_psi, 3 * cunfft->number * cunfft->spatial_expansion * sizeof(REAL));   
    cudaMalloc((void**) &gpuMalloc->gpu_g, cunfft->n_L * sizeof(Complex));
    cudaMalloc((void**) &gpuMalloc->gpu_n, 3 * sizeof(int));    
    cudaMalloc((void**) &gpuMalloc->gpu_f, 2 * cunfft->number * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_f_hat, 2 * cunfft->N_L * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv0, cunfft->N[0] * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv1, cunfft->N[1] * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv2, cunfft->N[2] * sizeof(REAL));
    cudaMalloc((void**) &gpuMalloc->gpu_N, 3 * sizeof(int)); 

    cudaMemcpy(gpuMalloc->gpu_n, cunfft->n, 3 * sizeof(int), cudaMemcpyHostToDevice);  
    cudaMemcpy(gpuMalloc->gpu_N, cunfft->N, 3 * sizeof(int), cudaMemcpyHostToDevice);  
    cudaMemcpy(gpuMalloc->gpu_c_phi_inv0, cunfft->c_phi_inv[0], cunfft->N[0] * sizeof(REAL), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_c_phi_inv1, cunfft->c_phi_inv[1], cunfft->N[1] * sizeof(REAL), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_c_phi_inv2, cunfft->c_phi_inv[2], cunfft->N[2] * sizeof(REAL), cudaMemcpyHostToDevice); 
}


extern "C" void cunfft_finalize(cunfft_plan* cunfft, gpu_malloc* gpuMalloc)
{
    free(cunfft->N);
    free(cunfft->n);
    free(cunfft->c_phi_inv);
    free(cunfft->position);
    free(cunfft->f);
    free(cunfft->f_hat);

    cudaFree(gpuMalloc->gpu_position);    
    cudaFree(gpuMalloc->gpu_g);
    cudaFree(gpuMalloc->gpu_f);
    cudaFree(gpuMalloc->gpu_f_hat);
    cudaFree(gpuMalloc->gpu_c_phi_inv0);
    cudaFree(gpuMalloc->gpu_c_phi_inv1);
    cudaFree(gpuMalloc->gpu_c_phi_inv2);
    cudaFree(gpuMalloc->gpu_N); 
    cudaFree(gpuMalloc->gpu_n);
    cudaFree(gpuMalloc->gpu_spatial_begin);
    cudaFree(gpuMalloc->gpu_psi);
}


/*
static __device__ inline REAL atomicAddReal(REAL* address, REAL value){
    REAL old = value;  
    REAL ret = atomicExch(address, 0.0f);
    REAL new_old = ret+old;
    while ((old = atomicExch(address, new_old))!=0.0f){
	new_old = atomicExch(address, 0.0f);
	new_old += old;
    }
    return ret;
}
*/


#ifdef SINGLE_PRECISION
static __device__ inline float atomicAddReal(float* address, float value)
{
    float old = value;  
    float ret = atomicExch(address, 0.0f);
    float new_old = ret+old;
    while ((old = atomicExch(address, new_old))!=0.0f)
	{
		new_old = atomicExch(address, 0.0f);
		new_old += old;
    }
    return ret;
}
#else 
static __device__ inline double atomicAddReal(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do
	{
        assumed = old;
        old = atomicCAS(address_as_ull,assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    }while (assumed != old);
    return __longlong_as_double(old);
}
#endif


__global__ static void spread(REAL* gpu_x, REAL* gpu_f, int* gpu_spatial_begin, REAL* gpu_psi, int* gpu_n, int expansion, int gpu_icon_num, CUFFTCOMPLEX* gpu_g, int precision, REAL box_length, REAL b)
{
    __shared__ int  n0, n1, n2; 
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;
    int i, i3, t0, t1, t2, l0, l1, l2, l_plain, t, u, l, lj;
    REAL phi_product, f_real, c;
    
    if(tid == 0)
         n0 = gpu_n[0]; n1 = gpu_n[1]; n2 = gpu_n[2];

    __syncthreads();
    i = bid * THREAD_NUM + tid;
    if(i < gpu_icon_num)
	{
        i3 = i * 3;

        for(t=0; t<3; t++)
		{
			c = gpu_x[i3+t] * gpu_n[t] / box_length;
			u = c;
			if(c < 0)
				u = u-1;
			u = u - precision;
			gpu_spatial_begin[i3+t] = u;
			
			for(l=u, lj=0; lj < expansion; l++, lj++)
				gpu_psi[(i3+t)*expansion+lj] = phi(gpu_x[i3+t]/box_length-((REAL)l)/gpu_n[t], gpu_n[t], b); 
		}

        f_real = gpu_f[2*i];
         
        for(t0=0,l0=gpu_spatial_begin[i3]; t0<(expansion); t0++,l0++)
	    for(t1=0,l1=gpu_spatial_begin[i3+1]; t1<(expansion); t1++,l1++)
	        for(t2=0,l2=gpu_spatial_begin[i3+2]; t2<(expansion); t2++,l2++)
			{             
				phi_product = gpu_psi[(i3)*(expansion)+t0]*gpu_psi[(i3+1)*(expansion)+t1]*gpu_psi[(i3+2)*(expansion)+t2];
                l_plain = (((l0+n0)%n0)*n1 + (l1+n1)%n1)*n2 + (l2+n2)%n2;                  
                atomicAddReal(&(gpu_g[l_plain].x), phi_product* f_real);                 
	        }
    }
}


__global__ static void scale(REAL* f_hat, REAL* c_phi_inv0, REAL* c_phi_inv1, REAL* c_phi_inv2, int* N, int* n, CUFFTCOMPLEX* g_hat)
{
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;

    int i, t0, t1, t2, k[3], ks[3], k_plain, ks_plain, N0=N[0],  N1 = N[1],  N2=N[2], n0=n[0], n1=n[1], n2=n[2];
    REAL c_phi_inv_k;

    i = bid * THREAD_NUM + tid;  

    if(i < N0*N1*N2)
	{  
        t0 = i / (N1*N2);
        t1 = i / N2 - t0 * N1;
        t2 = i % N2;
       
        if(t0 < N0/2)
		{
            k[0] = t0;
            ks[0] = N0/2 + t0;
        }
		else
		{
    	    k[0] = n0 - N0 + t0;
    	    ks[0] = t0 - N0/2;
        }

        if(t1 < N1/2)
		{
            k[1] = t1;
            ks[1] = N1/2 + t1;
        }
		else
		{
    	    k[1] = n1 - N1 + t1;
    	    ks[1] = t1 - N1/2;
        }

        if(t2 < N2/2)
		{
            k[2] = t2;
            ks[2] = N2/2 + t2;
        }
		else
		{
    	    k[2] = n2 - N2 + t2;
    	    ks[2] = t2 - N2/2;
        }

        c_phi_inv_k = c_phi_inv0[ks[0]]* c_phi_inv1[ks[1]]* c_phi_inv2[ks[2]];
        ks_plain = (ks[0]*N1 + ks[1])*N2 + ks[2];
        k_plain = (k[0]*n1 + k[1])*n2 + k[2];
        
        f_hat[2*ks_plain] = g_hat[k_plain].x * c_phi_inv_k;
        f_hat[2*ks_plain+1] = g_hat[k_plain].y * c_phi_inv_k;         
    }

}


__global__ static void subdivide(REAL* f_hat, REAL* c_phi_inv0, REAL* c_phi_inv1, REAL* c_phi_inv2, int* N, int* n, CUFFTCOMPLEX* g_hat)
{
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;
    int i, t0, t1, t2, k[3], ks[3], k_plain, ks_plain, N0=N[0],  N1=N[1],  N2=N[2], n0=n[0], n1=n[1], n2=n[2];
    REAL c_phi_inv_k;

    i = bid * THREAD_NUM + tid;  
    if(i < N0*N1*N2)
	{  
        t0 = i / (N1*N2);
        t1 = i / N2 - t0 * N1;
        t2 = i % N2;
       
        if(t0 < N0/2)
		{
            k[0] = t0;
            ks[0] = N0/2 + t0;
        }
		else
		{
    	    k[0] = n0 - N0 + t0;
    	    ks[0] = t0 - N0/2;
        }

        if(t1 < N1/2)
		{
            k[1] = t1;
            ks[1] = N1/2 + t1;
        }
		else
		{
    	    k[1] = n1 - N1 + t1;
    	    ks[1] = t1 - N1/2;
        }

        if(t2 < N2/2)
		{
            k[2] = t2;
            ks[2] = N2/2 + t2;
        }
		else
		{
    	    k[2] = n2 - N2 + t2;
    	    ks[2] = t2 - N2/2;
        }

        c_phi_inv_k = c_phi_inv0[ks[0]]* c_phi_inv1[ks[1]]* c_phi_inv2[ks[2]];
        ks_plain = (ks[0]*N1 + ks[1])*N2 + ks[2];
        k_plain = (k[0]*n1 + k[1])*n2 + k[2];

        g_hat[k_plain].x = f_hat[2*ks_plain]*c_phi_inv_k;
        g_hat[k_plain].y = f_hat[2*ks_plain+1]*c_phi_inv_k;
    }
}


__global__ static void interpolate(REAL box_length, REAL* gpu_x, int precision, REAL b, REAL* gpu_f, int* gpu_spatial_begin, REAL* gpu_psi, int* gpu_n, int expansion, int icon_num, CUFFTCOMPLEX* gpu_g)
{
    __shared__ int  n0, n1, n2; 
    const int tid = threadIdx.x;
    const int bid = blockIdx.x;
    int i, i3, t, t0, t1, t2, l, l0, l1, l2, lj, l_plain, u;
    REAL phi_product, c;
    
    if(tid == 0)
        n0 = gpu_n[0]; n1 = gpu_n[1]; n2 = gpu_n[2];
    __syncthreads();

    i = bid * THREAD_NUM + tid;  
    if(i < icon_num)
	{       
        i3 = i * 3;

        for(t=0; t<3; t++)
		{
			c = gpu_x[i3+t] * gpu_n[t] / box_length;
			u = c;
			if(c < 0)
				u = u-1;
			u = u - precision;
			gpu_spatial_begin[i3+t] = u;
			
			for(l=u, lj=0; lj < expansion; l++, lj++)
				gpu_psi[(i3+t)*expansion+lj] = phi(gpu_x[i3+t]/box_length-((REAL)l)/gpu_n[t], gpu_n[t], b); 
		}


		for(t0=0,l0=gpu_spatial_begin[i3]; t0<(expansion); t0++,l0++)
			for(t1=0,l1=gpu_spatial_begin[i3+1]; t1<(expansion); t1++,l1++)
				for(t2=0,l2=gpu_spatial_begin[i3+2]; t2<(expansion); t2++,l2++)
				{
					phi_product = gpu_psi[(i3)*(expansion)+t0]*gpu_psi[(i3+1)*(expansion)+t1]*gpu_psi[(i3+2)*(expansion)+t2];
                    l_plain = (((l0+n0)%n0)*n1 + (l1+n1)%n1)*n2 + (l2+n2)%n2; 		                        
                    gpu_f[2*i] += phi_product* gpu_g[l_plain].x;                    
                    gpu_f[2*i+1] += phi_product* gpu_g[l_plain].y;		                         			    
				}   
    }
}


extern "C" void cunfft_forward(cunfft_plan *cunfft, gpu_malloc *gpuMalloc)
{
    cufftHandle plan;
    cufftPlan3d(&plan, cunfft->n[0], cunfft->n[1], cunfft->n[2], CUFFT_TRANSFORM_TYPE);
    cudaMemcpy(gpuMalloc->gpu_position, cunfft->position, 3 * cunfft->number * sizeof(REAL), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_f, cunfft->f, 2 * cunfft->number * sizeof(REAL), cudaMemcpyHostToDevice);  
    cudaMemset(gpuMalloc->gpu_g, 0, cunfft->n_L * sizeof(Complex));  
    CUFFTCOMPLEX* gpu_g = (CUFFTCOMPLEX *)gpuMalloc->gpu_g; 

    spread<<<cunfft->blocks_M, THREAD_NUM>>>(gpuMalloc->gpu_position, gpuMalloc->gpu_f, gpuMalloc->gpu_spatial_begin, gpuMalloc->gpu_psi, gpuMalloc->gpu_n, cunfft->spatial_expansion, cunfft->number, gpu_g, cunfft->precision, cunfft->box_length, cunfft->b);

    CUFFTEXEC(plan, gpu_g, gpu_g, CUFFT_FORWARD);

    scale<<<cunfft->blocks_N_L, THREAD_NUM>>>(gpuMalloc->gpu_f_hat, gpuMalloc->gpu_c_phi_inv0, gpuMalloc->gpu_c_phi_inv1, gpuMalloc->gpu_c_phi_inv2, gpuMalloc->gpu_N, gpuMalloc->gpu_n, gpu_g);

    cudaMemcpy(cunfft->f_hat, gpuMalloc->gpu_f_hat, 2 * cunfft->N_L * sizeof(REAL), cudaMemcpyDeviceToHost);
    
    cufftDestroy(plan);
    return;
}      
    

extern "C" void cunfft_backward(cunfft_plan *cunfft, gpu_malloc *gpuMalloc)
{
    cufftHandle plan;
    cufftPlan3d(&plan, cunfft->n[0], cunfft->n[1], cunfft->n[2], CUFFT_TRANSFORM_TYPE);
    cudaMemcpy(gpuMalloc->gpu_position, cunfft->position, 3 * cunfft->number * sizeof(REAL), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_f_hat, cunfft->f_hat, 2 * cunfft->N_L * sizeof(REAL), cudaMemcpyHostToDevice); 
    cudaMemset(gpuMalloc->gpu_g, 0, cunfft->n_L * sizeof(Complex)); 
    CUFFTCOMPLEX* gpu_g_hat = (CUFFTCOMPLEX *)gpuMalloc->gpu_g;
 
    subdivide<<<cunfft->blocks_N_L, THREAD_NUM>>>(gpuMalloc->gpu_f_hat, gpuMalloc->gpu_c_phi_inv0, gpuMalloc->gpu_c_phi_inv1, gpuMalloc->gpu_c_phi_inv2, gpuMalloc->gpu_N, gpuMalloc->gpu_n, gpu_g_hat);   
     
    CUFFTEXEC(plan, gpu_g_hat, gpu_g_hat, CUFFT_INVERSE);
    
    cudaMemset(gpuMalloc->gpu_f, 0, 2 * cunfft->number * sizeof(REAL));  

    interpolate<<<cunfft->blocks_M, THREAD_NUM>>>(cunfft->box_length, gpuMalloc->gpu_position, cunfft->precision, cunfft->b, gpuMalloc->gpu_f, gpuMalloc->gpu_spatial_begin, gpuMalloc->gpu_psi, gpuMalloc->gpu_n, cunfft->spatial_expansion, cunfft->number, gpu_g_hat); 

    cudaMemcpy(cunfft->f, gpuMalloc->gpu_f, 2 * cunfft->number * sizeof(REAL), cudaMemcpyDeviceToHost);

    cufftDestroy(plan);
    return;
}
