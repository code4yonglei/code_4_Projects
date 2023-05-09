#define PI 3.1415926535897932384
#define THREAD_NUM 256

#ifdef SINGLE_PRECISION
    #define REAL float
    #define CUFFTCOMPLEX cufftComplex
    #define CUFFTEXEC cufftExecC2C
    #define CUFFT_TRANSFORM_TYPE CUFFT_C2C
    #define MPI_REAL_TYPE MPI_FLOAT
#else
    #define REAL double
    #define CUFFTCOMPLEX cufftDoubleComplex
    #define CUFFTEXEC cufftExecZ2Z
    #define CUFFT_TRANSFORM_TYPE CUFFT_Z2Z
    #define MPI_REAL_TYPE MPI_DOUBLE
#endif

typedef struct{
    REAL real;
    REAL image;
} Complex;

typedef struct{
    int number;
    REAL* position;
    REAL* f;
    int box_length;  

    int *N;
    int *n;
    int N_L;
    int n_L;
    int precision;        
    int spatial_expansion;    

    REAL sigma;    
    REAL b;    
    REAL **c_phi_inv;
    REAL *f_hat;

    int blocks_M;
    int blocks_N_L;
} cunfft_plan;


typedef struct{

    REAL *gpu_position;
    
    int *gpu_n;  
    int *gpu_N;  
    REAL *gpu_psi;
    int *gpu_spatial_begin;
    Complex* gpu_g;

    REAL *gpu_f;
    REAL *gpu_f_hat; 
    REAL *gpu_c_phi_inv0; 
    REAL *gpu_c_phi_inv1; 
    REAL *gpu_c_phi_inv2;
    
} gpu_malloc;

