#include <boost/python.hpp>

using namespace boost::python;

#include "ENUFForce.h"
#include <stdexcept>

#include <boost/bind.hpp>
using namespace boost;
using namespace std;


ENUFForce::ENUFForce(boost::shared_ptr<AllInfo> all_info,
                     boost::shared_ptr<NeighborList> nlist,
                     boost::shared_ptr<ParticleSet> group)
                    : Force(all_info), m_nlist(nlist), m_group(group),
                    m_first_run(true), m_params_set(false)
{

    m_box_changed = true;
    m_boxchange_connection = m_perf_conf->connectBoxChange(bind(&ENUFForce::slotBoxChanged, this));
    m_block_size = 256;
    m_ObjectName = "ENUFForce";
    if (m_perf_conf->isRoot())
        cout << "INFO : "<<m_ObjectName<<" has been created" << endl;
}


ENUFForce::~ENUFForce()
{
    m_boxchange_connection.disconnect();
    cuenuf_finalize(&m_recip, &m_gpuMalloc);
    cufftDestroy(plan);
}


void ENUFForce::setParams(Real alpha, Real sigma,
                        int precision, int Nx, int Ny, int Nz)
{
    m_params_set = true;
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_alpha = alpha;
    m_sigma = sigma;
    m_precision=precision;
    
    Real* h_charge = m_basic_info->getCharge()->getArray(location::host, access::read);
    m_q = 0.f; // get system charge
    m_q2 = 0.0;
    unsigned int Np = m_basic_info->getN();
    for(int i = 0; i < (int)Np; i++) {
        m_q += h_charge[i];
        m_q2 += h_charge[i]*h_charge[i];
    }

    if(fabs(m_q) > 0.0)
        cout << "Notice: system in not neutral, the net charge is " << m_q << endl;
}


Real phi_hut(int n, int k, Real b){
    Real phi_hut = exp(-(pow(M_PI*(k)/n,2.0) * b));
    return phi_hut;
}


void ENUFForce::cuenuf_init(unsigned int group_size,
                Real alpha, Real sigma, int precision,
                int Nx, int Ny, int Nz,
                recip_plan* recip, gpu_malloc* gpuMalloc,
                unsigned int block_size)
{

    recip->number = group_size;
    recip->alpha = alpha;

    recip->N = ToInt3(Nx, Ny, Nz);
    recip->kc = Nx;
    if(Ny>recip->kc)
        recip->kc = Ny;
    if(Nz>recip->kc)
        recip->kc = Nz;
    recip->kc /= 2;

    recip->sigma = sigma;
    recip->n = ToInt3(Nx * sigma, Ny * sigma, Nz * sigma);

    recip->precision = precision;
    recip->N_L = recip->N.x * recip->N.y * recip->N.z;
    recip->n_L = recip->n.x * recip->n.y * recip->n.z;
    
    recip->b = (2*sigma*precision) / ((2*sigma-1)*M_PI);        
    recip->spatial_expansion = 2*precision + 2;

    recip->c_phi_inv0 = (Real*)malloc(recip->N.x * sizeof(Real));
    recip->c_phi_inv1 = (Real*)malloc(recip->N.y * sizeof(Real));
    recip->c_phi_inv2 = (Real*)malloc(recip->N.z * sizeof(Real));
    recip->c_coeff = (Real*)malloc(recip->spatial_expansion * sizeof(Real));

    for(int ks=0; ks<recip->N.x; ks++)
        recip->c_phi_inv0[ks]= 1.0/(phi_hut(recip->n.x, ks-recip->N.x/2, recip->b));
    for(int ks=0; ks<recip->N.y; ks++)
        recip->c_phi_inv1[ks]= 1.0/(phi_hut(recip->n.y, ks-recip->N.y/2, recip->b));
    for(int ks=0; ks<recip->N.z; ks++)
        recip->c_phi_inv2[ks]= 1.0/(phi_hut(recip->n.z, ks-recip->N.z/2, recip->b));
    for(int i=0; i<recip->spatial_expansion; i++)
        recip->c_coeff[i]= exp(-(i-precision)*(i-precision)/recip->b);
    recip->block_size = block_size;
    recip->blocks_M = (recip->number + block_size - 1) / block_size;
    recip->blocks_N_L =  (recip->N_L + block_size - 1) / block_size;
    recip->uk = (Real*)malloc(recip->blocks_N_L * sizeof(Real));

    cudaMalloc((void**) &gpuMalloc->gpu_g, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMalloc((void**) &gpuMalloc->gpu_uk, recip->blocks_N_L * sizeof(Real));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv0, recip->N.x * sizeof(Real));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv1, recip->N.y * sizeof(Real));
    cudaMalloc((void**) &gpuMalloc->gpu_c_phi_inv2, recip->N.z * sizeof(Real));

    cudaMalloc((void**) &gpuMalloc->gpu_g_hatx, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMalloc((void**) &gpuMalloc->gpu_g_haty, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMalloc((void**) &gpuMalloc->gpu_g_hatz, recip->n_L * sizeof(CUFFTCOMPLEX)); 

    cudaMemcpy(gpuMalloc->gpu_c_phi_inv0, recip->c_phi_inv0,
        recip->N.x * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_c_phi_inv1, recip->c_phi_inv1,
        recip->N.y * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuMalloc->gpu_c_phi_inv2, recip->c_phi_inv2,
        recip->N.z * sizeof(Real), cudaMemcpyHostToDevice);

    cudaMemset(gpuMalloc->gpu_g_hatx, 0, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMemset(gpuMalloc->gpu_g_haty, 0, recip->n_L * sizeof(CUFFTCOMPLEX));
    cudaMemset(gpuMalloc->gpu_g_hatz, 0, recip->n_L * sizeof(CUFFTCOMPLEX)); 

    cudaMemset(gpuMalloc->gpu_uk, 0, recip->blocks_N_L * sizeof(Real));
}


void ENUFForce::cuenuf_finalize(recip_plan* recip, gpu_malloc* gpuMalloc)
{
    free(recip->c_phi_inv0);
    free(recip->c_phi_inv1);
    free(recip->c_phi_inv2);
    free(recip->c_coeff);
    free(recip->uk);

    cudaFree(gpuMalloc->gpu_uk);
    cudaFree(gpuMalloc->gpu_g);

    cudaFree(gpuMalloc->gpu_c_phi_inv0);
    cudaFree(gpuMalloc->gpu_c_phi_inv1);
    cudaFree(gpuMalloc->gpu_c_phi_inv2);

    cudaFree(gpuMalloc->gpu_g_hatx);
    cudaFree(gpuMalloc->gpu_g_haty);
    cudaFree(gpuMalloc->gpu_g_hatz);
}


void ENUFForce::computeForce(unsigned int timestep)
{

    if (!m_params_set)
    {
        cerr << endl << "setParams must be called prior to computeForces()" << endl;
        throw std::runtime_error("Error computing forces in ENUFForce");
    }

    unsigned int group_size = m_group->getNumMembers();
    if (group_size == 0)
        return;

    Real4* d_pos = m_basic_info->getPos()->getArray(location::device, access::read);
    Real* d_charge = m_basic_info->getCharge()->getArray(location::device, access::read);
    const BoxSize& box = m_basic_info->getGlobalBox();
    Real3 L = box.getL();
    Real volume = L.x * L.y * L.z;

    //begin ArrayHandle scope
    Real4* d_force = m_basic_info->getForce()->getArray(location::device, access::readwrite);
    unsigned int* d_group_members = m_group->getIdxGPUArray();
    
    if(m_first_run)
    {
        cuenuf_init(group_size, m_alpha, m_sigma, m_precision,
                    m_Nx, m_Ny, m_Nz, &m_recip, &m_gpuMalloc,
                    m_block_size);
        CHECK_CUDA_ERROR();
        cufftPlan3d(&plan, m_recip.n.x, m_recip.n.y, m_recip.n.z,
                    CUFFT_TRANSFORM_TYPE);
        m_first_run = false;
    }

    if(m_box_changed)
    {
        Real mxinv= Real(m_recip.n.x) / L.x;
        Real myinv= Real(m_recip.n.y) / L.y;
        Real mzinv= Real(m_recip.n.z) / L.z;
        m_recip.minv = ToReal3(mxinv, myinv, mzinv);
        m_box_changed = false;
    }

    // run all kernels in parallel on GPUs
    cuenuf(d_force, d_pos, d_charge, box, d_group_members,
            group_size, plan, &m_recip, &m_gpuMalloc);
    
    CHECK_CUDA_ERROR();

    // If there are exclusions, correct the long-range part of Electrostatic potential
    if( m_nlist->getExclusionsSet())
    {
        unsigned int* d_n_ex_idx = m_nlist->getNExIdx()->getArray(location::device, access::read);
        unsigned int* d_ex_list_idx= m_nlist->getExListIdx()->getArray(location::device, access::read);
        
        ForceLog force_log;
        force_log.virial = m_all_info->getLogFlags()[log_flag::virial];
        force_log.potential = m_all_info->getLogFlags()[log_flag::potential];
        force_log.virial_matrix = m_all_info->getLogFlags()[log_flag::virial_matrix]||m_all_info->getLogFlags()[log_flag::press_tensor];
        force_log.d_virial = m_basic_info->getVirial()->getArray(location::device, access::readwrite);
        force_log.d_virial_matrix = m_basic_info->getVirialMatrix()->getArray(location::device, access::readwrite);
        
        fix_exclusions2(d_force, force_log, d_pos, d_charge, box,
                d_n_ex_idx, d_ex_list_idx, m_nlist->getExListIndexer(),
                m_alpha, m_group->getIdxGPUArray(), group_size, m_block_size);
            
        CHECK_CUDA_ERROR();
    }

    bool cal_virial = m_all_info->getLogFlags()[log_flag::virial];
    bool cal_virial_matrix = m_all_info->getLogFlags()[log_flag::virial_matrix]||m_all_info->getLogFlags()[log_flag::press_tensor];
    bool cal_potential = m_all_info->getLogFlags()[log_flag::potential];
    
    if(cal_virial||cal_virial_matrix||cal_potential)
    {
        cudaMemcpy(m_recip.uk, m_gpuMalloc.gpu_uk,
                m_recip.blocks_N_L * sizeof(Real),
                cudaMemcpyDeviceToHost);
        
        m_recip.Uk = 0;
        for(int i = 0; i < m_recip.blocks_N_L; i++)
            m_recip.Uk += m_recip.uk[i];

        Real self_energy = m_q2*m_alpha/sqrt(M_PI);
        m_recip.Uk = m_recip.Uk*volume/(2*M_PI);

        // correction for self interactions of charged particles
        {
            Real4* h_force = m_basic_info->getForce()->getArray(location::host, access::readwrite);
            Real* h_virial = m_basic_info->getVirial()->getArray(location::host, access::readwrite);
            h_force[0].w += m_recip.Uk - self_energy;
            h_virial[0] += 0.0;
        }
    }
}


void export_ENUFForce()
{
    class_<ENUFForce, boost::shared_ptr<ENUFForce>, bases<Force>, boost::noncopyable >
        ("ENUFForce", init<boost::shared_ptr<AllInfo>, boost::shared_ptr<NeighborList>, boost::shared_ptr<ParticleSet> >()).def("setParams", &ENUFForce::setParams);
}

