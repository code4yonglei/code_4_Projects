#include "Force.h"
#include "ENUFForce.cuh"
#include "ParticleSet.h"

#ifndef __ENUFFORCE_H__
#define __ENUFFORCE_H__

class ENUFForce : public Force
{
    public:
        ENUFForce(boost::shared_ptr<AllInfo> all_info,
                boost::shared_ptr<NeighborList> nlist,
                boost::shared_ptr<ParticleSet> group);
        virtual ~ENUFForce();

        void setParams(Real alpha, Real sigma, int precision, 
                int Nx, int Ny, int Nz);

        void slotBoxChanged()
		{
            m_box_changed = true;
        }

    void cuenuf_init(unsigned int group_size, Real alpha, Real sigma,
            int precision, int Nx, int Ny, int Nz, recip_plan* recip,
            gpu_malloc* gpuMalloc, unsigned int block_size);
    void cuenuf_finalize(recip_plan* recip, gpu_malloc* gpuMalloc);
    protected:
        int m_Nx;             // Number of grid points in x direction
        int m_Ny;             // Number of grid points in y direction
        int m_Nz;             // Number of grid points in z direction
        int m_precision;      // Interpolation order
        Real m_alpha;         // Screening parameter for erfc(kappa*r)
        Real m_sigma;         // Real space cutoff
        Real m_q;             // Total system charge
        Real m_q2;            // Sum(q_i*q_i), where q_i is the charge of each particle
        bool m_box_changed;   // Set to true when box size has changed

        // Sum grid points for the calculation of pressure and energy (input)
        boost::shared_ptr<Array<Real> > m_energy_sum;

        // Sum grid points for virial
        boost::shared_ptr<Array<Real> > m_virial_sum;

        // Connection to ParticleData box size change signal
        boost::signals2::connection m_boxchange_connection;

        // Neighborlist for computation
        boost::shared_ptr<NeighborList> m_nlist;

        // Group to compute properties
        boost::shared_ptr<ParticleSet> m_group;

        int m_block_size;     // Block size to run calculation
        gpu_malloc m_gpuMalloc;
        recip_plan m_recip;
        cufftHandle plan;     // FFT calculation using GPU
        bool m_first_run;     // True if this is the first run
        bool m_params_set;    // Set to true when params are set

        virtual void computeForce(unsigned int timestep); // Compute forces
};

void export_ENUFForce();
#endif
