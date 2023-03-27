/*
GALAMOST - GPU-Accelerated Large-Scale Molecular Simulation Toolkit
COPYRIGHT
	GALAMOST Copyright (c) (2013) The group of Prof. Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of GALAMOST do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with GALAMOST are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	State Key Laboratory of Theoretical and Computational Chemistry,
	Institute of Theoretical Chemistry, Jilin University, China; 
	Prof. Zhong-Yuan Lu; 
	Email: luzhy@jlu.edu.cn
*/
//	Maintainer: Yong-Lei Wang, You-Liang Zhu

#include "EwaldForce.cuh"
Real4_tex_t pos_tex;
Real_tex_t charge_tex;

__global__ void gpu_compute_ewald_forces_kernel(Real4* d_force,
											ForceLog force_log, 
											Real4* d_pos,
											Real *d_charge,
											BoxSize box, 
											const unsigned int *d_n_neigh,
											const unsigned int *d_nlist, 
											Index2D nli, 
											Real *d_params, 
											int coeff_width,
											Real rcutsq,
											unsigned int *d_group_members,
											unsigned int group_size)
{
	extern __shared__ Real s_params[];
	for (unsigned int cur_offset = 0; cur_offset < coeff_width*coeff_width; cur_offset += blockDim.x)
	{
		if (cur_offset + threadIdx.x < coeff_width*coeff_width)
			s_params[cur_offset + threadIdx.x] = d_params[cur_offset + threadIdx.x];
	}
	__syncthreads();

    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size)
		return;	
    unsigned int idx = d_group_members[group_idx];	

	unsigned int n_neigh = d_n_neigh[idx];
	Real4 pos = texFetchReal4(d_pos, pos_tex, idx);
	Real qi = texFetchReal(d_charge, charge_tex, idx);

	Real4 force = d_force[idx];
	Real virial = Real(0.0);
	Real6 virial_matrix = ToReal6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	if(force_log.virial)
		virial = force_log.d_virial[idx];
	if(force_log.virial_matrix)
		virial_matrix = force_log.d_virial_matrix[idx];

    unsigned int cur_neigh = 0;
    unsigned int next_neigh = d_nlist[nli(idx, 0)];
	
	for (int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++)
	{
		cur_neigh = next_neigh;
		next_neigh = d_nlist[nli(idx, neigh_idx+1)];
		
		Real4 neigh_pos = texFetchReal4(d_pos, pos_tex, cur_neigh);
		Real qj = texFetchReal(d_charge, charge_tex, cur_neigh);		
		
		Real qiqj = qi*qj;

		Real dx = pos.x - neigh_pos.x;
		Real dy = pos.y - neigh_pos.y;
		Real dz = pos.z - neigh_pos.z;

		box.minDisImage(dx, dy, dz);
		Real rsq = dx*dx + dy*dy + dz*dz;

		if (rsq < rcutsq && qiqj != Real(0.0))
		{
			int typ_pair = __real_as_int(neigh_pos.w) * coeff_width + __real_as_int(pos.w);
			Real kappa = s_params[typ_pair];
			Real rinv = rsqrt_gala(rsq);
			Real r = Real(1.0) / rinv;
			Real r2inv = Real(1.0) / rsq;
			
			Real erfc_by_r_val = erfc_gala(kappa * r) * rinv;
			
			Real force_divr = qiqj * r2inv * (erfc_by_r_val + Real(2.0)*kappa*rsqrt_gala(Real(M_PI)) * exp_gala(-kappa*kappa* rsq));
			Real pair_eng = qiqj * erfc_by_r_val ;

			if(force_log.virial)
				virial += Real(1.0)/Real(6.0) * rsq * force_divr; 
			if(force_log.virial_matrix)
			{
				Real force_div2r = Real(0.5) * force_divr;
				virial_matrix.x +=  dx * dx * force_div2r;   // press_tensor_xx
				virial_matrix.y +=  dx * dy * force_div2r;   // press_tensor_xy
				virial_matrix.z +=  dx * dz * force_div2r;   // press_tensor_xz
				virial_matrix.w +=  dy * dy * force_div2r;   // press_tensor_yy
				virial_matrix.m +=  dy * dz * force_div2r;   // press_tensor_yz
				virial_matrix.n +=  dz * dz * force_div2r;   // press_tensor_zz			
			}			

			force.x += dx * force_divr;
			force.y += dy * force_divr;
			force.z += dz * force_divr;
		// energy is double counted: multiply by 0.5		
			force.w += pair_eng*Real(0.5);		
		}
	}

    d_force[idx] = force;
	if(force_log.virial)	
		force_log.d_virial[idx] = virial;
	if(force_log.virial_matrix)
		force_log.d_virial_matrix[idx] = virial_matrix;	
}


cudaError_t gpu_compute_ewald_forces(Real4* d_force, 
									ForceLog& force_log,
									Real4* d_pos,
									Real *d_charge,
									const BoxSize& box, 
									const unsigned int *d_n_neigh,
									const unsigned int *d_nlist,
									const Index2D& nli,
									Real *d_params, 
									int coeff_width,
									Real rcutsq,
									unsigned int *d_group_members,
									unsigned int group_size,
									unsigned int Ntot,
									unsigned int blocksize,
									unsigned int compute_capability)
{
    dim3 grid( (int)ceil((Real)group_size / (Real)blocksize), 1, 1);
    dim3 threads(blocksize, 1, 1);

	if (compute_capability < 350)
	{
		pos_tex.normalized = false;
		pos_tex.filterMode = cudaFilterModePoint;
		cudaError_t error = cudaBindTexture(0, pos_tex, d_pos, sizeof(Real4) * Ntot);
		if (error != cudaSuccess)
			return error;

		charge_tex.normalized = false;
		charge_tex.filterMode = cudaFilterModePoint;
		error = cudaBindTexture(0, charge_tex, d_charge, sizeof(Real) * Ntot);
		if (error != cudaSuccess)
			return error;
	}
	
	gpu_compute_ewald_forces_kernel<<< grid, threads, sizeof(Real)*coeff_width*coeff_width >>>(d_force,
												force_log,
												d_pos,
												d_charge,
												box, 
												d_n_neigh,
												d_nlist, 
												nli,		
												d_params, 
												coeff_width,
												rcutsq,
												d_group_members,
												group_size);
																							

    return cudaSuccess;
}
