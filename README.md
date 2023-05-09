# code_4_ENUF_project

## Computing electrostatic interactions in molecular simulations

**ENUF**, an abbreviation for the **E**wald summation method based on **N**on-**U**niform fast
**F**ourier transform (NFFT) technique, is a method proposed to calculate electrostatic interactions
between charged particles in molecular simulation systems.

The **ENUF** method was first implemented in the classic molecular dynamics simulation code, 
and thereafter in the dissipative particle dynamics (**DPD**) framework to calculate 
electrostatic interactions between charge density distributions at mesoscopic level. 
Selecting a set of optimal physical parameters, the ENUF method gives a good precision as desired 
and bears a computational complexity of *O*(*N*log*N*), which is comparable to that of 
PME (particle-mesh Ewald summation method) and PPPM (particle-particle particle-mesh 
Ewald summation method) methods.

The ENUF and ENUF-DPD methods were adopted to explore the dependence of conformational properties 
of polyelectrolytes on charge fraction, ion concentration and counterion valency of added salts, 
to investigate specific binding structures of dendrimers on bilayer membranes and the corresponding 
permeation mechanisms, and to study heterogeneous structures and dynamics in ionic liquids 
and how electrostatic interactions between charged particles affect these properties at extended 
spatiotemporal scales.

## 01_DPD_solvent_with_salt_FortranC
This folder contains the code (serial version) to simulate properties of ions in aqueous dilute solution.
The code is written using **Fortran** and **C**.
The [**FFTW**](https://www.fftw.org/) library should be installed before compilation of the source code.

## 02_DPD_solvent_with_salt_MPI_FortranC
This folder contains the paralled version of the code to simulate properties of ions in aqueous dilute solution. The code is written using Fortran and C and paralleled using **MPI**.

## 03_GALAMOST_CUDA_C++
This folder contains several files that have been implemented in the **GALAMOST** package.
The descriptions of these methods are available in Ref. [5].

## 04_hybrid_CUDA_MPI_2021CPCpaper
This folder contains files for the publication [7](https://www.sciencedirect.com/science/article/abs/pii/S0010465520302393).


## Relevant Publications
[1] Y.-L. Wang, A. Laaksonen*, Z.-Y. Lu*. [Implementation of non-uniform FFT based Ewald summation in dissipative particle dynamics method](https://www.sciencedirect.com/science/article/pii/S0021999112005542). J. Comput. Phys. 235, 666-682 (2013).

[2] Y.-L. Wang, F. Hedman, M. Porcu, F. Mocci, A. Laaksonen*. [Non-uniform FFT and its applications in particle simulations](https://www.scirp.org/journal/paperinformation.aspx?paperid=42807), Applied Mathematics 5, 520-541 (2014).

[3] S.-C. Yang, Y.-L. Wang, G.-S. Jiao, H.-J. Qian*, Z.-Y. Lu. [Accelerating electrostatic interaction calculations with graphical processing units based on new developments of Ewald method using non-uniform fast Fourier transform](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.24250), J. Comput.
Chem. 37, 378-387 (2016).

[4] S.-C. Yang*, Z.-Y. Lu, H.-J. Qian, Y.-L. Wang, J.-P. Han. [A hybrid Parallel architecture for electrostatic interactions in simulation of dissipative particle dynamics](https://www.sciencedirect.com/science/article/abs/pii/S0010465517302126). Comput. Phys. Commun. 220, 376-389 (2017).

[5] Y.-L. Wang∗, Y.-L. Zhu, Z.-Y. Lu, A. Laaksonen*. [Electrostatic interactions in soft particle systems: Mesoscale simulations of ionic liquids](https://pubs.rsc.org/en/content/articlehtml/2018/sm/c8sm00387d), Soft Matter 14, 4252-4267 (2018) (Highlighted on the inside back cover of Soft Matter).

[6] S.-C. Yang, B. Li, Y.-L. Zhu, A. Laaksonen, Y.-L. Wang*, [The ENUF method –- Ewald summation based on non-uniform fast Fourier transform: Implementation, parallelization, and application](https://onlinelibrary.wiley.com/doi/10.1002/jcc.26395), J. Comput. Chem. 41, 2316-2335 (2020).

[7] S.-C. Yang, Y.-L. Wang*. [A hybrid MPI-CUDA approach for nonequispaced discrete Fourier transformation](https://www.sciencedirect.com/science/article/abs/pii/S0010465520302393), Comput. Phys. Commun. 258, 107513 (2021).
