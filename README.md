# code_4_ENUF_project
Computing electrostatic interactions in molecular simulations

ENUF, an abbreviation for the Ewald summation method based on Non-Uniform fast Fourier transform (NFFT) technique, is a method proposed to calculate electrostatic interactions between charged particles in molecular simulation systems.

The ENUF method was first implemented in the classic molecular dynamics simulation code, and thereafter in the dissipative particle dynamics (DPD) framework to calculate electrostatic interactions between charge density distributions at mesoscopic level. Selecting a set of optimal physical parameters, the ENUF method gives a good precision as desired and bears a computational complexity of O(NlogN), which is comparable to that of PME and PPPM methods.

The ENUF and ENUF-DPD methods were adopted to explore the dependence of conformational properties of polyelectrolytes on charge fraction, ion concentration and counterion valency of added salts, to investigate specific binding structures of dendrimers on bilayer membranes and the corresponding permeation mechanisms, and to study heterogeneous structures and dynamics in ionic liquids and how electrostatic interactions between charged particles affect these properties at extended spatiotemporal scales.
