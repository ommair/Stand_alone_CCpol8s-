# Stand_alone_CCpol8s(2013) PES

Given code is written for the paper 

How well can polarization models of pairwise nonadditive forces describe liquid water?
Omololu Akin-Ojo and Krzysztof Szalewicz
Citation: J. Chem. Phys. 138, 024316 (2013); doi: 10.1063/1.4773821

Current version of the code is written by Ommair Ishaque, August 8, 2023 

It reads geometry file dimers_geo_xyz.dat which contains 3216 dimer configurations.
Position of all sites are required as input, see dimers_geo_xyz.dat

There are four options given for NB-induction model to run this code. You need to specify in the CONTROL file which option is required.

Below is the example CONTROL file for "damped_ind" with "physical_dipole" model.
 
ccpol8s_omo         # potential name (ccpol8s_omo)
user_config         # inintial configuration (fcc_config or user_config)
damped_ind          # damped_ind or undamped_ind
physical_dipole     # physical_dipole or mathematical_dipole

if you want to use "damped_ind" with mathemaitcal dipole then change the last line to "mathematical_dipole"

similarly for other possibilities.

In the OUTPUT file it prints vdw, Coulomb, NB-induction and CCpol8s (sum of latter three) in kcal/mol.

Finally use the make file to compile the code

make

example input with three dimers configurations

./execute < ccpol.inp  

 runninng water cluster with CCdpol8s_omo    calculation using damped_ind           physical_dipole

         dimer     VDW(kcal/mol)            COUL(kcal/mol)            IND(kcal/mol)      CCpol8s'+NB(ind)(kcal/mol)
           1   10.242648309741215       -7.5212509674435735       -1.3002060793115247        1.4211912629861168
           2   4.5076171437514336       -7.2480607737956166      -0.95243087077567701       -3.6928745008198600
           3   1.3453278177704959       -5.5222180677606874      -0.57889239123006853       -4.7557826412202600

or see OUTPUT file for fir energies

to clean use 

make clean

