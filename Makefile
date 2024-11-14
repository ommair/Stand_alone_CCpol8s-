# This is Makefile for creating executables

# This variable is for fortran compiler
FC = gfortran

# This variable will give compiler flags
FFLAGS = -O3 
FCFLAGS = -c 
FLAGS = -fbounds-check 

OBJS = common_arrays.o \
       utility.o       \
       molecular_sites.o \
       coulomb.o          \
       vdw_ccpol8s_omo.o    \
       nb_induction_model.o  \
       initialization.o    
    
       
          
EXECUTE = execute

$(EXECUTE): main.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $(EXECUTE) main.f90 $(OBJS)
 
# These lines produces the .mod and .o files
common_arrays.o: common_arrays.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) common_arrays.f90
utility.o: utility.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) utility.f90
molecular_sites.o: molecular_sites.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) molecular_sites.f90
initialization.o: initialization.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) initialization.f90
vdw_ccpol8s_omo.o: vdw_ccpol8s_omo.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol8s_omo.f90
coulomb.o: coulomb.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) coulomb.f90
nb_induction_model.o: nb_induction_model.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) nb_induction_model.f90
# This line is for clean up
clean: 
	rm -rf *.mod $(OBJS) *~ $(EXECUTE) OUTPUT* *.out *.xyz *.err *.dat
	
