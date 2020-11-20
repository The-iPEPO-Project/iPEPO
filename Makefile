HOSTNAME := $(shell hostname)

LAPACKLIB=-llapack -lblas
ARPACKLIB=-larpack

OPTFLAGS=-ffree-line-length-none -fbounds-check 

#-ffpe-trap=invalid,zero,overflow 
#-mcmodel=large
#DEBUG=-O3 -fopenmp

DEBUG=-O0 -ggdb -fbounds-check


# Which version of Arnoldi SVD to use (Arpack vs NAG)
ARNOLDI=lib_arnoldi_arpack_wrappers.o
NAGDIR=
NAGINC=
NAGLIB=

include Makefile.$(HOSTNAME)

LDFLAGS= $(DEBUG) $(LAPACKDIR) $(ARPACKDIR) $(ARCHLDFLAGS) $(NAGDIR)
FFLAGS= $(DEBUG) $(WARN) $(ARCHFFLAGS) $(OPTFLAGS) $(NAGINC)
LOADLIBES=$(LAPACKLIB) $(ARPACKLIB) $(NAGLIB)



### NB. to produce the makefile backslash: \ + Enter
TESTS = linsol_test  cl_ising_CTMRG  #mult_cl_ising_CTMRG

PROGRAMS_PEPS      = ising_ipeps        xyz_ipeps
PROGRAMS_IMAG_PEPS = ising_imag_ipeps   ising_imag_ipeps_FU  bhm_imag_ipeps  bhm_imag_ipeps_FU

PROGRAMS_MPS       = ising_imps xyz_imps kbhm_imps
PROGRAMS_IMAG_MPS  = ising_imag_imps bhm_imag_imps

PROGRAMS_TRAMTE    = ising_tramte_evecs ising_tramte_obs kbhm_tramte_evecs kbhm_tramte_obs


all: $(PROGRAMS_MPS) $(PROGRAMS_IMAG_MPS) $(PROGRAMS_TRAMTE) $(PROGRAMS_PEPS) $(PROGRAMS_IMAG_PEPS) $(TESTS) 
     


###### Utility ##########
UTIL = basic_functions.o \
       error_handling.o \
       array_utility.o \
       lib_lapack_wrappers.o \
       lib_arnoldi_utility.o \
       ${ARNOLDI} \
       TENS_helper.o \
       TENS_reshape.o \
       TENS_transpose.o \
       TENS_mult.o \
       TENS_mult_extension.o \
       matrix_functions.o \
       SVD_routines.o \
       utility.o \
       datafile_utility.o


###### Tensor Network operations ##########
MPS =  definitions_mps_mpo.o \
       definitions_peps_pepo.o \
       mps_peps_INOUT_files.o \
       iteration_helper.o \
       mps_mpo_utility.o \
       boundary_tensors.o \
       mps_mpo_algebra_finite.o \
       project_impo.o \
       mps_mpo_algebra_inf.o


PEPS = ${MPS} \
       project_ipeps.o \
       project_ipepo.o \
       peps_pepo_algebra.o


TRAMTE = power_method.o \
         tramte_network.o



###### TEBD ##########
TEBD_COMMON = simulation_parameters.o \
              propagator_utility.o


TEBD_PEPS = ${TEBD_COMMON} \
            observables_ipeps.o \
            tebd_callback_routines_ipeps.o \
            time_evolution_ipeps.o


TEBD_MPS = ${TEBD_COMMON} \
           observables_imps.o \
           tebd_callback_routines_imps.o \
           time_evolution_imps.o


TEBD_IMAG_MPS = ${TEBD_COMMON} \
                imaginary_observables_imps.o \
                imaginary_tebd_callback_routines_imps.o \
                imaginary_time_evolution_imps.o


TEBD_IMAG_PEPS = ${TEBD_COMMON} \
                 imaginary_observables_ipeps.o \
                 imaginary_tebd_callback_routines_ipeps.o \
                 imaginary_full_update_ipeps.o \
                 imaginary_time_evolution_ipeps.o


##### CTMRG #######
CTMRG = ctm_definitions.o \
        corner_transfer_matrix.o \
        TN_contractions.o 




############################################# Run TraMTE #########################################################

ising_tramte_evecs= $(UTIL) $(MPS) $(TEBD_MPS) $(TRAMTE) pauli_basis.o ising_propagator_pauli_basis.o run_ising_tramte_evecs.o
ising_tramte_evecs: $(ising_tramte_evecs)
	$(LD) $(LDFLAGS) $(ising_tramte_evecs) -o $@ $(LOADLIBES)


ising_tramte_obs= $(UTIL) $(MPS) $(TEBD_MPS) $(TRAMTE) pauli_basis.o ising_propagator_pauli_basis.o run_ising_tramte_obs.o
ising_tramte_obs: $(ising_tramte_obs)
	$(LD) $(LDFLAGS) $(ising_tramte_obs) -o $@ $(LOADLIBES)


kbhm_tramte_evecs= $(UTIL) $(MPS) $(TEBD_MPS) $(TRAMTE) gell_mann_basis.o sun_basis.o kbhm_propagator.o run_kbhm_tramte_evecs.o
kbhm_tramte_evecs: $(kbhm_tramte_evecs)
	$(LD) $(LDFLAGS) $(kbhm_tramte_evecs) -o $@ $(LOADLIBES)


kbhm_tramte_obs= $(UTIL) $(MPS) $(TEBD_MPS) $(TRAMTE) gell_mann_basis.o sun_basis.o kbhm_propagator.o run_kbhm_tramte_obs.o
kbhm_tramte_obs: $(kbhm_tramte_obs)
	$(LD) $(LDFLAGS) $(kbhm_tramte_obs) -o $@ $(LOADLIBES)


##################################################################################################################







########################################## Run iPEPS TEBD ###########################################################


ising_ipeps= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_PEPS) psi_pauli_basis.o ising_propagator_v2.o run_ising_ipeps.o
ising_ipeps: $(ising_ipeps)
	$(LD) $(LDFLAGS) $(ising_ipeps) -o $@ $(LOADLIBES)


xyz_ipeps= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_PEPS) psi_pauli_basis.o xyz_propagator.o run_xyz_ipeps.o
xyz_ipeps: $(xyz_ipeps)
	$(LD) $(LDFLAGS) $(xyz_ipeps) -o $@ $(LOADLIBES)


######################################################################################################################






################################################ Run iMPS TEBD #######################################################

ising_imps= $(UTIL) $(MPS) $(TEBD_MPS) pauli_basis.o ising_propagator_pauli_basis.o run_ising_imps.o
ising_imps: $(ising_imps)
	$(LD) $(LDFLAGS) $(ising_imps) -o $@ $(LOADLIBES)


xyz_imps= $(UTIL) $(MPS) $(TEBD_MPS) pauli_basis.o xyz_propagator_pauli_basis.o run_xyz_imps.o
xyz_imps: $(xyz_imps)
	$(LD) $(LDFLAGS) $(xyz_imps) -o $@ $(LOADLIBES)


kbhm_imps= $(UTIL) $(MPS) $(TEBD_MPS) gell_mann_basis.o sun_basis.o kbhm_propagator.o run_kbhm_imps.o
kbhm_imps: $(kbhm_imps)
	$(LD) $(LDFLAGS) $(kbhm_imps) -o $@ $(LOADLIBES)

######################################################################################################################






############################################# Run IMAGINARY iMPS TEBD ################################################

ising_imag_imps= $(UTIL) $(MPS) $(TEBD_IMAG_MPS) psi_pauli_basis.o ising_imaginary_propagator.o run_imaginary_ising_imps.o
ising_imag_imps: $(ising_imag_imps)
	$(LD) $(LDFLAGS) $(ising_imag_imps) -o $@ $(LOADLIBES)


bhm_imag_imps= $(UTIL) $(MPS) $(TEBD_IMAG_MPS) psi_sun_basis.o bhm_imaginary_propagator.o run_imaginary_bhm_imps.o
bhm_imag_imps: $(bhm_imag_imps)
	$(LD) $(LDFLAGS) $(bhm_imag_imps) -o $@ $(LOADLIBES)

#######################################################################################################################







###################################### Run IMAGINARY iPEPS TEBD ########################################################


bhm_imag_ipeps= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_IMAG_PEPS) psi_sun_basis.o bhm_imaginary_propagator.o run_imaginary_bhm_ipeps.o
bhm_imag_ipeps: $(bhm_imag_ipeps)
	$(LD) $(LDFLAGS) $(bhm_imag_ipeps) -o $@ $(LOADLIBES)


bhm_imag_ipeps_FU= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_IMAG_PEPS) psi_sun_basis.o bhm_imaginary_propagator.o run_imaginary_bhm_ipeps_FU.o
bhm_imag_ipeps_FU: $(bhm_imag_ipeps_FU)
	$(LD) $(LDFLAGS) $(bhm_imag_ipeps_FU) -o $@ $(LOADLIBES)


ising_imag_ipeps= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_IMAG_PEPS) psi_pauli_basis.o ising_imaginary_propagator.o run_imaginary_ising_ipeps.o
ising_imag_ipeps: $(ising_imag_ipeps)
	$(LD) $(LDFLAGS) $(ising_imag_ipeps) -o $@ $(LOADLIBES)


ising_imag_ipeps_FU= $(UTIL) $(PEPS) ${CTMRG} $(TEBD_IMAG_PEPS) psi_pauli_basis.o ising_imaginary_propagator.o run_imaginary_ising_ipeps_FU.o
ising_imag_ipeps_FU: $(ising_imag_ipeps_FU)
	$(LD) $(LDFLAGS) $(ising_imag_ipeps_FU) -o $@ $(LOADLIBES)

##########################################################################################################################





################################################## Testing ########################################################


cl_ising_CTMRG= $(UTIL) $(PEPS) ${CTMRG} classical_ising_2D.o test_classical_ising_CTMRG.o
cl_ising_CTMRG: $(cl_ising_CTMRG)
	$(LD) $(LDFLAGS) $(cl_ising_CTMRG) -o $@ $(LOADLIBES)


mult_cl_ising_CTMRG= $(UTIL) $(PEPS) ${CTMRG} classical_ising_2D.o test_mult_classical_ising_CTMRG.o
mult_cl_ising_CTMRG: $(mult_cl_ising_CTMRG)
	$(LD) $(LDFLAGS) $(mult_cl_ising_CTMRG) -o $@ $(LOADLIBES)


linsol_test= $(UTIL) $(PEPS) ${CTMRG} ${IMAG_FU_PEPS} test_linsol.o
linsol_test: $(linsol_test)
	$(LD) $(LDFLAGS) $(linsol_test) -o $@ $(LOADLIBES)

###################################################################################################################




print: 
	$(a2ps) $(patsubst %.o, %.f90, $(cca))

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.mod: %.o

clean:
	rm -f *.o *.mod $(PROGRAMS_TRAMTE) $(PROGRAMS_PEPS) $(PROGRAMS_IMAG_PEPS) $(PROGRAMS_MPS) $(PROGRAMS_IMAG_MPS) $(TESTS)

realclean: clean
	rm -f *~








