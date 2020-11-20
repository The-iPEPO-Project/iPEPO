#!/bin/bash

# Work out where we are, and where the code is
BASE=`pwd`
BIN=${BASE}/../ising_imag_ipeps

# Load variables for grid run
. grid-common.sh

# Read in parameters
. params-imaginary-ising-ipeps.sh

# Reload state from file
reload="N"  

setup_control() {
    RESDIR=${BASE}/imaginary_ising_g=${g}_eps=${eps}_chi=${chi}_eta=${eta}_epsENV=${epsENV}_chiENV=${chiENV}_etaENV=${etaENV}
    mkdir -p ${RESDIR}

    # Name of control file (within working dir)
    CONTROL=$RESDIR/CONTROL

    cat <<EOF > ${CONTROL}
$g, $delta
$N_steps_1, $N_steps_2, $N_steps_3, $N_steps_4, $N_steps_5
$dt_1, $dt_2, $dt_3, $dt_4, $dt_5
$test_interval_1, $test_interval_2, $test_interval_3, $test_interval_4, $test_interval_5
$min_sep, $max_sep
$eta, $chi, $eps
$etaENV, $chiENV, $epsENV
$xname, $x0, $dx, $nx
$reload
$s1, $s2
EOF
}

#"$reload_file"
#"$ctm_file"

#$s1, $s2


run_grid() {
    setup_control

    GridScript=${RESDIR}/GRID_SCRIPT

    job_name="ipeps-ising-imag"
    nproc=1

    write_grid_script
    
    qsub $GridScript
}


run_local() {
    setup_control

    cd ${RESDIR}

    
    export LD_LIBRARY_PATH=/opt/NAG/fll6a24dfl/lib/
    export GFORTRAN_UNBUFFERED_ALL="Y"
    $BIN < $CONTROL | tee log.txt

    cd ../
}


## Set params
g=2.0

## Setup parameter loop
xname="g"
x0=${g} 
dx=0.1
nx=1

## Run (start from scratch, without reloading)
run_local
     

