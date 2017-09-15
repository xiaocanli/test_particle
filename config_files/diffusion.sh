#!/bin/bash

conf=config_wlcs.dat
ch_deltay () {
    sed -i -e "s/$1/$2/" $conf
}

run_diffusion () {
    ch_deltay $1 $2
    mpirun -np 32 -npernode 2 -bysocket -bind-to-socket ./test_particle
    cp data/diff_coeffs.dat data/diff_coeffs_$3.dat
}

run_diffusion 0.5000000000 0.5000000000 05
run_diffusion 0.5000000000 0.6000000000 06
run_diffusion 0.6000000000 0.7000000000 07
run_diffusion 0.7000000000 0.8000000000 08
run_diffusion 0.8000000000 0.9000000000 09
ch_deltay 0.9000000000 0.5000000000
