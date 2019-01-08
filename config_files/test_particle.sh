#!/bin/bash

ana_path=/net/scratch3/xiaocanli/diffusion
pic_run=3D-Lx150-bg0.2-150ppc-2048KNL
run_path=$ana_path/$pic_run/
ana_config=$run_path/init.dat

ch_vel () {
    sed -i -e "s/\(Particle thermal speed (in light speed): \).*/\1$1/" $ana_config
}


ch_tframe () {
    sed -i -e "s/\(Initial time frame: \).*/\1$1/" $ana_config
}

mpi_size=64
fd_tinterval=1

run_test_particle () {
    cd $run_path
    ch_vel $1
    ch_tframe $2
    mkdir -p data
    srun -n $3 -N $4 -c $5 --cpu_bind=cores ./test_particle
    rm -rf data_${2}_${1}c
    mv data data_${2}_${1}c 
    ch_tframe 22170
    ch_vel 0.3
}

run_code () {
    for tframe in {22170..44340..2217}
    do
        echo "Time frame: $tframe" 
        for ptl_vel in 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
        do
            run_test_particle $ptl_vel $tframe $1 $2 $3
        done
    done
}

# ncores=64
# nnodes=32
# nthreads=18

# run_code $ncores $nnodes $nthreads

ncores=2
nnodes=1
nthreads=18

tframe=22170
ptl_vel=0.7 # in light speed
run_test_particle $ptl_vel $tframe $ncores $nnodes $nthreads
# ptl_vel=0.8 # in light speed
# run_test_particle $ptl_vel $tframe $ncores $nnodes $nthreads
# ptl_vel=0.9 # in light speed
# run_test_particle $ptl_vel $tframe $ncores $nnodes $nthreads
# ptl_vel=0.95 # in light speed
# run_test_particle $ptl_vel $tframe $ncores $nnodes $nthreads
