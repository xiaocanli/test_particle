#!/bin/bash

ana_path=/global/cscratch1/sd/xiaocan/test_particle
pic_run=qile_3d
run_path=$ana_path/$pic_run/
ana_config=$run_path/init_pic_hdf5.dat

ch_vel () {
    sed -i -e "s/\(Particle thermal speed (in light speed): \).*/\1$1/" $ana_config
}

ch_particle_domain () {
    sed -i -e "s/\(xmin_ptl = \).*/\1$1/" $ana_config
    sed -i -e "s/\(xmax_ptl = \).*/\1$2/" $ana_config
    # sed -i -e "s/\(ymin_ptl: \).*/\1$3/" $ana_config
    # sed -i -e "s/\(ymax_ptl: \).*/\1$4/" $ana_config
    # sed -i -e "s/\(zmin_ptl: \).*/\1$5/" $ana_config
    # sed -i -e "s/\(zmax_ptl: \).*/\1$6/" $ana_config
}

ch_tframe () {
    sed -i -e "s/\(Initial time frame: \).*/\1$1/" $ana_config
}

run_test_particle () {
    cd $run_path
    ch_vel $1
    ch_tframe $2
    ch_particle_domain $6 $7
    mkdir -p data
    srun -n $3 -N $4 -c $5 --cpu_bind=cores ./test_particle
    rm -rf data_${2}_${1}c_x$6_x$7
    mv data data_${2}_${1}c_x$6_x$7
    ch_tframe 22170
    ch_vel 0.3
}

run_code_velocity_scaling () {
    for tframe in {22170..44340..2217}
    do
        echo "Time frame: $tframe"
        for ptl_vel in 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95
        do
            run_test_particle $ptl_vel $tframe $1 $2 $3
        done
    done
}

mpi_size=1

# ncores=64
# nnodes=32
# nthreads=18
# run_code_velocity_scaling $ncores $nnodes $nthreads

ncores=1
nnodes=1
nthreads=1

tframe=22170
ptl_vel=0.7 # in light speed
# for xcenter in 100 220 290 440 510 645
for xcenter in 220 290 440 510 645
# for xcenter in 100
do
    xl=$(($xcenter-5))
    xr=$(($xcenter+5))
    run_test_particle $ptl_vel $tframe $ncores $nnodes $nthreads $xl $xr
done
