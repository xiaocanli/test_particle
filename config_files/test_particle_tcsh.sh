#!/bin/tcsh

setenv wlcs_conf init8.dat
setenv init init.dat
ch_nptl() {
    sed -i -e "s/\(Total number of particles: \).*/\1$1/" $init
}

# ch_mass () {
#     sed -i -e "s/\(Mass: \).*/\1$1/" $init
# }

# ch_charge () {
#     sed -i -e "s/\(Charge: \).*/\1$1/" $init
# }

# ch_tracking_time () {
#     sed -i -e "s/\(Total tracking time (s): \).*/\1$1/" $init
# }

# ch_nframes () {
#     sed -i -e "s/\(Number of diagnostic time frames: \).*/\1$1/" $init
# }

# ch_frequency () {
#     sed -i -e "s/$1/$2/" $wlcs_conf
# }

# ch_current () {
#     sed -i -e "s/$1/$2/" $wlcs_conf
# }

# ch_species () {
#     ch_nptl $1
#     ch_mass $2
#     ch_charge $3
# }

# ch_to_electron () {
#     ch_species 25000 0.000545 -1.0
# }

# run_test_particle () {
#     mpirun -np 2 -npernode 2 -bysocket -bind-to-socket ./test_particle
#     ch_current $1 $2
#     ch_tracking_time $3
#     cd data
#     cp diff_coeffs.dat diff_coeffs_$4.dat
#     cp espectrum.dat espectrum_$4.dat
#     cp espect_escape.dat espect_escape_$4.dat
#     cd ..
# }

# run_particle_sf1 () {
#     run_test_particle 0.25000000000000 1.25000000000000 $1 $7_025I0_0001Hz
#     run_test_particle 1.25000000000000 2.50000000000000 $2 $7_125I0_0001Hz
#     run_test_particle 2.50000000000000 25.0000000000000 $3 $7_25I0_0001Hz
#     run_test_particle 25.0000000000000 250.000000000000 $4 $7_250I0_0001Hz
#     run_test_particle 250.000000000000 2500.00000000000 $5 $7_2500I0_0001Hz
#     run_test_particle 2500.00000000000 0.25000000000000 $6 $7_25000I0_0001Hz
# }

# run_particle_sf2 () {
#     run_test_particle 0.25000000000000 1.25000000000000 $1 $7_025I0_001Hz
#     run_test_particle 1.25000000000000 2.50000000000000 $2 $7_125I0_001Hz
#     run_test_particle 2.50000000000000 25.0000000000000 $3 $7_25I0_001Hz
#     run_test_particle 25.0000000000000 250.000000000000 $4 $7_250I0_001Hz
#     run_test_particle 250.000000000000 2500.00000000000 $5 $7_2500I0_001Hz
#     run_test_particle 2500.00000000000 0.25000000000000 $6 $7_25000I0_001Hz
# }

# run_particle_sf3 () {
#     run_test_particle 0.25000000000000 1.25000000000000 $1 $7_025I0_01Hz
#     run_test_particle 1.25000000000000 2.50000000000000 $2 $7_125I0_01Hz
#     run_test_particle 2.50000000000000 25.0000000000000 $3 $7_25I0_01Hz
#     run_test_particle 25.0000000000000 250.000000000000 $4 $7_250I0_01Hz
#     run_test_particle 250.000000000000 2500.00000000000 $5 $7_2500I0_01Hz
#     run_test_particle 2500.00000000000 0.25000000000000 $6 $7_25000I0_01Hz
# }

# run_proton_mf () {
#     run_particle_sf1 0.8 0.4 0.04 0.004 0.0004 0.4 proton
#     ch_frequency 0.001000000 0.010000000
#     ch_nframes 500
#     run_particle_sf2 0.08 0.04 0.004 0.0004 0.00004 0.04 proton
#     ch_frequency 0.010000000 0.100000000
#     ch_nframes 200
#     run_particle_sf3 0.008 0.004 0.0004 0.00004 0.000004 4.0 proton
#     ch_frequency 0.100000000 0.001000000
#     ch_nframes 1000
# }

# run_electron_mf () {
#     run_particle_sf1 0.02 0.01 0.001 0.0001 0.00001 0.01 electron
#     ch_frequency 0.001000000 0.010000000
#     ch_nframes 500
#     run_particle_sf2 0.002 0.001 0.0001 0.00001 0.000001 0.001 electron
#     ch_frequency 0.010000000 0.100000000
#     ch_nframes 200
#     run_particle_sf3 0.0002 0.0001 0.00001 0.000001 0.0000001 0.1 electron
#     ch_frequency 0.100000000 0.001000000
#     ch_nframes 1000
# }

# ch_tracking_time 4.0
# ch_species 16 1.0 1.0
# run_proton_mf
# ch_species 16 0.000545 -1.0
# ch_tracking_time 0.1
# run_electron_mf
