#!/bin/bash

ana_path=pic_analysis
ana_config=config_files/analysis_config.dat
conf=config_files/conf.dat
ch_inductive () {
    sed -i -e "s/\(inductive = \).*/\1$1/" $ana_config
}

ch_htt () {
    sed -i -e "s/\(httx = \).*,\(.*\)/\1$1,\2/" $conf
    sed -i -e "s/\(htty = \).*,\(.*\)/\1$2,\2/" $conf
    sed -i -e "s/\(httz = \).*,\(.*\)/\1$3,\2/" $conf
}

run_dissipation () {
    cd $1
    ch_inductive 1
    ch_htt 64 1 2
    git pull
    git checkout pfields
    cd build
    make
    make install
    cd ..
    mpirun -np 128 ./dissipation
}

fpath=/net/scratch2/xiaocanli/mime25-sigma01-beta02-200-100/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/mime25-sigma033-beta006-200-100/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/guofan/sigma1-mime25-beta001/$ana_path
run_dissipation $fpath
cp data/* /net/scratch2/xiaocan/sigma1-mime25-beta001/$ana_path/data
fpath=/net/scratch2/xiaocanli/mime25-guide0-beta0007-200-100/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/sigma1-mime100-beta001-mustang/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/mime25-guide0-beta001-200-100/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/mime25-guide0-beta001-200-100-sigma033/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/mime25-sigma1-beta002-200-100-noperturb/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/xiaocanli/mime25-sigma1-beta002-guide1-200-100/$ana_path
run_dissipation $fpath
fpath=/net/scratch2/guofan/sigma1-mime25-beta001-track-3/$ana_path
run_dissipation $fpath
