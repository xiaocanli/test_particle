"""
Plotting procedures for reconnection
"""
import argparse
import collections
import gc
import math
import os
import os.path
import struct
import subprocess
import sys

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style as style
import multiprocessing
import numpy as np
from joblib import Parallel, delayed
from matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage.filters import median_filter, gaussian_filter

import palettable
from shell_functions import mkdir_p

style.use(['seaborn-white', 'seaborn-paper'])
# rc('font', **{'family': 'serif', 'serif': ["Times", "Palatino", "serif"]})
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc("font", family="Times New Roman")
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

colors_Set1_9 = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
colors_Dark2_8 = palettable.colorbrewer.qualitative.Dark2_8.mpl_colors
colors_Paired_12 = palettable.colorbrewer.qualitative.Paired_12.mpl_colors
colors_Tableau_10 = palettable.tableau.Tableau_10.mpl_colors
colors_GreenOrange_6 = palettable.tableau.GreenOrange_6.mpl_colors
colors_Bold_10 = palettable.cartocolors.qualitative.Bold_10.mpl_colors

font = {
    'family': 'serif',
    #'color'  : 'darkred',
    'color': 'black',
    'weight': 'normal',
    'size': 24,
    }

def plot_spectrum(run_dir):
    """ plotting particle energy spectrum.
    """
    fname = run_dir + 'data/espectrum.dat'
    with open(fname, 'r') as f:
        data = np.genfromtxt(f, delimiter='')
    dimx, dimy = data.shape

    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig = plt.figure(figsize=[7, 5])
    ax1 = fig.add_axes([xs, ys, w1, h1])
    for j in range(1, dimx, 20):
        ax1.loglog(data[0,:],data[j,:],linewidth=2)

    ax1.set_xlabel(r'$\gamma - 1$', fontdict=font, fontsize=20)
    ax1.set_ylabel(r'$f(\gamma - 1)$', fontdict=font, fontsize=20)
    ax1.tick_params(labelsize=16)
    ax1.set_xlim([1E-4, 1E1])
    ax1.grid(True)
    
    plt.show()


def get_variable_value(variable_name, current_line, content, split_symbol='='):
    """
    Get the value of one variable from the content of the information file.

    Args:
        variable_name: the variable name.
        current_line: current line number.
        content: the content of the information file.
    Returns:
        variable_value: the value of the variable.
        line_number: current line number after the operations.
    """
    line_number = current_line
    while not variable_name in content[line_number]:
        line_number += 1
    single_line = content[line_number]
    line_splits = single_line.split(split_symbol)
    variable_value = float(line_splits[1])
    return (variable_value, line_number)


def get_num_steps(run_dir):
    """
    """
    fname = run_dir + 'init.dat'
    with open(fname) as f:
        content = f.readlines()
    nlines = len(content)
    current_line = 0
    ttot, current_line = get_variable_value('Total tracking time',
                                            current_line, content,
                                            split_symbol=':')
    dt, current_line = get_variable_value('Tracking time step',
                                          current_line, content,
                                          split_symbol=':')
    nstep = int(ttot / dt)
    return nstep


def plot_particle_trajectory(run_dir):
    """
    """
    fname = run_dir + 'data/particle_diagnostics.h5'
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    ntraj = get_num_steps(run_dir) / 10 + 1
    if ntraj > 1E6:
        ntraj = 1E6
    nptl = sz / ntraj
    for iptl in range(10):
        px = np.array(dset_ptl['x'][iptl*ntraj:(iptl+1)*ntraj])
        py = np.array(dset_ptl['y'][iptl*ntraj:(iptl+1)*ntraj])
        pz = np.array(dset_ptl['z'][iptl*ntraj:(iptl+1)*ntraj])
        ux = np.array(dset_ptl['ux'][iptl*ntraj:(iptl+1)*ntraj])
        uy = np.array(dset_ptl['uy'][iptl*ntraj:(iptl+1)*ntraj])
        uz = np.array(dset_ptl['uz'][iptl*ntraj:(iptl+1)*ntraj])
        t = np.array(dset_ptl['t'][iptl*ntraj:(iptl+1)*ntraj-1])
        Ex = np.array(dset_emf['Ex'][iptl*ntraj:(iptl+1)*ntraj])
        Ey = np.array(dset_emf['Ey'][iptl*ntraj:(iptl+1)*ntraj])
        Ez = np.array(dset_emf['Ez'][iptl*ntraj:(iptl+1)*ntraj])
        Bx = np.array(dset_emf['Bx'][iptl*ntraj:(iptl+1)*ntraj])
        By = np.array(dset_emf['By'][iptl*ntraj:(iptl+1)*ntraj])
        Bz = np.array(dset_emf['Bz'][iptl*ntraj:(iptl+1)*ntraj])
        gamma = np.sqrt(1.0 + ux*ux + uy*uy + uz*uz)
        ene = gamma - 1.0
        xs, ys = 0.15, 0.15
        w1, h1 = 0.8, 0.8
        fig = plt.figure(figsize=[7, 5])
        ax1 = fig.add_axes([xs, ys, w1, h1], projection='3d')
        ax1.plot(px, py, pz)
        ax1.tick_params(labelsize=16)
        plt.show()
    # file.close


def transfer_to_h5part(run_dir):
    """Transfer current HDF5 file to H5Part format
    
    All particles at the same time step are stored in the same time step

    Args:
        run_dir: simulation directory
    """
    fname = run_dir + 'data/particle_diagnostics.h5'
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    nstep = get_num_steps(run_dir) / 10 + 1
    nptl = sz / nstep
    fname_out = run_dir + 'data/particle_diagnostics.h5part'
    toffset = 100
    with h5py.File(fname_out, 'w') as fh_out:
        for tf in range(0, nstep, toffset):
            print("Time frame: %d" % tf)
            x = np.array(dset_ptl['x'][tf::nstep])
            y = np.array(dset_ptl['y'][tf::nstep])
            z = np.array(dset_ptl['z'][tf::nstep])
            ux = np.array(dset_ptl['ux'][tf::nstep])
            uy = np.array(dset_ptl['uy'][tf::nstep])
            uz = np.array(dset_ptl['uz'][tf::nstep])
            gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
            t = np.array(dset_ptl['t'][tf::nstep])
            Ex = np.array(dset_emf['Ex'][tf::nstep])
            Ey = np.array(dset_emf['Ey'][tf::nstep])
            Ez = np.array(dset_emf['Ez'][tf::nstep])
            Bx = np.array(dset_emf['Bx'][tf::nstep])
            By = np.array(dset_emf['By'][tf::nstep])
            Bz = np.array(dset_emf['Bz'][tf::nstep])
            grp = fh_out.create_group('Step#' + str(tf//toffset))
            grp.create_dataset('x', (nptl, ), data=x)
            grp.create_dataset('y', (nptl, ), data=y)
            grp.create_dataset('z', (nptl, ), data=z)
            grp.create_dataset('ux', (nptl, ), data=ux)
            grp.create_dataset('uy', (nptl, ), data=uy)
            grp.create_dataset('uz', (nptl, ), data=uz)
            grp.create_dataset('gamma', (nptl, ), data=gamma)
            grp.create_dataset('t', (nptl, ), data=t)
            grp.create_dataset('Ex', (nptl, ), data=Ex)
            grp.create_dataset('Ey', (nptl, ), data=Ey)
            grp.create_dataset('Ez', (nptl, ), data=Ez)
            grp.create_dataset('Bx', (nptl, ), data=Bx)
            grp.create_dataset('By', (nptl, ), data=By)
            grp.create_dataset('Bz', (nptl, ), data=Bz)


def plot_diffusion_coefficients(run_dir, run_name):
    """
    """
    fname = run_dir + 'data/diff_coeffs.dat'
    fdata = np.genfromtxt(fname, skip_header=1)
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fdir = '../img/' + run_name + '/'
    mkdir_p(fdir)
    def plot_one_diff(index, label1, img_name):
        fig = plt.figure(figsize=[7, 5])
        ax1 = fig.add_axes([xs, ys, w1, h1])
        nptl = fdata[:, -1]
        ax1.plot(fdata[:, index] / nptl, linewidth=2, label=label1, color='k')
        ax1.tick_params(labelsize=16)
        leg = ax1.legend(loc=2, prop={'size': 20}, ncol=1,
                         shadow=False, fancybox=False, frameon=False)
        fname = fdir + img_name + '.eps'
        fig.savefig(fname)
    plot_one_diff(0, r'$D_{rr}$', 'drr')
    plot_one_diff(1, r'$D_{xx}$', 'dxx')
    plot_one_diff(2, r'$D_{yy}$', 'dyy')
    plot_one_diff(3, r'$D_{zz}$', 'dzz')
    plot_one_diff(4, r'$D_{pp}$', 'dpp')
    plot_one_diff(5, r'$D_{\mu\mu}$', 'duu')
    plot_one_diff(6, r'$D_{aa}$', 'daa')
    plot_one_diff(7, r'$D_{ee}$', 'dee')
    plt.show()


if __name__ == "__main__":
    run_dir = '../bin/3d_reconnection/'
    run_name = '3d_reconnection'
    # plot_spectrum(run_dir)
    # plot_particle_trajectory(run_dir)
    plot_diffusion_coefficients(run_dir, run_name)
    # transfer_to_h5part(run_dir)
