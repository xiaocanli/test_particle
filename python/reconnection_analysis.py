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
    ax1.set_xlim([1E-4, 1E0])
    ax1.grid(True)
    
    plt.show()


def plot_particle_trajectory(run_dir):
    """
    """
    fname = run_dir + 'data/particle_diagnostics.h5'
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    ntraj = 2000
    nptl = sz / ntraj
    # subplot positions
    xs = 0.15
    xe = 0.8
    yint = 0.2
    for iptl in range(10):
        ye = 0.95
        px = np.array(dset_ptl['x'][iptl*ntraj+1:(iptl+1)*ntraj])
        py = np.array(dset_ptl['y'][iptl*ntraj+1:(iptl+1)*ntraj])
        pz = np.array(dset_ptl['z'][iptl*ntraj+1:(iptl+1)*ntraj])
        ux = np.array(dset_ptl['ux'][iptl*ntraj+1:(iptl+1)*ntraj])
        uy = np.array(dset_ptl['uy'][iptl*ntraj+1:(iptl+1)*ntraj])
        uz = np.array(dset_ptl['uz'][iptl*ntraj+1:(iptl+1)*ntraj])
        t = np.array(dset_ptl['t'][iptl*ntraj+1:(iptl+1)*ntraj-1])
        Ex = np.array(dset_emf['Ex'][iptl*ntraj+1:(iptl+1)*ntraj])
        Ey = np.array(dset_emf['Ey'][iptl*ntraj+1:(iptl+1)*ntraj])
        Ez = np.array(dset_emf['Ez'][iptl*ntraj+1:(iptl+1)*ntraj])
        Bx = np.array(dset_emf['Bx'][iptl*ntraj+1:(iptl+1)*ntraj])
        By = np.array(dset_emf['By'][iptl*ntraj+1:(iptl+1)*ntraj])
        Bz = np.array(dset_emf['Bz'][iptl*ntraj+1:(iptl+1)*ntraj])
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


def plot_diffusion_coefficients(run_dir):
    """
    """
    fname = run_dir + 'data/diff_coeffs.dat'
    fdata = np.genfromtxt(fname)
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig = plt.figure(figsize=[7, 5])
    ax1 = fig.add_axes([xs, ys, w1, h1])
    ax1.plot(fdata[1:, 1])
    ax1.tick_params(labelsize=16)
    plt.show()


if __name__ == "__main__":
    run_dir = '../bin/3d_reconnection/'
    # plot_spectrum(run_dir)
    # plot_particle_trajectory(run_dir)
    plot_diffusion_coefficients(run_dir)
