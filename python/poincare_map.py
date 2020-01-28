"""
Plotting procedures for Poincaré Map
"""
import argparse
import math
import multiprocessing
import subprocess

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
import palettable
import pandas as pd
from evtk.hl import pointsToVTK
from joblib import Parallel, delayed
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage.filters import gaussian_filter, median_filter

from shell_functions import mkdir_p

plt.style.use("seaborn-deep")
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = \
[r"\usepackage{amsmath, bm}",
 r"\DeclareMathAlphabet{\mathsfit}{\encodingdefault}{\sfdefault}{m}{sl}",
 r"\SetMathAlphabet{\mathsfit}{bold}{\encodingdefault}{\sfdefault}{bx}{sl}",
 r"\newcommand{\tensorsym}[1]{\bm{\mathsfit{#1}}}"]
COLORS = palettable.colorbrewer.qualitative.Set1_9.mpl_colors


def plot_poincare_map(plot_config):
    """
    """
    run_dir = plot_config['run_dir']
    run_name = plot_config['run_name']
    tframe = plot_config['tframe']
    fname = run_dir + 'data/poincare_map_' + str(tframe) + '.h5'
    fh = h5py.File(fname, 'r')
    grp = fh['poincare_map']
    dset = grp['xpos']
    sz = dset.shape
    xpos = np.zeros(sz, dtype=dset.dtype)
    ypos = np.zeros(sz, dtype=dset.dtype)
    zpos = np.zeros(sz, dtype=dset.dtype)
    dset.read_direct(xpos)
    dset = grp['ypos']
    dset.read_direct(ypos)
    dset = grp['zpos']
    dset.read_direct(zpos)
    fh.close()
    fig = plt.figure(figsize=[18, 12])
    rect = [0.10, 0.80, 0.8, 0.18]
    nh, nv = 3, 1
    COLORS = palettable.tableau.Tableau_10.mpl_colors
    for i in range(0, nh):
        ax = fig.add_axes(rect)
        rect[1] -= rect[3] + 0.01
        ax.set_prop_cycle('color', COLORS)
        for j in range(0, nv):
            point = j * nh + i
            p1 = ax.scatter(xpos[point, 0], zpos[point, 0], s=50,
                            color=COLORS[j])
            ax.scatter(xpos[point, :], zpos[point, :], s=2,
                       color=COLORS[j])
        # ax.set_xlim([0, 750])
        # ax.set_ylim([-120, 120])
    plt.show()


def poincare_map_hui(plot_config):
    """
    """
    run_dir = plot_config['run_dir']
    run_name = plot_config['run_name']
    tframe = plot_config['tframe']
    fname = run_dir + 'data/poincare_map_' + str(tframe) + '.h5'
    fh = h5py.File(fname, 'r')
    grp = fh['poincare_map']
    dset = grp['xpos']
    sz = dset.shape
    xpos = np.zeros(sz, dtype=dset.dtype)
    ypos = np.zeros(sz, dtype=dset.dtype)
    zpos = np.zeros(sz, dtype=dset.dtype)
    dset.read_direct(xpos)
    dset = grp['ypos']
    dset.read_direct(ypos)
    dset = grp['zpos']
    dset.read_direct(zpos)
    fh.close()
    nh, nv = 101, 1
    i = 2
    xmin, xmax = 0, 1000
    COLORS = palettable.tableau.Tableau_10.mpl_colors
    # fig = plt.figure(figsize=[14, 5])
    fig = plt.figure(figsize=[10, 5])
    rect = [0.08, 0.11, 0.9, 0.86]
    ax = fig.add_axes(rect)
    ax.set_prop_cycle('color', COLORS)
    for j in range(0, nv):
        for i in range(0, nh):
            point = j * nh + i
            # rect[1] -= rect[3] + 0.01
            color = plt.cm.jet((point + 0.5)/float(nh*nv), 1)
            cond = xpos[point, :] < xmax
            ax.scatter(xpos[point, cond], zpos[point, cond], marker='o', s=1,
                       color=COLORS[0])
            # p1 = ax.scatter(xpos[point, 0], zpos[point, 0], s=50,
            #                 color=COLORS[1])
            # # ax.tick_params(axis='x', labelbottom=False)
            # # ax.tick_params(axis='y', labelleft=False)
    ax.tick_params(bottom=True, top=False, left=True, right=True)
    ax.tick_params(axis='x', which='minor', direction='in', top=True)
    ax.tick_params(axis='x', which='major', direction='in', top=True)
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    ax.set_xlim([0, 1000])
    ax.set_ylim([-250, 250])
    # ax.set_xlim([0, 750])
    # ax.set_ylim([-125, 125])
    ax.set_xlabel(r'$x/d_e$', fontsize=16)
    ax.set_ylabel(r'$z/d_e$', fontsize=16, labelpad=0)
    ax.tick_params(labelsize=12)

    # fdir = '../img/poincare_map/tframe_' + str(tframe) + '/'
    # mkdir_p(fdir)
    # fname = fdir + 'poincare_map_' + str(point) + '.jpg'
    fdir = '../img/poincare_map/'
    mkdir_p(fdir)
    fname = fdir + 'poincare_map_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=400)
    # plt.close()
    plt.show()


def get_cmd_args():
    """Get command line arguments
    """
    default_run = 'trinity_run'
    default_run_dir = ('/net/scratch3/xiaocanli/diffusion/' + default_run + '/')
    parser = argparse.ArgumentParser(description='Analysis for Poincare map')
    parser.add_argument('--run_name', action="store",
                        default=default_run, help='PIC/MHD run name')
    parser.add_argument('--run_dir', action="store",
                        default=default_run_dir, help='PIC/MHD run directory')
    parser.add_argument('--tinterval', action="store", default='2217', type=int,
                        help='Time interval')
    parser.add_argument('--tframe', action="store", default='10', type=int,
                        help='Time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to use a time loop to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='40', type=int,
                        help='ending time frame')
    parser.add_argument('--poincare_map', action="store_true", default=False,
                        help='whether to plot Poincare map')
    parser.add_argument('--poincare_map_hui', action="store_true", default=False,
                        help='whether to plot Poincare map for Hui')
    return parser.parse_args()


def analysis_single_frames(plot_config, args):
    """Analysis for multiple time frames
    """
    tframe = args.tframe
    if args.poincare_map:
        plot_poincare_map(plot_config)
    if args.poincare_map_hui:
        poincare_map_hui(plot_config)


def process_input(plot_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    pass


def analysis_multi_frames(plot_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tstart"], plot_config["tend"] + 1)
    if args.time_loop:
        pass
    else:
        ncores = multiprocessing.cpu_count()
        ncores = 8
        Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, args, tframe)
                                for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    plot_config = {}
    plot_config["run_name"] = args.run_name
    plot_config["run_dir"] = args.run_dir
    plot_config["tframe"] = args.tframe
    plot_config["tstart"] = args.tstart
    plot_config["tend"] = args.tend
    plot_config["tinterval"] = args.tinterval
    if args.multi_frames:
        analysis_multi_frames(plot_config, args)
    else:
        analysis_single_frames(plot_config, args)


if __name__ == "__main__":
    main()
