"""
Plotting procedures for reconnection
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

FONT = {'family': 'serif',
        'color': 'black',
        'weight': 'normal',
        'size': 24}

def plot_spectrum(plot_config):
    """ plotting particle energy spectrum.
    """
    run_dir = plot_config["run_dir"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/espectrum.dat')
    fdata = np.genfromtxt(fname)
    nframes, nbins = fdata.shape

    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig = plt.figure(figsize=[7, 5])
    ax1 = fig.add_axes([xs, ys, w1, h1])
    for tframe in range(1, nframes, 50):
        ax1.loglog(fdata[0, :], fdata[tframe, :], linewidth=2)

    ax1.set_xlabel(r'$\gamma - 1$', fontdict=FONT, fontsize=20)
    ax1.set_ylabel(r'$f(\gamma - 1)$', fontdict=FONT, fontsize=20)
    ax1.tick_params(labelsize=16)
    ax1.set_xlim([1E-2, 2E1])
    ax1.set_ylim([1E0, 1E6])
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
    """Get total number of time steps for particle tracking
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
    print("Total number of steps: %d" % nstep)
    return int(nstep)


def get_traj_tinterval(run_dir):
    """Get time interval for particle trajectory diagnostics
    """
    fname = run_dir + 'init.dat'
    with open(fname) as f:
        content = f.readlines()
    nlines = len(content)
    current_line = 0
    msg = 'Time step interval for trajectory diagnostics'
    tinterval_traj, current_line = get_variable_value(msg, current_line,
                                                      content, split_symbol=':')
    print("Time interval for trajectory diagnostics: %d" % tinterval_traj)
    return int(tinterval_traj)


def get_domain_sizes(run_dir):
    """Get domain sizes
    """
    fname = run_dir + 'init.dat'
    with open(fname) as f:
        content = f.readlines()
    nlines = len(content)
    current_line = 0
    xmin, current_line = get_variable_value('xmin = ', current_line,
                                            content, split_symbol='=')
    xmax, current_line = get_variable_value('xmax = ', current_line,
                                            content, split_symbol='=')
    ymin, current_line = get_variable_value('ymin = ', current_line,
                                            content, split_symbol='=')
    ymax, current_line = get_variable_value('ymax = ', current_line,
                                            content, split_symbol='=')
    zmin, current_line = get_variable_value('zmin = ', current_line,
                                            content, split_symbol='=')
    zmax, current_line = get_variable_value('zmax = ', current_line,
                                            content, split_symbol='=')
    lx = xmax - xmin
    ly = ymax - ymin
    lz = zmax - zmin
    return (lx, ly, lz)


def adjust_pos(pos, length):
    """Adjust position for periodic boundary conditions.

    Args:
        pos: the position along one axis
        length: the box size along that axis
    """
    crossings = []
    offsets = []
    offset = 0
    nt, = pos.shape
    pos_b = np.zeros(nt)
    pos_b = np.copy(pos)
    for i in range(nt - 1):
        if (pos[i] - pos[i + 1] > 0.1 * length):
            crossings.append(i)
            offset += length
            offsets.append(offset)
        if (pos[i] - pos[i + 1] < -0.1 * length):
            crossings.append(i)
            offset -= length
            offsets.append(offset)
    nc = len(crossings)
    if nc > 0:
        crossings = np.asarray(crossings)
        offsets = np.asarray(offsets)
        for i in range(nc - 1):
            pos_b[crossings[i] + 1:crossings[i + 1] + 1] += offsets[i]
        pos_b[crossings[nc - 1] + 1:] += offsets[nc - 1]
    return pos_b


def plot_particle_trajectory(plot_config):
    """
    """
    run_dir = plot_config["run_dir"]
    run_name = plot_config["run_name"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    tinterval_traj = get_traj_tinterval(run_dir)
    nsteps_tot = get_num_steps(run_dir)
    if nsteps_tot > 1E6:
        nsteps_tot = int(1E6)
    ntraj = nsteps_tot // tinterval_traj + 1

    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/particle_diagnostics.h5')
    with h5py.File(fname,'r') as fh:
        group = fh['/particles_fields']
        dset_ptl = group['particles']
        dset_emf = group['fields']
        sz, = dset_ptl.shape
        nptl = sz // ntraj
        x = np.array(dset_ptl['x']).reshape((nptl, -1))
        y = np.array(dset_ptl['y']).reshape((nptl, -1))
        z = np.array(dset_ptl['z']).reshape((nptl, -1))
        ux = np.array(dset_ptl['ux']).reshape((nptl, -1))
        uy = np.array(dset_ptl['uy']).reshape((nptl, -1))
        uz = np.array(dset_ptl['uz']).reshape((nptl, -1))
        t = np.array(dset_ptl['t']).reshape((nptl, -1))
        Ex = np.array(dset_emf['Ex']).reshape((nptl, -1))
        Ey = np.array(dset_emf['Ey']).reshape((nptl, -1))
        Ez = np.array(dset_emf['Ez']).reshape((nptl, -1))
        Bx = np.array(dset_emf['Bx']).reshape((nptl, -1))
        By = np.array(dset_emf['By']).reshape((nptl, -1))
        Bz = np.array(dset_emf['Bz']).reshape((nptl, -1))
    lx, ly, lz = 750, 375, 312.5
    dxx = np.zeros(ntraj)
    nptl_select = 0
    for iptl in range(nptl):
        px = x[iptl, :]
        px = adjust_pos(px, 750)
        sign_switch = np.where((np.diff(np.sign(np.gradient(px))) != 0) * 1 == 1)
        sign_switch = np.squeeze(sign_switch)
        sign_switch = np.insert(sign_switch, 0, 0)
        sign_switch = np.append(sign_switch, ntraj - 1)
        interval = np.diff(sign_switch)
        sz, = interval.shape
        if sz > 0:
            max_interval = np.max(interval)
        else:
            max_interval = ntraj
        # if (np.max(px) - np.min(px)) < 500:
        if max_interval < ntraj // 4:
            nptl_select += 1
            dxx += (px - px[0])**2
    print("Number of selected particles: %d" % nptl_select)
    dxx /= np.linspace(1, 10000, ntraj) * nptl_select * 2
    plt.plot(dxx)
    plt.show()
    fdir1 = '../img/' + run_name + '/trajectory/'
    fdir2 = '../img/' + run_name + '/trajectory_scattering/'
    fdir3 = '../img/' + run_name + '/trajectory_streaming/'
    mkdir_p(fdir1)
    mkdir_p(fdir2)
    mkdir_p(fdir3)
    # for iptl in range(nptl):
    #     print("Particle: %d" % iptl)
    #     px = adjust_pos(x[iptl, :], 750)
    #     py = adjust_pos(y[iptl, :], 375)
    #     pz = adjust_pos(z[iptl, :], 312.5)
    #     pt = t[iptl, :]

    #     # fig = plt.figure(figsize=[14, 5])
    #     # rect = [0.07, 0.15, 0.4, 0.8]
    #     # ax1 = fig.add_axes(rect, projection='3d')
    #     # ax1.plot(px, py, pz)
    #     # ax1.tick_params(labelsize=12)
    #     # ax1.set_xlabel(r'$x/d_e$', fontdict=FONT, fontsize=16)
    #     # ax1.set_ylabel(r'$y/d_e$', fontdict=FONT, fontsize=16)
    #     # ax1.set_zlabel(r'$z/d_e$', fontdict=FONT, fontsize=16)

    #     # rect[0] += rect[2] + 0.1
    #     # ax2 = fig.add_axes(rect)
    #     # ax2.plot(px, pz)
    #     # ax2.tick_params(labelsize=12)
    #     # ax2.set_xlabel(r'$x/d_e$', fontdict=FONT, fontsize=16)
    #     # ax2.set_ylabel(r'$z/d_e$', fontdict=FONT, fontsize=16)

    #     fig = plt.figure(figsize=[7, 5])
    #     rect = [0.15, 0.15, 0.8, 0.8]
    #     ax1 = fig.add_axes(rect)
    #     ax1.plot(pt, px)
    #     ax1.set_xlim([pt.min(), pt.max()])
    #     ax1.tick_params(labelsize=12)
    #     ax1.set_xlabel(r'$t\omega_{pe}$', fontdict=FONT, fontsize=16)
    #     ax1.set_ylabel(r'$x/d_e$', fontdict=FONT, fontsize=16)
    #     sign_switch = np.where((np.diff(np.sign(np.gradient(px))) != 0) * 1 == 1)
    #     sign_switch = np.squeeze(sign_switch)
    #     sign_switch = np.insert(sign_switch, 0, 0)
    #     sign_switch = np.append(sign_switch, ntraj - 1)
    #     interval = np.diff(sign_switch)
    #     sz, = interval.shape
    #     if sz > 0:
    #         max_interval = np.max(interval)
    #     else:
    #         max_interval = ntraj
    #     # if (np.max(px) - np.min(px)) < 2000:
    #     if max_interval < ntraj // 4:
    #         fname = fdir2 + 'traj_' + str(iptl) + '.pdf'
    #     else:
    #         fname = fdir3 + 'traj_' + str(iptl) + '.pdf'
    #     fig.savefig(fname)
    #     plt.close()
    # #     plt.show()


def transfer_to_h5part(plot_config):
    """Transfer current HDF5 file to H5Part format

    All particles at the same time step are stored in the same time step
    """
    run_dir = plot_config["run_dir"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/particle_diagnostics.h5')
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    tinterval_traj = get_traj_tinterval(run_dir)
    nsteps_tot = get_num_steps(run_dir)
    if nsteps_tot > 1E6:
        nsteps_tot = int(1E6)
    ntraj = nsteps_tot // tinterval_traj + 1
    nptl = sz / ntraj
    fname_out = (run_dir + 'data_' + str(tindex) + '_' +
                 str(ptl_vel) + 'c/particle_diagnostics.h5part')
    print(fname_out)
    toffset = 10
    with h5py.File(fname_out, 'w') as fh_out:
        for tf in range(0, ntraj, toffset):
            print("Time frame: %d" % tf)
            x = np.array(dset_ptl['x'][tf::ntraj])
            y = np.array(dset_ptl['y'][tf::ntraj])
            z = np.array(dset_ptl['z'][tf::ntraj])
            ux = np.array(dset_ptl['ux'][tf::ntraj])
            uy = np.array(dset_ptl['uy'][tf::ntraj])
            uz = np.array(dset_ptl['uz'][tf::ntraj])
            gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
            t = np.array(dset_ptl['t'][tf::ntraj])
            Ex = np.array(dset_emf['Ex'][tf::ntraj])
            Ey = np.array(dset_emf['Ey'][tf::ntraj])
            Ez = np.array(dset_emf['Ez'][tf::ntraj])
            Bx = np.array(dset_emf['Bx'][tf::ntraj])
            By = np.array(dset_emf['By'][tf::ntraj])
            Bz = np.array(dset_emf['Bz'][tf::ntraj])
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


def transfer_to_csv(plot_config):
    """Transfer current HDF5 file to CSV

    Each CSV file contains one trajectory

    """
    run_dir = plot_config["run_dir"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/particle_diagnostics.h5')
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    tinterval_traj = get_traj_tinterval(run_dir)
    nsteps_tot = get_num_steps(run_dir)
    if nsteps_tot > 1E6:
        nsteps_tot = int(1E6)
    ntraj = nsteps_tot // tinterval_traj + 1
    nptl = sz / ntraj
    fdir = run_dir + 'data_' + str(tindex) + '_' + str(ptl_vel) + 'c/'
    fdir += 'traj_csv/'
    mkdir_p(fdir)
    pdata = np.zeros([14, ntraj])
    # for iptl in range(nptl):
    for iptl in range(2):
        print(iptl)
        ps, pt = ntraj * iptl, ntraj * (iptl + 1)
        pdata[0] = np.array(dset_ptl['x'][ps:pt])
        pdata[1] = np.array(dset_ptl['y'][ps:pt])
        pdata[2] = np.array(dset_ptl['z'][ps:pt])
        pdata[3] = np.array(dset_ptl['ux'][ps:pt])
        pdata[4] = np.array(dset_ptl['uy'][ps:pt])
        pdata[5] = np.array(dset_ptl['uz'][ps:pt])
        pdata[6] = np.sqrt(1.0 + np.sum(pdata[3:6]**2, axis=0))
        pdata[7] = np.array(dset_ptl['t'][ps:pt])
        pdata[8] = np.array(dset_emf['Ex'][ps:pt])
        pdata[9] = np.array(dset_emf['Ey'][ps:pt])
        pdata[10] = np.array(dset_emf['Ez'][ps:pt])
        pdata[11] = np.array(dset_emf['Bx'][ps:pt])
        pdata[12] = np.array(dset_emf['By'][ps:pt])
        pdata[13] = np.array(dset_emf['Bz'][ps:pt])
        fname = fdir + 'traj_' + str(iptl) + '.csv'
        # np.savetxt(fname, pdata.T, delimiter=",",
        #            header="x,y,z,ux,uy,uz,gamma,t,Ex,Ey,Ez,Bx,By,Bz")
        df = pd.DataFrame(pdata.T)
        df.to_csv(fname, mode='w', index=True,
                  header=["x", "y", "z", "ux", "uy", "uz", "gamma", "t",
                          "Ex", "Ey", "Ez", "Bx", "By", "Bz"])


def get_crossings(pos, length):
    """Get the crossing points along one axis

    Args:
        pos: the position along one axis
        length: the box size along that axis
    """
    crossings = []
    offsets = []
    offset = 0
    nt, = pos.shape
    pos_b = np.zeros(nt)
    pos_b = np.copy(pos)
    for i in range(nt-1):
        if (pos[i]-pos[i+1] > 0.1*length):
            crossings.append(i)
            offset += length
            offsets.append(offset)
        if (pos[i]-pos[i+1] < -0.1*length):
            crossings.append(i)
            offset -= length
            offsets.append(offset)
    nc = len(crossings)
    if nc > 0:
        crossings = np.asarray(crossings)
        offsets = np.asarray(offsets)
    return (nc, crossings, offsets)


def adjust_pos(pos, length):
    """Adjust position for periodic boundary conditions.

    Args:
        pos: the position along one axis
        length: the box size along that axis
    """
    nc, crossings, offsets = get_crossings(pos, length)
    if nc > 0:
        for i in range(nc - 1):
            pos_b[crossings[i] + 1:crossings[i + 1] + 1] += offsets[i]
        pos_b[crossings[nc - 1] + 1:] += offsets[nc - 1]
    return pos_b


def piecewise_trajectory(plot_config):
    """Save piecewise trajectory
    """
    iptl = plot_config['iptl']
    tint = plot_config['tint']
    run_dir = plot_config["run_dir"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/particle_diagnostics.h5')
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    tinterval_traj = get_traj_tinterval(run_dir)
    nsteps_tot = get_num_steps(run_dir)
    if nsteps_tot > 1E6:
        nsteps_tot = int(1E6)
    nframes = nsteps_tot // tinterval_traj + 1
    nptl = sz / nframes
    lx, ly, lz = get_domain_sizes(run_dir)
    pdata = np.zeros([6, nframes])
    ps, pt = nframes * iptl, nframes * (iptl + 1)
    pdata[0] = np.array(dset_ptl['x'][ps:pt])
    pdata[1] = np.array(dset_ptl['y'][ps:pt])
    pdata[2] = np.array(dset_ptl['z'][ps:pt])
    pdata[3] = np.array(dset_ptl['ux'][ps:pt])
    pdata[4] = np.array(dset_ptl['uy'][ps:pt])
    pdata[5] = np.array(dset_ptl['uz'][ps:pt])
    fdir = run_dir + 'data_' + str(tindex) + '_' + str(ptl_vel) + 'c/'
    fdir += 'piecewise_trajectory/ptl_' + str(iptl) + '/'
    mkdir_p(fdir)
    nc_x, cross_x, offsets_x = get_crossings(pdata[0], lx)
    nc_y, cross_y, offsets_y = get_crossings(pdata[1], ly)
    nc_z, cross_z, offsets_z = get_crossings(pdata[2], lz)
    crossings = np.unique(np.sort(np.hstack((cross_x, cross_y, cross_z))))
    ncross, = crossings.shape
    icross1 = 0
    icross2 = 0
    for tframe in range(nframes):
        if tframe > crossings[icross2]:
            icross1 = icross2
            if icross2 < ncross - 1:
                icross2 += 1
            else:
                icross2 = ncross - 1
        neighbors = [0, crossings[icross1],
                     tframe - tint, tframe, tframe + tint,
                     crossings[icross2], nframes]
        neighbors = sorted(list(set(neighbors)))
        tindex = neighbors.index(tframe)
        if tframe == 0:
            ts = 0
            te = int(neighbors[tindex + 1])
        elif tframe == crossings[icross2]:
            ts = int(neighbors[tindex - 1])
            te = tframe + 1
        else:
            ts = int(neighbors[tindex - 1])
            te = int(neighbors[tindex + 1])
        fname = fdir + "tframe_" + str(tframe)
        pointsToVTK(fname, pdata[0, ts+1:te], pdata[1, ts+1:te], pdata[2, ts+1:te],
                    data = {"ux": pdata[3, ts+1:te],
                            "uy": pdata[4, ts+1:te],
                            "uz": pdata[5, ts+1:te]})


def diffusion_coefficients_multi(plot_config):
    """
    """
    run_dir = plot_config["run_dir"]
    run_name = plot_config["run_name"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    vels = np.linspace(0.3, 0.9, 7)
    vels = np.append(vels, 0.95)
    gammas = 1.0 / np.sqrt(1 - vels**2)
    moms = vels * gammas
    tframes = np.linspace(10, 20, 11, dtype=int)
    # tframes = np.linspace(10, 10, 1, dtype=int)
    tindices = tframes * plot_config["tinterval"]
    dxxs = np.zeros((len(tframes), len(vels)))

    fdir = '../img/' + run_name + '/'
    mkdir_p(fdir)
    for itindex, tindex in enumerate(tindices):
        fig = plt.figure(figsize=[7, 2.5])
        rect = [0.08, 0.18, 0.39, 0.76]
        xlen = fig.get_figwidth()
        ylen = fig.get_figheight()
        ax1 = fig.add_axes(rect)
        rect[0] += rect[2] + 0.1
        ax2 = fig.add_axes(rect)
        for ivel, vel in enumerate(vels[::-1]):
            fname = (run_dir + 'data_' + str(tindex) + '_' +
                     str(vel) + 'c/diff_coeffs.dat')
            fdata = np.genfromtxt(fname, skip_header=1)
            nptl = fdata[:, -1]
            kappa = fdata[:, 1] / nptl
            nframes, = kappa.shape
            tfs = np.arange(nframes)
            tstart = nframes//4
            # coeff = np.polyfit(tfs[tstart:], kappa[tstart:], 4)
            # p = np.poly1d(coeff)
            # kappa_fit = p(tfs)
            # kappa_fit -= kappa_fit[0]
            coeff = np.polyfit(tfs[tstart:], np.log(kappa[tstart:]), 3)
            p = np.poly1d(coeff)
            kappa_fit = np.exp(p(tfs))
            kappa_fit -= kappa_fit[0]
            kappa_adjust = kappa - kappa_fit
            label1 = r'$' + str(vel) + 'c$'
            if ivel == 4:
                p1, = ax1.plot(tfs, kappa, linewidth=1, linestyle='--', color='k')
                p2, = ax1.plot(tfs, kappa_adjust, linewidth=1,
                               linestyle='-', color=p1.get_color())
                p3, = ax1.plot(tfs, kappa_fit, linewidth=1,
                               linestyle='-.', color=p1.get_color())
                vel_str = r'$' + str(vel) + 'c$'
                ax1.text(0.02, 0.95, vel_str, color='k', fontsize=10,
                         bbox=dict(facecolor='none', alpha=1.0,
                                   edgecolor='none', pad=10.0),
                         horizontalalignment='left', verticalalignment='top',
                         transform=ax1.transAxes)
                xlim = ax1.get_xlim()
                ylim = ax1.get_ylim()
                lx = xlim[1] - xlim[0]
                ly = ylim[1] - ylim[0]
                t0 = 7000
                xextent = 2000 / lx
                yextent = (kappa[t0+1000] - kappa[t0-1000]) / ly
                angle = math.atan((yextent * ylen)/(xextent * xlen))
                angle *= 180 / math.pi
                angle += 10
                yoffset = (kappa[t0] - ylim[0]) / (ylim[1] - ylim[0]) - 0.03
                ax1.text(0.6, yoffset, 'all particles', color='k',
                         fontsize=10, rotation=angle,
                         bbox=dict(facecolor='none', alpha=1.0,
                                   edgecolor='none', pad=10.0),
                         horizontalalignment='left', verticalalignment='bottom',
                         transform=ax1.transAxes)
                yoffset = (kappa_fit[t0] - ylim[0]) / (ylim[1] - ylim[0]) - 0.03
                ax1.text(0.6, yoffset, 'streaming particles', color='k',
                         fontsize=10, rotation=angle,
                         bbox=dict(facecolor='none', alpha=1.0,
                                   edgecolor='none', pad=10.0),
                         horizontalalignment='left', verticalalignment='bottom',
                         transform=ax1.transAxes)
                yoffset = (kappa_adjust[t0] - ylim[0]) / (ylim[1] - ylim[0]) + 0.01
                ax1.text(0.8, yoffset, 'rest', color='k', fontsize=10,
                         bbox=dict(facecolor='none', alpha=1.0,
                                   edgecolor='none', pad=10.0),
                         horizontalalignment='left', verticalalignment='bottom',
                         transform=ax1.transAxes)
            p2, = ax2.plot(tfs, kappa_adjust, linewidth=1, label=label1,
                           linestyle='-', color='k')
            kappa_mean = np.mean(kappa_adjust[nframes//5:])
            dxxs[itindex, ivel] = kappa_mean
            ylim = ax2.get_ylim()
            k0 = (kappa_adjust[9000] - ylim[0]) / (ylim[1] - ylim[0])
            vel_str = r'$' + str(vel) + 'c$'
            ax2.text(0.98, k0-0.01, vel_str, color='k', fontsize=10,
                     bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
                     horizontalalignment='right', verticalalignment='top',
                     transform=ax2.transAxes)
        ax1.set_xlim([tfs.min(), tfs.max()])
        ax1.tick_params(labelsize=10)
        ax1.tick_params(axis='x', which='minor', direction='in', top='on')
        ax1.tick_params(axis='x', which='major', direction='in', top='on')
        ax1.tick_params(axis='y', which='minor', direction='in', left='on')
        ax1.tick_params(axis='y', which='major', direction='in')
        ax1.set_xlabel(r'$t\omega_{pe}$', fontdict=FONT, fontsize=12)
        # ax1.set_ylabel(r'$d_e^2\omega_{pe}$', fontdict=FONT, fontsize=12)
        ax1.set_ylabel(r'$D_{xx}/(v_Ad_i)$', fontdict=FONT, fontsize=12)

        ax2.set_xlim([tfs.min(), tfs.max()])
        ax2.tick_params(labelsize=10)
        ax2.tick_params(axis='x', which='minor', direction='in', top='on')
        ax2.tick_params(axis='x', which='major', direction='in', top='on')
        ax2.tick_params(axis='y', which='minor', direction='in', left='on')
        ax2.tick_params(axis='y', which='major', direction='in')
        ax2.set_xlabel(r'$t\omega_{pe}$', fontdict=FONT, fontsize=12)
        # ax2.set_ylabel(r'$d_e^2\omega_{pe}$', fontdict=FONT, fontsize=12)
        ax2.set_ylabel(r'$D_{xx}/(v_Ad_i)$', fontdict=FONT, fontsize=12)

        fname = fdir + "dxx_" + str(tindex) + ".pdf"
        fig.savefig(fname)
        plt.close()
        # plt.show()

    fdir = '../data/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'dxxs.dat'
    dxxs.tofile(fname)


def plot_dxx(plot_config):
    """
    """
    run_dir = plot_config["run_dir"]
    run_name = plot_config["run_name"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    vels = np.linspace(0.3, 0.9, 7)
    vels = np.append(vels, 0.95)
    vels = vels[::-1]
    gammas = 1.0 / np.sqrt(1 - vels**2)
    moms = vels * gammas
    vthe = 0.1 * math.sqrt(2)
    gamma0 = 1.0 / math.sqrt(1 - vthe**2)
    p0 = gamma0 * vthe
    moms /= p0
    tframes = np.linspace(10, 20, 11, dtype=int)
    tindices = tframes * plot_config["tinterval"]

    fdir = '../data/' + run_name + '/'
    fname = fdir + 'dxxs.dat'
    dxxs = np.fromfile(fname)
    dxxs = dxxs.reshape((len(tframes), len(vels)))

    dxx_mean = np.mean(dxxs, axis=0)
    dxx_std = np.std(dxxs, axis=0)

    dxx_min = np.min(dxxs, axis=0)
    dxx_max = np.max(dxxs, axis=0)

    fig = plt.figure(figsize=[3.25, 2.5])
    rect = [0.16, 0.18, 0.8, 0.77]
    ax1 = fig.add_axes(rect)

    p1, = ax1.plot(moms, dxx_mean)
    ax1.fill_between(moms, dxx_min, dxx_max, color=p1.get_color(),
                     alpha=0.5)
    ax1.tick_params(labelsize=10)
    ax1.tick_params(axis='x', which='minor', direction='in', top='on')
    ax1.tick_params(axis='x', which='major', direction='in', top='on')
    ax1.tick_params(axis='y', which='minor', direction='in', left='on')
    ax1.tick_params(axis='y', which='major', direction='in')
    ax1.set_xlabel(r'$p/p_0$', fontdict=FONT, fontsize=12)
    # ax1.set_ylabel(r'$d_e^2\omega_{pe}$', fontdict=FONT, fontsize=12)
    ax1.set_ylabel(r'$D_{xx}/(v_Ad_i)$', fontdict=FONT, fontsize=12)

    fdir = '../img/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'dxxs.pdf'
    fig.savefig(fname)

    plt.show()


def plot_diffusion_coefficients(plot_config):
    """
    """
    run_dir = plot_config["run_dir"]
    run_name = plot_config["run_name"]
    tindex = plot_config["tframe"] * plot_config["tinterval"]
    ptl_vel = plot_config["ptl_vel"]
    fname = (run_dir + 'data_' + str(tindex) + '_' +
             str(ptl_vel) + 'c/diff_coeffs.dat')
    fdata = np.genfromtxt(fname, skip_header=1)
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fdir = '../img/' + run_name + '/'
    mkdir_p(fdir)
    def plot_one_diff(index, label1, img_name):
        fig = plt.figure(figsize=[7, 5])
        ax1 = fig.add_axes([xs, ys, w1, h1])
        nptl = fdata[:, -1]
        kappa = fdata[:, index] / nptl
        nframes, = kappa.shape
        tfs = np.arange(nframes)
        tstart1, tend1 = nframes//4, 6*nframes//10
        coeff = np.polyfit(tfs[tstart1:tend1], kappa[tstart1:tend1], 1)
        p = np.poly1d(coeff)
        kappa_fit1 = p(tfs)
        tstart2 = tend1
        coeff = np.polyfit(tfs[tstart2:], kappa[tstart2:], 1)
        p = np.poly1d(coeff)
        kappa_fit2 = p(tfs)
        kappa_fit1 -= kappa_fit1[0]
        kappa_fit2 -= kappa_fit2[0]
        kappa_fit2 += kappa_fit1[tstart2] - kappa_fit2[tstart2]
        kappa_adjust = np.copy(kappa)
        kappa_adjust[:tend1] -= kappa_fit1[:tend1]
        kappa_adjust[tstart2:] -= kappa_fit2[tstart2:]
        # ax1.plot(tfs, kappa, linewidth=2, label=label1, color='k')
        # ax1.plot(tfs, kappa_adjust, linewidth=2, label=label1, color='k')
        # coeff = np.polyfit(tfs[tstart1:], kappa[tstart1:], 3)
        # p = np.poly1d(coeff)
        # kappa_fit1 = p(tfs)
        coeff = np.polyfit(tfs[tstart1:], np.log(kappa[tstart1:]), 1)
        p = np.poly1d(coeff)
        kappa_fit1 = np.exp(p(tfs))
        ax1.plot(tfs, np.log(kappa), linewidth=2, label=label1, color='k')
        # ax1.plot(tfs, kappa_fit1, linewidth=2, label=label1, color='k')
        # ax1.plot(tfs, kappa_fit1 - kappa_fit1[0], linewidth=2,
        #          label=label1, color='k')
        # ax1.plot(tfs, kappa-kappa_fit1, linewidth=2, label=label1, color='k')
        # grad_kappa = np.gradient(kappa)
        # # grad_kappa = gaussian_filter(grad_kappa, 9)
        # grad_grad_kappa = np.gradient(grad_kappa)
        # # grad_grad_kappa = gaussian_filter(grad_grad_kappa, 7)
        # ax1.plot(grad_grad_kappa[2000:], linewidth=2,
        #          label=label1, color='k')
        ax1.tick_params(labelsize=16)
        ax1.tick_params(axis='x', which='minor', direction='in', top='on')
        ax1.tick_params(axis='x', which='major', direction='in', top='on')
        ax1.tick_params(axis='y', which='minor', direction='in', left='on')
        ax1.tick_params(axis='y', which='major', direction='in')
        leg = ax1.legend(loc=2, prop={'size': 20}, ncol=1,
                         shadow=False, fancybox=False, frameon=False)
        fname = fdir + img_name + '.eps'
        fig.savefig(fname)
    # plot_one_diff(0, r'$D_{rr}$', 'drr')
    plot_one_diff(1, r'$D_{xx}$', 'dxx')
    # plot_one_diff(2, r'$D_{yy}$', 'dyy')
    # plot_one_diff(3, r'$D_{zz}$', 'dzz')
    # plot_one_diff(4, r'$D_{pp}$', 'dpp')
    # plot_one_diff(5, r'$D_{\mu\mu}$', 'duu')
    # plot_one_diff(6, r'$D_{aa}$', 'daa')
    # plot_one_diff(7, r'$D_{ee}$', 'dee')
    plt.show()


def get_cmd_args():
    """Get command line arguments
    """
    default_run = '3D-Lx150-bg0.2-150ppc-2048KNL'
    default_run_dir = ('/net/scratch3/xiaocanli/diffusion/' + default_run + '/')
    parser = argparse.ArgumentParser(description='Analysis for test-particle simulations')
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
    parser.add_argument('--ptl_vel', action="store", default='0.3', type=float,
                        help='Particle speed in light speed')
    parser.add_argument('--iptl', action="store", default='1', type=int,
                        help='particle id')
    parser.add_argument('--tint', action="store", default='20', type=int,
                        help='Number of steps before and after current step for piecewise trajectory')
    parser.add_argument('--plot_spectrum', action="store_true", default=False,
                        help="whether to plot particle energy spectrum")
    parser.add_argument('--ptl_traj', action="store_true", default=False,
                        help="whether to plot particle trajectory")
    parser.add_argument('--diff_coeff', action="store_true", default=False,
                        help="whether to plot diffusion coefficients")
    parser.add_argument('--multi_vel', action="store_true", default=False,
                        help="whether to plot for multiple velocities")
    parser.add_argument('--tran_hdf5', action="store_true", default=False,
                        help="whether to transfer trajectory to HDF5 format")
    parser.add_argument('--tran_csv', action="store_true", default=False,
                        help="whether to transfer trajectory to CSV format")
    parser.add_argument('--plot_dxx', action="store_true", default=False,
                        help="whether to plot dxx")
    parser.add_argument('--piecewise_traj', action="store_true", default=False,
                        help="whether to get piecewise trajectory")
    return parser.parse_args()


def analysis_single_frames(plot_config, args):
    """Analysis for multiple time frames
    """
    tframe = args.tframe
    if args.plot_spectrum:
        plot_spectrum(plot_config)
    elif args.ptl_traj:
        plot_particle_trajectory(plot_config)
    elif args.diff_coeff:
        if args.multi_vel:
            diffusion_coefficients_multi(plot_config)
        else:
            plot_diffusion_coefficients(plot_config)
    elif args.tran_hdf5:
        transfer_to_h5part(plot_config)
    elif args.tran_csv:
        transfer_to_csv(plot_config)
    elif args.plot_dxx:
        plot_dxx(plot_config)
    elif args.piecewise_traj:
        piecewise_trajectory(plot_config)


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
    plot_config["tinterval"] = args.tinterval
    plot_config["ptl_vel"] = args.ptl_vel
    plot_config["iptl"] = args.iptl
    plot_config["tint"] = args.tint
    if args.multi_frames:
        analysis_multi_frames(plot_config, args)
    else:
        analysis_single_frames(plot_config, args)


if __name__ == "__main__":
    main()
