"""
Analysis procedures for particle energy spectrum.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import scipy
from scipy.special import ellipk, ellipe
import math
import os.path
import struct
import collections
import palettable
import h5py

colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {'family' : 'serif',
        #'color'  : 'darkred',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24,
        }


def mag_field_loop():
    rho_loop = 1
    len0 = 10.0
    ntot = 256
    xmax = 3
    norm = 0.2
    x = np.linspace(-xmax, xmax, ntot)
    z = np.linspace(-xmax, xmax, ntot)
    xv, zv = np.meshgrid(x, z)
    theta = np.arctan(xv/zv)
    theta[np.where(theta < 0)] += math.pi
    stheta = np.sin(theta)
    ctheta = np.cos(theta)
    r = np.sqrt(xv**2 + zv**2)
    k = np.sqrt(4 * rho_loop * r * stheta /
            (rho_loop**2 + r**2 + 2*rho_loop*r*stheta))
    ellk = scipy.special.ellipk(k)
    elle = scipy.special.ellipe(k)
    rtot2 = rho_loop**2 + r**2
    rho1 = np.sqrt(rtot2 + 2*rho_loop*r*stheta)
    tmp1 = rtot2 / (rtot2 - 2*rho_loop*r*stheta)
    Brho = ctheta * (-ellk + tmp1*elle) / (stheta * rho1)
    Bz = (ellk - (r**2 - rho_loop**2)*elle/(rtot2 - 2*rho_loop*r*stheta)) / rho1
    Btot = np.sqrt(Brho**2 + Bz**2)
    Aphi = 2 * rho_loop * ((2-k**2)*ellk - 2*elle) / rho1 / k**2
    Bz *= norm
    Brho *= norm
    Aphi *= norm

    fig = plt.figure(figsize=[10, 5])
    width = 0.39
    height = 0.78
    xs = 0.07
    ys = 0.95 - height
    ax = fig.add_axes([xs, ys, width, height])
    dmax = 2
    im = ax.imshow(Brho, cmap=plt.cm.seismic, vmin=-dmax, vmax=dmax,
            extent=[-xmax, xmax, -xmax, xmax])
    ax.set_xlabel(r'$\rho/a$', fontdict=font, fontsize=20)
    ax.set_ylabel(r'$z/a$', fontdict=font, fontsize=20)
    ax.tick_params(labelsize=16)
    levels = np.linspace(np.max(Aphi)*0.2, np.min(Aphi), 20)
    # ax.contour(x, z, Aphi, colors='black', linewidths=0.5, levels=levels)
    ax.text(0.95, 0.85, r'$B_\rho$', color='k', fontsize=24,
            bbox=dict(facecolor='none', alpha=1.0,
                      edgecolor='none', pad=10.0),
            horizontalalignment='right', verticalalignment='bottom',
            transform = ax.transAxes)
    xs += width + 0.03
    ax1 = fig.add_axes([xs, ys, width, height])
    im1 = ax1.imshow(Bz, cmap=plt.cm.seismic, vmin=-dmax, vmax=dmax,
            extent=[-xmax, xmax, -xmax, xmax])
    ax1.tick_params(axis='y', labelleft='off')
    ax1.set_xlabel(r'$\rho/a$', fontdict=font, fontsize=20)
    ax1.tick_params(labelsize=16)
    xs += width + 0.01
    cax = fig.add_axes([xs, ys, 0.02, height])
    cbar = fig.colorbar(im1, cax=cax)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.set_ylabel(r'$[I/a]$', fontsize=20)
    levels = np.linspace(np.max(Aphi)*0.2, np.min(Aphi), 20)
    ax1.text(0.95, 0.85, r'$B_z$', color='k', fontsize=24,
            bbox=dict(facecolor='none', alpha=1.0,
                      edgecolor='none', pad=10.0),
            horizontalalignment='right', verticalalignment='bottom',
            transform = ax1.transAxes)
    # ax1.contour(x, z, Aphi, colors='black', linewidths=0.5, levels=levels)
    # strm = ax0.streamplot(xv, zv, Brho, Bz, color=Brho, linewidth=2,
    #         cmap=plt.cm.autumn)
    # fig0.colorbar(strm.lines)

    # fig = plt.figure(figsize=[7, 5])
    # width = 0.69
    # height = 0.8
    # xs = 0.16
    # ys = 0.95 - height
    # ax = fig.add_axes([xs, ys, width, height])
    # ax.plot(Brho[0,:])
    fname = 'img/emf_loop.jpg'
    fig.savefig(fname, dpi=300)

    plt.show()

def emf_wire():
    """Electromagnetic fields of wire current.

    The fields are in unit of omega I / c
    """
    omega = 1.0
    ntot = 256
    rho = np.linspace(0.01, 10, ntot)
    x = omega * rho
    j0 = scipy.special.j0(x)
    j1 = scipy.special.j1(x)
    y0 = scipy.special.y0(x)
    y1 = scipy.special.y1(x)
    t = math.pi / 4
    dt = math.pi / 60
    emf0 = 0.1 * math.pi
    for i in range(30, 31):
        print i
        t = i * dt
        omega_t = omega * t
        Ez = -(y0*math.sin(omega_t) + j0*math.cos(omega_t)) * emf0
        Bphi = (-y1*math.cos(omega_t) + j1*math.sin(omega_t)) * emf0

        fig = plt.figure(figsize=[7, 5])
        width = 0.8
        height = 0.8
        xs = 0.16
        ys = 0.95 - height
        ax = fig.add_axes([xs, ys, width, height])
        ax.plot(rho, Bphi, linewidth=2, color=colors[0], label=r'$B_\phi$')
        ax.plot(rho, Ez, linewidth=2, color=colors[1], label=r'$E_z$')
        ax.plot([0, np.max(rho)], [0,0], color='k', linestyle='--')
        ax.set_ylim([-2, 2])
        ax.set_xlabel(r'$x$', fontdict=font, fontsize=24)
        ax.set_ylabel(r'Fields [$\omega I c^{-1}$]', fontdict=font, fontsize=24)
        ax.tick_params(labelsize=20)
        leg = ax.legend(loc=4, prop={'size':24}, ncol=1,
                shadow=False, fancybox=False, frameon=False)
        tname = r"$\omega t = $" + r"{ %.2f}" % omega_t
        ax.text(0.95, 0.85, tname, color='k', fontsize=24,
                bbox=dict(facecolor='none', alpha=1.0,
                          edgecolor='none', pad=10.0),
                horizontalalignment='right', verticalalignment='bottom',
                transform = ax.transAxes)
        fname = 'img/eb_wire_' + str(i).zfill(3) + '.jpg'
        fig.savefig(fname, dpi=300)
        fname = 'img/eb_wire_' + str(i).zfill(3) + '.eps'
        fig.savefig(fname)
        plt.close()

    # plt.show()


def plot_spectrum():
    """
    Plotting the energy spectrum.
    """
    f = open('../data/espectrum.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    f.close
    dimx,dimy = data.shape
    ene_log = data[0,:];
    fig = plt.figure(figsize=[7, 5])
    width = 0.69
    height = 0.8
    xs = 0.16
    ys = 0.95 - height
    ax = fig.add_axes([xs, ys, width, height])
    for j in range(1, dimx-1):
        color = plt.cm.jet(j/float(dimx), 1)
        ax.loglog(ene_log, data[j,:], color=color, linewidth=2)
    pindex = -1;
    fpower_log = ene_log**pindex;
    ax.loglog(ene_log, fpower_log*1E2, color='k', linestyle='--', linewidth=2)
    
    ax.set_xlim([1E-5, 10]);
    ax.set_ylim([1E1, 1E8]);
    ax.set_xlabel('Energy (MeV)', fontdict=font)
    ax.set_ylabel('Flux (#/MeV)', fontdict=font)
    ax.tick_params(labelsize=16)
    # plt.grid(True)
    # plt.savefig('ene_spectrum.eps')
    
    plt.show()


def plot_diff_coeffs():
    """
    Plotting diffusion coefficients.
    """
    f = open('../data/diff_coeffs.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    f.close
    dimx,dimy = data.shape
    drr = data[:,0]
    dpp = data[:,1]
    duu = data[:,2]
    nptl = data[:,3]

    fig = plt.figure(figsize=[7, 5])
    width = 0.69
    height = 0.8
    xs = 0.16
    ys = 0.95 - height
    ax = fig.add_axes([xs, ys, width, height])
    ax.plot(drr/nptl, color='r', linewidth=2)
    ax.set_ylim([0.02, 0.025])
    
    ax.tick_params(labelsize=16)
    
    plt.show()


def plot_diff_coeffs_multi():
    """
    Plotting diffusion coefficients for multiple runs.
    """
    fig = plt.figure(figsize=[7, 5])
    width = 0.69
    height = 0.8
    xs = 0.16
    ys = 0.95 - height
    ax = fig.add_axes([xs, ys, width, height])
    for i in range(5, 10):
        fname = '../data/diff_coeffs_' + str(i).zfill(2) + '.dat'
        f = open(fname, 'r')
        data = np.genfromtxt(f, delimiter='')
        f.close
        dimx,dimy = data.shape
        drr = data[:,0]
        dpp = data[:,1]
        duu = data[:,2]
        nptl = data[:,3]

        ax.plot(drr/nptl, linewidth=2)
    ax.set_ylim([0.003, 0.009])
    
    ax.tick_params(labelsize=16)
    
    plt.show()



def plot_trajectory(run_name):
    """
    Plotting particle trajectories
    """
    fname = '../data/particle_diagnostics_' + run_name + '.h5'
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    print sz
    ntraj = 100000
    nptl = sz / ntraj
    print("Total number of particles: ", nptl)
    if not os.path.isdir('img_traj/'):
        os.makedirs('img_traj/')
    odir = 'img_traj/' + run_name + '/'
    if not os.path.isdir(odir):
        os.makedirs(odir)
    for iptl in range(nptl):
        px = np.array(dset_ptl['x'][iptl*ntraj+1:(iptl+1)*ntraj-1])
        py = np.array(dset_ptl['y'][iptl*ntraj+1:(iptl+1)*ntraj-1])
        pz = np.array(dset_ptl['z'][iptl*ntraj+1:(iptl+1)*ntraj-1])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # ax = fig.gca()
        ax.plot(px, py, pz, color=colors[1])
        # title = 'Trajectory of particle ' + str(iptl).zfill(3)
        # Title = ax.set_title(title)
        xLabel = ax.set_xlabel('$x$', fontsize=24)
        yLabel = ax.set_ylabel('$y$', fontsize=24)
        zLabel = ax.set_zlabel('$z$', fontsize=24)
        ax.tick_params(labelsize=16)
        plt.tight_layout()
        fname = odir + 'ptl' + str(iptl).zfill(3) + '_traj.png'
        plt.savefig(fname)
        plt.close()
        # plt.show()
    file.close

if __name__ == "__main__":
    mag_field_loop()
    # emf_wire()
    # plot_spectrum()
    # plot_diff_coeffs()
    # plot_diff_coeffs_multi()
    # run_names = ['wire', 'loop_1MK', 'loop_10MK', '100MK', '1wlcs_symmetric',
    #         '1wlcs_IL1_y001', '1wlcs_IL1_y01', '8wlcs']
    # for rname in run_names:
    #     plot_trajectory(rname)
