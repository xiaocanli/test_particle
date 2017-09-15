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
from evtk.hl import pointsToVTK 

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
        # tname = r"$\omega t = $" + r"{ %.2f}" % omega_t
        tname = r"$\omega t = \pi/2$"
        ax.text(0.95, 0.85, tname, color='k', fontsize=24,
                bbox=dict(facecolor='none', alpha=1.0,
                          edgecolor='none', pad=10.0),
                horizontalalignment='right', verticalalignment='bottom',
                transform = ax.transAxes)
        fname = 'img/eb_wire_' + str(i).zfill(3) + '.jpg'
        fig.savefig(fname, dpi=300)
        fname = 'img/eb_wire_' + str(i).zfill(3) + '.eps'
        fig.savefig(fname)
        # plt.close()

    plt.show()


def get_nptl(ene, f):
    """Get particle number from energy spectrum.
    """
    nbins, = ene.shape
    emin = ene[0] * 2.0
    emax = ene[-1]
    logemin = math.log10(emin)
    logemax = math.log10(emax)
    logde = (logemax - logemin) / (nbins - 2.0)
    einterval = np.zeros(nbins)
    einterval[0] = emin
    for i in range(1, nbins - 1):
        einterval[i] = 10**logemin * (10**(logde*i)-10**(logde*(i-1)))
    einterval[nbins-1] = 10**logemin * (10**(logde*nbins) -
            10**(logde*(nbins-1)))

    return np.sum(f*einterval)


def plot_spectrum(**kwargs):
    """
    Plotting the energy spectrum.
    """
    species = kwargs['species'] if 'species' in kwargs else 'proton'
    normI = kwargs['normI'] if 'normI' in kwargs else '025'
    omega = kwargs['omega'] if 'omega' in kwargs else '0001'
    fname_in = '../data/espectrum_' + species + '_' + normI + \
               'I0_' + omega + 'Hz.dat'
    fname_out = '../data/espect_escape_' + species + '_' + normI + \
               'I0_' + omega + 'Hz.dat'
    fin = open(fname_in, 'r')
    fout = open(fname_out, 'r')
    data_in = np.genfromtxt(fin, delimiter='')
    data_out = np.genfromtxt(fout, delimiter='')
    fin.close()
    fout.close()
    dimx, dimy = data_in.shape
    ene_log = data_in[0,:]
    ns = kwargs['nstep'] if 'nstep' in kwargs else 20
    nt = kwargs['nt'] if 'nt' in kwargs else dimx - 1
    tseries = np.zeros(ns+1)
    dt = nt**(1.0/ns)
    tseries[0] = 1
    for i in range(1, ns+1):
        tseries[i] = tseries[i-1] * dt
    tseries.astype(int)
    fig = plt.figure(figsize=[7, 5])
    width = 0.73
    height = 0.8
    xs = 0.13
    ys = 0.95 - height
    ax = fig.add_axes([xs, ys, width, height])
    norm = mpl.colors.Normalize(vmin=0, vmax=ns)
    c_m = mpl.cm.rainbow
    s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    for i in range(ns+1):
        ct = int(tseries[i])
        color = c_m(i/float(ns+1), 1)
        ax.loglog(ene_log, data_in[ct,:], color=color, linewidth=2)
    n0 = get_nptl(ene_log, data_in[1, :])
    n1 = get_nptl(ene_log, data_in[-1, :])
    print "Fraction of particles remains: ", n1 / n0
    fname_in = '../data/diff_coeffs_' + species + '_' + normI + \
               'I0_' + omega + 'Hz.dat'
    data = np.genfromtxt(fname_in, delimiter='')
    nptl = data[:, -1]
    fraction = r'%.1f' % (nptl[-1] * 100.0 / nptl[0] - 0.1)
    espect_out = np.sum(data_out[1:,:], axis=0)
    ax.loglog(ene_log, espect_out, color='k',
            linestyle='-', linewidth=3, label='Escape spectrum')
    xs += width + 0.02
    cax = fig.add_axes([xs, ys, 0.03, height])
    cbar = plt.colorbar(s_m, cax=cax)
    cbar.ax.tick_params(labelsize=16)
    lname = str(ns) + r'$\log t / \log t_\text{max}$'
    cbar.ax.set_ylabel(lname, fontsize=20)
    pindex = kwargs['pindex'] if 'pindex' in kwargs else -0.8
    pnorm = kwargs['pnorm'] if 'pnorm' in kwargs else 1E7
    fpower_log = ene_log**pindex
    ax.loglog(ene_log, fpower_log*pnorm, color='k', linestyle='--', linewidth=2)
    text_pos = kwargs['text_pos'] if 'text_pos' in kwargs else [0.95, 0.75]
    tname = r'$\sim E^{' + str(pindex) + '}$'
    ax.text(text_pos[0], text_pos[1], tname, color='k', fontsize=20,
            bbox=dict(facecolor='none', alpha=1.0,
                      edgecolor='none', pad=10.0),
            horizontalalignment='right', verticalalignment='bottom',
            transform = ax.transAxes)
    tname1 = r'$F_\text{in} = ' + fraction + '\%$'
    ax.text(0.98, 0.9, tname1, color='k', fontsize=20,
            bbox=dict(facecolor='none', alpha=1.0,
                      edgecolor='none', pad=10.0),
            horizontalalignment='right', verticalalignment='bottom',
            transform = ax.transAxes)
    leg = ax.legend(loc=3, prop={'size':20}, ncol=1,
            shadow=False, fancybox=False, frameon=False)

    xlims = kwargs['xlims'] if 'xlims' in kwargs else [1E-5, 10]
    ylims = kwargs['ylims'] if 'ylims' in kwargs else [1E1, 1E10]
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel('$E$ [MeV]', fontdict=font, fontsize=20)
    ax.set_ylabel('$f(E)$ [MeV$^{-1}$]', fontdict=font, fontsize=20)
    ax.tick_params(labelsize=16)
    fname = 'espect_' + species + '_' + normI + '_' + omega + '.eps'
    plt.savefig('img_spectrum/' + fname)

    plt.close()
    # plt.show()


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
    # fname = '../data/particle_diagnostics_' + run_name + '.h5'
    fname = '../data/particle_diagnostics.h5'
    file = h5py.File(fname,'r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    # ntraj = 12500
    # nptl = sz / ntraj
    nptl = 128
    ntraj = sz / nptl
    print("Total number of particles: ", nptl)
    if not os.path.isdir('img_traj/'):
        os.makedirs('img_traj/')
    odir = 'img_traj/' + run_name + '/'
    if not os.path.isdir(odir):
        os.makedirs(odir)
    for iptl in range(10):
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
        # fname = odir + 'ptl' + str(iptl).zfill(3) + '_traj.png'
        # plt.savefig(fname)
        # plt.close()
        plt.show()
    file.close


def plot_spectrum_multi():
    """Plot multiple particle energy spectrum
    """
    plot_spectrum()
    kwargs = {'species': 'proton', 'normI': '125', 'omega': '0001',
              'pindex': -0.8, 'pnorm': 1E7}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25', 'omega': '0001',
            'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 20]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '250', 'omega': '0001',
              'pindex': -0.55, 'pnorm': 1E7, 'xlims': [1E-5, 30]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '2500', 'omega': '0001',
              'pindex': -0.55, 'pnorm': 1E7, 'xlims': [1E-5, 30]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25000', 'omega': '0001',
              'pindex': -0.55, 'pnorm': 1E7, 'xlims': [1E-5, 30]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '025', 'omega': '001',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 1E2],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '125', 'omega': '001',
              'pindex': -0.5, 'pnorm': 1E7, 'xlims': [1E-5, 1E2],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25', 'omega': '001',
              'pindex': -0.5, 'pnorm': 1E7, 'xlims': [1E-5, 100],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '250', 'omega': '001',
              'pindex': -0.5, 'pnorm': 1E7, 'xlims': [1E-5, 200],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '2500', 'omega': '001',
              'pindex': -0.5, 'pnorm': 1E7, 'xlims': [1E-5, 300],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25000', 'omega': '001',
              'pindex': -0.5, 'pnorm': 1E7, 'xlims': [1E-5, 100],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '025', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 1E2],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '125', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 1E2],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E6, 'xlims': [1E-5, 100],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '250', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 200],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '2500', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 300],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'proton', 'normI': '25000', 'omega': '01',
              'pindex': -0.7, 'pnorm': 1E7, 'xlims': [1E-5, 100],
              'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '025', 'omega': '0001',
              'pindex': -0.75, 'pnorm': 1E5, 'xlims': [1E-5, 3E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '125', 'omega': '0001',
              'pindex': -0.7, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '25', 'omega': '0001',
              'pindex': -0.7, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '250', 'omega': '0001',
              'pindex': -0.6, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '2500', 'omega': '0001',
              'pindex': -0.6, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '025', 'omega': '001',
              'pindex': -0.75, 'pnorm': 1E5, 'xlims': [1E-5, 3E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '125', 'omega': '001',
              'pindex': -0.7, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)
    kwargs = {'species': 'electron', 'normI': '25', 'omega': '001',
              'pindex': -0.7, 'pnorm': 2E5, 'xlims': [1E-5, 5E0],
              'ylims': [10, 2E8], 'text_pos': [0.95, 0.65]}
    plot_spectrum(**kwargs)


def particle_trajectory_h5part():
    """Transfer particle trajectory to H5Part format
    """
    fname = '../data/particle_diagnostics.h5'
    fh = h5py.File(fname,'r')
    group = fh['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    nptl = 128
    ntraj = sz / nptl
    print("Total number of particles: ", nptl)
    x = np.array(dset_ptl['x'])
    y = np.array(dset_ptl['y'])
    z = np.array(dset_ptl['z'])
    q = np.arange(nptl)
    fname_out = '../data/particle_diagnostics.h5part'
    with h5py.File(fname_out, 'w') as fh_out:
        stride = 81
        for ct in range(1, ntraj, stride):
            print ct
            grp = fh_out.create_group('Step#'+str((ct-1)/stride))
            px = x[ct::ntraj]
            py = y[ct::ntraj]
            pz = z[ct::ntraj]
            grp.create_dataset('q', (nptl, ), data=q)
            grp.create_dataset('X', (nptl, ), data=px)
            grp.create_dataset('Y', (nptl, ), data=py)
            grp.create_dataset('Z', (nptl, ), data=pz)
    fh.close


def particle_trajectory_vtk():
    """Transfer particle trajectory to vtk format
    """
    fname = '../data/particle_diagnostics.h5'
    fh = h5py.File(fname,'r')
    group = fh['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    nptl = 128
    ntraj = sz / nptl
    print("Total number of particles: ", nptl)
    x = np.array(dset_ptl['x'])
    y = np.array(dset_ptl['y'])
    z = np.array(dset_ptl['z'])
    stride = 9
    q = np.arange(ntraj / stride)
    fname_out = '../data/particle_diagnostics.vtu'
    for iptl in range(nptl):
        px = np.array(x[iptl*ntraj+1:(iptl+1)*ntraj-1:stride])
        py = np.array(y[iptl*ntraj+1:(iptl+1)*ntraj-1:stride])
        pz = np.array(z[iptl*ntraj+1:(iptl+1)*ntraj-1:stride])
        fname = "./points" + str(iptl)
        pointsToVTK(fname, px, py, pz, data = {"q" : q})
    fh.close


if __name__ == "__main__":
    # mag_field_loop()
    # emf_wire()
    # plot_spectrum_multi()
    # plot_diff_coeffs()
    # plot_diff_coeffs_multi()
    # run_names = ['wire', 'loop_1MK', 'loop_10MK', '100MK', '1wlcs_symmetric',
    #         '1wlcs_IL1_y001', '1wlcs_IL1_y01', '8wlcs']
    # for rname in run_names:
    # plot_trajectory('')
    particle_trajectory_h5part()
    # particle_trajectory_vtk()
