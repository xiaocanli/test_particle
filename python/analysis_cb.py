"""
Plotting procedures for force-free field case.
"""
import h5py
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def PlotSpectrum():
    """
    Plotting the energy spectrum.
    """
    font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 16,
            }

    f = open('espectrum.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    f.close
    dimx,dimy = data.shape
    for j in range(1, dimy, 19):
        plt.loglog(data[0,:],data[j,:],linewidth=2)
    
    plt.title('Energy spectrum', fontdict=font)
    plt.xlabel('Energy (MeV)', fontdict=font)
    plt.ylabel('Flux (#/MeV)', fontdict=font)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig('ene_spectrum.eps')
    
    plt.show()

def PlotTrajectory():
    """
    Plotting particle trajectories
    """
    file = h5py.File('particle_diagnostics.h5','r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    sz, = dset_ptl.shape
    ntraj = 100000
    nptl = sz / ntraj
    for iptl in range(nptl):
        px = np.array(dset_ptl['x'][iptl*ntraj:(iptl+1)*ntraj-1])
        py = np.array(dset_ptl['y'][iptl*ntraj:(iptl+1)*ntraj-1])
        pz = np.array(dset_ptl['z'][iptl*ntraj:(iptl+1)*ntraj-1])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(px, py, pz)
        title = 'Trajectory of particle ' + str(iptl).zfill(3)
        Title = ax.set_title(title)
        xLabel = ax.set_xlabel('$x$', fontsize=24)
        yLabel = ax.set_ylabel('$y$', fontsize=24)
        zLabel = ax.set_zlabel('$z$', fontsize=24)
        fname = 'image_trajectory/ptl' + str(iptl).zfill(3) + '_traj.png'
        plt.savefig(fname)
        plt.close(fig)
    file.close

def PlotPtlFields():
    """
    Plotting particles information and fields components in
    the save plot.
    """
    file = h5py.File('particle_diagnostics.h5','r')
    group = file['/particles_fields']
    dset_ptl = group['particles']
    dset_emf = group['fields']
    sz, = dset_ptl.shape
    ntraj = 100000
    rest_ene_proton = 938.2724046
    # subplot positions
    xs = 0.15
    xe = 0.8
    yint = 0.2
    rest_ene = rest_ene_proton
    nptl = sz / ntraj
    for iptl in range(nptl):
        ye = 0.95
        px = np.array(dset_ptl['x'][iptl*ntraj:(iptl+1)*ntraj-1])
        py = np.array(dset_ptl['y'][iptl*ntraj:(iptl+1)*ntraj-1])
        pz = np.array(dset_ptl['z'][iptl*ntraj:(iptl+1)*ntraj-1])
        vx = np.array(dset_ptl['vx'][iptl*ntraj:(iptl+1)*ntraj-1])
        vy = np.array(dset_ptl['vy'][iptl*ntraj:(iptl+1)*ntraj-1])
        vz = np.array(dset_ptl['vz'][iptl*ntraj:(iptl+1)*ntraj-1])
        t = np.array(dset_ptl['t'][iptl*ntraj:(iptl+1)*ntraj-1])
        Ex = np.array(dset_emf['Ex'][iptl*ntraj:(iptl+1)*ntraj-1])
        Ey = np.array(dset_emf['Ey'][iptl*ntraj:(iptl+1)*ntraj-1])
        Ez = np.array(dset_emf['Ez'][iptl*ntraj:(iptl+1)*ntraj-1])
        Bx = np.array(dset_emf['Bx'][iptl*ntraj:(iptl+1)*ntraj-1])
        By = np.array(dset_emf['By'][iptl*ntraj:(iptl+1)*ntraj-1])
        Bz = np.array(dset_emf['Bz'][iptl*ntraj:(iptl+1)*ntraj-1])
        beta = np.sqrt(vx*vx+vy*vy+vz*vz)
        gamma = 1.0/np.sqrt(1.0-beta*beta)
        ene = (gamma-1.0)*rest_ene
        fig = plt.figure()
        ax = fig.gca(projection='3d')
#        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
        ax1 = plt.subplot(411)
        x_plt, = ax1.plot(t[0:ntraj:100], px[0:ntraj:100], 'r', label='$x$')
        y_plt, = ax1.plot(t[0:ntraj:100], py[0:ntraj:100], 'g', label='$y$')
        z_plt, = ax1.plot(t[0:ntraj:100], pz[0:ntraj:100], 'b', label='$z$')
        nbins = len(ax1.get_yticklabels())
        ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins/2, prune='upper'))
        ax1.set_position([xs, ye-yint, xe-xs, yint-0.01])
        ax1.legend(handles=[x_plt, y_plt, z_plt], loc='center left', 
                fontsize=16, bbox_to_anchor=(1, 0.5), fancybox=True,
                shadow=False)
        ax1.axes.get_xaxis().set_visible(False)
        plt.tick_params(labelsize=16)

        ax2 = plt.subplot(412)
        vx_plt, = ax2.plot(t[0:ntraj:100], vx[0:ntraj:100], 
                'r', label='$v_x$')
        vy_plt, = ax2.plot(t[0:ntraj:100], vy[0:ntraj:100],
                'g', label='$v_y$')
        vz_plt, = ax2.plot(t[0:ntraj:100], vz[0:ntraj:100],
                'b', label='$v_z$')
        #ax2.plot(t[0:ntraj:100], ene[0:ntraj:100], 'k')
        vmax = max(vx.max(), vy.max(), vz.max())
        vmin = min(vx.min(), vy.min(), vz.min())
        ax2.set_ylim([vmin, vmax])
        nbins = len(ax2.get_yticklabels())
        ye = ye - yint
        ax2.set_position([xs, ye-yint, xe-xs, yint-0.01])
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins/2, prune='upper'))
        ax2.legend(handles=[vx_plt, vy_plt, vz_plt], loc='center left', 
                fontsize=16, bbox_to_anchor=(1, 0.5), fancybox=True,
                shadow=False)
        ax2.axes.get_xaxis().set_visible(False)
        plt.tick_params(labelsize=16)

        ax3 = plt.subplot(413)
        Ex_plt, = ax3.plot(t[0:ntraj:100], Ex[0:ntraj:100], 
                'r', label='$E_x$')
        Ey_plt, = ax3.plot(t[0:ntraj:100], Ey[0:ntraj:100], 
                'g', label='$E_y$')
        Ez_plt, = ax3.plot(t[0:ntraj:100], Ez[0:ntraj:100], 
                'b', label='$E_z$')
        ye = ye - yint
        ax3.set_position([xs, ye-yint, xe-xs, yint-0.01])
        nbins = len(ax3.get_yticklabels())
        ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins/2, prune='upper'))
        ax3.legend(handles=[Ex_plt, Ey_plt, Ez_plt], loc='center left', 
                fontsize=16, bbox_to_anchor=(1, 0.5), fancybox=True,
                shadow=False)
        ax3.axes.get_xaxis().set_visible(False)
        plt.tick_params(labelsize=16)

        ax4 = plt.subplot(414)
        Bx_plt, = ax4.plot(t[0:ntraj:100], Bx[0:ntraj:100], 
                'r', label='$B_x$')
        By_plt, = ax4.plot(t[0:ntraj:100], By[0:ntraj:100], 
                'g', label='$B_y$')
        Bz_plt, = ax4.plot(t[0:ntraj:100], Bz[0:ntraj:100], 
                'b', label='$B_z$')
        ax4.set_xlabel('$t (sec)$', fontsize=24)
        ye = ye - yint
        ax4.set_position([xs, ye-yint, xe-xs, yint-0.01])
        nbins = len(ax4.get_yticklabels())
        ax4.yaxis.set_major_locator(MaxNLocator(nbins=nbins/2, prune='upper'))
        ax4.legend(handles=[Bx_plt, By_plt, Bz_plt], loc='center left', 
                fontsize=16, bbox_to_anchor=(1, 0.5), fancybox=True,
                shadow=False)
        plt.tick_params(labelsize=16)
        fname = 'image_ptl_emf/ptl' + str(iptl).zfill(3) + '_emf.eps'
        plt.savefig(fname)
        plt.close(fig)
        #plt.show()
    file.close

if __name__ == "__main__":
    PlotSpectrum()
