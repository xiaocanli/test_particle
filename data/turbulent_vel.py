"""
Plotting procedures for force-free field case.
"""
import h5py
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def PlotVelocityDistribution():
    """
    Plotting the velocity distributions.
    """
    font = {'family' : 'serif',
            #'color'  : 'darkred',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 24,
            }
    # ux
    f = open('ux_distribution.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    vel = data[:,0]
    fux = data[:,1]
    f.close
    # uy
    f = open('uy_distribution.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    fuy = data[:,1]
    f.close
    # uz
    f = open('uz_distribution.dat', 'r')
    data = np.genfromtxt(f, delimiter='')
    fuz = data[:,1]
    f.close
    
    fux_max = np.amax(fux)
    fuy_max = np.amax(fuy)
    fuz_max = np.amax(fuz)
    fmax = max(fux_max, fuy_max, fuz_max)
    fux = fux / fmax
    fuy = fuy / fmax
    fuz = fuz / fmax
    
    dimx,dimy = data.shape
    p1, = plt.plot(vel,fux,linewidth=2, label='$f(u_x)$', color='r')
    p2, = plt.plot(vel,fuy,linewidth=2, label='$f(u_y)$', color='g')
    p3, = plt.plot(vel,fuz,linewidth=2, label='$f(u_z)$', color='b')
    
    plt.xlim([-1,1])
    #plt.yscale('log')
    plt.title('Velocity distribution', fontdict=font)
    plt.xlabel('$u$', fontdict=font)
    plt.ylabel('$f(u)$', fontdict=font)
    plt.tick_params(labelsize=16)
    legend = plt.legend(fontsize=24)
    plt.tight_layout()
    plt.grid(True)
    
    plt.savefig('vel_dist.eps')
    plt.show()

def PlotContourVelocity():
    """
    Contour plot for velocity fields.
    """

    # Open the file, group and dataset using default properties.
    fid = h5py.h5f.open('u4d.h5')
    group = h5py.h5g.open(fid, '/u4d')
    dset = h5py.h5d.open(group, 'ux4d')
    dimsf = dset.shape
    nt, nz, ny, nx = dimsf
    print dimsf
    print dset.dtype
    space = dset.get_space()
    count = (1, 1, ny, nx)
    offset = (1, nz/2, 0, 0)
    stride = (1, 1, 1, 1)
    block = (1, 1, 1, 1)
    space.select_hyperslab(offset, count)
    # Using the data using previously selected hyperslab.
    rdata = np.zeros((1, 1, ny, nx), dtype=np.float64)
    #dset.read(h5py.h5s.ALL, space, rdata) 

    del dset
    del space
    del group
    del fid

if __name__ == "__main__":
    PlotVelocityDistribution()
    #PlotContourVelocity()
