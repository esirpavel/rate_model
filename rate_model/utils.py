import matplotlib.pylab as pl
import numpy as np

def plot_simul(net, figsize=(10, 8),*, figure=None, cax=None):
    if figure:
        pass
    else:
        pl.figure(figsize=figsize)
    net_mesh = pl.pcolormesh(net.tm, np.degrees(net.pos), net.ActRE.T)
    pl.xlim((0, net.sim_time))
    pl.ylim((-180, 180))
    pl.xlabel('Time (s)')
    pl.ylabel(r'$\theta$')
    if cax:
        pl.colorbar(net_mesh, cax=cax)
    else:
        pl.colorbar()

def check_diff(phi2, phi1):
    dphi = phi2 - phi1
    if abs(dphi) < np.pi:
        return dphi
    else:
        if phi2 > 0:
            return dphi - 2*np.pi
        else:
            return dphi + 2*np.pi
