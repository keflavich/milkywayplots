import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib
from matplotlib.patches import Circle,Wedge
import numpy as np

from  mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot

from mpl_toolkits.axisartist import SubplotHost
try:
    from mpl_toolkits.axisartist import ParasiteAxesAuxTrans
except ImportError:
    from mpl_toolkits.axes_grid1.parasite_axes import ParasiteAxes

import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D


def make_polar_axis(figure):
    """
    Generate a polar axis.

    Examples
    --------
    >>> from pylab import *
    >>> f = figure()
    >>> ax1,ax2 = make_polar_axis(f)
    >>> f.add_subplot(ax1)
    """
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()

    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # 20, 20 : number of sampling points along x, y direction
    extreme_finder = angle_helper.ExtremeFinderCycle(40, 40,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(0, np.inf),
                                                     )

    grid_locator1 = angle_helper.LocatorDMS(12)
    # Find a grid values appropriate for the coordinate (degree,
    # minute, second).

    tick_formatter1 = angle_helper.FormatterDMS()
    # And also uses an appropriate formatter.  Note that,the
    # acceptable Locator and Formatter class is a bit different than
    # that of mpl's, and you cannot directly use mpl's Locator and
    # Formatter here (but may be possible in the future).

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1,
                                        #tick_formatter2=matplotlib.ticker.FuncFormatter(lambda x: x * kpc_per_pix)
                                        )


    ax1 = SubplotHost(figure, 1, 1, 1, grid_helper=grid_helper, axisbg='#333333')

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)

    # let right axis shows ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks=0
    # let bottom axis shows ticklabels for 2nd coordinate (radius)
    ax1.axis["bottom"].get_helper().nth_coord_ticks=1

    try:
        ax2 = ParasiteAxesAuxTrans(ax1, tr, "equal")
    except NameError:
        ax2 = ParasiteAxes(parent_axes=ax1, aux_transform=tr, viewlim_mode='equal')

    ax1.parasites.append(ax2)

    return ax1,ax2

class MilkywayPlot(object):

    def __init__(self, figure, galrad=15.0, earthrad=8.5, npix=1024):

        self.kpc_per_pix = galrad/(npix/2.0)
        self.center = npix/2. - 0.5
        self.earthxpos = earthrad / self.kpc_per_pix
        self.earthrad = earthrad
        self.galrad = galrad
        self.npix = npix

        self.center = self.npix/2. - 0.5

        self.figure = figure

        self.setup_milkyway_plot()

    def setup_milkyway_plot(self, zorder=2):

        ax1,ax2 = make_polar_axis(self.figure)
        self.figure.add_subplot(ax1)

        ax1.set_aspect(1.)
        ax1.set_xlim(-1*self.earthxpos*self.kpc_per_pix,(self.npix-self.earthxpos)*self.kpc_per_pix)
        ax1.set_ylim(-1*(self.npix/2.)*self.kpc_per_pix,(self.npix/2.)*self.kpc_per_pix)
        ax1.set_xlim(-1*self.earthxpos*self.kpc_per_pix,(self.npix-self.earthxpos)*self.kpc_per_pix+5)
        ax1.set_ylim(-1,(self.npix/2.)*self.kpc_per_pix)

        #dropboxdeletedthis OKim = OKzone.astype('int')+OKgal.astype('int')
        #ax1.imshow(OKim[(npix/2.)-250:,:],extent=[(-1*earthxpos)*kpc_per_pix,(npix-earthxpos)*kpc_per_pix,-250*kpc_per_pix,(npix/2.)*kpc_per_pix],cmap=matplotlib.cm.gray)
        #plt.plot(linspace(-pi,pi,1000),
        ax1.add_artist(Circle([self.earthrad,0],self.galrad,facecolor='gray',edgecolor='none',zorder=zorder))
        #ax1.add_artist(Wedge([0,0],maxdist,0,90,facecolor='white',edgecolor='none'))

        self.ax1,self.ax2 = ax1,ax2

    def mark_completeness_zone(self, maxdist=17.5, mindist=5, low_longitude_cutoff=6, zorder=3):
        yy,xx = np.indices([self.npix,self.npix])
        rr = np.sqrt((xx-self.center)**2+(yy-self.center)**2)*self.kpc_per_pix
        theta = np.arctan2(xx - self.center, yy - self.center)

        x1 = 0
        y1 = (self.galrad**2 - self.earthrad**2)**0.5
        #x2 = (galrad**2-distance_cut**2-earthrad**2) / (-2*distance_cut*earthrad) * distance_cut
        #y2 = (maxdist**2 - x2**2)**0.5
        anglecut = low_longitude_cutoff/180.0 * np.pi
        x3 = np.cos(anglecut) * maxdist
        y3 = np.sin(anglecut) * maxdist
        if self.galrad+self.earthrad > maxdist:
            th2 = np.arccos(-1*(self.galrad**2-maxdist**2-self.earthrad**2)/(2*maxdist*self.earthrad))
            x2 = np.cos(th2)*maxdist
            y2 = np.sin(th2)*maxdist
        else:
            # daaaang geometry
            x2 = x3 = (2*self.earthrad + ((2*self.earthrad)**2-4*(np.sin(anglecut)**2+1)*(self.earthrad**2-self.galrad**2))**0.5) / (2*(np.sin(anglecut)**2+1))
            y2 = y3 = np.sin(anglecut)*(x2)
            maxdist = (x3**2+y3**2)**0.5

        x4 = 0
        y4 = mindist
        x5 = np.cos(anglecut) * mindist
        y5 = np.sin(anglecut) * mindist
        xx1 = np.concatenate([np.linspace(0,x2,1000), np.zeros(1000)])
        xx2 = np.concatenate([np.zeros(1000), np.linspace(x2,x3,1000)])
        xx = xx1+xx2
        obs_curve = (self.galrad**2-(xx-self.earthrad)**2)**0.5 *(xx<=x2)  + (maxdist**2-xx**2)**0.5 * (xx>x2)
        low_curve = (mindist**2-(xx*(xx<=x5))**2)**0.5 * (xx<=x5) + xx*np.sin(anglecut) * (xx>x5)
        self.ax1.fill_between(xx,low_curve,obs_curve,facecolor='white',edgecolor='none',zorder=zorder)
        #print "Vertices: " + "\nVertices: ".join(["%f,%f" % (x,y) for (x,y) in zip((x1,x2,x3,x4,x5),(y1,y2,y3,y4,y5))])

        xxclose = np.linspace(0,x5,1000)
        high_curve2 = (mindist**2-(xxclose)**2)**0.5 
        low_curve2  = xxclose*np.sin(anglecut)
        self.ax1.fill_between(xxclose, low_curve2, high_curve2, facecolor='#BBBBBB', edgecolor='none', zorder=zorder+1)

        self.ax1.axis["lon"] = axis = self.ax1.new_floating_axis(1, maxdist)
        axis.axes.set_zorder(20)
        axis.axis.set_zorder(20)
        axis.set_zorder(20)

    def mark_earth(self,zorder=50):
        self.ax2.plot(0,0,'o',color='gold',markersize=10,markeredgecolor='gold',markerfacecolor='none',zorder=zorder,alpha=1,markeredgewidth=2)
        self.ax2.plot(0,0,'.',color='gold',markersize=3,markeredgecolor='gold',markerfacecolor='gold',zorder=zorder,alpha=1,markeredgewidth=2)

    def mark_gc(self,zorder=50):
        self.ax2.plot(0,8.5,'o',color=(0.3,1,0.3,0.5),markersize=30,markeredgecolor=(0.3,1,0.3,0.5),zorder=zorder)
        self.ax2.text(0,8.5,'GC',horizontalalignment='center',verticalalignment='center',zorder=zorder+1)

    def label(self):
        # this results in weird duplication ax1.set_xlabel("Heliocentric Distance (kpc)")
        self.ax1.axis['bottom'].label.set_text("Heliocentric Distance (kpc)")
        self.ax1.set_ylabel("Heliocentric Distance (kpc)")
        # label.set_text http://matplotlib.sourceforge.net/examples/axes_grid/demo_floating_axes.html
        self.ax1.axis['top'].label.set_text("Galactic Longitude")
        self.ax2.axis['top'].label.set_text("Galactic Longitude")
        self.ax1.axis['right'].label.set_text("Galactic Longitude")
        self.ax2.axis['right'].label.set_text("Galactic Longitude")
        self.ax1.axis['top'].label.set_visible(True)
        self.ax2.axis['top'].label.set_visible(True)
        self.ax1.axis['right'].label.set_visible(True)
        self.ax2.axis['right'].label.set_visible(True)


    def add_grid(self):
        self.ax1.grid(True,zorder=11)

    def mark_corotation(self):
        self.ax1.add_artist(Circle([8.5,0],8.5,facecolor='none',edgecolor='r',linestyle='dotted',zorder=5))

    def do_everything(self, **kwargs):
        self.mark_earth()
        self.mark_gc()
        self.label()
        self.add_grid()
        self.mark_corotation()
        self.mark_completeness_zone(**kwargs)

    def plot_tangents(self, npts=100, zorder=105):
        import kdist
        lon = np.linspace(0,90,npts)
        lat = np.zeros(npts)
        velocity = np.ones(npts) * 1000 # just pick a velocity that's too high
        tanpt = kdist.kdist(lon,lat,velocity,silent=True)/1000.

        self.ax2.plot(lon,tanpt,'--',color='m', zorder=zorder)

    def plot_londist(self, lon, dist, zorder=105, **kwargs):
        self.ax2.plot(lon, dist, zorder=zorder, **kwargs)

    def scatter_londist(self, lon, dist, zorder=105, **kwargs):
        self.ax2.scatter(lon, dist, zorder=zorder, **kwargs)
