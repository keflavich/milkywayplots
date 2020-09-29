import numpy as np
import PIL

import matplotlib.pyplot as plt

from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from mpl_toolkits.axisartist import SubplotHost, ParasiteAxesAuxTrans

import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D

# annotated version
# http://upload.wikimedia.org/wikipedia/commons/8/89/236084main_MilkyWay-full-annotated.jpg

import astropy.utils.data as aud
import os

def getfile(url,fn=None):
    if fn is None:
        fn = os.path.split(url)[-1]
    if not os.path.exists(fn):
        import requests
        response = requests.get(url)
        response.raise_for_status()

        with open(fn,'wb') as of:
            of.write(response.content)


def get_image(mw_img_url="http://upload.wikimedia.org/wikipedia/commons/0/09/Milky_Way_2005.jpg"):
    """
    Download an image - default is non-annotated Robert Hurt image

    Annotated version here:
    http://upload.wikimedia.org/wikipedia/commons/8/89/236084main_MilkyWay-full-annotated.jpg
    """
    getfile(mw_img_url)

def make_mw_plot(fig=None, mw_img_name="Milky_Way_2005.jpg", solar_rad=8.5,
                 fignum=5):
    """
    Generate a "Milky Way" plot with Robert Hurt's Milky Way illustration as
    the background.

    .. TODO:
        Figure out how to fix the axis labels.  They don't work now!

    Parameters
    ----------
    fig : matplotlib.figure instance
        If you want to start with a figure instance, can specify it
    mw_img_name: str
        The name of the image on disk
    solar_rad : float
        The assumed Galactocentric orbital radius of the sun
    fignum : int
        If Figure not specified, use this figure number
    """

    # load image
    mw = np.array(PIL.Image.open(mw_img_name))[:,::-1]

    # set some constants
    npix = mw.shape[0] # must be symmetric
    # Galactic Center in middle of image
    gc_loc = [x/2 for x in mw.shape]

    # Sun is at 0.691 (maybe really 0.7?) length of image
    sun_loc = mw.shape[0]/2,int(mw.shape[1]*0.691)
    # determine scaling
    kpc_per_pix = solar_rad / (sun_loc[1]-gc_loc[1])
    boxsize = npix*kpc_per_pix

    # most of the code below is taken from:
    # http://matplotlib.sourceforge.net/examples/axes_grid/demo_curvelinear_grid.html
    # and http://matplotlib.sourceforge.net/examples/axes_grid/demo_floating_axis.html

    if fig is None:
        fig = plt.figure(fignum)
    plt.clf()

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    # this defines the polar coordinate system @ Galactic center
    tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()

    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # grid helper stuff, I think (grid is off by default)
    # This may not apply to the image *at all*, but would if you
    # used the grid
    # 40, 40 : number of sampling points along x, y direction
    extreme_finder = angle_helper.ExtremeFinderCycle(40, 40,
                                                     lon_cycle = 360,
                                                     lat_cycle = None,
                                                     lon_minmax = None,
                                                     lat_minmax = (0, np.inf),
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


    ax = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)
    fig.add_subplot(ax)
    # ax.transData is still a (rectlinear) pixel coordinate. Only the
    # grids are done in galactocentric coordinate.

    # show the image
    ax.imshow(mw,extent=[-boxsize/2,boxsize/2,-boxsize/2,boxsize/2])

    ax_pixgrid = ax.twin() # simple twin will give you a twin axes,
                           # but with normal grids.

    # to draw heliocentric grids, it is best to update the grid_helper
    # with new transform.

    # need to rotate by -90 deg to get into the standard convention
    tr_helio = Affine2D().scale(np.pi/180., 1.).translate(-np.pi/2.,0) + \
               PolarAxes.PolarTransform() + \
               Affine2D().translate(0,solar_rad)
    # Note that the transform is from the heliocentric coordinate to
    # the pixel coordinate of ax (i.e., ax.transData).

    ax.get_grid_helper().update_grid_finder(aux_trans=tr_helio)

    # Now we defina parasite axes with galactocentric & heliocentric
    # coordinates.

    # A parasite axes with given transform
    gc_polar = ParasiteAxesAuxTrans(ax, tr, "equal")
    ax.parasites.append(gc_polar)
    # note that ax2.transData == tr + galactocentric_axis.transData
    # Anthing you draw in ax2 will match the ticks and grids of galactocentric_axis.

    hc_polar = ParasiteAxesAuxTrans(ax, tr_helio, "equal")
    ax.parasites.append(hc_polar)


    return ax, ax_pixgrid, gc_polar, hc_polar, tr, tr_helio

if __name__=="__main__":
    get_image()

    ax, ax_pixgrid, gcp, hcp, tr_gal, tr_helio = make_mw_plot()
    ax.grid(color="w")
    ax_pixgrid.grid(color="r")
    ax_pixgrid.axis[:].major_ticklabels.set_color("r")
    ax.set_autoscale_on(False)

    # hcp.grid(color='white')
    hcp.plot(0,0,'o',color='gold',markersize=10,markeredgecolor='gold',markerfacecolor='none',zorder=50,alpha=1,markeredgewidth=2)
    hcp.plot(0,0,'.',color='gold',markersize=3,markeredgecolor='gold',markerfacecolor='gold',zorder=50,alpha=1,markeredgewidth=2)
    ax.plot(0,0,'o',color=(0.3,1,0.3,0.5),markersize=30,markeredgecolor=(0.3,1,0.3,0.5))
    ax.text(0,0,'GC',horizontalalignment='center',verticalalignment='center')

    gcp.plot(-30,3,'s')
    hcp.plot(20,3,'^')
    hcp.plot(45,6,'h')
    hcp.plot(-45,6,'o')
    hcp.plot(np.linspace(0,90),np.ones(50)*3)
    hcp.plot(np.linspace(90,180),np.ones(50)*4)

    plt.show()
