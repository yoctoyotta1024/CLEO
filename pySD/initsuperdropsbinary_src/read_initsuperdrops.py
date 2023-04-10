import numpy as np
import matplotlib.pyplot as plt

from .create_initsuperdrops import initSDsinputsdict, ManyInitAttrs
from ..readbinary import readbinary
from ..gbxboundariesbinary_src.read_gbxboundaries import get_gbxvols_from_gridfile

def get_superdroplet_attributes(configfile, constsfile, initSDsfile):
    ''' get gridbox boundaries from binary file and 
    re-dimensionalise usign COORD0 const from constsfile '''

    inputs = initSDsinputsdict(configfile, constsfile)
    
    attrs = read_dimless_superdrops_binary(initSDsfile)

    # re-dimensionalise SD attributes
    attrs.radius = attrs.radius * inputs["R0"]
    attrs.m_sol = attrs.m_sol * inputs["MASS0"]
    attrs.coord3 = attrs.coord3 * inputs["COORD0"]
    attrs.coord1 = attrs.coord1 * inputs["COORD0"]
    attrs.coord2 = attrs.coord2 * inputs["COORD0"]

    return attrs


def read_dimless_superdrops_binary(filename):
    ''' return dimenionsless gbx boundaries by reading binary file'''

    datatypes = [np.uintc, np.uint, np.double, np.double]
    datatypes += [np.double]*3
    data, ndata_pervar = readbinary(filename)

    idxs = [0,0,0,0,0,0] # indexs for division of data list between each variable
    for n in range(1, len(ndata_pervar)):
        idxs[n-1] = np.sum(ndata_pervar[:n])

    attrs = ManyInitAttrs()
    attrs.sd_gbxindex = np.asarray(data[:idxs[0]], dtype=datatypes[0])
    attrs.eps = np.asarray(data[idxs[0]:idxs[1]], dtype=datatypes[1])
    attrs.radius = np.asarray(data[idxs[1]:idxs[2]], dtype=datatypes[2])
    attrs.m_sol = np.asarray(data[idxs[2]:idxs[3]], dtype=datatypes[3])
    attrs.coord3 = np.asarray(data[idxs[3]:idxs[4]], dtype=datatypes[4])
    attrs.coord1 = np.asarray(data[idxs[4]:idxs[5]], dtype=datatypes[5])
    attrs.coord2 = np.asarray(data[idxs[5]:], dtype=datatypes[6])

    print("attribute shapes: ", attrs.sd_gbxindex.shape, attrs.eps.shape,
          attrs.radius.shape, attrs.m_sol.shape, attrs.coord3.shape,
          attrs.coord1.shape, attrs.coord2.shape)
    
    return attrs


def plot_initdistribs(configfile, constsfile, initSDsfile,
                      gridfile, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    gbxvols = get_gbxvols_from_gridfile(gridfile, constsfile=constsfile)
    attrs = get_superdroplet_attributes(configfile,constsfile, initSDsfile)

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))
    axs = axs.flatten()

    # create nbins evenly spaced in log10(r)
    nbins = 100
    minr, maxr = np.min(attrs.radius)/10, np.max(attrs.radius)*10
    hedgs = np.linspace(np.log10(minr), np.log10(maxr),
                        nbins+1)  # edges to lnr bins

    unique_idxs = np.unique(attrs.sd_gbxindex)
    for j, idx in enumerate(unique_idxs):
        vol = gbxvols[j]
        i2plt = np.where(attrs.sd_gbxindex == idx)
        l0 = plot_radiusdistrib(axs[0], hedgs, 
                                attrs.radius[i2plt], attrs.eps[i2plt])

        l1 = plot_numconcdistrib(axs[1], hedgs, attrs.eps[i2plt],
                                 attrs.radius[i2plt], vol)

        l3 = plot_masssolutedistrib(axs[3], hedgs, attrs.eps[i2plt],
                                    attrs.radius[i2plt], attrs.m_sol[i2plt],
                                    vol)
        
        if attrs.coord3 != []:
            l2 = plot_coord3distrib(axs[2], hedgs, attrs.coord3[i2plt],
                                    attrs.radius[i2plt])

    fig.tight_layout()
    if savefig:
        fig.savefig(binpath+"/initdistribs.png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    plt.show()


def log10r_frequency_distribution(radius, hedgs, wghts):
    ''' get distribution of data with weights 'wghts' against 
    log10(r). Uses np.histogram to get frequency of a particular
    value of data that falls in each bin (with each bin defined
    by it's edges 'hedgs'). Return distirbution alongside the radius
    bin centers and widths in [m]'''

    if type(wghts) != np.ndarray:
        wghts = np.full(np.shape(radius), wghts)

    hist, hedgs = np.histogram(np.log10(radius), bins=hedgs,
                               weights=wghts, density=None)

    # convert [m] to [micron]
    hedgs = (10**(hedgs))*1e6
    # radius bin widths [micron]
    hwdths = hedgs[1:] - hedgs[:-1]
    # radius bin centres [micron]
    hcens = (hedgs[1:]+hedgs[:-1])/2

    return hist, hedgs, hwdths, hcens


def plot_radiusdistrib(ax, hedgs, radius, eps):
    ''' get and plotthe superdroplet radius in each log10(r)
    bin and as a scatter on a twinx axis with their multiplicities'''

    l1 = ax.scatter(radius*1e6, eps, zorder=1,
                    label="multiplicities")

    ax2 = ax.twinx()
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radius, hedgs, 1)
    l2 = ax2.step(hcens, hist, where='mid', alpha=0.8, zorder=0,
                  color="grey", label="number distribution")

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_yscale("log")

    ax.set_ylabel("superdroplet multiplicity")
    ax2.set_ylabel("superdroplet number distribution")

    if not ax.get_legend():
        ax.legend(loc="lower left")
        ax2.legend(loc="lower right")

    return [l1, l2]


def plot_numconcdistrib(ax, hedgs, eps, radius, vol):
    ''' get and plot frequency of real droplets in each log10(r) bin '''

    wghts = eps / vol / 1e6  # [cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(
                                            radius, hedgs, wghts)

    line = ax.step(hcens, hist, label="binned distribution", where='mid')
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("real droplet number concentration / cm$^{-3}$")
    
    if not ax.get_legend():
        ax.legend(loc="lower left")

    return line


def plot_masssolutedistrib(ax, hedgs, eps, radius, m_sol, vol):
    ''' get and plot frequency of real droplets in each log10(r) bin '''

    wghts = m_sol*eps/vol * 1000 / 1e6  # [g cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(
        radius, hedgs, wghts)

    line = ax.step(hcens, hist, where='mid')
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("solute mass per unit volume / g cm$^{-3}$")

    return line


def plot_coord3distrib(ax, hedgs, coord3, radius):

    line = None
    if any(coord3):
        line = ax.scatter(radius*1e6, coord3)

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("superdroplet coord3 / m")

    return line