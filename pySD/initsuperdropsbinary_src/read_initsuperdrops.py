import numpy as np
import matplotlib.pyplot as plt

from .create_initsuperdrops import initSDsinputsdict
from ..readbinary import readbinary


def get_superdroplet_attributes(configfile, constsfile, initSDsfile):
    ''' get gridbox boundaries from binary file and 
    re-dimensionalise usign COORD0 const from constsfile '''

    inputs = initSDsinputsdict(configfile, constsfile)
    
    epss, radii, m_sols, coord3s = read_dimless_superdrops_binary(initSDsfile)

    radii = radii * inputs["R0"]
    m_sols = m_sols * inputs["MASS0"]
    coord3s = coord3s * inputs["COORD0"]

    return epss, radii, m_sols, coord3s


def read_dimless_superdrops_binary(filename):
    ''' return dimenionsless gbx boundaries by reading binary file'''

    data, ndata_pervar = readbinary(filename)

    idxs = []
    for n in range(1, len(ndata_pervar)):
        # indexs for division of data list between each variable
        idxs.append(np.sum(ndata_pervar[:n]))

    epss = np.asarray(data[:idxs[0]], dtype=np.uint)
    radii = np.asarray(data[idxs[0]:idxs[1]], dtype=np.double)
    m_sols = np.asarray(data[idxs[1]:idxs[2]], dtype=np.double)
    coord3s = np.asarray(data[idxs[2]:], dtype=np.double)

    print("attribute shapes: ", epss.shape, radii.shape, m_sols.shape, coord3s.shape)
    
    return epss, radii, m_sols, coord3s


def plot_initdistribs(configfile, constsfile, initSDsfile,
                      vol, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    epss, radii, m_sols, coord3s = get_superdroplet_attributes(configfile,
                                                               constsfile,
                                                               initSDsfile)

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))
    axs = axs.flatten()

    # create nbins evenly spaced in log10(r)
    nbins = 100
    minr, maxr = np.min(radii)/10, np.max(radii)*10
    hedgs = np.linspace(np.log10(minr), np.log10(maxr),
                        nbins+1)  # edges to lnr bins

    l0 = plot_radiidistrib(axs[0], hedgs, radii, epss)

    l1 = plot_numconcdistrib(axs[1], hedgs, epss, radii, vol)

    l3 = plot_masssolutedistrib(axs[3], hedgs, epss, radii, m_sols, vol)

    l2 = plot_coord3distrib(axs[2], hedgs, coord3s, radii)

    fig.tight_layout()
    if savefig:
        fig.savefig(binpath+"/initdistribs.png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    plt.show()


def log10r_frequency_distribution(radii, hedgs, wghts):
    ''' get distribution of data with weights 'wghts' against 
    log10(r). Uses np.histogram to get frequency of a particular
    value of data that falls in each bin (with each bin defined
    by it's edges 'hedgs'). Return distirbution alongside the radius
    bin centers and widths in [m]'''

    if type(wghts) != np.ndarray:
        wghts = np.full(np.shape(radii), wghts)

    hist, hedgs = np.histogram(np.log10(radii), bins=hedgs,
                               weights=wghts, density=None)

    # convert [m] to [micron]
    hedgs = (10**(hedgs))*1e6
    # radius bin widths [micron]
    hwdths = hedgs[1:] - hedgs[:-1]
    # radius bin centres [micron]
    hcens = (hedgs[1:]+hedgs[:-1])/2

    return hist, hedgs, hwdths, hcens


def plot_radiidistrib(ax, hedgs, radii, epss):
    ''' get and plotthe superdroplet radii in each log10(r)
    bin and as a scatter on a twinx axis with their multiplicities'''

    l1 = ax.scatter(radii*1e6, epss, zorder=1,
                    color="purple", label="multiplicities")

    ax2 = ax.twinx()
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radii, hedgs, 1)
    l2 = ax2.step(hcens, hist, where='mid', alpha=0.8, zorder=0,
                  color="grey", label="number distribution")

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_yscale("log")

    ax.set_ylabel("superdroplet multiplicity")
    ax2.set_ylabel("superdroplet number distribution")

    ax.legend(loc="lower left")
    ax2.legend(loc="lower right")

    return [l1, l2]


def plot_numconcdistrib(ax, hedgs, epss, radii, vol):
    ''' get and plot frequency of real droplets in each log10(r) bin '''

    wghts = epss / vol / 1e6  # [cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(
        radii, hedgs, wghts)

    line = ax.bar(hcens, hist, hwdths, color="teal",
                  label="binned distribution")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("real droplet number concentration / cm$^{-3}$")
    ax.legend(loc="lower left")

    return line


def plot_masssolutedistrib(ax, hedgs, epss, radii, m_sols, vol):
    ''' get and plot frequency of real droplets in each log10(r) bin '''

    wghts = m_sols*epss/vol * 1000 / 1e6  # [g cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(
        radii, hedgs, wghts)

    line = ax.bar(hcens, hist, hwdths, color="teal")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("solute mass per unit volume / g cm$^{-3}$")

    return line


def plot_coord3distrib(ax, hedgs, coord3s, radii):

    line = None
    if any(coord3s):
        line = ax.scatter(radii*1e6, coord3s, c="purple")

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("superdroplet coord3 / m")

    return line
