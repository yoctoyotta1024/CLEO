import numpy as np
import matplotlib.pyplot as plt

from .create_gbxboundaries import get_COORD0_from_constsfile
from ..readbinary import readbinary


def get_gridboxboundaries(gridfile, COORD0=False, constsfile=""):
    ''' get gridbox boundaries from binary file and 
    re-dimensionalise usign COORD0 const from constsfile '''

    if not COORD0:
        COORD0 = get_COORD0_from_constsfile(constsfile)
    
    gbxbounds =  read_dimless_gbxboundaries_binary(gridfile, COORD0) 

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds)

    return zhalf, xhalf, yhalf

def get_domainvol_from_gridfile(gridfile, COORD0=False, constsfile=""):
    ''' get total domain volume from binary file '''

    zhalf, xhalf, yhalf = get_gridboxboundaries(gridfile, COORD0=COORD0,
                                               constsfile=constsfile)
    
    return calc_domainvol(zhalf, xhalf, yhalf)

def read_dimless_gbxboundaries_binary(filename, COORD0=False):
    ''' return dictionary for gbx indicies to gbx boundaries by
    reading binary file. Return dimensionless version if COORD0
    not give (=False). '''

    data, ndata_pervar = readbinary(filename)
    datatypes = [np.uintc, np.double]
    
    idxs = np.asarray(data[:ndata_pervar[0]], dtype=datatypes[0])
    ngridboxes = len(idxs)
    boundsdata = np.asarray(data[ndata_pervar[0]:], dtype=datatypes[1])
    boundsdata = np.reshape(boundsdata, [ngridboxes, len(boundsdata)//ngridboxes])
    
    if COORD0:
        boundsdata = boundsdata * COORD0
    
    gbxbounds = {idxs[i]: boundsdata[i] for i in range(ngridboxes)}
    
    return gbxbounds

def halfcoords_from_gbxbounds(gbxbounds):
    ''' returns half coords of gbx boundaries in lists obtained
     from gbxbounds dictionary '''

    boundsdata = np.asarray(list(gbxbounds.values()))
    
    zhalf = np.unique(np.sort(boundsdata[:,0]))
    zhalf = np.append(zhalf, np.amax(boundsdata[:,1]))

    xhalf = np.unique(np.sort(boundsdata[:,2]))
    xhalf = np.append(xhalf, np.amax(boundsdata[:,3]))

    yhalf = np.unique(np.sort(boundsdata[:,4]))
    yhalf = np.append(yhalf, np.amax(boundsdata[:,5]))

    print("zhalf: ", zhalf)
    print("xhalf: ", xhalf)
    print("yhalf: ", yhalf)

    return zhalf, xhalf, yhalf

def plot_gridboxboundaries(constsfile, gridfile, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    zhalf, xhalf, yhalf = get_gridboxboundaries(gridfile, 
                                                constsfile=constsfile)

    halfs = [zhalf, xhalf, yhalf]
    deltas = []
    fulls = []
    for half in halfs:
        full, delta = get_fullcell_and_cellspacing(half)
        deltas.append(delta)
        fulls.append(full)

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))

    for i, crd in enumerate(["z", "x", "y"]):
        xlims = [0.8*np.amin(deltas[i])/1000, np.amax(deltas[i])/1000*1.2]
        ylims = [np.amin(halfs[i])/1000, np.amax(halfs[i])/1000]
        axs[i].scatter(deltas[i]/1000, fulls[i]/1000,
                       color="k", label="centres")
        axs[i].hlines(halfs[i]/1000, xlims[0], xlims[1],
                      color="grey", alpha=0.8, linewidth=0.8, label="boundaries")
        axs[i].set_xlim(xlims)
        axs[i].set_ylim(ylims)
        axs[i].set_xlabel("gridbox spacing, \u0394 "+crd+" /km")
        axs[i].set_ylabel("gridbox centres, "+crd+"f /km")
        axs[i].legend()

    fig.tight_layout()
    if savefig:
        fig.savefig(binpath+"/gridboxboundaries.png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    plt.show()


def get_fullcell_and_cellspacing(halfcell):

    fullcell = (halfcell[1:]+halfcell[:-1])/2
    cellwidth = abs(halfcell[1:]-halfcell[:-1])

    return fullcell, cellwidth


def calc_domainvol(zhalf, xhalf, yhalf):

    widths = []
    for half in [zhalf, xhalf, yhalf]:
        widths.append(np.amax(half) - np.amin(half))

    domainvol = np.prod(widths)

    return domainvol


def calc_gridboxvols(zhalf, xhalf, yhalf):

    widths = []
    for half in [xhalf, yhalf]:
        widths.append(np.amax(half) - np.amin(half))

    area = np.prod(widths)
    zwidths = abs(zhalf[1:] - zhalf[:-1])
    gridboxvols = zwidths * area

    return gridboxvols


def calc_domaininfo(zhalf, xhalf, yhalf):

    domainvol = calc_domainvol(zhalf, xhalf, yhalf)

    gridboxvols = calc_gridboxvols(zhalf, xhalf, yhalf)
    num_gridboxes = len(gridboxvols)

    if num_gridboxes == 1:
        SDnspace = 0
    else:
        SDnspace = 1

    return domainvol, gridboxvols, num_gridboxes, SDnspace


def print_domain_info(constsfile, gridfile):
    ''' create values from constants file & config file
    required as inputs to create initial 
    superdroplet conditions '''

    zhalf, xhalf, yhalf = get_gridboxboundaries(gridfile,
                                                constsfile=constsfile)

    domainvol, gridboxvols, num_gridboxes, SDnspace = calc_domaininfo(
        zhalf, xhalf, yhalf)

    zwdths = abs(zhalf[1:]-zhalf[:-1])
    xwdths = abs(xhalf[1:]-xhalf[:-1])
    ywdths = abs(yhalf[1:]-yhalf[:-1])

    ztot = abs(np.amax(zhalf) - np.amin(zhalf))
    xtot = abs(np.amax(xhalf) - np.amin(xhalf))
    ytot = abs(np.amax(yhalf) - np.amin(yhalf))

    print("\n------ DOMAIN / GRIDBOXES INFO ------")
    print("------------ "+str(SDnspace)+"-D MODEL ------------")
    print("domain dimensions: ({:3g}x{:3g}x{:3g})m^3".format(ztot, xtot, ytot))
    print("domain no. gridboxes: "+str(len(zwdths)) +
          "x"+str(len(xwdths))+"x"+str(len(ywdths)))
    print("domain (upper,lower) z limits: ({:3g},{:3g})m".format(
        np.amax(zhalf), np.amin(zhalf)))
    print("domain (upper,lower) x limits: ({:3g},{:3g})m".format(
        np.amax(xhalf), np.amin(xhalf)))
    print("domain (upper,lower) y limits: ({:3g},{:3g})m".format(
        np.amax(yhalf), np.amin(yhalf)))
    print("avg gridbox z spacing: {:3g} m".format(np.mean(zwdths)))
    print("avg gridbox x spacing: {:3g} m".format(np.mean(xwdths)))
    print("avg gridbox y spacing: {:3g} m".format(np.mean(ywdths)))
    print("avg gridbox volume: {:3g}".format(np.mean(gridboxvols))+" m^3")
    print("total domain volume: {:3g} m^3".format(domainvol))
    print("total no. gridboxes:", num_gridboxes)
    print("------------------------------------\n")
