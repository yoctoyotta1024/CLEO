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

def get_gbxvols_from_gridfile(gridfile, COORD0=False, constsfile=""):
    ''' get total domain volume from binary file '''

    if not COORD0:
        COORD0 = get_COORD0_from_constsfile(constsfile)
    
    gbxbounds =  read_dimless_gbxboundaries_binary(gridfile, COORD0) 

    return calc_gridboxvols(gbxbounds)

def fullcell_fromhalfcoords(zhalf, xhalf, yhalf):

    zfull = get_fullcell_and_cellspacing(zhalf)[0] 
    xfull = get_fullcell_and_cellspacing(xhalf)[0] 
    yfull = get_fullcell_and_cellspacing(yhalf)[0] 
    
    return zfull, xfull, yfull

def fullcoords_forallgridboxes(gbxbounds, ndims):

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds)
    zfull, xfull, yfull = fullcell_fromhalfcoords(zhalf, xhalf, yhalf)

    zfullcoords = np.tile(zfull, int(ndims[1]*ndims[2])) # zfull of every gridbox in order of gbxindex
    xfullcoords = np.tile(np.repeat(xfull, ndims[0]), int(ndims[2]))
    yfullcoords = np.repeat(yfull, ndims[0]*ndims[1])

    return zfullcoords, xfullcoords, yfullcoords

def halfcoords_forallgridboxes(gbxbounds, ndims):

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds)

    zhalfcoords = np.repeat(zhalf, 2)[1:-1]
    xhalfcoords = np.repeat(xhalf, 2)[1:-1]
    yhalfcoords = np.repeat(yhalf, 2)[1:-1]

    zhalfcoords = np.tile(zhalfcoords, int(ndims[1]*ndims[2])) # zhalf of every gridbox in order of gbxindex
    xhalfcoords = np.tile(np.repeat(xhalfcoords, ndims[0]), int(ndims[2]))
    yhalfcoords = np.repeat(yhalfcoords, ndims[0]*ndims[1])

    return zhalfcoords, xhalfcoords, yhalfcoords

def allgbxfullcoords_fromgridfile(gridfile, COORD0=False):

    gbxbounds, ndims = read_dimless_gbxboundaries_binary(gridfile,
                                                        COORD0=COORD0,
                                                        return_ndims=True)
    return fullcoords_forallgridboxes(gbxbounds, ndims)

def read_dimless_gbxboundaries_binary(filename, COORD0=False,
                                      return_ndims=False):
    ''' return dictionary for gbx indicies to gbx boundaries by
    reading binary file. Return dimensionless version if COORD0
    not give (=False). '''

    data, ndata_pervar = readbinary(filename)
    datatypes = [np.uint, np.uintc, np.double]

    ll = [0,0] # indexs for division of data list between each variable
    for n in range(1, len(ndata_pervar)):
        ll[n-1] = np.sum(ndata_pervar[:n])
    
    ndims = np.asarray(data[:ll[0]], dtype=datatypes[0])
    gbxidxs = np.asarray(data[ll[0]:ll[1]], dtype=datatypes[1]) 

    ngridboxes = int(np.prod(ndims))
    if len(gbxidxs) != ngridboxes:
        err = "number of gridbox indexes not consistent with (z,x,y) dims"
        raise ValueError(err)
    
    boundsdata = np.asarray(data[ll[1]:], dtype=datatypes[2])
    boundsdata = np.reshape(boundsdata, [ngridboxes, len(boundsdata)//ngridboxes])
     
    if COORD0:
        boundsdata = boundsdata * COORD0
    
    gbxbounds = {gbxidxs[i]: boundsdata[i] for i in range(ngridboxes)}
    
    if return_ndims:
        return gbxbounds, ndims
    else:
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
        fig.savefig(binpath+"gridboxboundaries.png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"gridboxboundaries.png")
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


def calc_gridboxvols(gbxbounds):

    gbxvols = []
    for gbxindex, bounds in gbxbounds.items():
        zwidth = bounds[1] - bounds[0]
        xwidth = bounds[3] - bounds[2]
        ywidth = bounds[5] - bounds[4]
        
        gbxvols.append(zwidth * xwidth * ywidth )
    
    return gbxvols


def domaininfo(gbxbounds):

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds) 
    domainvol = calc_domainvol(zhalf, xhalf, yhalf)

    gridboxvols = calc_gridboxvols(gbxbounds)
    num_gridboxes = len(gridboxvols)

    return domainvol, gridboxvols, num_gridboxes

def grid_dimensions(gbxbounds):

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds) 
    
    if len(zhalf) == 2:
        griddims = 0
    elif len(xhalf) == 2:
        griddims = 1
    elif len(yhalf) == 2:
        griddims = 2
    else:
        griddims = 3
    
    extents, spacings = [], []
    for half in [zhalf, xhalf, yhalf]:
        extents.append([np.amin(half), np.amax(half)])
        spacings.append(abs(half[1:]-half[:-1]))
    
    return extents, spacings, griddims

def print_domain_info(constsfile, gridfile):
    ''' create values from constants file & config file
    required as inputs to create initial 
    superdroplet conditions '''
    
    COORD0 = get_COORD0_from_constsfile(constsfile)
    gbxbounds =  read_dimless_gbxboundaries_binary(gridfile, COORD0) 
    
    domainvol, gridboxvols, num_gridboxes = domaininfo(gbxbounds)
    xtns, spacings, griddims = grid_dimensions(gbxbounds) 
    ztot = abs(xtns[0][0] - xtns[0][1])
    xtot = abs(xtns[1][0] - xtns[1][1])
    ytot = abs(xtns[2][0] - xtns[2][1])

    inforstr = "\n------ DOMAIN / GRIDBOXES INFO ------\n"+\
    "------------- "+str(griddims)+"-D MODEL -------------\n"+\
    "domain dimensions: ({:3g}x{:3g}x{:3g})m^3\n".format(ztot, xtot, ytot)+\
    "domain no. gridboxes: "+str(len(spacings[0])) +\
          "x"+str(len(spacings[1]))+"x"+str(len(spacings[2]))+"\n"+\
    "domain z limits: ({:3g},{:3g})m\n".format(np.amin(xtns[0]),np.amax(xtns[0]))+\
    "domain x limits: ({:3g}, {:3g})m\n".format(np.amin(xtns[1]),np.amax(xtns[1]))+\
    "domain y limits: ({:3g}, {:3g})m\n".format(np.amin(xtns[2]), np.amax(xtns[2]))+\
    "mean gridbox z spacing: {:3g} m\n".format(np.mean(spacings[0]))+\
    "mean gridbox x spacing: {:3g} m\n".format(np.mean(spacings[1]))+\
    "mean gridbox y spacing: {:3g} m\n".format(np.mean(spacings[2]))+\
    "mean gridbox volume: {:3g}".format(np.mean(gridboxvols))+" m^3\n"+\
    "total domain volume: {:3g} m^3\n".format(domainvol)+\
    "total no. gridboxes: "+str(num_gridboxes)+\
    "\n------------------------------------\n"
    print(inforstr)
