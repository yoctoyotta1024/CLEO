"""
----- CLEO -----
File: sdtracing.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
!!!NOTE!!! DEPRECIATED SINCE V0.40.0, REPLACED BY SUPERDROPS MODULE !!!!!!
functions to extract attribute data for specifc superdroplets based on the sdIds
e.g. for tracing their trajectories
"""


import numpy as np
import awkward as ak
import random
import warnings

warnings.warn(
    "sdtracing module is depreciated since CLEO v0.40.0,\nreplaced by superdrops"
)


def attr_for_superdroplet(sddata, Id, attr):
    """selects attribute from sddata belonging
    to superdroplet with identitiy 'Id'
    at every output time"""

    bools = ak.Array(sddata.sdId == Id)  # True/False id is found in sdId at each time
    attr4Id = sddata[attr][bools]  # attribute where sdId = Id

    num = ak.num(
        attr4Id
    )  # at each time, number of positions where sdId is found (should be 0 or 1)
    if any(num[num != 1]):
        errmsg = (
            "attribute timeseries has times when more"
            + " than one sdId==Id. num should be list of either 1 or 0"
            + " (for Id found in sddata at given time or not)"
        )
        raise ValueError(errmsg)

    attr4Id = ak.where(
        num == 0, ak.Array([[np.nan]]), attr4Id
    )  # replace empty values with np.nan

    return ak.flatten(attr4Id, axis=1)  # remove excess dimension


def attributes_for1superdroplet(sddata, Id, attrs):
    """selects attributes in 'attrs' from sddata
    belonging to superdroplet with identitiy 'Id'
    at every output time"""

    attrs4Id = {}
    for attr in attrs:
        attrs4Id[attr] = attr_for_superdroplet(sddata, Id, attr)

    return attrs4Id


def attribute_for_superdroplets_sample(
    sddata, attr, ndrops2sample=0, minid=0, maxid=0, ids=[]
):
    """returns 2D array with dimensions [time, SD]
    containing attribute data over time for a sample of
    superdroplets. Sample is either for superdroplets with
    specific Ids in 'ids' list, or sample of 'ndrops2sample'
    randomly selected superdrops with Ids in the range
    [minid, maxid]"""

    if np.any(ids):
        sample = ids
    else:  # ids == []
        population = list(range(minid, maxid, 1))
        if ndrops2sample == 0:
            ndrops2sample = maxid
        sample = random.sample(population, ndrops2sample)

    ndrops_attr = []
    for id in sample:
        attr4Id = attr_for_superdroplet(sddata, id, attr)
        ndrops_attr.append(attr4Id)

    return np.asarray(ndrops_attr).T


def attr_at_times(attrdata, time, times2sel):
    """selects attribute (for all superdroplets)
    at times closest to 'times2sel'"""

    inds = []  # list containing indexes of times closest to times2sel
    for t in times2sel:
        inds.append(np.argmin(abs(time - t)))

    return attrdata[inds]


def attributes_at_times(sddata, time, times2sel, attrs2sel):
    """selects attributes at given times from
    sddata (for all superdroplets in sddata)"""

    selected_data = {}  # dict containting selected attributes at selected times

    for attr in attrs2sel:
        selattr_data = attr_at_times(sddata[attr], time, times2sel)
        selected_data[attr] = selattr_data

    return selected_data


def attrs_for_superdroplets_sample(
    sddata, attrs, ndrops2sample=0, minid=0, maxid=0, ids=[]
):
    """returns dictionary of 2D arrays (with dimensions [time, SD])
    for each attribute in 'attrs' list for a sample of
    superdroplets. Sample is either for superdroplets with
    specific Ids in 'ids' list, or sample of 'ndrops2sample'
    randomly selected superdrops with Ids in the range
    [minid, maxid]"""

    if np.any(ids):
        sample = ids
    else:  # ids == []
        population = list(range(minid, maxid, 1))
        if ndrops2sample == 0:
            ndrops2sample = maxid
        sample = random.sample(population, ndrops2sample)

    data = {}
    for a in attrs:
        data[a] = attribute_for_superdroplets_sample(sddata, a, ids=sample)

    return data
