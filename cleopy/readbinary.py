"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: readbinary.py
Project: cleopy
Created Date: Tuesday 26th August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""


import numpy as np
import struct

from .writebinary import DataTypeCodes, MetadataPerVariable


def readbinary(filename, isprint=True):
    """return list of vectors containing dimenionsless
    data read from binary file"""

    print("Reading binary file:\n " + str(filename))

    nvars, metabytes, metapervar = read_metadata(filename, isprint=isprint)

    data, ndata_pervar = read_data(filename, nvars, metabytes, metapervar)

    return data, ndata_pervar


def read_metadata(filename, isprint=True):
    """read global metadata and return the variable specific metadata
    in a 2D array with each row being the metadata for a different variable"""

    dtc = DataTypeCodes()

    metabytes, gblmeta_bytes, nvars, mpv_bytes = first4uints(filename, dtc)

    gblmetadata, metapervar = get_metadatapervar(
        filename, dtc, gblmeta_bytes, nvars, metabytes
    )

    if isprint:
        print("Metadata: \n", "'" + gblmetadata + "'")

    return metapervar.shape[0], metabytes, metapervar


def first4uints(filename, dtc):
    with open(filename, mode="rb") as binaryfile:
        bytes4uints = dtc.format_size("IIII")
        uints = struct.unpack("<IIII", binaryfile.read(bytes4uints))

    metabytes, gblmeta_bytes, nvars, metapervar_bytes = uints

    return metabytes, gblmeta_bytes, nvars, metapervar_bytes


def get_metadatapervar(filename, dtc, gblmeta_bytes, nvars, metabytes):
    mpv_format = MetadataPerVariable().mpv_format

    nchars = int(gblmeta_bytes / dtc.dtype2bytesize(type("c")))
    format = "<" + "c" * nchars + mpv_format * nvars

    skipbytes = dtc.format_size("IIII")
    nbytes2read = metabytes - skipbytes
    with open(filename, mode="rb") as binaryfile:
        binaryfile.seek(skipbytes)
        mdata = struct.unpack(format, binaryfile.read(nbytes2read))
        metastr = mdata[:nchars]
        mpervar = mdata[nchars:]

    gblmetadata = "".join([str(i.decode()) for i in metastr])
    metapervar = np.reshape(mpervar, (nvars, int(len(mpervar) / nvars)))

    return gblmetadata, metapervar


def read_data(filename, nvars, metabytes, metapervar):
    dataformat, ndata_pervar = get_dataformat(nvars, metapervar)
    with open(filename, mode="rb") as binaryfile:
        binaryfile.seek(metabytes)
        data = struct.unpack(dataformat, binaryfile.read())

    if np.sum(ndata_pervar) != len(data):
        err = (
            str(len(data))
            + " datapoints read in total is incorrect given"
            + " metadata states no. datapoints per variable is "
            + str(ndata_pervar)
        )
        raise ValueError(err)

    return data, ndata_pervar


def get_dataformat(nvars, metapervar):
    dtypes = metapervar[:, 3]  # struct code for datatype of each variable
    ndata_pervar = np.asarray(
        metapervar[:, 2], dtype=np.uintc
    )  # no. datapoints of each variable

    dataformat = "<"
    for n in range(nvars):
        format = [dtypes[n]] * ndata_pervar[
            n
        ]  # list of binary encoded characters for struct format
        format = "".join(
            [str(f.decode()) for f in format]
        )  # convert to one continous decoded string
        dataformat += format

    return dataformat, ndata_pervar
