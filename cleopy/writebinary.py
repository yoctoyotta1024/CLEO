"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: writebinary.py
Project: cleopy
Created Date: Tuesday 7th May 2024
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


def writebinary(filename, data, ndata, datatypes, units, scale_factors, metastr):
    """'data' is 1D array containing continuous list of variables.
    The number of data points of each variable is give by it's index in
    'ndata', likewise it's datatype, unit and scale_factor are in
    those lists. 'data' is written to binary file with this metadata
    beforehand and a global metadata string explaingnig how to interpret
    the file"""

    check_validinputs(data, ndata, datatypes, units, scale_factors)

    nvars = np.uintc(len(ndata))

    metamaker = CreateMetadataForBinaryArray(
        nvars, ndata, datatypes, units, scale_factors, metastr
    )
    metadata, metaformat = metamaker.get_metadata()

    dataformat = get_dataformat(nvars, ndata, datatypes)

    array2write = metadata + data
    format = metaformat + dataformat

    print("Writing gridbox boundaries binary file to:\n " + str(filename))
    s = struct.pack(format, *array2write)
    f = open(filename, "wb")
    f.write(s)
    f.close()


def check_validinputs(data, ndata, datatypes, units, scale_factors):
    """' check that each variable's scale_factor is a double, unit is
    a binary encoded single character and that the dataype given matches
    the type of the data"""

    if any([type(s) != np.double for s in scale_factors]):
        raise ValueError("type of scale_factors must be C type double")

    for u in units:
        if not isinstance(u, bytes) or np.size(u) != 1:
            raise ValueError("type of units is not binary C type char")

    i = 0
    for j, n in enumerate(ndata):
        if any([type(d) != datatypes[j] for d in data[i : i + n]]):
            err = (
                "stated datatype "
                + str(datatypes[j])
                + " doesn't match type(data) "
                + str(type(data[n]))
            )
            raise ValueError(err)
        i += n


def get_dataformat(nvars, ndata, datatypes):
    dtc = DataTypeCodes()

    dataformat = ""
    for n in range(nvars):
        nvar = np.uintc(ndata[n])
        dataformat += dtc.d2f[datatypes[n]] * nvar

    return dataformat


class DataTypeCodes:
    def __init__(self):
        self.d2f = {
            # dict for converting dtype to struct formatting code
            np.uintc: "I",
            np.double: "d",
            np.uint: "Q",
            type("c"): "c",
        }

        self.d2binaryf = {
            # dict for converting dtype to binary encoded struct formatting
            np.uintc: b"I",
            np.double: b"d",
            np.uint: b"Q",
            type("c"): b"c",
        }

    def dtype2bytesize(self, datatype):
        """returns C type unsigned int for the size of the
        datatype given int's format when stored in binary using
        python's struct module"""

        dcode = self.d2f[datatype]

        return np.uintc(struct.calcsize(dcode))

    def format_size(self, format):
        """returns size in bytes of data stored in a given format
        using python's struct module"""

        bytesize = 0
        for c in format:
            bytesize += struct.calcsize(c)

        return bytesize


class MetadataPerVariable(DataTypeCodes):
    def __init__(self):
        super(MetadataPerVariable, self).__init__()

        self.mpv_format = "IIIccd"
        self.mpv_bytesize = self.format_size(self.mpv_format)

    def varmetadata(self, datap0, ndata, datatype, unit, scale_factor):
        bytespos0 = np.uintc(datap0)  # position of first datapoint of var
        varsz = self.dtype2bytesize(datatype)  # size in bytes of 1 datapoint of var
        nvar = np.uintc(ndata)  # number of datapoints of var
        vartype = self.d2binaryf[datatype]  # binary char symbolising datatype of var
        varunts = unit  # binary char symbolising units of var when * sclae_factor
        varsf = np.double(scale_factor)  # double for scale_factor constant

        metapervar = [bytespos0, varsz, nvar, vartype, varunts, varsf]

        varbytes = varsz * nvar

        return metapervar, varbytes


class CreateMetadataForBinaryArray(MetadataPerVariable):
    def __init__(self, nvars, ndata, datatypes, units, scale_factors, metastr):
        super(CreateMetadataForBinaryArray, self).__init__()

        self.nvars = nvars
        self.ndata = ndata
        self.datatypes = datatypes
        self.units = units
        self.scale_factors = scale_factors

        self.gblmetastr = (
            "4 unsigned ints before this metadata string are"
            + " [1. position of first byte of data (after all the metadata),"
            + " 2. no. bytes of (this) global metadata string, 3. no. bytes"
            + " per variable specific metadata, 4. no. of variables in data]."
            + " After this global metadata string comes variable specific"
            + " metadata. For each variable, this is 3 unsigned ints, 2 chars"
            + " and then a double; it states: [1. position of first databyte,"
            + " 2. size (in bytes) of one datapoint, 3. no. of datapoints,"
            + " 4. char to indicate python struct type, 5. char to indicate"
            + " the units once multiplied by, 6. the scale factor]. "
            + metastr
        )

    def get_metadata(self):
        gblmeta, gblmeta_format, gblmeta_bytes = self.metastr_to_chars(self.gblmetastr)

        datap0 = (
            self.format_size("IIII") + gblmeta_bytes + self.mpv_bytesize * self.nvars
        )
        metapvars, metapvars_format, metapv_bytes = self.variables_metadata(datap0)

        metaformat = "<IIII" + gblmeta_format + metapvars_format
        metadata = [datap0, gblmeta_bytes, self.nvars, metapv_bytes]
        metadata += gblmeta + metapvars

        return metadata, metaformat

    def metastr_to_chars(self, metadatastr):
        """returns metadata string as list of characters encoded
        as binary bytes alongside the corresponding format interpretation
        and total size of the metadata characters (in bytes)"""

        metachars = [m.encode() for m in [*metadatastr]]

        # interpret this metadata as c type characters
        metaformat = "c" * len(metachars)

        metabytes = struct.calcsize("c") * len(metachars)

        return metachars, metaformat, metabytes

    def variables_metadata(self, datap0):
        metapervars = []
        metapervars_format = ""

        # datap0 is position in bytes of the first datapoint of a variable in file
        for n in range(self.nvars):
            metapvar, varbytes = self.varmetadata(
                datap0,
                self.ndata[n],
                self.datatypes[n],
                self.units[n],
                self.scale_factors[n],
            )

            metapervars.extend(metapvar)
            metapervars_format += self.mpv_format

            datap0 += varbytes

        return metapervars, metapervars_format, self.mpv_bytesize
