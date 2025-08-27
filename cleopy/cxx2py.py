"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: cxx2py.py
Project: cleopy
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
function(s) for converting c++ values into python ones
"""


def print_dict_statement(filename, dictname, dict):
    """print key, value pairs for dictionary created by reading filename"""

    print("\n---- " + dictname + " from ", filename, "-----")
    for c in dict:
        print(c, "=", dict[c])
    print("---------------------------------------------\n")


def remove_excess_line(line):
    """removes white spaces and
    and comments from a line"""

    line = line.strip()
    line = line.replace(" ", "")
    if "#" in line:
        line = line[: line.find("#")]

    return line


def line_with_assignment(line):
    """find all the lines in a file
    which have form [...] = [...] ;
    ie. lines which assign variables. return
    list of these lines truncatd at ;"""

    line = line.replace(" ", "")
    if "double" in line or "int" in line:
        ind1 = line.find("=")
        ind2 = line.find(";")
        if ind2 != -1 and ind1 != -1 and line[0] != "/":
            goodline = line[: ind2 + 1]
            return goodline

    else:
        return False


def where_typename_ends(line):
    """finds index where c++ name for
    variable type ends for a const double,
    double, int or const int."""

    if line[:3] == "int":
        x = "int"
    elif line[:6] == "double":
        x = "double"
    elif line[:11] == "constdouble":
        x = "constdouble"
    elif line[:8] == "constint":
        x = "constint"
    elif line[:12] == "constexprint":
        x = "constexprint"
    elif line[:15] == "constexprdouble":
        x = "constexprdouble"
    elif line[:20] == "constexprunsignedint":
        x = "constexprunsignedint"
    elif line[:17] == "constexpruint64_t":
        x = "constexpruint64_t"
    else:
        raise Exception("c++ type for variable not understood")
    return len(x)


def read_cxxconsts_into_floats(filename):
    """returns dictionary of value: float from
    (const) doubles and (const) ints
    assigned in a c++ file. Also returns
    dictionary of notfloats for values that
    couldn't be converted"""

    floats = {}
    notfloats = {}
    with open(filename) as file:
        rlines = []
        filelines = file.readlines()
        for line in filelines:
            goodline = line_with_assignment(line)
            if goodline:
                rlines.append(goodline)

        for line in rlines:
            x = where_typename_ends(line)
            name = line[x : line.find("=")]
            value = line[line.find("=") + 1 : line.find(";")]

            try:
                floats[name] = float(value)
            except ValueError:
                notfloats[name] = value

    return floats


def derive_more_floats(consts):
    """return mconsts dictionary containing
    some derived key,values from values in
    consts dictionary"""

    mconsts = {
        "COORD0": consts["TIME0"]
        * consts["W0"],  # characteristic coordinate grid scale [m]
        "RGAS_DRY": consts["RGAS_UNIV"]
        / consts["MR_DRY"],  # specific gas constant for dry air [J/Kg/K]
        "RGAS_V": consts["RGAS_UNIV"]
        / consts["MR_WATER"],  # specific gas constant for water [J/Kg/K]
        "CP0": consts["CP_DRY"],  # characteristic Heat capacity [J/Kg/K]
        "Mr_ratio": consts["MR_WATER"] / consts["MR_DRY"],
    }

    mconsts["RHO0"] = consts["P0"] / (
        mconsts["CP0"] * consts["TEMP0"]
    )  # characteristic density [Kg/m^3]
    mconsts["MASS0"] = (consts["R0"] ** 3) * mconsts["RHO0"]  # characteristic mass [Kg]

    return mconsts
