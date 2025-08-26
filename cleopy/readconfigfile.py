"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: readconfigfile.py
Project: cleopy
Created Date: Wednesday 17th April 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

from ruamel.yaml import YAML


def extract_floats(node, floats, notfloats):
    """Function to recursively searches YAML node
    and add any value that is convertible to a float to a dictionary
    with its key"""

    if isinstance(node, dict):
        for key, value in node.items():
            if isinstance(value, (dict, list)):
                extract_floats(value, floats, notfloats)
            else:
                try:
                    floats[key] = float(value)
                except ValueError:
                    notfloats[key] = value
    elif isinstance(node, list):
        for item in node:
            extract_floats(item, floats, notfloats)

    return floats, notfloats


def read_configparams_into_floats(filename):
    """returns dictionary of {key, values} pairs from
    keys from a config yaml file which can be assigned float values.
    Also obtains dictionary of notfloats for values that couldn't be converted
    and converts nspacedims to integer if found in floats."""

    yaml = YAML(typ="safe")
    with open(filename, "r") as file:
        config = yaml.load(file)

    floats, notfloats = {}, {}
    floats, notfloats = extract_floats(config, floats, notfloats)

    try:
        floats["nspacedims"] = int(floats["nspacedims"])  # no spatial coords to SDs
    except KeyError as e:
        print("Warning, ignoring error", e)
        pass

    return floats
