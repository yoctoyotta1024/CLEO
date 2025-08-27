"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: editconfigfile.py
Project: cleopy
Created Date: Wednesday 17th January 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

from ruamel.yaml import YAML


def update_param(node, param, new_value):
    """Function to recursively searches for 'param' key in YAML node
    and updates it's value to with 'new_value' when found"""

    if isinstance(node, dict):
        if param in node:
            node[param] = new_value  # update value
            return True
        else:
            for key, val in node.items():
                is_success = update_param(val, param, new_value)
                if is_success:
                    return True
    elif isinstance(node, list):
        for item in node:
            is_success = update_param(item, param, new_value)
            if is_success:
                return True
    return False


def edit_config_params(filename, params2change):
    """rewrites config YAML file with key,value pairs listed in params2change updated to new values
    whilst preserving original YAML file's formatting and comments etc."""

    yaml = YAML()

    # Load the YAML file
    with open(filename, "r") as file:
        data = yaml.load(file)

    # Update the parameters from the YAML file
    for param, new_value in params2change.items():
        is_success = update_param(data, param, new_value)
        if not is_success:
            errmsg = param + " could not be updated to new value: " + str(new_value)
            raise ValueError(errmsg)

    # Overwrite the YAML file
    with open(filename, "w") as file:
        yaml.dump(data, file)
