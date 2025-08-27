"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: test_publish.py
Project: tests
Created Date: Tuesday 26th August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Check that basic features of cleopy work, test catches cases where
e.g. files are missing so the import doesn't work.
It is recommended to check that e.g. assets are included.
"""

from cleopy import (
    cxx2py,
    editconfigfile,
    geninitconds,
    readbinary,
    readconfigfile,
    writebinary,
)
import cleopy.gbxboundariesbinary_src as gbxboundariesbinary_src
import cleopy.initsuperdropsbinary_src as initsuperdropsbinary_src
import cleopy.sdmout_src as sdmout_src
import cleopy.thermobinary_src as thermobinary_src


def test_import():
    """Check that the imports works."""
    assert cxx2py is not None
    assert editconfigfile is not None
    assert geninitconds is not None
    assert readbinary is not None
    assert readconfigfile is not None
    assert writebinary is not None
    assert gbxboundariesbinary_src is not None
    assert initsuperdropsbinary_src is not None
    assert sdmout_src is not None
    assert thermobinary_src is not None
