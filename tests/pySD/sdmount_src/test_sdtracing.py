import numpy as np
import awkward as ak
import pytest

from pySD.sdmout_src.sdtracing import (
    get_awkward_shape,
    ak_digitize_2D,
    ak_digitize_3D,
)


@pytest.mark.parametrize(
    "a, should",
    [
        (ak.Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), [10]),
        (ak.Array([[1, 2, 3], [4, 5, 6]]), [2, 3]),
        (ak.Array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]), [2, 2, 2]),
        (ak.Array([[1, 2, 3], [4, 5, 6, 7]]), [2, np.nan]),
    ],
)
def test_get_awkward_shape(a, should):
    print(get_awkward_shape(a))
    assert get_awkward_shape(a) == should


@pytest.mark.parametrize(
    "x, bins, should",
    [
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]]),
            np.array([0, 3, 6, 100]),
            ak.Array([[1, 1, 2, 0, 4], [2, 2, 3, 4, 3, 4]]),
        ),
    ],
)
def test_ak_digitize_2D(x, bins, should):
    """Tests returns of the ak_digitize_2D function"""
    x_binned = ak_digitize_2D(x, bins)
    # check for equality
    assert ak.sum(x_binned != should) == 0
    # check for same shape
    assert get_awkward_shape(x_binned) == get_awkward_shape(should)
    # check for same counts
    assert ak.count(x_binned) == ak.count(should)


@pytest.mark.parametrize(
    "x, bins, exception",
    [
        (
            ak.Array([[[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]], [None]]),
            np.array([0, 3, 6, 100]),
            ValueError,
        ),
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]]),
            np.array([[0], [3]]),
            ValueError,
        ),
    ],
)
def test_array_to_bin_index_exception(x, bins, exception):
    """Tests for exceptions in the ak_digitize_2D function"""
    with pytest.raises(exception):
        ak_digitize_2D(x, bins)


@pytest.mark.parametrize(
    "x, bins, should",
    [
        (
            ak.Array([[[1, 2, 3], [4, 5]], [[1], [2]]]),
            np.array([0, 2, 4, 6]),
            ak.Array([[[1, 2, 2], [3, 3]], [[1], [2]]]),
        ),
    ],
)
def test_ak_digitize_3D(x, bins, should):
    """Tests returns of the ak_digitize_3D function"""
    x_binned = ak_digitize_3D(x, bins)
    # check for equality
    assert ak.sum(x_binned != should) == 0
    # check for same shape
    assert get_awkward_shape(x_binned) == get_awkward_shape(should)
    # check for same counts
    assert ak.count(x_binned) == ak.count(should)


@pytest.mark.parametrize(
    "x, bins, exception",
    [
        # The first test is not a valid case, because only values in the last axis are allowed
        (
            ak.Array([[[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]], [None, 1]]),
            np.array([0, 3, 6, 100]),
            ValueError,
        ),
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None, 1]]),
            np.array([[0], [3]]),
            ValueError,
        ),
    ],
)
def test_ak_digitize_3D_exception(x, bins, exception):
    """Tests for exceptions in the ak_digitize_3D function"""
    with pytest.raises(exception):
        ak_digitize_3D(x, bins)
