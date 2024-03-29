#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Taken from https://github.com/edm1/rank-based-INT on 2021-06-15

import sys
import os
import numpy as np
import pandas as pd
import scipy.stats as ss


def rank_INT(series, c=3.0 / 8, stochastic=True):
    """ Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.

        Args:
            param1 (pandas.Series):   Series of values to transform
            param2 (Optional[float]): Constand parameter (Bloms constant)
            param3 (Optional[bool]):  Whether to randomise rank of ties
        
        Returns:
            pandas.Series
    """

    # Check input
    assert (isinstance(series, pd.Series))
    assert (isinstance(c, float))
    assert (isinstance(stochastic, bool))

    print("1 - series\n", series)
    # Set seed
    np.random.seed(123)

    # Take original series indexes
    orig_idx = series.index
    print("2 - orig_idx\n", orig_idx)
    # Drop NaNs
    series = series.loc[~pd.isnull(series)]
    print("3 - series after drop\n", series)

    # Get ranks
    if stochastic == True:
        # Shuffle by index
        series = series.loc[np.random.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence why we randomised the series)
        rank = ss.rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)
    print("4 - rank\n", rank)

    # Convert rank to normal distribution
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))

    return transformed.reindex(index=orig_idx)  # Tweaked by Antton 2022-11-11


def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2 * c + 1)
    return ss.norm.ppf(x)


def test():
    # Test
    s = pd.Series([2, 1, 1, np.nan, 4, 3], index=["a", "b", "c", "d", "e", "f"])
    res = rank_INT(s, stochastic=True)
    print(res)

    return 0


if __name__ == '__main__':
    test()
