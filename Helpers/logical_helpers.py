"""
Boolean Helper Functions for Element-Wise Ops on >2 Arrays.
"""

import numpy as np

def xlogical_and(arr):
    """Element-Wise .AND. Operation Over All Arrays."""
    for iarr in range(len(arr)):
        if iarr == 0:
            tmp = arr[0]
        else:
            tmp = np.logical_and(tmp, arr[iarr])
    return tmp
