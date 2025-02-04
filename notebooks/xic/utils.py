import numpy as np
import numba as nb


@nb.njit
def correlation_coefficient(x: np.ndarray, ys: np.ndarray):
    """
    Calculate the correlation coefficient between x and each y in ys.
    Returns a numpy array of the same length as ys.

    Args:
        x: numpy array of shape (n,)
        ys: numpy array of shape (m, n)

    Returns:
        numpy array of shape (m,) where elements are correlation coefficients.
        Returns 0 for cases where either x or y has zero variance.
    """
    n = len(x)
    # Calculate means
    mx = x.mean()
    # Calculate mean for each y array manually since axis parameter isn't supported
    m = len(ys)
    my = np.zeros(m)
    for i in range(m):
        my[i] = np.sum(ys[i]) / n
    
    # Initialize array for results
    result = np.zeros(m)

    var_x = np.sum((x - mx) * (x - mx)) / n
    
    # Calculate correlation coefficient for each y in ys
    for i in range(m):
        # Calculate covariance and variances
        cov = np.sum((x - mx) * (ys[i] - my[i])) / n
        var_y = np.sum((ys[i] - my[i]) * (ys[i] - my[i])) / n
        
        # Handle zero variance cases
        if var_x == 0 or var_y == 0:
            result[i] = 0
        else:
            result[i] = cov / np.sqrt(var_x * var_y)
    
    return result

@nb.njit
def normalize_profiles(intensity_slice, center_dilations=2):
    """
    Calculate normalized intensity profiles from dense array.
    
    Args:
        dense: numpy array where first dimension represents different measurements,
              and subsequent dimensions represent mz and rt
        center_dilations: number of points to consider around center for normalization
    
    Returns:
        numpy array of normalized intensity profiles where values > 0 in center
    """
    
    center_idx = intensity_slice.shape[1] // 2

    # Calculate mean manually instead of using axis parameter
    center_intensity = np.ones((intensity_slice.shape[0], 1))
    intensity_mask = np.zeros(intensity_slice.shape[0], dtype=np.bool_)

    for i in range(intensity_slice.shape[0]):
        window = intensity_slice[i, center_idx-center_dilations:center_idx+center_dilations]
        center_intensity[i, 0] = np.sum(window) / window.shape[0]
        intensity_mask[i] = center_intensity[i, 0] > 0

    # filter for > 0
    center_intensity_divisor = center_intensity[intensity_mask]
    intensity_slice = intensity_slice[intensity_mask]

    # intensity slice normalized by max intensity
    center_intensity_normalized = intensity_slice / center_intensity_divisor
    return center_intensity_normalized

@nb.njit
def median_axis(array, axis=0):
    """Calculate the median along a specified axis.
    
    Args:
        array: Input array
        axis: Axis along which to calculate median (default 0)
    
    Returns:
        Array of medians
    """
    if axis == 0:
        result = np.zeros(array.shape[1])
        for i in range(array.shape[1]):
            result[i] = np.median(array[:,i])
    else:  # axis == 1
        result = np.zeros(array.shape[0])
        for i in range(array.shape[0]):
            result[i] = np.median(array[i,:])
    
    return result