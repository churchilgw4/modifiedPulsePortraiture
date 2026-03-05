import pygeostat
import numpy as np
def gaussian_test(points, log_w):
    skew_phi = pygeostat.statistics.utils.weighted_skew(points.T[0], wts=np.exp(log_w))
    skew_DM = pygeostat.statistics.utils.weighted_skew(points.T[1], wts=np.exp(log_w))
    kurtosis_phi = pygeostat.statistics.utils.weighted_kurtosis(points.T[0], wts=np.exp(log_w))
    kurtosis_DM = pygeostat.statistics.utils.weighted_kurtosis(points.T[1], wts=np.exp(log_w))
    cov = pygeostat.statistics.utils.weighted_covariance(points.T[0], points.T[1], wt=np.exp(log_w))
    gt = {'skew_phi': skew_phi, 'skew_DM': skew_DM, 'kurtosis_phi': kurtosis_phi, 'kurtosis_DM': kurtosis_DM, 'covariance': cov}
    return gt