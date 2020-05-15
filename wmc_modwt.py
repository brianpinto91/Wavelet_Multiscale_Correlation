'''******************************************************************************
*Programming language: Python
*Filename: wmc_modwt.py
*Description: This module is used to calculate wavelet multiscale correlations between time series with regular time interval.
*Institute: DLR / ILT
*Author: Brian Pinto
*Version: 1.0
*Date: 21.10.2019
*The script has been tested with the following version:
    Python: 3.7.3
    Spyder: 3.3.6
    numpy: 1.16.4
    pywt: 1.0.3
'''

import numpy as np
import pywt #PyWavelets is open source wavelet transform software for Python.
import math

FILTER = 'db1' #Type of wavelet filter to be used. Refer the manual for the package PyWavelets https://github.com/pistonly/modwtpy

def circular_convolve_d(h_t, v_j_1, j):
    '''
    This function is imported from https://github.com/pistonly/modwtpy
    It calculates the wavelet coefficients for a particular time scale using the given high pass filter and the low passed time series.
    
    Args:
        *h_t (1D np array): A vector representation of high pass filter.
        
        *v_j_1 (1D np array): The time series after previous (j-1) th decomposition. It has the same length as that of the time series. 
        
        *j (int): Represents j th time scale 
    
    Returns:
        *w_j (1D np array): Wavelet coefficients at j-1 th time scale 
    '''
    N = len(v_j_1)
    L = len(h_t)
    w_j = np.zeros(N)
    l = np.arange(L)
    for t in range(N):
        index = np.mod(t - 2 ** (j - 1) * l, N)
        v_p = np.array([v_j_1[ind] for ind in index])
        w_j[t] = (np.array(h_t) * v_p).sum()
    return w_j


def modwt_1D(series, level):
    '''
    This function is imported from https://github.com/pistonly/modwtpy
    It calculates the wavelet coefficients for all the time scales upto the passed 'level' parameter. 
    
    Args:
        *series (1D np array): A time series as a numpy array.
        
        *level (int): Maximum number of time scale for which wavelet coefficients are required. (Python indices starting from 0 are used for convinience. So 0 means 1st time scale in wavelet theory)
    
    Returns:
        *wavecoeff (nD np array): rows represent the wavelet at different time scales. The columns are of the same lenght as that of the input time series. 
    '''
    wavelet = pywt.Wavelet(FILTER)
    h = wavelet.dec_hi
    g = wavelet.dec_lo
    h_t = -1 * np.array(h) / np.sqrt(2) #Here -1 is added to be consistent with matlab filters
    g_t = np.array(g) / np.sqrt(2)
    wavecoeff = []
    v_j_1 = series
    for j in range(1,level+1,1):
        w = circular_convolve_d(h_t, v_j_1, j)
        v_j_1 = circular_convolve_d(g_t, v_j_1, j)
        wavecoeff.append(w)
    #wavecoeff.append(v_j_1) #we do not need the scaling coefficients
    return np.vstack(wavecoeff)


def var_j(wavCoeff_j, j, M):
    '''
    This function calculates the variance of a time series at a particular time scale using the passed wavelet coefficeints 
    The calculations are based on the formulas provided in the thesis appendix chapter for wavelets. (Note that the python indexing starts from 0. So the formulas are accordingly adopted from the theoretical definitions)
    The function is called from within the function: wmcMatrixForScale().
    
    Args:
        *wavCoeff_j (1D np array): wavelet coefficients at j th scale with the same length as that of the time series.
        
        *j (int): the jth time scale for which wavelet coefficients are required.
        
        *M (int): Length of the wavelet filter used.
        
    Returns:
        *var_j (float): Variance at jth scale of the time series passed in the form of its jth scale wavelet coefficients.
    '''
    T = wavCoeff_j.shape[0]
    M_j = (2**(j) - 1) * (M - 1) + 1
    T_j = T - M_j + 1
    var = 0
    for i in range (M_j-1, T):
        var = var + wavCoeff_j[i]**2
    return var / T_j
    


def cov_j(wavCoeff_X, wavCoeff_Y, scale_j, M):
    '''
    This function calculates the co-variance of two time series at a particular time scale using their passed wavelet coefficeints 
    The calculations are based on the formulas provided in the thesis appendix chapter for wavelets. (Note that the python indexing starts from 0. So the formulas are accordingly adopted from the theoretical definitions).
    The function is called from within the function: wmcMatrixForScale().
    
    Args:
        *wavCoeff_X (nD np array): wavelet coefficients of X at j th scale with the same length as that of the time series X.
        
        *wavCoeff_X (nD np array): wavelet coefficients of Y at j th scale with the same length as that of the time series Y.
        
        *j (int): the jth time scale for which wavelet coefficients are required.
        
        *M (int): Length of the wavelet filter used.
        
    Returns:
        *cov_j (float): Co-variance at jth scale between the time series passed in the form of their jth scale wavelet coefficients.
    '''
    T = wavCoeff_X.shape[0]   #length of the time series 
    M_j = (2**(scale_j) - 1) * (M - 1) + 1
    T_j = T - M_j + 1
    var = 0
    for i in range (M_j-1, T):
        var = var + (wavCoeff_X[i] * wavCoeff_Y[i])
    return var / T_j

def wmcMatrixForScale(wav_coeff, j):
    '''
    This function calculates the wavelet correlation matrix at a particular scale indicated by the passed parameters. 
    The matrix will be a square symmetrical matrix. It dimension will be equal to the number of timeseries considered in the analysis.
    The function is called from within the function: wmcMatrixScales().
    
    Args:
        *wav_coeff (2D np array): Each row represents the wavelet coefficents at jth scale for each time series.
            
        *j (int): The jth time scale for which wavelet coefficients are required.
        
    Returns:
        *wmc_scale_mat (n by n numpy array): Wavelet correlation matrix at a particular scale.
    '''
    wf_l = len(pywt.Wavelet(FILTER).dec_hi)
    wmc_scale_mat = []
    var_vec = []
    
    for x in range(0,wav_coeff.shape[0],1):
        wc_x = wav_coeff[x,:]
        var_x = var_j(wc_x, j + 1, wf_l)
        var_vec.append(var_x)
    var_vec = np.array(var_vec).reshape(1,-1)
    denom = np.sqrt(var_vec.T.dot(var_vec))    
    covar_mat = []
    for x in range(0,wav_coeff.shape[0],1):
        row_covar = []
        for y in range(0,wav_coeff.shape[0],1): 
            wc_x = wav_coeff[x,:]
            wc_y = wav_coeff[y,:]
            covar_xy = cov_j(wc_x, wc_y, j + 1, wf_l)
            row_covar.append(covar_xy)
        covar_mat.append(np.array(row_covar))
    numer = np.array(covar_mat)
    wmc_scale_mat = np.divide(numer, denom, where=abs(denom)>0)
    np.fill_diagonal(wmc_scale_mat, 1)
    return wmc_scale_mat

def wmcMatrixScales(tseries, scales=-1):
    '''
    This function calculates the wavelet correlation matrix at all the requred scales indicated by the passed parameters. 
    The matrix will be a square symmetrical matrix. It dimension will be equal to the number of timeseries considered in the analysis.
    The function is called from within the python script: Yearly_analysis.py
    
    Args:
        *tseries (2D np array): Each row represents the time bins and columns represent different time series
            
        *scales (int): The number of scales at which the correlations are to be calculated. Ex: Scales 3 provides analysis of 1st, 2nd and 3rd time scale. Default is -1 which means calculate at all scales
        
    Returns:
        *wmc_mat (j by n by n numpy array): Wavelet correlation matrix at all the required scales upto j.
    '''
    lengthTseries = tseries.shape[0] 
    numItems = tseries.shape[1] #Number of timeseries
    if scales == -1: #if no maximum scale is passed, calculate the correlations at all possible scales
        scales = int(math.log(lengthTseries,2)) #analysis is possible only for time scales in dyadic of 2
    wav_coeff = []
    for n in range(numItems):
        wav_coeff.append(modwt_1D(tseries[:,n], scales))
    wav_coeff = np.array(wav_coeff)
    wmc_mat = [] #create a list to store the correlation matrix at each scale
    for j in range(1,scales+1,1):
        wmc_mat.append(wmcMatrixForScale(wav_coeff[:,j-1,:], j-1))
    wmc_mat = np.array(wmc_mat)    
    return np.array(wmc_mat)


if __name__ == '__main__':
    #test this module using two sine waves
    t_sam = 0.01
    t = np.arange(0,10.01,t_sam)
    x = np.sin(2*math.pi*1*t).reshape(-1,1)
    y = -np.sin(2*math.pi*1*t).reshape(-1,1)
    n = len(x) 
    tseries = np.concatenate((x,y),axis=1)
    #wavelet multiscale correlation
    wmc_mat = wmcMatrixScales(tseries,1)
   