import numpy as np
from math import sqrt

def gaussian(x, x_0, A, sigma):
    '''
    Parameters
    ----------
    A (float) = area under curve
    sigma (float) = scale parameter
    x_0 (float) = location parameter
    x (float) = array of number(s)
    
    Returns
    -------
    Probability of x given normal distribution with parameters above as an array
    '''
    return A/sigma/sqrt(2*np.pi)*np.exp(-np.power(x-x_0,2)/(2.0*sigma**2))
