import numpy as np
from math import sqrt

def gaussian(x, x_0, A, sigma):
    '''
    Parameters
    ----------
    A = area under curve
    sigma = scale parameter
    x_0 = location parameter
    x = array of number(s)
    
    Returns
    -------
    Probability of x given normal distribution with parameters above as an array
    '''
    return A/sigma/sqrt(2*np.pi)*np.exp(-np.power(x-x_0,2)/(2.0*sigma**2))

import csv
import matplotlib.pyplot as plt
import os

def spectrum_grapher(filename, c='k', alpha=1, size=1):
    '''
    This function accepts a CSV file containing absorbance data and plots wavelength versus absorbance.
    IMPORTANT: This script must be in the same directory as the folder containing the spectral data. 
    Ex// This jupyter notebook is in in the same folder as 'CTABSOS-Patman3'
    
    Parameters
    ----------
    filename = Name of csv file (with .csv at the end)
    
    Returns
    -------
    UV/Vis spectrum of data (matplotlib figure)
    '''
    # Begin by getting current working directory to find the filenames containing the spectral data
    cwd = os.getcwd()
    
    # Open the filename passed into the function
    # I copied the 'with' statement from a previous lecture
    with open(cwd + '/CTABSOS-Patman3/{}'.format(filename), 'r', newline = '') as csvfile:
        
        # Create a list of wavelengths and absorbances to be plotted later
        reader = csv.reader(csvfile, delimiter = ',', quotechar = '|')
        wavelengths = []; absorbances = []
        for row in reader:
            # x is wavelength, y is absorbance
            x = float(row[0]); y=float(row[1])
            wavelengths.append(x); absorbances.append(y)
        array = np.vstack((wavelengths, absorbances))
    
    # Plotting
    plt.scatter(array[0], array[1], s=size, color=c, label=filename.strip('.csv'), alpha=alpha)
        
    # Title the axes and create a legend
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absorbance')
    return

import csv
import matplotlib.pyplot as plt
import os

def spectrum_data(filename):
    '''
    This function accepts a CSV file containing absorbance data and plots wavelength versus absorbance.
    IMPORTANT: This script must be in the same directory as the folder containing the spectral data. 
    Ex// This jupyter notebook is in in the same folder as 'CTABSOS-Patman3'
    
    Parameters
    ----------
    filename = Name of csv file (with .csv at the end)
    
    Returns
    -------
    spectrum_data[0] = array of wavelengths 
    spectrum_data[1] = array of absorbances
    spectrum_data[2] = wavelength corresponding to maximum absorbance 
    spectrum_data[3] = maximum absorbance 
    '''
    # Begin by getting current working directory to find the filenames containing the spectral data
    cwd = os.getcwd()
    
    # Open the filename passed into the function
    # I copied the 'with' statement from a previous lecture
    with open(cwd + '/CTABSOS-Patman3/{}'.format(filename), 'r', newline = '') as csvfile:
        
        # Create a list of wavelengths and absorbances to be plotted later
        reader = csv.reader(csvfile, delimiter = ',', quotechar = '|')
        wavelengths = []; absorbances = []
        for row in reader:
            # x is wavelength, y is absorbance
            x = float(row[0]); y=float(row[1])
            wavelengths.append(x); absorbances.append(y)
    return wavelengths, absorbances, wavelengths[absorbances.index(max(absorbances))], max(absorbances)

%matplotlib inline

# Get current working directory
cwd = os.getcwd()

# Get list of all spectral filenames in directory to iterate through and plot
files = os.listdir(cwd + '/CTABSOS-Patman3/')

# Create lists for each of the different types of data
# test_standards = the 3 test standards
# standards = the 3 spectral standards
# samples = the 13 sample spectra
test_standards = [name for name in files if 'TestStd' in name]
standards = [name for name in files if 'Std' in name and name not in test_standards]
samples = [name for name in files if 'Sample' in name]
samples.sort()

# Plotting
for name in samples:
    spectrum_grapher(name, c='c', alpha=0.5, size=2)
spectrum_grapher('VesicleStd.csv', c='r', size=2)
spectrum_grapher('MicelleStd.csv', c='b', size=2)
spectrum_grapher('WaterStd.csv', c='g', size=2)

plt.legend()

%matplotlib inline
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
import warnings

# Get current working directory
cwd = os.getcwd()

# Get list of all spectral filenames in directory to iterate through and plot
files = os.listdir(cwd + '/CTABSOS-Patman3/')

# Here we'll make create fit parameters for each of the sample data using the Gaussiun function defined above.
# Firstly I make some guesses at the parameters of the function to help curve_fit make an estimate. 
best_fit_data = {}
for f in files:
    # This filter is added because we don't really care about the covariance parameters since
    # all of the parameters we're interested in are measured independently of one another.
    warnings.filterwarnings('ignore', message='Covariance of the parameters could not be estimated')
    dt = spectrum_data(f)
    
    # The estimated center of the data is located at the wavelength corresponding to the maximum absorbance.
    # This wavelength has already been determined by the function spectrum_data
    est_center = dt[2]
    
    # The estimed standard deviation (or spread) of the data is found using an incredibly rough formula,
    # one-fourth of the range of the data
    est_std = (max(dt[1])-min(dt[1]))/4
    
    # The estimated area is found by integrating under the curve. Effectively, this area is used to scale
    # the fitted curves.
    est_area = trapezoid(dt[1], dt[0])
    
    # All of these parameters are stored in a dictionary
    best_fit_data[f] = curve_fit(gaussian, dt[0], dt[1], p0=[est_center, est_area, est_std])
    
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance')

# Create lists for each of the different types of data and graph them
# test_standards = the 3 test standards
# standards = the 3 spectral standards
# samples = the 13 sample spectra
test_standards = ['TestStd1.csv', 'TestStd2.csv', 'TestStd3.csv']
for d in test_standards:
    x0, A, std = best_fit_data[d][0]
    plt.plot(spectrum_data(d)[0], 
             gaussian(spectrum_data(d)[0], x0, A, std), 
             c='b', lw=2, label=d.strip('.csv'))
    
standards = ['VesicleStd.csv', 'MicelleStd.csv', 'WaterStd.csv']
for d in standards:
    x0, A, std = best_fit_data[d][0]
    plt.plot(spectrum_data(d)[0], 
             gaussian(spectrum_data(d)[0], x0, A, std), 
             c='r', lw=2, label=d.strip('.csv'))

samples = [name for name in files if 'Sample' in name]
for d in samples:
    x0, A, std = best_fit_data[d][0]
    plt.plot(spectrum_data(d)[0], 
             gaussian(spectrum_data(d)[0], x0, A, std), 
             c='c', lw=2, alpha=0.4)
plt.legend()

%matplotlib inline
from scipy.optimize import linprog

# First establish the fit parameters we found above to refer back to later
micelle_fit = best_fit_data['MicelleStd.csv']
vesicle_fit = best_fit_data['VesicleStd.csv']
water_fit = best_fit_data['WaterStd.csv']

# Now we'll start to build the system of equations which will define the linear program
# Building the constraint library
constraints = {}
for f in samples:
    temp = best_fit_data[f]
    constraints[f] = [temp[0][0], temp[0][1], temp[0][2]]
 
# Creating a list of the constraints put in place by the areas under the curves of the 
# Micelles, Vesicles, and Water
area_constraints = [micelle_fit[0][1], vesicle_fit[0][1], water_fit[0][1]]

# Bringing together the lists we defined above into larger arrays (the matrices depicted above)
# X_sys represents the equality constraints, X_inq represents inequality constraints
A_sys = [[1, 1, 1], area_constraints] 
A_inq = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
b_inq = [0, 0, 0]
convex_combinations = {}
for f in samples:
    temp = best_fit_data[f]
    b_sys = [1, temp[0][1]]
    convex_combinations[f] = linprog(c=[1, -1, -1], A_eq=A_sys, A_ub=A_inq, b_eq=b_sys, b_ub=b_inq)
    
# Printing the data collected using linprog
for i in samples:
    print('{2} Sample Numbers\nMicelle = {0:0.4f}\nVesicle = {1:0.4f}\n'.format(
        convex_combinations[i].x[0], convex_combinations[i].x[1], i.strip('.csv')))

for i in samples:
        prop_m, prop_v = convex_combinations[i].x[0], convex_combinations[i].x[1]
        tot = prop_m + prop_v
        print('{2} Sample Numbers\nMicelle = {0:0.4f}\nVesicle = {1:0.4f}\n'.format(
            prop_m/tot, prop_v/tot, i.strip('.csv')))
        plt.figure()
        spectrum_grapher(i, c='k', alpha=1)
        plt.plot(spectrum_data('MicelleStd.csv')[0], gaussian(spectrum_data(i)[0],
                                           best_fit_data['MicelleStd.csv'][0][0],
                                           prop_m/tot * best_fit_data[i][0][1],
                                           best_fit_data['MicelleStd.csv'][0][2]), 
                                           c='r', alpha=0.6, label='Micelle')
        plt.plot(spectrum_data('VesicleStd.csv')[0], gaussian(spectrum_data(i)[0],
                                           best_fit_data['VesicleStd.csv'][0][0],
                                           prop_v/tot * best_fit_data[i][0][1],
                                           best_fit_data['VesicleStd.csv'][0][2]),
                                           c='b', alpha=0.6, label='Vesicle')
        plt.legend()
