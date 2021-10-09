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
