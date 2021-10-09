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
