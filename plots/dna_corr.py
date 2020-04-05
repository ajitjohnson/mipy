# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:54:30 2020
@author: Ajit Johnson Nirmal
DNA correlation over t-cycif cycles plot to assess tissue loss
"""

# Import Packages
import sys
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing as mp


# Function to calculate  correlation
def corr(reference_dna, testing_dna):
    correlation = np.corrcoef(reference_dna.flat, testing_dna.flat)[0, 1]
    return correlation

# Main Function
def dna_corr(image_filepath, output_filepath, file_name=None, channels_per_cycle=4, reference_dna_channel=0):
    
    """   

    Parameters
    ----------
    image_filepath : string
        Path to the image (OME-TIF).
    output_filepath : string
        Path to the location where the plot and CSV file will be saved.
    file_name : string, optional
        Custom file name to save the results. The default is None.
    channels_per_cycle : int, optional
        Number of channels per t-cycif cycle. This assumes that the first channel of every cycle is always
        stained for DNA. The default is 4.
    reference_dna_channel : int, optional
        The DNA channel that will serve as the reference to which all other DNA channels will be compared.
        By default the first DNA channel is considered as reference channel. The default is 0.

    Returns
    -------
    correlation : PDF and CSV
        A plot showing the correlation coefficient between DNA will be saved along with a
        CSV file containing the raw values.
        
    Example
    -------
    image_filepath = '/Users/aj/Desktop/image-1/registration/1.ome.tif'
    output_filepath = '/Users/aj/Desktop/image-1/plots'
    # Run the function
    dna_corr(image_filepath, output_filepath, file_name = 'image_1')

    """
    
    # read image
    image = tiff.imread(image_filepath)
    # get number of cycles, using channels per cycle
    num_cycle = image.shape[0] // channels_per_cycle
    # reference DNA intensity
    reference_dna = image[reference_dna_channel]

    # Run the correlation function
    test_dna = image[[i * channels_per_cycle for i in list(range(num_cycle))]]
    pool = mp.Pool(mp.cpu_count())
    final_corr = pool.starmap(corr, [(reference_dna, testing_dna) for testing_dna in test_dna]) # Apply function
    pool.close()
    
    ## Generate Figure
    fig, ax = plt.subplots()
    plt.plot(final_corr, 'o-')
    plt.ylim(0, 1)
    plt.xticks(np.arange(0, num_cycle, 1.0))
    plt.title("DNA correlation")
    plt.xlabel('Cycle Number')
    ax.set_xticklabels([i + 1 for i in list(range(num_cycle))])
    plt.ylabel('Correlation Coefficient')
    if file_name is None:
        plt.savefig(output_filepath + "/dna_correlation_plot.pdf")
    else:
        plt.savefig(output_filepath + "/" + file_name + "_dna_correlation_plot.pdf")

    ## Generate the excel sheet and save it
    corrcoef_df = pd.DataFrame(final_corr)
    corrcoef_df = corrcoef_df.rename(columns={0: "correlation coefficient"})
    corrcoef_df.insert(0, 'Cycle', range(1, 1 + len(corrcoef_df)))
    if file_name is None:
        corrcoef_df.to_csv(output_filepath + "/dna_correlation.csv", index=False)
    else:
        corrcoef_df.to_csv(output_filepath + "/" + file_name + "_dna_correlation.csv", index=False)

# Run the wrapping function
dna_corr(image_filepath = sys.argv[1], output_filepath = sys.argv[2])


