Jupyter matlab script, parameter files and Jupyter notebooks to analyse and plot microscopy data with the ImAnalaysis pipeline and Python.
Folders contain files that are named and contain data approximately as stated below:

scriptImAnatoExcel.m
--------------------
Matlab scripts that runs the ImAnalysis pipeline of the ELfLab to segement cells using a Unet and detecting dots using Wavelet. The Parameter files are written for even and odd Positions, as the imaging is done on one chip witth two lanes that contain a different strain each. Each parameter file is adjused for the Region of Interest (ROI) of each strain.

ParameterFile.txt
-----------------
The Parameter file with the specifications of ROI and wavelet detection parameters.

CheckDotDetectionParams.ipynb
-----------------------------
Script that is used to process gaussian fits to fluorescent dots obtained during microscopy and detected with the image processing pipeline of the ElfLab. The pipeline output contains the dots detected in each postition, the time point that the image has been recorded, the intensity of the dot as well as the sigma of the gaussian fit. It also contains the number of cells segmented in each position to calculate the average number of redorded dots per cell.

PlotResults.ipynb
-----------------------------
The repo also contains a jupyter notebook for plotting of the detected spots and processing data obtained in different experiments.

DotDetectionFunctions.py
-------------------------
Contains the commonly used functions for analysis of the microfluidics/microscopy experiment. 
