# Neural-Spike-Sorter
Program that detects and sorts various subcomponents of inputted neural signals - coded in MATLAB

This program's code is divided into seven blocks: 

* Filtering - uses 200Hz high pass FIR filter OR a custom 'smart filtering' filter
* Neural spike Detection - detects by either approximate energy maximization OR by findpeaks
* Spike alignment - centering spikes based on center of energy OR based on peak location
* Feature extraction - using dimensionality reduction (PCA)
* Clustering - using k-means clustering, by generic k-means clustering or by a personal modification I made
* Classification - creates cluster templates, counts number of spikes in each cluster, and visualizes. 
* Analysis

Note that the project report is attatched in 'Report.pptx'
