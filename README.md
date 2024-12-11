# Neuropixels-footprint
get footprint metric from all Neuropixels probe

## Description
This repository calculates footprint metric from Neuropixels mean waveforms, a metric to evaluate waveform spead from peak site.
Footprint is defined as the radius from peak site, at which the mean waveform peak voltage drops below 30 μV.
Units with footpirnt <20 μm are mostly to be axons. Details in [Ye & Shelton et al, 2023](https://www.biorxiv.org/content/10.1101/2023.08.23.554527v3).

Input: 
mean waveform: ncluster x nchan(384) x tSampleN(82)
chanMap to indicate site positions in : xcoords (384x1), ycoords (384x1) 

Output:
footprint in μm.

Workflow:
(1) calculate waveform amplitude at each site
(2) 2d interpolate waveform amplitude into high resolution matrix (at 1 μV), using site positions in chanMap
(3) extract amplitude samples at different radius from peak amplitude site
(4) average across amplitude samples at different radius, and find radius when mean amplitudes drops to 30 μV. 

![examples](https://github.com/zhiwen10/Neuropixels-footprint/blob/main/examples.png)

