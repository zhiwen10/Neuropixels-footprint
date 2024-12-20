# Neuropixels-footprint
Get footprint metric from all Neuropixels probe

## Description
This repository calculates footprint metric from Neuropixels mean waveforms, a metric to evaluate waveform spead from peak site.
Footprint is defined as the radius from peak site, at which the mean waveform peak voltage drops below 30 μV.
Units with footpirnt <20 μm are most likely axons. For details, check out [Ye & Shelton et al, 2023](https://www.biorxiv.org/content/10.1101/2023.08.23.554527v3).

![examples](https://github.com/zhiwen10/Neuropixels-footprint/blob/main/examples2.png)

### Input: 

Mean waveforms: ncluster x nchan(384) x tSampleN(82)

ChanMap to indicate site positions in : xcoords (384x1), ycoords (384x1) 

### Output:

Footprint in μm.

### Workflow:

(1) calculate waveform amplitude at each site

(2) transform waveform amplitude to high resolution matrix (at 1 μV) by 2d interpolate, using site positions in chanMap

(3) extract amplitude samples at different radius from peak amplitude site

(4) average across amplitude samples at different radius, and find radius where mean amplitudes drops to 30 μV. 

### Example data:
We include example mean waveforms in the "data" folder (NP1.0, NP2.0 and NP-UHD data in *.mat and *.npy format).

The corresponding channel maps are in the "chanMaps" folder.   

To extract mean waveforms from your own Neuropixels data *.bin files, we provide a Matlab script "data_preprocessing.m". 

## python version

Instructions
```bash
# locate to folders with code
cd root/Neuropixels-footprint/python_version
# Create new virtual environment
conda env create --file requirements.yaml
# Activate the environment
conda activate footprint
# open jupyter notebook and navitage to footprint_all.ipynb
jupyter notebook
```
