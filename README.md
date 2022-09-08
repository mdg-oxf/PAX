# PAX
This repository contains Matlab code and example array data for demonstrating co-registration of PAM and B-mode ultrasound images using the PAX method, as described in "PAX (passive-active crossing) Method for Sub-Millimeter Co-registration of Passive Acoustic Mapping and B-mode Images"
Michael D. Gray and Constantin C. Coussios, in review (as of 08Sep2022) at IEEE UFFC.

This demonstration uses simulated backscatter data from a point source on the central axis of a C52v curvilinear array.
Signals were propagated through a two-layer model consisting of a fast (muscle/liver) ambient medium and a peripheral layer of slow (fat) material.

The top level script is PAX_curvilinear_array_demo_08Sep2022.m, and it operates on the curvilinear array data in: PAX_curvilinear_array_demo_data_08Sep2022.mat

The script calls the following additional functions:

  teaImg_08Sep2022.m: creates a Time Exposure Acoustics image using conventional passive beamforming
  
  delayfm.m:          applys explicit delays/leads to a matrix of time domain signal data
  
  intersections.m:    finds the intersection of two curves in 2D space (written by Douglas M. Schwarz, dmschwarz@ieee.org)
                      This function was obtained from the matlab file exchange:
                      https://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections?s_tid=ta_fx_results
