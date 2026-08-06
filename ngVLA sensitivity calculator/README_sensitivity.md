ngVLA sensitivity calculator written in Python by Dr. Viviana Rosero.


This version uses  T'sys (i.e., the antenna system temperature corrected to the top of the atmosphere --see ngVLA memo #96) and aperture efficiency values from Wes Grammer's front-end cascade model v15 2024.01.07. This version also accounts for "rms" confusion using the cumulative 1.4 GHZ source counts from Matthews et al. (2021) and provided by Jim Condon (private communication). The “rms” confusion is added in quadrature with the theoretical thermal rms for both the point source and brightness sensitivity. Additionally, the allowed synthesized beam sizes for the requested subarray are now bounded to always provide physically achievable results.  Note that in line mode, this script calculates the sensitivity and brightness temperature at the requested frequency without accounting for variation over the bandwidth, which is assumed to be small and representative of a single channel. In continuum mode, the sensitivity and brightness temperature are averaged over the entire receiver bandwidth for bands 1-5 and the results are presented at the band center. For continuum mode at band 6, the calculations are always averaged over 20 GHz of instantaneous correlated bandwidth and the results are reported at the center of this bandwidth. When a desired frequency is provided such that part of the 20 GHz would fall outside the band 6 frequency limits, the bandwidth is shifted just enough to lie completely within these limits.


Requires:   scipy version >= 1.0.0

It should work in Python versions 2 and 3 provided the above
requirement is satisfied.  It will not work in CASA version 5 because
the scipy version is too old, but will work in CASA version 6.

It can be run from the command line or the individual functions
can be imported and used directly.

See the command line help for further instructions, i.e.,

./ngVLA_sensitivity_calculator.py -h
