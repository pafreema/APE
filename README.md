This set of files is to complement and expand on the Analytical Protostellar Environment (APE) code by Dr. Pierre Marchand.

Information about APE can be found in Marchand et al. (2025) or at https://bitbucket.org/pmarchan/ape_code/src/master/

The workflow I use for APE and how to use these files is outlined in ape_workflow.sh.

If you would like to make your own linelist files for molecules not included here, use linelists.ipynb. These are needed for step 6, making synthetic observations.
This will need the installation of CASSIS LTE Python, from Dr. Sandrine Bottinelli, found here https://gitlab.in2p3.fr/sandrine.bottinelli/cassis-lte-python

I used APE for ngVLA and ALMA synthetic observations and initial results (moment 0 maps and spectra) are found in the Results folder for several tests combining these parameters, with a 1 hour integration: 
- a 0.2 Msun, 2 Msun, and 8 Msun cloud
- at Class 0 and Class 1 stages
- for CH3OH, CH3CN, CH3CHO, CH3OCH3, HCOOCH3, NH2CHO (transitions in ALMA band 3/ngVLA band 6)
- ALMA configurations 6 and 10 and ngVLA configurations core and core+spiral
- sources located in Taurus and Perseus (to test DEC and distance)
