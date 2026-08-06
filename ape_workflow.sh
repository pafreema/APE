## THIS IS A SET OF STEPS TO RUN APE, WITHOUT THE EXTENDED DESCRIPTIONS OF THE USER MANUAL
## YOU WILL NEED APE, NAUTILUS, RADMC3D AND IMAGER INSTALLED
## https://bitbucket.org/pmarchan/ape_code/src/master/
## https://gitlab.in2p3.fr/LAB/astrochem-tools/pnautilus
## https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/
## https://imager.oasu.u-bordeaux.fr/
## THESE STEPS ALSO CONTAIN EXTRA STEPS BASED ON P FREEMAN'S WORKFLOW
## THE STEPS ARE NUMBERED BASED ON THE APE USER MANUAL

## Set up a folder in your APE installation directory for your run
## Ensure the following are in your folder:
## 1. APE files and folders: parameters.nml, scripts, lib, ape executable, Nautilus_files
## 2. Nautilus files and folders: nautilus_scripts, Nautilus executables: nautilus_code, nautilus_outputs, nautilus_trace_major, nautilus_major_reactions (in your Nautilus installation, from Nautilus/build/src folder)
## 2b. add Nautilus executables to ./Nautilus_files/ folder
## 3. Extra files from https://drive.google.com/drive/folders/1ILehbs3lgz96K2mNJyh1qHaKnXLuw2Cp?usp=sharing
## 3a. do_all_casa.py, plot_temp_dens.py, plot_synthetic_image.py, plot_synthetic_image_casa.py copied into your ./scripts/ folder
## 3b. additional linelist files copied into your ./lib/ folder (you can make your own linelist files from the linelists.ipynb in the same google drive; you will need the CASSIS LTE PYTHON package from https://gitlab.in2p3.fr/sandrine.bottinelli/cassis-lte-python)
## 3c. ngVLA_sensitivity_calculator.py, subarray_data.pkl, receiver_data.pkl
## 3d. RevF_array_configurations - ngVLA configuration files

## for plotting data from fits files, you will need spectral-cube https://spectral-cube.readthedocs.io/en/latest/installing.html
## for CASA simulated observations, you will need casatools and casatasks https://casadocs.readthedocs.io/en/stable/api/casatools.html https://casadocs.readthedocs.io/en/stable/api/casatasks.html

################################################################


## STEP 3.1 - MAKE AN APE SNAPSHOT, SHOWING YOUR PHYSICAL SETUP
## set your parameters in parameters.nml (mass, time, grid size, etc) and ensure the following parameters are set:
## radmc_output=.true.
## use_radmc_temp=.false.
## particle=.false.
## grid_of_part=.false.
open ./example/parameters.nml

## run ape
./ape parameters.nml

## STEP 4.3 - RUN RADMC3D RADIATIVE TRANSFER CODE TO FIND TEMPERATURE EVOLUTION. THIS WILL CREATE FOLDERS EACH CONTAINING A SNAPSHOT IN TIME
## set your parameters: start time, end time, number of runs. this will create snapshots at evenly spaced times based on these parameters.
open ./scripts/setup_radmc3d.sh

./scripts/setup_radmc3d.sh

## enter number of CPUs to use
open ./scripts/do_all_radmc_snapshots.sh

## WARNING this may take hours depending on the number of snapshots and CPUs you have
cd simu_radmc3d/setup_name
../../scripts/do_all_radmc_snapshots.sh

## produce temperature files
../../scripts/gather_average_temperature.sh

cd ../..

## STEP 3.3 - RUN APE 'GRID OF PARTICLES MODE'. THIS WILL CREATE FOLDERS FOR EACH PARTICLE IN YOUR GRID
## keep your physical parameters in parameters.nml the same (mass, time, grid size, etc) noting that the grid size will determine how many particles you have
## ensure the following parameters are set:
## use_radmc_temp=.true.
## particle=.true.
## reverse=.true.
## grid_of_part=.true.
open ./parameters.nml

## if you have compiled the code with mpif90, run the following with N the number of CPUs
#mpirun -np N ./ape parameters.nml
## else
./ape parameters.nml

## plot temperature and density results from run_00*, this is saved in map_dens_temp in this folder
python3 ./scripts/plot_temp_dens.py 00*

## ADDITIONAL STEP BEFORE RUNNING STEP 5.3 - RUN NAUTILUS OUTSIDE OF APE, MODELLING A PRESTELLAR STAGE TO GET SPECIES ABUNDANCES FOR THE PROTOSTELLAR STAGE
## change directory to your Nautilus installation folder
cd ../../Nautilus

cp -r ./example_simulation ./prestellar_phase

cd ./prestellar_phase

## set your physical parameters for the prestellar stage
open parameters.in

## change initial abundances of species (I use the same as Marchand et al (2025) and Ruaud et al (2018))
open abundances.in

../build/src/nautilus_code

## copy abundances.tmp into your Nautilus_files folder in APE
cp abundances.tmp ape_path/Nautilus_files/

## move back to your APE run folder
cd ../../ape_path/

## replace the old abundances file with the new
mv Nautilus_files/abundances.tmp Nautilus_files/abundances.in

## STEP 5.3 - RUN NAUTILUS GAS-GRAIN CODE ON ALL PARTICLES
## Make sure that the Nautilus parameter files, your chemical network, and the nautilus_code and nautilus_outputs executables are present in the Nautilus_files directory
## set the number of CPUs you want to use and ensure line 12 matches the name of your nautilus or nautilus_code executable
## uncomment 'rm -rf *.out' in order to save space and remove output files as you go
## add ./nautilus_outputs below line 12 if you want to create /ab and /struct files, which will allow you to look at individual particle evolutions (this will also take up lots of space)
open ./scripts/do_all_nautilus.sh

## edit the Nautilus parameters if you need
open ./Nautilus_files/parameters.in

## WARNING this step might take hours or days depending on your grid size and computational
./scripts/do_all_nautilus.sh

## edit this script:
## line 10: f=open("../output_ape_sph.dat","r")
## line 46:     file1="../../../particles/part_"+ichar+"/abundances.tmp"
open ./scripts/generate_numberdens_file.py

## edit this script:
## line 15: f=open("../output_ape_sph.dat","r")
open ./scripts/plot_abundance_map.py

## make number density and abundance files, where SSSSS is the name of your molecule
cd simu_radmc3d/setup_name/run_00*
python3 ../../../scripts/generate_numberdens_file.py SSSSS

## plot abundance maps
python3 ../../../scripts/plot_abundance_map.py SSSSS

## STEP 6 - MAKE SYNTHETIC OBSERVATIONS IN IMAGER

## ensure that dust_temperature.dat, number_density_SSSSS.inp are already in this folder
## you will need to copy the following files into this folder:
cp ../../../scripts/make_synthetic_observation.py .
cp ../../../lib/*.inp .

## edit the following for your observational parameters
## go to https://almascience.eso.org/proposing/sensitivity-calculator and find the noise level based on your ALMA configuration and observation parameters
## I change these names where the number is the ALMA configuration (important as the plotting and CASA script below rely on this naming system)
## SimulationDirectorySuffix="_6"
## ImageName="Image_"+species+"_6"
## FitsImage="Image_"+species+"_6.fits"
open make_synthetic_observation.py

## watch for error: if it doesn't succeed the first time you may get multiple folders with *_6_2, *_6_3, ...; edit the successful folder name back to *_6

## WARNING this can take time, and will make a new folder XXXXX_YYYYYGHz_ZZZZZ has been created, with XXXXX the species, YYYYY the average frequency of your spectral window in GHz, and ZZZZZ a user-choosen suffix
python3 make_synthetic_observation.py

cd XXXXX_YYYYYGHz_ZZZZZ

## launch IMAGER
imager @script_XXXXX.ima

## exit IMAGER

## plot image—run is the number of the snapshot, molecule is your species name, transition is the frequency value (no units, same as in the folder name), ext is the ALMA configuration used as the extension in your files (see above)
## or change the script so it calls to the right path
cd ..
python3 ../../../scripts/plot_synthetic_image.py run molecule transition ext

## STEP 6 ADDITION TO APE MANUAL - MAKE SYNTHETIC OBSERVATIONS IN CASA
## copy ngVLA files into run_00* folder
cp ../../../ngVLA_sensitivity_calculator.py ../../../subarray_data.pkl ../../../receiver_data.pkl .
cp -r ../../../RevF_array_configurations .

## this assumes you have already made synthetic observations for ALMA above as it uses the model image. Currently I match a model image for ALMA configuration 6 will ngVLA core configurations, and ALMA configuration 10 with ngVLA spiral configuration.
cp ../../../scripts/do_all_casa.py .

## set file paths, observational and cleaning parameters
open do_all_casa.py

## open CASA
casa

## in CASA run script
execfile('do_all_casa.py')

## exit CASA

## plot image—run is the number of the snapshot, molecule is your species name, transition is the frequency value (no units), ext is the ngVLA configuration used as the extension in your files
python3 ../../../scripts/plot_synthetic_image_casa.py run molecule transition ext

cd ../../..

####################################################################


## OPTIONAL ADDITIONAL STEPS TO LOOK AT INDIVIDUAL PARTICLE EVOLUTION

## STEP 3.2 - RUN APE IN PARTICLE MODE
## some of this assumes you have already run STEP 3.3 - GRID OF PARTICLE MODE above
## set z_ini and x_ini based on particle locations in ./particles/grid_of_part.dat, make note of which particle number it is
## ensure the following parameters are set:
## use_radmc_temp=.true.
## particle=.true.
## reverse=.true.
## grid_of_part=.false.
open parameters.nml

./ape parameters.nml

python3 ./scripts/plot_traj.py

## STEP 5.2 - RUN NAUTILUS ON INDIVIDUAL PARTICLE
## change into the folder for the particle number you chose above
cd particles/part_0****

## I copy the previous outputs into this folder for organization. NOTE if you have already run grid of particle mode you will already have a structure_evolution.dat folder, they should be identical copies.
cp ../../part_traj.dat ../../dens_temp_* ../../structure_evolution.dat .
cp ../../nautilus_scripts/* ../../Nautilus_files/* .

## run Nautilus
./nautilus_code

## create abundance files for all species
./nautilus_outputs

python3 plot_abundances.py species=XXXXX,YYYYY,ZZZZZ

## find important formation and destruction reactions, it will return a plot and some files, compare the IDs in the plot to those of the *.reaction files
./nautilus_trace_major
python3 trace_species.py species=XXXXX

## find important formation and destruction reactions at a specific time, follow inputs on screen
./nautilus_major_reactions

## You may want to delete the large number of *.out files that Nautilus creates. Note if you do this, you will NOT be able to use ./nautilus_trace_reactions or ./nautilus_major_reactions. You will still be able to use plot_abundances.py.
rm -rf *.out
