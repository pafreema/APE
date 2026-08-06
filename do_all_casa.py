import os
import sys
import numpy as np
import time
from casatools import synthesisutils
from casatasks.private import simutil
from astropy.io import fits

## PARAMETERS AND SET UP

path='/Users/pfreeman/Documents/postdoc/ape_code/example_copy/simu_radmc3d/cloud_2msol/run_003/'
molecule = 'CH3OH'
transition = 107.01
ngvla_config = 'core'
alma_config = '6'
tsource = 3600 # in seconds
tint = 100 # in seconds
direction = 'J2000 04:39:53 +26:03:09' # Taurus # Position of the source that we want to observe

## Set the name of the configuration file for the ngVLA Main subarray:
conf_file = 'RevF_array_configurations/ngvla-revF.{}.cfg'.format(ngvla_config) #'ngvla-main-revC.cfg'
conf_dir = path
conf_path = os.path.join(conf_dir,conf_file)
sys.path.append(path)

#################################################
## Only edit below here if you want to change the simobserve (lines 31-46) or tclean (lines 113-125) parameters

## STEP 1 - SIM OBSERVE

model = conf_dir+'{}_{}GHz_{}/Image_{}_{}_0.fits'.format(molecule,str(transition),alma_config,molecule,alma_config)
simobserve(project = '{}_{}GHz_{}'.format(molecule,transition,ngvla_config),
                    skymodel = model,
                    setpointings = True,
                    indirection = direction,
                    incenter = '',
                    inwidth = '',
                    integration = str(tint)+'s',
                    obsmode = 'int',
                    antennalist = conf_path, 
                    refdate = '2014/05/01',
                    hourangle = 'transit', 
                    totaltime = str(tsource)+'s', #1hrs
                    thermalnoise = '',
                    overwrite=True,
                    graphics = 'none',
                    verbose = True)

################################################
## STEP 2 - ADD SIM NOISE

## Use simutil to read the .cfg file:
u = simutil.simutil()
xx,yy,zz,diam,padnames,_,telescope,posobs = u.readantenna(conf_dir+conf_file)

## for the noise
## ensure ngVLA_senstitvity calculator is in the correct directory
from ngVLA_sensitivity_calculator import *
tsource = tsource
delta_f = 141e-6 # in GHz
clight = 2.99792458e8 # in m/s
delta_v = delta_f*clight/transition # in m/s

a_line,b_line,c_line = sigma_ps_fn(ngvla_config,  freq=int(transition),type_cal = 'line', theta=-1, t_int=tsource, delta_v=delta_v, verbose=True )
sigma_NA = a_line[1]*1e-6   #in Jy/beam  #original in uJy/beam
print ("line noise = {} Jy/beam".format(sigma_NA))

nantennas = len(xx)
nbaselines = nantennas*(nantennas-1)/2
print ('number of baselines: ', nbaselines)
npol = 2.0
nchan = 1.0   #noise per channel
tint = tint
nintegrations = tsource/tint
noise_input = sigma_NA * np.sqrt(nchan* npol * nbaselines * nintegrations ) # jy

sigma_simple ='{}Jy'.format(noise_input)
print ("the input noise is {}".format(sigma_simple))

## Create a copy of the noise-free MS:
os.system('cp -r {0}_{1}GHz_{2}/{0}_{1}GHz_{2}.ngvla-revF.{2}.ms {0}_{1}GHz_{2}/{0}_{1}GHz_{2}_noisy.ms'.format(molecule,str(transition),ngvla_config))

## Open the MS we want to add noise to with the sm tool:
sm.openfromms('{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_noisy.ms'.format(molecule,str(transition),ngvla_config))

## Set the noise level using the simplenoise parameter estimated in the section on Estimating the Scaling Parameter for Adding Thermal Noise:
sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)

## Add noise to the 'DATA' column (and the 'CORRECTED_DATA' column if present):
sm.corrupt()

## Close the sm tool:
sm.done()

################################################
## STEP 3 - DO TCLEAN

t0 = time.time()

model = conf_dir+'{}_{}GHz_{}/Image_{}_{}_0.fits'.format(molecule,str(transition),alma_config,molecule,alma_config)

cellsize = imhead(model,'summary')['incr'][0]*206265
ra_shape = imhead(model,'summary')['shape'][0]
su = synthesisutils()
imsize = su.getOptimumSize(ra_shape)

def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = '{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_{3}_noisy.ms'.format(molecule,str(transition),ngvla_config,n)
    else:
        imagename = imagename = '{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_n{3}_noisy.ms'.format(molecule,str(transition),ngvla_config,n)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile,
            datacolumn='', imagename=imagename, spw='0',
            imsize = [imsize],startmodel='',
            cell=str(cellsize)+'arcsec', specmode='cube',
            gridder='standard',
            deconvolver='hogbom',scales=[0,8,24],
            weighting='briggs',robust=r, uvtaper=str(4*cellsize)+'arcsec',
            threshold = '100.0uJy',niter=n, interactive=False)
    return n

msfile = ['{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_noisy.ms'.format(molecule,str(transition),ngvla_config)]
n =  1000 #100000
r = 0.5

run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minutes'.format(time.time()-t0, (time.time()-t0)/60))

ia.open(conf_dir+'{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_n{3}_noisy.ms.image'.format(molecule,str(transition),ngvla_config,n))
beam=ia.restoringbeam(channel=0,polarization=0)

exportfits(imagename='{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_n{3}_noisy.ms.image'.format(molecule,str(transition),ngvla_config,n), fitsimage='{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_n{3}_noisy.fits'.format(molecule,str(transition),ngvla_config,n))

## add the beam parameters to the FITS header
with fits.open('{0}_{1}GHz_{2}/{0}_{1}GHz_{2}_clean_n{3}_noisy.fits'.format(molecule,str(transition),ngvla_config,n), mode="update") as filehandle:
    filehandle[0].header["BMAJ"] = beam['major']['value']/3600
    filehandle[0].header["BMIN"] = beam['minor']['value']/3600
    filehandle[0].header["BPA"] = beam['positionangle']['value']
