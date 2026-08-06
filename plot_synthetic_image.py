# Script to plot the synthetic images
 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs import WCS
from astropy import units as u
from astropy.visualization.wcsaxes import add_beam
import sys


setup_name='cloud_2msol'
run=sys.argv[1:][0]
molecule=sys.argv[1:][1]
transition=sys.argv[1:][2]+'GHz'
ext=sys.argv[1:][3]
datfile =fits.open('./'+molecule+'_'+transition+'_'+ext+'/baseline_synthetic_Image_'+molecule+'_'+ext+'_0.fits')

data = SpectralCube.read(datfile)
data_vel = data.with_spectral_unit(u.km/u.s, velocity_convention='radio')
wcs = WCS(datfile[0].header)
hdr = datfile[0].header
mom0 = data_vel.moment(order=0)

lo=int(np.round(data.shape[1]/2)-np.round(data.shape[1]/24))
hi=int(np.round(data.shape[1]/2)+np.round(data.shape[1]/24))

fig = plt.figure(figsize=(12,8), facecolor='w')
ax = fig.add_subplot(1,1,1, projection = mom0[lo:hi,lo:hi].wcs)
im = ax.imshow(mom0[lo:hi,lo:hi].value, cmap='inferno')
ax.set_xlabel('\n RA', fontsize=16, color='k')
ax.set_ylabel('DEC', fontsize=16, color='k')
cbar = plt.colorbar(im)
im.set_clim(0, np.nanmax(mom0.value))
cbar.set_label('\n Integrated Intensity (Jy/beam*km/s)', size=16)
add_beam(ax, header=hdr, corner='bottom left', facecolor='w')
#plt.show()
plt.savefig('./baseline_synthetic_Image_'+molecule+'_'+transition+'_'+ext+'_'+run+'_mom0.png')
plt.close()

data_sub=data_vel.subcube(xlo=lo, xhi=hi, ylo=lo, yhi=hi)
spectrum=data_sub.mean(axis=(1,2))

fig = plt.figure(figsize=(12,8), facecolor='w')
ax = fig.add_subplot(1,1,1)
im = ax.step(data_sub.spectral_axis, spectrum)
ax.set_xlabel('\n Velocity (km/s)', fontsize=16, color='k')
ax.set_ylabel('Intensity (Jy/beam)', fontsize=16, color='k')
#plt.show()
plt.savefig('./baseline_synthetic_Image_'+molecule+'_'+transition+'_'+ext+'_'+run+'_spectrum.png')
plt.close()
