# Script to plot the density in a snapshot
 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys

#######################################################

# To choose
# d = density
# t = temperature
# vr = Radial velocity
# vphi = azimuthal velocity
# vx = x-velocity
# vz = z-velocity
# Av = Extinction
# l = location (1=envelope,2=disk,4=outflow)

variable="d"
rmaxplot_au=200
val_min=1e5
val_max=1e14

#######################################################

rad_path="simu_radmc3d/cloud_2msol/"
run="run_"+sys.argv[1:][0]

# Read data
f=open("output_ape_sph.dat","r")
data=f.read().split("\n")
data.pop(-1) # The last line is empty
f.close()

# Read parameters
loggrid=int(data[0].split(" -- ")[0])
rad_min=float(data[1].split(" -- ")[0])
rad_max=float(data[2].split(" -- ")[0])
nrad=int(data[3].split(" -- ")[0])
theta_min=float(data[4].split(" -- ")[0])
theta_max=float(data[5].split(" -- ")[0])
ntheta=int(data[6].split(" -- ")[0])
nvar=13

var_index= {
    "d": [3,3.83e-24,True,"Density (cm$^{-3}$)"],
    "t": [4,1.0,False,"Temperature (K)"],
    "vr": [5,1e5,False,"V_r (km s$^{-1}$)"],
    "vphi": [7,1e5,False,"V_phi (km s$^{-1}$)"],
    "vx": [8,1e5,False,"V_x (km s$^{-1}$)"],
    "vz": [9,1e5,False,"V_z (km s$^{-1}$)"],
    "Av": [10,1.0,True,"Av"],
    "l": [12,1.0,False,"Location"]
    }

ivar=var_index[variable][0]
var_norm=var_index[variable][1]
logplot=var_index[variable][2]
label=var_index[variable][3]
if ivar==12:
  val_min=1
  val_max=4
  
f=open(rad_path+run+"/gas_temperature.dat","r")
data_t=f.read().split("\n")
data_t.pop(-1)
f.close()

# Build grid
theta=np.linspace(np.pi/2.0-theta_min,np.pi/2.0-theta_max,ntheta)
if loggrid:
  r=np.logspace(np.log10(rad_min),np.log10(rad_max),nrad)
else:
 r=np.linspace(rad_min,rad_max,nrad)
R,THETA=np.meshgrid(r,theta)
TEMP=np.zeros((ntheta,nrad))
DENS=np.zeros((ntheta,nrad))
sphere=np.zeros((nrad*ntheta,nvar))
sphere_t=np.zeros(nrad*ntheta)
sphere_d=np.zeros((nrad*ntheta,nvar))
print(sphere)
# Put data in array
for i, dat in enumerate(data[10:]):
  cell=[float(x) for x in dat.split(" ") if x.strip()]
  sphere_d[i,0:nvar]=cell#[0:nvar-1]
for i,dat in enumerate(data_t):
  sphere_t[i]=float(dat)

# Read density
for i in range(ntheta):
    TEMP[:][i]=sphere_t[i*nrad:(i+1)*nrad]
    DENS[:][i]=sphere_d[i:nrad*ntheta:ntheta,ivar]/var_norm

# Plot things
# plt.figure(figsize=(11,6))
# ax1=plt.subplot(121,projection='polar')
plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(12,9))
ax1=plt.subplot(111,projection='polar')
ax1.set_thetamin(0)
ax1.set_thetamax(90)
ax1.grid(False)
if logplot:
  im=ax1.pcolormesh(THETA,R,DENS,cmap='inferno',norm=matplotlib.colors.LogNorm(vmin=val_min,vmax=val_max),shading="auto")
else:
  im=ax1.pcolormesh(THETA,R,DENS,cmap='inferno',vmin=val_min,vmax=val_max,shading="auto")
plt.axis([np.pi/2.0-theta_max,np.pi/2.0-theta_min,0,rmaxplot_au])
cc=ax1.contour(THETA,R,TEMP,[10,20,30,50,75,100,300],colors="y")
ax1.clabel(cc, inline=True,fmt='%d')
cbar=plt.colorbar(im,ax=ax1,fraction=0.046,pad=0.06)
cbar.set_label(label,rotation=270,labelpad=20)
plt.gca().set_rmin(0.0)

plt.savefig("Map_temp_dens.png")
plt.show()
