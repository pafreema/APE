#!/usr/bin/env python
'''ngVLA subarray sensitivity and key performance metric calculator 
V. Rosero '''



import sys
import pickle
import argparse        
import numpy as np
from scipy.interpolate import PPoly,BSpline,CubicSpline, interp1d



## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
## 
## Function to calculate sensitivity of a subarray
##
## input parameters:
## 
##     subarray ## name of array (string): sba, core, spiral, mid, long, main, 
##                 main+long, spiral+mid, spiral+core mid+long
## 
##     freq     ## frequency in GHz (float): e.g., 17
##         
##     theta    ## resolution in arcsec (float): e.g., 0.5
##                 theta=-1 will calculate native resolution (natural, no taper)
## 
##     t_obs    ## on-source time in hours (float): e.g., 1.  [optional: default=1.0]
## 
##     delta_v  ## channel width in m/s (float): e.g., 2e3.  [optional: default=10e3]
##
## 
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def calculate_sensitivity( subarray, freq=-1, theta=-1, t_obs=1., delta_v=10e3, verbose=True, calc_minmax_theta = False ):
    
    if subarray == 'sba':
        D = 6.            
        if theta != -1:
            print('\n***WARNING*** The calculator does not support tapering of the sba. Values will be reported for natural resolution.')
            theta=-1
    else:
        D = 18.

    if calc_minmax_theta:
        print(f'\ncalculating minimum and maximum resolution for subarray {subarray} and frequency {freq} GHz')    
        x = sigma_ps_fn( subarray=subarray, freq = freq, D = D, type_cal = 'calc_minmax_theta', delta_v=delta_v, t_int= t_obs*3600., theta= theta, verbose=verbose)
        return

    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, theta= theta, delta_v=delta_v,  verbose=False)    

    type_cal = 'continuum'
    print('\n\n {0} calculation:'.format(type_cal))
    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, type_cal = type_cal, delta_v=delta_v, t_int= t_obs*3600., theta= theta, verbose=verbose)
    band_key, theta2, delta_freq, eta_w, freq2 = c[0], c[1], c[2], c[3], a[0]
    str1 = '{0} array, {4:.1f} hours at frequency {1} GHz ({2}) with resolution= {3} mas  (eta_w = {5:.2f})'.format(subarray, freq2, band_key, theta2*1e3, t_obs, eta_w)
    print('\n'+str1+'\n'+'-'*len(str1))
    print('continuum point source sensitivity: {0} uJy/beam'.format(a[1] ))
    print('continuum brightness sensitivity: {0} K'.format(a[2] ))

    type_cal = 'line'
    print('\n\n {0} calculation:'.format(type_cal))
    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, type_cal = type_cal, delta_v=delta_v, t_int= t_obs*3600., theta= theta, verbose=verbose)
    band_key, theta2, delta_freq, eta_w, freq2 = c[0], c[1], c[2], c[3], a[0]
    str1 = '{0} array, {4:.1f} hours at frequency {1} GHz ({2}) with resolution= {3} mas  (eta_w = {5:.2f})'.format(subarray, freq2, band_key, theta2*1e3, t_obs, eta_w)
    print('\n'+str1+'\n'+'-'*len(str1))
    print('line point source sensitivity in a {0} km/s ({1:.3f} kHz) channel: {2}uJy/beam'.format(delta_v/1e3, delta_freq/1e3, a[1] ))
    print('line brightness sensitivity in a {0} km/s ({1:.3f} kHz) channel: {2}K'.format(delta_v/1e3, delta_freq/1e3, a[2] ))



# nu in GHz, theta in arcsec.  min theta ~ 678 mas.  returns rms confusion in uJy
# Data from Jim Condon. 
# The cumulative counts of sources at 1.4 GHz from Matthews et al. (2021) are listed below:
def calc_confusion_condon2( nu, theta ):

   logSJy = np.array([ 1.85,  1.75,  1.65,  1.55,  1.45,  1.35,  1.25,  1.15,  1.05,
        0.95,  0.85,  0.75,  0.65,  0.55,  0.45,  0.35,  0.25,  0.15,
        0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75,
       -0.85, -0.95, -1.05, -1.15, -1.25, -1.35, -1.45, -1.55, -1.65,
       -1.75, -1.85, -1.95, -2.05, -2.15, -2.25, -2.35, -2.45, -2.55,
       -2.65, -2.75, -2.85, -2.95, -3.05, -3.15, -3.25, -3.35, -3.45,
       -3.55, -3.65, -3.75, -3.85, -3.95, -4.05, -4.15, -4.25, -4.35,
       -4.45, -4.55, -4.65, -4.75, -4.85, -4.95, -5.05, -5.15, -5.25,
       -5.35, -5.45, -5.55, -5.65, -5.75, -5.85, -5.95, -6.05, -6.15,
       -6.25, -6.35, -6.45, -6.55, -6.65, -6.75, -6.85, -6.95, -7.05,
       -7.15, -7.25, -7.35, -7.45, -7.55, -7.65, -7.75, -7.85, -7.95,
       -8.05, -8.15, -8.25, -8.35, -8.45, -8.55, -8.65, -8.75, -8.85,
       -8.95, -9.05])

   Ntotsr = np.array([3.71e-02, 9.17e-02, 1.73e-01, 2.93e-01, 4.72e-01, 7.40e-01,
       1.14e+00, 1.74e+00, 2.63e+00, 3.98e+00, 5.99e+00, 9.00e+00,
       1.35e+01, 2.02e+01, 3.01e+01, 4.49e+01, 6.67e+01, 9.88e+01,
       1.46e+02, 2.15e+02, 3.14e+02, 4.57e+02, 6.62e+02, 9.49e+02,
       1.35e+03, 1.90e+03, 2.64e+03, 3.62e+03, 4.90e+03, 6.54e+03,
       8.61e+03, 1.12e+04, 1.43e+04, 1.81e+04, 2.26e+04, 2.79e+04,
       3.42e+04, 4.14e+04, 4.98e+04, 5.95e+04, 7.07e+04, 8.35e+04,
       9.82e+04, 1.15e+05, 1.34e+05, 1.57e+05, 1.83e+05, 2.13e+05,
       2.50e+05, 2.94e+05, 3.48e+05, 4.17e+05, 5.06e+05, 6.25e+05,
       7.89e+05, 1.02e+06, 1.34e+06, 1.79e+06, 2.43e+06, 3.34e+06,
       4.60e+06, 6.34e+06, 8.70e+06, 1.18e+07, 1.60e+07, 2.13e+07,
       2.80e+07, 3.63e+07, 4.65e+07, 5.88e+07, 7.32e+07, 9.00e+07,
       1.09e+08, 1.31e+08, 1.56e+08, 1.84e+08, 2.14e+08, 2.48e+08,
       2.84e+08, 3.23e+08, 3.66e+08, 4.11e+08, 4.60e+08, 5.11e+08,
       5.66e+08, 6.24e+08, 6.85e+08, 7.50e+08, 8.18e+08, 8.89e+08,
       9.64e+08, 1.04e+09, 1.12e+09, 1.21e+09, 1.30e+09, 1.39e+09,
       1.49e+09, 1.59e+09, 1.70e+09, 1.81e+09, 1.93e+09, 2.05e+09,
       2.18e+09, 2.31e+09, 2.45e+09, 2.60e+09, 2.75e+09, 2.92e+09,
       3.09e+09, 3.27e+09])

   Omega = 1. / (25 * Ntotsr)
   theta_fwhm = np.sqrt(4* np.log(2) * Omega / np.pi) * np.rad2deg(1.)*3600. # arcsec
   confusion = 1e6*(10**logSJy)/5.  # uJy   scalar = (nu/1.4)**-0.7 
   scalar = (nu/1.4)**-0.7 

   my_interp = interp1d( theta_fwhm, confusion*scalar, bounds_error=False, fill_value=0 )
   sigma_conf = my_interp(theta)

   return sigma_conf

                       
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
##                 
## Function to calculate the Point Source  Sensitivity  (sigma_rms)
## and the Surface Brightness Sensitivity  (sigma_T)
##
## using the equations and values from ngVLA memo 17     
##              
## c = 2.9979246e8 ## [m/s] speed of light
## k_B = 1.38e-23  ## [J/K]  Boltzmann Constant
## eta_c = 0.98    ## correlator efficiency 
## n_pol = 2.      ## number of polarizations
## eta_Q =  0.9625 ## Digitizer quantization efficiency
## D = 18.         ## [m] Diameter of the antenna 
## N_ant           ## Total number of antennas in subarray (from table)
## T_s             ## [K] System Temperature value below the atmosphere (from table)
## Tprime_s        ## [K] System Temperature  value to the top of the atmosphere  (from table; used for calculations below)
## t_int           ## [s] integration (or observation) time
## A               ## [m^2] Geometric area of single antenna
## eta_A           ## aperture (antenna) efficiency (from table)
## SEFD            ## [Jy] System Equivalent flux Density  of an antenna
## delta_v         ## [m/s] velocity resolution
## delta_nu        ## [Hz] correlated bandwidth  (continuum: from table, line: calculated from delta_v)
## theta           ## [arcsec] resolution (clean beam size)
## b_max           ## [m] physical length of longest baseline in subarray
##
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def sigma_ps_fn( subarray = 'core', freq = -1, type_cal = 'continuum', \
                t_int = 3600., N_ant = 214., k_B = 1.38e-23, delta_v = 10e3, theta = -1, verbose=False, \
                eta_c = 0.98, eta_Q =  0.9625, n_pol = 2., D = 18., c = 2.9979246e8 ):

    if subarray == 'sba':
        D = 6.            
        if theta != -1:
            print('\n***WARNING*** The calculator does not support tapering of the sba. Values will be reported for natural resolution.')
            theta=-1
    else:
        D = 18.
 
    
    if isinstance(freq,str):  nu = float( freq.rstrip('GHz') )
    else:  nu = freq

    if not ((theta == -1) or (theta > 0)):
        print ('theta {0} must be either -1 or >0'.format(theta))
        sys.exit(0)
        
    try:
        with open( 'receiver_data.pkl','rb') as in1:
            receiver_data = pickle.load(in1)                           
    except:
        with open( 'receiver_data.pkl','rb') as in1:
            receiver_data = pickle.load(in1,encoding='latin1')                           


    try:
        with open( 'subarray_data.pkl','rb') as in1:
            subarray_parameter_data = pickle.load(in1)                           
    except:
        with open( 'subarray_data.pkl','rb') as in1:
            subarray_parameter_data = pickle.load(in1,encoding='latin1')                                   
    
    try:
        subarray_parameters = subarray_parameter_data[subarray]
    except:
        print ('subarray {0} must be one of:\nsba, core, spiral, mid, long, main, main+long, spiral+mid, spiral+core, \nmid+long'.format(subarray))
        sys.exit(0)



    # find band from frequency
    band_keys = receiver_data.keys()
    valid_band_keys = []
    for key in band_keys:
        band_freqs = receiver_data[key]['freq']             
        if (nu < np.max(band_freqs)) and (nu >= np.min(band_freqs)):
            valid_band_keys.append(key)


    if len(valid_band_keys) >= 1:
        if verbose: print('frequency {0} GHz found in band(s): {1}'.format(nu, valid_band_keys))
    if len(valid_band_keys) > 1:
        band_key = valid_band_keys[-1]
        if verbose: print('selecting band {0} based on lower Tprimesys'.format( band_key))
    elif len(valid_band_keys) == 0:
        print('frequency {0} GHz not found within any receiver bands'.format(nu))
        sys.exit(0)
    else:
        band_key = valid_band_keys[0]        
            
        
    max_bw = receiver_data[band_key]['max_bw']     
    band_freqs = receiver_data[band_key]['freq'] 
    fC = receiver_data[band_key]['freq_center'] 
    all_T_s = receiver_data[band_key]['tSys']
    all_Tprime_s = receiver_data[band_key]['tprimeSys'] 
    all_eta_A = receiver_data[band_key]['antEff'] 

    fL, fH = np.min(band_freqs), np.max(band_freqs)
    cs_T_s = CubicSpline(band_freqs, all_T_s)
    cs_Tprime_s = CubicSpline(band_freqs, all_Tprime_s)
    cs_eta_A = CubicSpline(band_freqs, all_eta_A) 
    
    b_max = subarray_parameters['b_max']
    b_min = subarray_parameters['b_min']
    N_ant = subarray_parameters['N_ant']
    
    if subarray_parameters['spline_type'] == 'cubic':
        cs = PPoly(*subarray_parameters['spline_params'])
    elif subarray_parameters['spline_type'] == 'univariate':
        cs = BSpline(*subarray_parameters['spline_params'])

    nu_line = nu

    if type_cal == 'line':
        T_s = cs_T_s(nu)
        Tprime_s = cs_Tprime_s(nu)
        eta_A = cs_eta_A(nu)
        
    else:
        if band_key != 'BAND_6':
            x_freqs = np.linspace(fL,fH,100)
            nu = np.mean([fL,fH])
            if verbose: print('using entire continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,fL,fH))

        else:
            if (nu-10 < fL): 
                if verbose: print('continuum bandwidth extends beyond receiver edges, shifting center frequency from {0} to {1}'.format(nu, fL+10))
                nu = fL + 10
                if verbose: print('using continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,nu-10,nu+10))
            elif (nu+10 > fH): 
                if verbose: print('continuum bandwidth extends beyond receiver edges, shifting center frequency from {0} to {1}'.format(nu, fH-10))
                nu = fH - 10
                if verbose: print('using entire continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,nu-10,nu+10))
            x_freqs = np.linspace(nu-10,nu+10,100) 

        T_s = np.mean(cs_T_s(x_freqs))
        Tprime_s = np.mean(cs_Tprime_s(x_freqs))
        eta_A = np.mean(cs_eta_A(x_freqs))

    if verbose: print('printing interpolated  Tsys (below atm): {0} K at frequency {2} GHz (band average: {1})'.format(T_s,receiver_data[band_key]['Tsys'],nu))    
    if verbose: print('using interpolated  Tprimesys (top atm): {0} K at frequency {2} GHz (band average: {1})'.format(Tprime_s,receiver_data[band_key]['Tprimesys'],nu))
    if verbose: print('using interpolated eta_A: {0} at frequency {2} GHz (band average: {1})'.format(eta_A,receiver_data[band_key]['eta_A'],nu))        

    
    theta_min_30, theta_max_30 = np.min(10**cs.x), np.max(10**cs.x)
    if type_cal == 'calc_minmax_theta':
        print(f'\ncontinuum (nu= {nu:.3f} GHz):')
        print(f'theta_min={theta_min_30* 30./nu /1e3:.3f} arcsec, theta_max={theta_max_30* 30./nu/1e3:.3f} arcsec')
        print(f'\nline (nu= {nu_line:.3f} GHz):')
        print(f'theta_min={theta_min_30* 30./nu_line /1e3:.3f} arcsec, theta_max={theta_max_30* 30./nu_line/1e3:.3f} arcsec')

#        return (theta_min, theta_max)
        return cs


    if theta == -1: 
        eta_w = 1.
        theta = subarray_parameters['theta_nat_30_GHz'] * 30./nu /1e3 #arcsec
        if verbose: print('using native resolution {0} mas at frequency {1} GHz'.format(theta*1e3,nu))        

    else:

        theta_min_nu, theta_max_nu = theta_min_30* 30./nu /1e3, theta_max_30* 30./nu /1e3
        theta_min_nu_line, theta_max_nu_line = theta_min_30* 30./nu_line /1e3, theta_max_30* 30./nu_line /1e3

        theta_OK = True
        if (theta > theta_max_nu) or (theta < theta_min_nu) or (theta > theta_max_nu_line) or (theta < theta_min_nu_line):
            print(f'\ntheta is not compatible with the {subarray} subarray')
            theta_OK = False

        if (theta > theta_max_nu) or (theta < theta_min_nu):
            print(f'\ncontinuum (nu= {nu:.3f} GHz):')
            print(f'theta_min={theta_min_30* 30./nu /1e3:.3f} arcsec, theta_max={theta_max_30* 30./nu/1e3:.3f} arcsec')

        if (theta > theta_max_nu_line) or (theta < theta_min_nu_line):
            print(f'\nline (nu= {nu_line:.3f} GHz):')
            print(f'theta_min={theta_min_30* 30./nu_line /1e3:.3f} arcsec, theta_max={theta_max_30* 30./nu_line/1e3:.3f} arcsec')

        if not theta_OK:
            sys.exit(0)


        eta_w = float(np.log10( theta*1e3*nu/30. ) )
        
        if verbose: print('using inefficiency factor eta_w: {0} at frequency {1} GHz and resolution {2} mas)'.format(eta_w,nu,theta*1e3))


    



    delta_nu = nu * delta_v/c * 1e9  ## line width in Hz    
    A = np.pi * (D/2.)**2  ## [m**2]
    Field_view = 1.02*(206265/60.)*c/(D*nu*1e9)   ## arcmin  Uniform illumination
    total_eff_A= eta_A * N_ant * A
    Res_max_base = (206265.)*c*1e3/(nu*1e9*b_max)  ## marcsec
    LAS = (206265.)*c/(nu*1e9*b_min)  ## arcsec
    
    
    SEFD = (2* k_B * Tprime_s/ (eta_Q * eta_A * A))/1e-26   ## [Jy]

    avg_SEFD = SEFD * (receiver_data[band_key]['Tprimesys'] / Tprime_s) * (receiver_data[band_key]['eta_A'] / eta_A)

    if verbose: print('calculated SEFD: {0} Jy at frequency {1} GHz (band average: {2})'.format(SEFD,nu,avg_SEFD))        


    data_to_print_performance = (nu, fL, fH,   max_bw, \
                                 Field_view, eta_A, \
                                 total_eff_A/1e3, Tprime_s,\
                                 SEFD, Res_max_base, LAS)

    band_info = (band_key, theta, delta_nu, eta_w)
    
    if type_cal == 'continuum':
        delta_nu = max_bw * 1e9  ## GHz
    elif type_cal == 'line':        
        delta_nu = nu * delta_v/c * 1e9  ## line width in GHz


    sigma_ps = (SEFD / (eta_c * np.sqrt(n_pol* delta_nu * t_int * N_ant*(N_ant-1))))/1e-6  ##uJy        
    sigma_rms = eta_w * sigma_ps
        
    if verbose: print('calculated image thermal noise: {0} uJy/beam at frequency {1} GHz'.format(sigma_rms,nu))        

    sigma_conf = calc_confusion_condon2( nu, theta)

    if verbose: print('calculated confusion RMS: {0} uJy/beam at frequency {1} GHz'.format(sigma_conf,nu))        

    # adding confusion to sigma_T         
    sigma_T = 1.216 * (sigma_rms/(nu**2 * theta**2))   ## K when sigma_ps in uJy
    
    if verbose: print('calculated image brightness sensitivity: {0} K at frequency {1} GHz'.format(sigma_T,nu))        

    sigma_T_conf = 1.216 * (sigma_conf/(nu**2 * theta**2))   ## K when sigma_ps in uJy
    
    if verbose: print('calculated brightness confusion RMS: {0} K at frequency {1} GHz'.format(sigma_T_conf,nu))        


    sigma_rms = (sigma_rms**2 + sigma_conf**2)**.5
    sigma_T = (sigma_T**2 + sigma_T_conf**2)**.5

    data_to_print = (nu,   sigma_rms,  sigma_T, sigma_conf)

    
    return data_to_print, data_to_print_performance, band_info




## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
##                 
## Code to support command line execution
##    
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subarray_help = '''name of array (string): sba, core, spiral, mid, long, main, 
                main+long, spiral+mid, spiral+core, mid+long'''

freq_help = ''' frequency in GHz (float): e.g., 24.5'''

theta_help = '''resolution in arcsec (float or -1): e.g., 0.5. 
                theta=-1 will calculate native resolution (natural, no taper)'''

t_obs_help = '''on-source time in hours (float): e.g., 4.  Default: 1'''

delta_v_help = '''channel width in m/s (float): e.g., 2e3.  Default: 10e3'''

calc_minmax_theta_help = '''for mode theta=-1, also return the minimum and maximum resolution in arcsec for the given frequency and subarray'''

    
if __name__ == '__main__':     

    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('subarray', type=str, help=subarray_help)
    parser.add_argument('frequency', type=float, help=freq_help)
    parser.add_argument('theta', type=float, help=theta_help)
    parser.add_argument('--t_obs', type=float, default=1.0, help=t_obs_help)
    parser.add_argument('--delta_v',type=float, default=10e3, help=delta_v_help)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('--calc_minmax_theta', action="store_true", help=calc_minmax_theta_help)

    args = parser.parse_args()    
    print(__doc__)
    calculate_sensitivity( args.subarray, args.frequency, args.theta, args.t_obs, args.delta_v, args.verbose, args.calc_minmax_theta )    
    


#
