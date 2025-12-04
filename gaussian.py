from constants import * 
import numpy as np 

def gaussian_kernel_lumo_density_new(central_wavelength_nm,beta_fwhm,wavelength_range_nm):
    
    lam_fwhm_nm = central_wavelength_nm * beta_fwhm
    #convert to standard deviation.
    lam_sigma_nm = lam_fwhm_nm/2.355
    
    expon = (central_wavelength_nm - wavelength_range_nm)/(lam_sigma_nm)
    norm_factor = 1.0 / (SQRT_TWOPI * lam_sigma_nm *10.0) #factor of ten to put in per angstrom.

    return norm_factor * np.exp (-np.power(expon,2)*0.5) 

def gaussian_kernel_lumo_density(central_wavelength_nm,beta,wavelength_range_nm):
    #returns a gaussian Luminosity density (in units of erg s^{-1} ang^{-1}) with unit luminosity.
    #i.e: an integral over this gaussian returns 1 erg/s

    #so, multiplying this Kernel by the desired lumoninosity gaurantees an integal over the line gives desired lumo,
        #print(central_wavelength_nm)
    expon = (central_wavelength_nm - wavelength_range_nm)/central_wavelength_nm
        #print(expon[0],central_wavelength_nm)

        #if any(expon  == np.nan):
            #print(central_wavelength_nm)

    
    #for this simple case, i assume we integrate over the whole Gaussian. 

    #lumo in ergs/s/ang
    norm_factor = 1.0 / (SQRT_PI * beta * central_wavelength_nm*10)
    
    expon/= beta 

    return norm_factor * np.exp (-np.power(expon,2) ) 

def trap(wl,lumo):
    integral = 0
    for ii in range(0,len(wl)-1,1):
        integral += 0.5 * (wl[ii+1]-wl[ii]) * (lumo[ii]+lumo[ii+1])
    return integral