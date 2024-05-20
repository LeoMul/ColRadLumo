path_to_colradpy = 'ColRadPy/'
import sys
sys.path.append(path_to_colradpy)

from colradpy import colradpy 
from time import process_time 

import numpy as np

#masses
NUCLEON_MASS_KG = 1.67493e-27
SOLAR_MASS_KG = 1.989e+30 
SQRT_TWOPI = 2.5066282746310002
#Planck's constant times c in cgs e.g: erg cm s 
HC_CGS = 6.63e-34 * 3e8 * 1e7 * 1e2



class requested_lines:
    def __init__(self,wavelengths,avalues,pecs,pec_levels,lumo_ph,lumo_erg,csf_labels,angular_momenta,energies):
        self.wl_vac_nm = wavelengths
        self.avalues     = avalues

        self.pec_levels = pec_levels
        self.pec        = pecs
        self.lumo_ph    = lumo_ph
        self.lumo_erg   = lumo_erg

        self.csfs_labels = csf_labels
        self.angular_momenta = angular_momenta
        self.energy_levels_cm = energies

        self.header = 'wlvac(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        self.string_format = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        self.rule = 136*'-'
        strings = [] 

        for jj in range(0,len(wavelengths)):
            lumo = self.lumo_ph[jj]
            flux = self.lumo_erg[jj]
            upper = self.pec_levels[jj][0]+1
            lower = self.pec_levels[jj][1]+1
                    
            upper_j   = self.angular_momenta[upper-1]
            lower_j   = self.angular_momenta[lower-1]
            upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
            lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)

            if len(lower_csf) < 14:
                lower_csf = (14-len(lower_csf))*' ' + lower_csf
            if len(upper_csf) < 14:
                upper_csf = (14-len(upper_csf))*' ' + upper_csf

            upper_wav = self.energy_levels_cm[upper-1]
            lower_wav = self.energy_levels_cm[lower-1]
            avalue = self.avalues[upper-1,lower-1]
            wavel = self.wl_vac_nm[jj]
            this_pec = self.pec[jj][0]
            strings.append(self.string_format.format(wavel,
                                                               lower,
                                                               upper,
                                                               lower_wav,
                                                               lower_csf,
                                                               upper_wav,
                                                               upper_csf,
                                                               avalue,
                                                               this_pec,
                                                               lumo,
                                                               flux))
            self.strings = strings

    def display(self):
        print(self.rule)
        print(self.header)
        print(self.rule)

        for string in self.strings:
            print(string)
        print(self.rule)


def combine_requested_lines(list_of_requested_lines:list[requested_lines]) -> requested_lines:
    #possibly to combine the stronget lines from many datasets...
    #for the table output would need to keep track of a number of things..


    return 


class colradlumo_calc:
    def __init__(self,adf04_path:str,density:float,temp:float,atomic_mass_number:float):
        self.adf04_path = adf04_path
        self.density = density
        self.temp = temp 
        self.atomic_mass_number = atomic_mass_number
        self.atomic_mass_kg = atomic_mass_number * NUCLEON_MASS_KG

        met = np.array([0])#never change this unless you know what you're doing...
        density = np.array([density])
        temp = np.array([temp])
        t = process_time()

        #this is the default option, i 
        #just coded it this way as a failsafe.
        norm_pops_for_pecs = False

        colradpy_run = colradpy(adf04_path,met,temp,density,use_ionization=False,suppliment_with_ecip=False,use_recombination_three_body=False,use_recombination=False,rate_interp_col='cubic',default_pop_norm=norm_pops_for_pecs)
        colradpy_run.make_electron_excitation_rates()
        colradpy_run.make_ioniz_from_reduced_ionizrates()
        colradpy_run.suppliment_with_ecip()
        colradpy_run.populate_cr_matrix()
        colradpy_run.solve_quasi_static() 

        self.colradpy_class = colradpy_run

        self.pec_levels =  colradpy_run.data['processed']['pec_levels']

        t = process_time() - t 
        print('ColRadPy cpu time (sec) - {:7.2f}'.format(t))
        
        pec =  colradpy_run.data['processed']['pecs'][:,met,0,0]
        
        #don't ever change this unless you want to give me a headache
        #if we didn't already normalise pops by the paritition function (which is the default option)
        #then do it.
        if not norm_pops_for_pecs:
            pops = colradpy_run.data['processed']['pops_no_norm'][:,met,0,0]
            #ground =1 , so not in array. add it on.
            sum_pops = 1 + np.sum(pops)
            pec /= sum_pops

        wl_vac_nm = colradpy_run.data['processed']['wave_vac']
        wl_air_nm = colradpy_run.data['processed']['wave_air']
        self.pec = pec 

        self.wl_vac_nm  = wl_vac_nm 
        self.wl_vac_ang = wl_vac_nm*10
        self.wl_air_nm  = wl_air_nm 
        self.wl_air_ang = wl_air_nm*10        

        num_ions_in_a_solar_mass = SOLAR_MASS_KG / self.atomic_mass_kg

        #we already normalised the pec
        self.luminosity_photons_per_solar_mass = (pec * density  * num_ions_in_a_solar_mass).flatten()

        #CONVERT  NM TO CM
        self.photon_energies_ergs = HC_CGS / (wl_vac_nm*1e-7) .flatten()

        self.luminosity_ergs_per_solar_mass = self.luminosity_photons_per_solar_mass * self.photon_energies_ergs

        #atomic data
        self.avalues = colradpy_run.data['cr_matrix']['A_ji']
        self.energy_levels_cm = colradpy_run.data['atomic']['energy']
        self.csfs_labels = colradpy_run.data['atomic']['config']
        #self.terms = colradpy_run.data['atomic']['nist_conf_form']
        self.angular_momenta = colradpy_run.data['atomic']['w']
        #self.orbital_ang = 
        self.statistical_weights = 2*colradpy_run.data['atomic']['w'] + 1

        self.csf_id = colradpy_run.data['atomic']['id_groups']

        self.scaled_lumo_photos = [] 
        self.scaled_lumo_ergs = []

    def predict_mass_for_requested_lumo_wl(self,wavelength_requested_nm,spectral_line_requested_ergs_s):
        print('Requesting L = {:7.5E} ph/s for spectral line λ = {:11.4f} nm'.format(spectral_line_requested_ergs_s,wavelength_requested_nm) )
        
        index = np.argmin(np.abs(self.wl_vac_nm - wavelength_requested_nm))
        required_mass = spectral_line_requested_ergs_s / self.luminosity_ergs_per_solar_mass[index]
        print('Closest wavelength found: λ  = {:11.4f} nm'.format(self.wl_vac_nm[index]))
        print('Luminosity in one solar mass: {:7.5E} ergs/s , {:7.5E} ph/s'.format(self.luminosity_ergs_per_solar_mass[index],self.luminosity_photons_per_solar_mass[index]))
        print('Require ion-mass of {:11.4f} M_solar for requested luminosity. '.format(required_mass))



    def scale_lumo_by_ion_mass(self,mass_of_ion_solar_units:float):
        self.scaled_lumo_photos = self.luminosity_photons_per_solar_mass * mass_of_ion_solar_units
        self.scaled_lumo_ergs = self.luminosity_ergs_per_solar_mass * mass_of_ion_solar_units

    def select_strongest_n_lines(self,n_select:int) -> requested_lines:
        #arguments for requested_lines_class
        #self,wavelengths,avalues,pecs,pec_levels,lumo_ph,lumo_erg,csf_labels,angular_momenta,energies

        #selecting max n
        arguments =  np.argpartition(self.scaled_lumo_ergs, -n_select)[-n_select:]

        #sorting from strongest to weakest of this subset
        arguments = arguments[np.argsort(self.scaled_lumo_ergs[arguments])][::-1]

        #if the user did not scale by a mass yet, just set it to one solar mass
        if len(self.scaled_lumo_ergs) == 0:
            print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
            self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)

        #putting all this in the class. 
        requestlines = requested_lines(self.wl_vac_nm[arguments],
                                       self.avalues,
                                       self.pec[arguments],
                                       self.pec_levels[arguments],
                                       self.scaled_lumo_photos[arguments],
                                       self.scaled_lumo_ergs[arguments],
                                       self.csfs_labels,
                                       self.angular_momenta,
                                       self.energy_levels_cm)
        return requestlines


    def display_requested_lines_array(self,requested_lines:float,wl_tol:float):
        header = 'wlvac(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        
        #header_string = ' {:8},     {:2} - {:2}, {:14},  {:1,}, {:14}, {:10}, {:10}, {:10}, {:10}, {:10}'
        #header_string = header_string.format(
        #                     'wlvac(nm)',
        #                     'i',
        #                     'j',
        #                     'E_j (cm-1)',
        #                     'Level j',
        #                     'E_i (cm-1)'                                                          
        #                     'Level i', 
        #                     'Level i', 
        #                     'pec cm^3/s',    
        #                     'L (10^50 ph/s)',  
        #                     'L (10^38 erg/s)'                   
        #                     )

        string = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        #
        ##if the user hasn't scaled the lumo yet:
        if len(self.scaled_lumo_ergs) == 0:
            print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
            self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)

        strings_to_be_printed = []

        for ii in range(0,len(requested_lines)):
            pec_ind = np.where(np.abs(self.wl_vac_nm - requested_lines[ii]) < wl_tol)

            if len(pec_ind[0]) > 0:
                if len(pec_ind[0]) > 1:
                    print(72*'-')
                    print('duplicate wavelength for requested lune at {} found - displaying all.'.format(requested_lines[ii]))
                    print(72*'-')

                for jj in range(0,len(pec_ind[0])):

                    pec_ind_this = pec_ind[0][jj]
                    lumo = self.scaled_lumo_photos[pec_ind_this]
                    flux = self.scaled_lumo_ergs[pec_ind_this]
                    upper = self.pec_levels[pec_ind_this][0]+1
                    lower = self.pec_levels[pec_ind_this][1]+1
                    
                    upper_j   = self.angular_momenta[upper-1]
                    lower_j   = self.angular_momenta[lower-1]
                    upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
                    lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)

                    if len(lower_csf) < 14:
                        lower_csf = (14-len(lower_csf))*' ' + lower_csf
                    if len(upper_csf) < 14:
                        upper_csf = (14-len(upper_csf))*' ' + upper_csf

                    upper_wav = self.energy_levels_cm[upper-1]
                    lower_wav = self.energy_levels_cm[lower-1]
                    avalue = self.avalues[upper-1,lower-1]
                    wavel = self.wl_vac_nm[pec_ind_this]
                    this_pec = self.pec[pec_ind_this][0]
                    strings_to_be_printed.append(string.format(wavel,
                                                               lower,
                                                               upper,
                                                               lower_wav,
                                                               lower_csf,
                                                               upper_wav,
                                                               upper_csf,
                                                               avalue,
                                                               this_pec,
                                                               lumo,
                                                               flux))
            else:
                print("wavelength {:11.4f} not found - skipping".format(requested_lines[ii])) 
        
        
        print(136*'-')

        print(header)
        print(136*'-')

        for string in strings_to_be_printed:
            print(string)
        print(136*'-')

    def line_broadening_lumo_density(self,velocity_frac_speed_light,wavelength_array):
        broadened_spec = np.zeros_like(wavelength_array) 
        for (index,wavelength) in enumerate(self.wl_air_nm.flatten()):
            broadened_spec += self.scaled_lumo_ergs[index] * gaussian_kernel_lumo_density(wavelength,velocity_frac_speed_light,wavelength_array)
        

        return broadened_spec

def gaussian_kernel_lumo_density(central_wavelength_nm,beta,wavelength_range_nm):
    #returns a gaussian Luminosity density (in units of erg s^{-1} ang^{-1}) with unit luminosity.
    #i.e: an integral over this gaussian returns 1 erg/s

    #so, multiplying this Kernel by the desired lumoninosity gaurantees an integal over the line gives desired lumo,
    expon = (central_wavelength_nm - wavelength_range_nm)/central_wavelength_nm
    
    #for this simple case, i assume we integrate over the whole Gaussian. 

    #lumo in ergs/s/ang
    norm_factor = 1.0 / (SQRT_TWOPI * beta * central_wavelength_nm*10)
    
    expon/= np.sqrt(2)*beta 

    return norm_factor * np.exp (-np.power(expon,2) ) 

def trap(wl,lumo):
    integral = 0
    for ii in range(0,len(wl)-1,1):
        integral += 0.5 * (wl[ii+1]-wl[ii]) * (lumo[ii]+lumo[ii+1])
    return integral