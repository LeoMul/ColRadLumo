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
SQRT_PI = 1.7724538509055159

#Planck's constant times c in cgs e.g: erg cm s 
HC_CGS = 6.63e-34 * 3e8 * 1e7 * 1e2

C_CGS = 3.00E10
ELECTRON_MASS_G = 9.11E-28
ELECTRIC_CHARGE_CGS_ESU = 4.8e-10
SOBOLEV_CONST = np.pi * ELECTRIC_CHARGE_CGS_ESU * ELECTRIC_CHARGE_CGS_ESU\
                /(ELECTRON_MASS_G * C_CGS)


OSCILLATOR_CONST = 6.6702E15

class requested_lines:
    def __init__(
                 self,
                 wavelengths,
                 wavelengths_air,
                 avalues,
                 pecs,
                 pec_levels,
                 lumo_ph,
                 lumo_erg,
                 csf_labels,
                 angular_momenta,
                 energies,
                 mass=0.0,
                 temp=0.0,
                 dens=0.0,
                 elementcode='',
                 beta = 0.0
                 ):

        self.wl_vac_nm = wavelengths
        self.avalues     = avalues

        self.pec_levels = pec_levels
        self.pec        = pecs
        self.lumo_ph    = lumo_ph
        self.lumo_erg   = lumo_erg

        self.csfs_labels = csf_labels
        self.angular_momenta = angular_momenta
        self.energy_levels_cm = energies
        self.wl_value = []
        self.wl_values_air = []
        self.temp = temp 
        self.mass = mass 
        self.dens = dens
        self.element_code = elementcode

        self.header = 'wlvac(nm), wlair(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        self.string_format = ' {:12.4f}, {:12.4f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        
        self.header_with_element_code = ' ElemCod, wlvac(nm), wlair(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        self.string_format_with_element_code = '{:8}, {:12.4f}, {:12.4f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        
        
        if beta != 0:
            self.string_format = self.string_format + ', {:8.2f}'
        self.rule = 136*'-'
        strings = [] 
        strings_with_codes = []
        self.labels =[]
        for jj in range(0,len(wavelengths)):
            lumo = self.lumo_ph[jj]
            flux = self.lumo_erg[jj]
            upper = self.pec_levels[jj][0]+1
            lower = self.pec_levels[jj][1]+1
                    
            upper_j   = self.angular_momenta[upper-1]
            lower_j   = self.angular_momenta[lower-1]
            upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
            lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)

            if len(lower_csf) < 16:
                lower_csf = (16-len(lower_csf))*' ' + lower_csf
            if len(upper_csf) < 16:
                upper_csf = (16-len(upper_csf))*' ' + upper_csf

            self.labels.append(r'{} $\to$ {}'.format(upper_csf,lower_csf))

            upper_wav = self.energy_levels_cm[upper-1]
            lower_wav = self.energy_levels_cm[lower-1]
            avalue = self.avalues[upper-1,lower-1]
            wavel = self.wl_vac_nm[jj]
            self.wl_value.append(wavel)
            self.wl_values_air.append(wavelengths_air[jj])
            this_pec = self.pec[jj]
            this_fwhm = 2.355 * beta * wavel / np.sqrt(2)
            if beta == 0 :
                strings.append(self.string_format.format(wavel,
                                                        wavelengths_air[jj],
                                                        lower,
                                                        upper,
                                                        lower_wav,
                                                        lower_csf,
                                                        upper_wav,
                                                        upper_csf,
                                                        avalue,
                                                        this_pec,
                                                        lumo,
                                                        flux
                                                        )
                                                        )
            else:
                strings.append(self.string_format.format(wavel,
                                                        wavelengths_air[jj],
                                                        lower,
                                                        upper,
                                                        lower_wav,
                                                        lower_csf,
                                                        upper_wav,
                                                        upper_csf,
                                                        avalue,
                                                        this_pec,
                                                        lumo,
                                                        flux,
                                                        this_fwhm
                                                        )
                                                        )
            strings_with_codes.append(self.string_format_with_element_code.format(self.element_code,
                                                     wavel,
                                                     wavelengths_air[jj],
                                                     lower,
                                                     upper,
                                                     lower_wav,
                                                     lower_csf,
                                                     upper_wav,
                                                     upper_csf,
                                                     avalue,
                                                     this_pec,
                                                     lumo,
                                                     flux
                                                     )
                                                     )

            
            self.strings = strings
            self.strins_with_codes = strings_with_codes
    def display(self):
        print(self.rule)
        print(self.header)
        print(self.rule)

        for string in self.strings:
            print(string)
        print(self.rule)

    def write(self,filename='',write_head = True):
        #temp_K = temp_eV*11600
        
        if filename == '': 
            filename = 'line_lumo_pec_{}_{:5.3e}Msun_{:5.0e}eV_{:5.3e}cm-3.dat'.format(
                self.element_code,
                self.mass,
                self.temp,
                self.dens)
        


        f=  open(filename,'w')
        #f.write('T_e = {:5.0f} K\n'.format(temp_K))
        #f.write('n_e = {:5.3e} cm^(-3)\n'.format(dens))
        #f.write('M_Te2+ = {:5.3e} M_sun\n'.format(mass))

        #f.write(self.rule+'\n')
        if write_head:
            f.write(self.header_with_element_code+'\n')
        #f.write(self.rule+'\n')

        for string in self.strins_with_codes:
            f.write(string+'\n')
        #f.write(self.rule+'\n')
        f.close()




def combine_requested_lines(list_of_requested_lines:list[requested_lines]) -> requested_lines:
    #possibly to combine the stronget lines from many datasets...
    #for the table output would need to keep track of a number of things..


    return 


class colradlumo_calc:
    def __init__(self,adf04_path:str,density,temp,atomic_mass_number:float):


        self.adf04_path = adf04_path
        self.density = density
        self.temp = temp 
        #print(density)
        self.atomic_mass_number = atomic_mass_number
        self.atomic_mass_kg = atomic_mass_number * NUCLEON_MASS_KG
        met = np.array([0])#never change this unless you know what you're doing...
        #density = np.array([density])
        #temp = np.array([temp])
        t = process_time()

        #this is the default option, i 
        #just coded it this way as a failsafe.
        norm_pops_for_pecs = False

        num_temps = len(temp)
        num_dens = len(density)
        


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
        
        pec =  colradpy_run.data['processed']['pecs'][:,0,:,:]

        pec[pec<0] = 0.0 #1e-30
        #don't ever change this unless you want to give me a headache
        #if we didn't already normalise pops by the paritition function (which is the default option)
        #then do it.
        #I'm essentially overriding anywhere where the user can mess with the code.

        
        pops = colradpy_run.data['processed']['pops_no_norm'][:,0,:,:]
        
        #need to do a sum here over the right axis
        sum_pops = 1.0 + np.sum(pops,axis=0)

        self.sum_pops = sum_pops
        self.pops_normed = pops / sum_pops
        if not norm_pops_for_pecs:
            #ground =1 , so not in array. add it on.
            pec /= sum_pops



        wl_vac_nm = colradpy_run.data['processed']['wave_vac']
        wl_air_nm = colradpy_run.data['processed']['wave_air']
        self.pec = pec 
        #print(pec[pec<0])

        self.wl_vac_nm  = wl_vac_nm 
        self.wl_vac_ang = wl_vac_nm*10
        self.wl_air_nm  = wl_air_nm 
        self.wl_air_ang = wl_air_nm*10        
        self.num_dens = num_dens
        self.num_temps = num_temps
        self.num_wl = len(wl_air_nm)
        self.escape_prob = np.ones_like(self.wl_vac_nm)

        num_ions_in_a_solar_mass = SOLAR_MASS_KG / self.atomic_mass_kg

        self.num_ions_in_a_solar_mass = num_ions_in_a_solar_mass

        lumo_per_ion = np.zeros([len(wl_air_nm),num_temps,num_dens])
        #print('dim check:')
        #print(np.shape(pec))
        #print(np.shape(lumo_per_ion))

        for ii in range(0,num_dens):
            lumo_per_ion[:,:,ii] = pec[:,:,ii] * density[ii]

        lumo_per_ion = lumo_per_ion * num_ions_in_a_solar_mass


        #we already normalised the pec
        self.luminosity_photons_per_solar_mass = lumo_per_ion.copy()

        #CONVERT  NM TO CM
        self.photon_energies_ergs = HC_CGS / (wl_vac_nm*1e-7) #.flatten()
        self.luminosity_ergs_per_solar_mass = self.luminosity_photons_per_solar_mass.copy()

        print(np.shape(self.photon_energies_ergs))
        print(self.photon_energies_ergs[0])
        for ii in range(0,len(wl_air_nm)):
            self.luminosity_ergs_per_solar_mass[ii,:,:] = self.luminosity_photons_per_solar_mass[ii,:,:] * self.photon_energies_ergs[ii]

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

    def predict_mass_for_requested_lumo_wl(self,wavelength_requested_nm,spectral_line_requested_ergs_s,temp,density):
        print('Requesting L = {:7.5E} ph/s for spectral line λ = {:11.4f} nm'.format(spectral_line_requested_ergs_s,wavelength_requested_nm) )
        
        index = np.argmin(np.abs(self.wl_vac_nm - wavelength_requested_nm))
        
        temp_index = np.argmin(np.abs(self.temp - temp))
        dens_index = np.argmin(np.abs(self.density - density))
        print('using temp: {:10.2f}eV ~ {:10.2f} K'.format(self.temp[temp_index],self.temp[temp_index]*11600))
        print('using dens: {:10.2e} cm-3'.format(self.density[dens_index]))

        required_mass = spectral_line_requested_ergs_s / self.luminosity_ergs_per_solar_mass[index,temp_index,dens_index]
        print('Closest wavelength found: λ  = {:11.4f} nm'.format(self.wl_vac_nm[index]))
        print('Luminosity in one solar mass: {:7.5E} ergs/s , {:7.5E} ph/s'.format(self.luminosity_ergs_per_solar_mass[index,temp_index,dens_index],self.luminosity_photons_per_solar_mass[index,temp_index,dens_index]))
        print('Require ion-mass of {:11.4f} M_solar for requested luminosity. '.format(required_mass))



    def scale_lumo_by_ion_mass(self,mass_of_ion_solar_units):

        self.scaled_lumo_photos = np.zeros([
            self.num_wl,
            self.num_temps,
            self.num_dens,
            len(mass_of_ion_solar_units)
        ])
        self.mass = mass_of_ion_solar_units
        self.scaled_lumo_ergs = self.scaled_lumo_photos.copy()

        for ii in range(0,len(mass_of_ion_solar_units)):
            self.scaled_lumo_photos[:,:,:,ii] = self.luminosity_photons_per_solar_mass * mass_of_ion_solar_units[ii]
            self.scaled_lumo_ergs[:,:,:,ii] = self.luminosity_ergs_per_solar_mass * mass_of_ion_solar_units[ii]

    def select_strongest_n_lines(self,n_select:int,temperature,density,mass,beta=0.0) -> requested_lines:
        #arguments for requested_lines_class
        #self,wavelengths,avalues,pecs,pec_levels,lumo_ph,lumo_erg,csf_labels,angular_momenta,energies

        #selecting max n
        #print(np.shape(self.scaled_lumo_ergs))

        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        mass_index = np.argmin(np.abs(self.mass - mass))
        arguments =  np.argpartition(self.scaled_lumo_ergs[:,temp_index,density_index,mass_index].flatten(), -n_select)[-n_select:]
        #sorting from strongest to weakest of this subset
        print('Using parameters: ')
        print('Te = {} eV'.format(self.temp[temp_index]))
        print('Ne = {} cm-3'.format(self.density[density_index]))
        print('M  = {} Msun'.format(self.mass[mass_index]))
        array_to_be_sorted = self.scaled_lumo_ergs[:,temp_index,density_index].flatten()
        array_to_be_sorted = self.wl_air_nm[:].flatten()

        arguments = arguments[np.argsort(array_to_be_sorted[arguments])]#[::-1]
        #if the user did not scale by a mass yet, just set it to one solar mass
        if len(self.scaled_lumo_ergs) == 0:
            print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
            self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)

        #putting all this in the class. 

        element_code = self.colradpy_class.data['atomic']['element'].replace(' ', '')+str(self.colradpy_class.data['atomic']['charge_state'])+'+'

        requestlines = requested_lines(self.wl_vac_nm[arguments],
                                       self.wl_air_nm[arguments],
                                       self.avalues,
                                       self.pec[arguments,temp_index,density_index].flatten(),
                                       self.pec_levels[arguments],
                                       self.scaled_lumo_photos[arguments,temp_index,density_index,mass_index].flatten(),
                                       self.scaled_lumo_ergs[arguments,temp_index,density_index,mass_index].flatten(),
                                       self.csfs_labels,
                                       self.angular_momenta,
                                       self.energy_levels_cm,
                                       self.mass[mass_index],
                                       self.temp[temp_index],
                                       self.density[density_index],
                                       element_code,
                                       beta
                                       )
        return requestlines

    def optical_depth(self,temperature,density,time_exp_days,velocity_c,mass_solarunits,ion_string='',printing=False):
        
        #select the densities closest to the user selectrion.
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))

        #for getting the populations
        populations = np.zeros(len(self.pops_normed[:,temp_index,density_index])+1)

        #ColRadPy's pops_no_norm as ground=1 and not included in the arry, include this and normalise.
        populations[1:] = self.colradpy_class.data['processed']['pops_no_norm'][:,0,temp_index,density_index]#self.pops_normed[:,temp_index,density_index]
        populations[0] = 1.0 #/ self.sum_pops[temp_index,density_index]
        populations/= np.sum(populations)
        
        #Setting the Sob density array
        num_lines = len(self.pec_levels)
        optical_depth = np.ones(num_lines)* SOBOLEV_CONST#Sobolev const = e^2 / (m_e c)
        optical_depth_test = np.zeros(num_lines)
        time_exp_cgs = time_exp_days * 24.0 * 3600.0
        volume = 4.0 * np.pi * (C_CGS* time_exp_cgs*velocity_c)**3 / 3.0
        num_ions = mass_solarunits * self.num_ions_in_a_solar_mass
        eps = 1.0 / (4.0 * np.pi)

        header = '      wlvac(nm),  transition,     E_j (cm-1),         Level j,      Pop(j)     E_i (cm-1),        Level i,      Pop(j),  A_ij(s^-1),   SobDepth'


        if printing:
            print('Requested density: {:10.2e}'.format(density))
            print('Using density:     {:10.2e}'.format(self.density[density_index] ))
            print('Requested temp:    {:10.2f}'.format(temperature))
            print('Using temp:        {:10.2f}'.format(self.temp[temp_index] ))
            print('Volume:            {:10.2e}'.format(volume))
            print('Num ions:          {:10.2e}'.format(num_ions))
            print('Ion Density:       {:10.2e}'.format(num_ions/volume))
        fvalue_array = np.zeros(num_lines)
        strings = []
        string_format = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:10.2E}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        #print(header)
        for ii in range(0, num_lines):
            
            upper = self.pec_levels[ii][0]+1
            lower = self.pec_levels[ii][1]+1
            correction = 1.0 - self.statistical_weights[lower-1]*populations[upper-1]\
                                /self.statistical_weights[upper-1]*populations[lower-1]
            
            #Convert A -value to f value
            f_value = self.avalues[upper-1,lower-1] * self.statistical_weights[upper-1]\
                        /self.statistical_weights[lower-1] * self.wl_vac_ang[ii]**2
            f_value /= OSCILLATOR_CONST
            fvalue_array[ii] = f_value

           # print(f_value)
            #calculate optical depth - put wavelength in MKS as so is the Sob factor... time_exp in secs?

            factor = f_value * self.wl_vac_nm[ii] * 1e-7 * time_exp_cgs * populations[lower-1] * num_ions * correction / volume
            optical_depth[ii] *= factor
            optical_depth_test[ii] = (self.wl_air_nm[ii]*1e-7)**3 * eps*self.avalues[upper-1,lower-1] * 0.5 * populations[lower-1]*(self.statistical_weights[upper-1]/self.statistical_weights[lower-1] - populations[upper-1]/populations[lower-1])*time_exp_cgs *num_ions/volume
            #print(optical_depth[ii],optical_depth_test[ii])
            upper_j   = self.angular_momenta[upper-1]
            lower_j   = self.angular_momenta[lower-1]
            upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
            lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)
            if len(lower_csf) < 15:
                lower_csf = (15-len(lower_csf))*' ' + lower_csf
            if len(upper_csf) < 15:
                upper_csf = (15-len(upper_csf))*' ' + upper_csf
            upper_wav = self.energy_levels_cm[upper-1]
            lower_wav = self.energy_levels_cm[lower-1]
            wavel = self.wl_vac_nm[ii]

            strings.append(string_format.format(          wavel,
                                                               lower,
                                                               upper,
                                                               lower_wav,
                                                               lower_csf,
                                                               populations[lower-1],
                                                               upper_wav,
                                                               upper_csf,
                                                               populations[upper-1],
                                                               self.avalues[upper-1,lower-1],
                                                               optical_depth[ii],
                                                               (1.0 - np.exp(-optical_depth[ii])) / optical_depth[ii]
                                                        ))
        print(header)
        for ii in range(0,len(optical_depth)):
                upper = self.pec_levels[ii][0]
                lower = self.pec_levels[ii][1]
                print('{:6}'.format(ion_string),strings[ii])

                #if ii ==0 :
                #    #print(fvalue_array[ii])
                #if self.avalues[upper,lower] > 1e4:
                ##print(lower,populations[lower-1])
                #    #print(fvalue_array[ii])
                #    print('{:6}'.format(ion_string),strings[ii])
                #    break
        self.escape_prob = (1.0 - np.exp(-optical_depth)) / optical_depth
        self.sobolov = optical_depth

    #def display_requested_lines_array(self,requested_lines:float,wl_tol:float):
    #    header = 'wlvac(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
    #    
    #    #header_string = ' {:8},     {:2} - {:2}, {:14},  {:1,}, {:14}, {:10}, {:10}, {:10}, {:10}, {:10}'
    #    #header_string = header_string.format(
    #    #                     'wlvac(nm)',
    #    #                     'i',
    #    #                     'j',
    #    #                     'E_j (cm-1)',
    #    #                     'Level j',
    #    #                     'E_i (cm-1)'                                                          
    #    #                     'Level i', 
    #    #                     'Level i', 
    #    #                     'pec cm^3/s',    
    #    #                     'L (10^50 ph/s)',  
    #    #                     'L (10^38 erg/s)'                   
    #    #                     )
#
    #    string = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
    #    #
    #    ##if the user hasn't scaled the lumo yet:
    #    if len(self.scaled_lumo_ergs) == 0:
    #        print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
    #        self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)
#
    #    strings_to_be_printed = []
#
    #    for ii in range(0,len(requested_lines)):
    #        pec_ind = np.where(np.abs(self.wl_vac_nm - requested_lines[ii]) < wl_tol)
#
    #        if len(pec_ind[0]) > 0:
    #            if len(pec_ind[0]) > 1:
    #                print(72*'-')
    #                print('duplicate wavelength for requested lune at {} found - displaying all.'.format(requested_lines[ii]))
    #                print(72*'-')
#
    #            for jj in range(0,len(pec_ind[0])):
#
    #                pec_ind_this = pec_ind[0][jj]
    #                lumo = self.scaled_lumo_photos[pec_ind_this]
    #                flux = self.scaled_lumo_ergs[pec_ind_this]
    #                upper = self.pec_levels[pec_ind_this][0]+1
    #                lower = self.pec_levels[pec_ind_this][1]+1
    #                
    #                upper_j   = self.angular_momenta[upper-1]
    #                lower_j   = self.angular_momenta[lower-1]
    #                upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
    #                lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)
#
    #                if len(lower_csf) < 14:
    #                    lower_csf = (14-len(lower_csf))*' ' + lower_csf
    #                if len(upper_csf) < 14:
    #                    upper_csf = (14-len(upper_csf))*' ' + upper_csf
#
    #                upper_wav = self.energy_levels_cm[upper-1]
    #                lower_wav = self.energy_levels_cm[lower-1]
    #                avalue = self.avalues[upper-1,lower-1]
    #                wavel = self.wl_vac_nm[pec_ind_this]
    #                this_pec = self.pec[pec_ind_this][0]
    #                strings_to_be_printed.append(string.format(wavel,
    #                                                           lower,
    #                                                           upper,
    #                                                           lower_wav,
    #                                                           lower_csf,
    #                                                           upper_wav,
    #                                                           upper_csf,
    #                                                           avalue,
    #                                                           this_pec,
    #                                                           lumo,
    #                                                           flux))
    #        else:
    #            print("wavelength {:11.4f} not found - skipping".format(requested_lines[ii])) 
    #    
    #    
    #    print(136*'-')
#
    #    print(header)
    #    print(136*'-')
#
    #    for string in strings_to_be_printed:
    #        print(string)
    #    print(136*'-')


    #def write_requested_lines_array(self,requested_lines:float,wl_tol:float,file:str):
    #    header = 'wlvac(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
    #    
    #    #header_string = ' {:8},     {:2} - {:2}, {:14},  {:1,}, {:14}, {:10}, {:10}, {:10}, {:10}, {:10}'
    #    #header_string = header_string.format(
    #    #                     'wlvac(nm)',
    #    #                     'i',
    #    #                     'j',
    #    #                     'E_j (cm-1)',
    #    #                     'Level j',
    #    #                     'E_i (cm-1)'                                                          
    #    #                     'Level i', 
    #    #                     'Level i', 
    #    #                     'pec cm^3/s',    
    #    #                     'L (10^50 ph/s)',  
    #    #                     'L (10^38 erg/s)'                   
    #    #                     )
#
    #    string = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
    #    #
    #    f = open(file,'w')
    #    ##if the user hasn't scaled the lumo yet:
    #    if len(self.scaled_lumo_ergs) == 0:
    #        print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
    #        self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)
#
    #    strings_to_be_printed = []
#
    #    for ii in range(0,len(requested_lines)):
    #        pec_ind = np.where(np.abs(self.wl_vac_nm - requested_lines[ii]) < wl_tol)
#
    #        if len(pec_ind[0]) > 0:
    #            if len(pec_ind[0]) > 1:
    #                print(72*'-')
    #                print('duplicate wavelength for requested lune at {} found - displaying all.'.format(requested_lines[ii]))
    #                print(72*'-')
#
    #            for jj in range(0,len(pec_ind[0])):
#
    #                pec_ind_this = pec_ind[0][jj]
    #                lumo = self.scaled_lumo_photos[pec_ind_this]
    #                flux = self.scaled_lumo_ergs[pec_ind_this]
    #                upper = self.pec_levels[pec_ind_this][0]+1
    #                lower = self.pec_levels[pec_ind_this][1]+1
    #                
    #                upper_j   = self.angular_momenta[upper-1]
    #                lower_j   = self.angular_momenta[lower-1]
    #                upper_csf = self.csfs_labels[upper-1]+'{:5}'.format(upper_j)
    #                lower_csf = self.csfs_labels[lower-1]+'{:5}'.format(lower_j)
#
    #                if len(lower_csf) < 14:
    #                    lower_csf = (14-len(lower_csf))*' ' + lower_csf
    #                if len(upper_csf) < 14:
    #                    upper_csf = (14-len(upper_csf))*' ' + upper_csf
#
    #                upper_wav = self.energy_levels_cm[upper-1]
    #                lower_wav = self.energy_levels_cm[lower-1]
    #                avalue = self.avalues[upper-1,lower-1]
    #                wavel = self.wl_vac_nm[pec_ind_this]
    #                this_pec = self.pec[pec_ind_this][0]
    #                strings_to_be_printed.append(string.format(wavel,
    #                                                           lower,
    #                                                           upper,
    #                                                           lower_wav,
    #                                                           lower_csf,
    #                                                           upper_wav,
    #                                                           upper_csf,
    #                                                           avalue,
    #                                                           this_pec,
    #                                                           lumo,
    #                                                           flux))
    #        else:
    #            f.write("wavelength {:11.4f} not found - skipping\n".format(requested_lines[ii])) 
    #    
    #    
    #    f.write(136*'-'+'\n')
#
    #    f.write(header+'\n')
    #    f.write(136*'-'+'\n')
#
    #    for string in strings_to_be_printed:
    #        f.write(string+'\n')
    #    f.write(136*'-'+'\n')


    def line_broadening_lumo_density(self,velocity_frac_speed_light,wavelength_array,temperature,density):
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        #print(density_index,temp_index)
        broadened_spec = np.zeros_like(wavelength_array) 
        for (index,wavelength) in enumerate(self.wl_air_nm.flatten()):
            #print(index)
            if (wavelength != 0.0) and (wavelength<1000*wavelength_array[-1]):
                broadened_spec += self.scaled_lumo_ergs[index,temp_index,density_index] * gaussian_kernel_lumo_density(wavelength,velocity_frac_speed_light,wavelength_array)
                #print(self.scaled_lumo_ergs[0],wavelength)

        return broadened_spec
    
    def line_broadening_lumo_density_fwhm(self,velocity_fwhm_cunits,wavelength_array,temperature,density):
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        #print(density_index,temp_index)
        broadened_spec = np.zeros_like(wavelength_array) 
        for (index,wavelength) in enumerate(self.wl_air_nm.flatten()):
            #print(index)
            if (wavelength != 0.0) and (wavelength<1000*wavelength_array[-1]):
                broadened_spec += self.scaled_lumo_ergs[index,temp_index,density_index] * gaussian_kernel_lumo_density_new(wavelength,velocity_fwhm_cunits,wavelength_array)
                #print(self.scaled_lumo_ergs[0],wavelength)

        return broadened_spec

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