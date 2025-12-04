path_to_colradpy = 'ColRadPy/'
import sys
sys.path.append(path_to_colradpy)
#imports
from atomic_masses import *
from colradpy import colradpy 
from time import process_time 
import numpy as np
from constants import * 
from requested_lines import * 
from gaussian import * 

class colradlumo_calc:
    def __init__(self,adf04_path:str,
                 density,
                 temp,
                 atomic_mass_number:float = 0.0,
                 use_Opacity=False):


        self.adf04_path = adf04_path
        self.density = density
        self.temp = temp 
        self.getElementCode()
        
        atomic_number = getAtomicNumber(self.element)
        print(atomic_number)
        
        if atomic_mass_number == 0.0:
            self.atomic_mass_number = amu[atomic_number]
            atomic_mass_number = amu[atomic_number]
            print(atomic_number,self.atomic_mass_number)
        else:
            self.atomic_mass_number = atomic_mass_number
            
        self.atomic_mass_kg = self.atomic_mass_number * NUCLEON_MASS_KG

        
        met = np.array([0])#never change this unless you know what you're doing...
        #density = np.array([density])
        #temp = np.array([temp])
        t = process_time()

        #this is the default option, i 
        #just coded it this way as a failsafe.
        norm_pops_for_pecs = False
        self.norm_pops_for_pecs = norm_pops_for_pecs
        num_temps = len(temp)
        num_dens = len(density)
        


        colradpy_run = colradpy(adf04_path,met,temp,density,use_ionization=False,suppliment_with_ecip=False,use_recombination_three_body=False,use_recombination=False,rate_interp_col='cubic',default_pop_norm=norm_pops_for_pecs)
        colradpy_run.make_electron_excitation_rates()
        colradpy_run.make_ioniz_from_reduced_ionizrates()
        colradpy_run.suppliment_with_ecip()
        colradpy_run.populate_cr_matrix()
        colradpy_run.solve_quasi_static() 
        
        #is this redundant? refactor necessary.
        #self.colradpy_class = colradpy_run
        self.colradpy_run = colradpy_run

        
        
        self.pec_levels =  colradpy_run.data['processed']['pec_levels']

        t = process_time() - t 
        print('ColRadPy cpu time (sec) - {:7.2f}'.format(t))
        

        #don't ever change this unless you want to give me a headache
        #if we didn't already normalise pops by the paritition function (which is the default option)
        #then do it.
        #I'm essentially overriding anywhere where the user can mess with the code.
        self.normalisePops()
        
        wl_vac_nm = colradpy_run.data['processed']['wave_vac']
        wl_air_nm = colradpy_run.data['processed']['wave_air']
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
        #CONVERT  NM TO CM
        self.photon_energies_ergs = HC_CGS / (wl_vac_nm*1e-7) #.flatten()

        self.pdiff=0.0
        self.init_lumo()        

        #atomic data
        self.avalues = colradpy_run.data['cr_matrix']['A_ji']
        #print('check: ',np.shape(self.avalues))
        self.energy_levels_cm = colradpy_run.data['atomic']['energy']
        
        self.wl_vac_ang_matrix = np.zeros_like(self.avalues)
        nlev = len(self.energy_levels_cm)
        for ii in range(0,nlev):
            for jj in range(0,nlev):
                if (ii != jj):
                    self.wl_vac_ang_matrix[ii,jj] =abs ( 1e8 / (self.energy_levels_cm[ii] - self.energy_levels_cm[jj]))
        
        
        
        self.csfs_labels = colradpy_run.data['atomic']['config']
        #self.terms = colradpy_run.data['atomic']['nist_conf_form']
        self.angular_momenta = colradpy_run.data['atomic']['w']
        #self.orbital_ang = 
        self.statistical_weights = 2*colradpy_run.data['atomic']['w'] + 1

        self.csf_id = colradpy_run.data['atomic']['id_groups']

        self.scaled_lumo_photos = [] 
        self.scaled_lumo_ergs = []
        
        aval_save = colradpy_run.data['rates']['a_val']
        self.aval_save = aval_save
        self.avalue_save = self.avalues.copy()
        
        self.converged = False
        self.suma_old = np.sum(aval_save)
        self.suma_new = 0.0

        
        
        
        iter = 0 
        self.tracker  = np.zeros(4)
        self.esc_old = np.ones_like(aval_save)
        
                

    def init_lumo(self):
        lumo_per_ion = np.zeros([len(self.wl_air_nm),self.num_temps,self.num_dens])
        #print('INIT LUMO',self.pec[0,0,0])
        for ii in range(0,self.num_dens):
            lumo_per_ion[:,:,ii] = self.pec[:,:,ii] * self.density[ii]
        lumo_per_ion = lumo_per_ion * self.num_ions_in_a_solar_mass


        #we already normalised the pec
        self.luminosity_photons_per_solar_mass = lumo_per_ion.copy()
        self.luminosity_ergs_per_solar_mass = self.luminosity_photons_per_solar_mass.copy()

        for ii in range(0,len(self.wl_air_nm)):
            self.luminosity_ergs_per_solar_mass[ii,:,:] = self.luminosity_photons_per_solar_mass[ii,:,:] * self.photon_energies_ergs[ii]

    def getElementCode(self):
        #gets elements from the adf04 file
        try:
            f = open(self.adf04_path,'r')
            firstLine = f.readline()
            f.close()

            #print(firstLine)

            #firstLine = firstLine.replace(' ','')

            element = firstLine[0:2]
            charge  = firstLine[3:5]

            self.element = element
            self.charge_plus  = int(charge)
            self.charge_roman = int_to_roman(int(charge)+1)
            print(element,self.charge_roman)

        except:
            print('failed to get element code - debug please :)')
        


    def normalisePops(self):
        pec =  self.colradpy_run.data['processed']['pecs'][:,0,:,:]
        pec[pec<0] = 0.0 #1e-30
        
        pops = self.colradpy_run.data['processed']['pops_no_norm'][:,0,:,:]
        #need to do a sum here over the right axis
        sum_pops = 1.0 + np.sum(pops,axis=0)
        self.sum_pops = sum_pops
        self.pops_normed = pops / sum_pops
        if not self.norm_pops_for_pecs:
            #ground =1 , so not in array. add it on.
            pec /= sum_pops
            
        self.pec = pec 
        #print('check lpm',pec[0,0,0])
                
    def convergeOpacity(self,
                        desired_temp,
                        desired_density,
                        time_exp_days,
                        velocity_c,
                        mass_solar,
                        oscillator_breaker=False,
                        debug_printing = False,
                        fraction_override = 0.0,
                        relaxation_steps=9999):
        ns =  min(relaxation_steps,MAX_OPACITY_ITER)
        self.escape_prob = np.ones_like(self.aval_save)
        self.allowed = self.aval_save > 1e2
        self.forbidden = np.invert(self.allowed)
        for ii in range(0,ns):
            
            self.opacityIteration(
                         desired_temp,
                         desired_density,
                         time_exp_days,
                         velocity_c,
                         mass_solar,
                         ii,
                         oscillator_breaker,
                         debug_printing,
                         fraction_override
                         )

            if self.converged:
                print('Exiting  opacity self consitency run at iteration {:3}.'.format(ii))
                print('Population diff: {:8.2e}'.format(self.pdiff))
                print('Avalue diff: {:8.2e}'.format(self.avalue_sum_check))

                break
        
        if (ii == ns-1):
            print('warning - opacity calculation may not be converged for {} {}'.format(self.element,self.charge_roman))
        self.init_lumo()
        self.scale_lumo_by_ion_mass(np.array([mass_solar]))
        
    def opacityIteration(self,
                         desired_temp,
                         desired_density,
                         time_exp_days,
                         velocity_c,
                         mass_solar,
                         iter,
                         oscillator_breaker,
                         debug_printing = False,
                         fraction_override=0.0
                         ):
        self.escape_prob_old = self.escape_prob
        mean1  = np.mean(self.escape_prob_old)
        self.optical_depth(desired_temp,
                           desired_density,
                           time_exp_days,
                           velocity_c,
                           mass_solar,
                           printing=False,
                           fraction_override=fraction_override)
        p1 = self.pops_normed
        self.escape_prob = 0.5 * (self.escape_prob + self.escape_prob_old) #take average for stability
        print(72*'-')
        form1= '{:10}{:10}{:10}{:10}{:10}'
        print("Iter:      ",form1.format('Min','Max','Average','Median','StDev'))
        form = '{:10.3E}{:10.3E}{:10.3E}{:10.3E}{:10.3E}'
        print('Full:     ',form.format(np.min(self.escape_prob),np.max(self.escape_prob),np.mean(self.escape_prob),np.median(self.escape_prob),np.std(self.escape_prob)))
        print('Allowed:  ',form.format(np.min(self.escape_prob[self.allowed]),np.max(self.escape_prob[self.allowed]),np.mean(self.escape_prob[self.allowed]),np.median(self.escape_prob[self.allowed]),np.std(self.escape_prob[self.allowed])))
        print('Forbidden:',form.format(np.min(self.escape_prob[self.forbidden]),np.max(self.escape_prob[self.forbidden]),np.mean(self.escape_prob[self.forbidden]),np.median(self.escape_prob[self.forbidden]),np.std(self.escape_prob[self.forbidden])))
        self.beta_allowed_mean   = np.mean(self.escape_prob[self.allowed])
        self.beta_allowed_geomean = np.exp(np.log(self.escape_prob[self.allowed]).mean())
        self.beta_allowed_median = np.median(self.escape_prob[self.allowed])
        self.beta_allowed_stdev  = np.std(self.escape_prob[self.allowed])

        self.colradpy_run.data['rates']['a_val'] = self.aval_save * self.escape_prob
        self.colradpy_run.populate_cr_matrix()
        self.colradpy_run.solve_quasi_static() 
        self.normalisePops()
        p2 = self.pops_normed 
        self.pdiff = np.sum(np.power(p1-p2,2))
        
        self.avalues = self.colradpy_run.data['cr_matrix']['A_ji']
        sumavalues = np.sum(self.avalues)
        self.suma_new = sumavalues
        
        
        if debug_printing:
            print('population check: ',self.pdiff)
            print('{:5}, {:10.2e}'.format(iter,np.sum(self.avalues)))
        
        tracker = self.tracker
        tracker[iter%4] = sumavalues


        self.esc_old = self.escape_prob
        self.avalue_sum_check = abs( self.suma_old / self.suma_new -1.0)
        if ( self.avalue_sum_check < 0.001):
            self.converged = True
        
        #if (self.pdiff < 1e-3):
        #    self.converged = True
        
        self.suma_old = self.suma_new

        return 0

                
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

        return index,required_mass


    def scale_lumo_by_ion_mass(self,mass_of_ion_solar_units):
        #print(mass_of_ion_solar_units)
        self.scaled_lumo_photos = np.zeros([
            self.num_wl,
            self.num_temps,
            self.num_dens,
            len(mass_of_ion_solar_units)
        ])
        #print('hello',np.shape(self.scaled_lumo_photos))
        self.mass = mass_of_ion_solar_units
        self.scaled_lumo_ergs = self.scaled_lumo_photos.copy()

        for ii in range(0,len(mass_of_ion_solar_units)):
            self.scaled_lumo_photos[:,:,:,ii] = self.luminosity_photons_per_solar_mass * mass_of_ion_solar_units[ii]
            self.scaled_lumo_ergs[:,:,:,ii] = self.luminosity_ergs_per_solar_mass * mass_of_ion_solar_units[ii]

    def get_param_indices(self,temperature,density,mass):
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        mass_index = np.argmin(np.abs(self.mass - mass))
        return density_index,temp_index,mass_index 
    
    def return_lines(self,arguments,temp_index,density_index,mass_index,beta = 0.0): 
        element_code = self.colradpy_run.data['atomic']['element'].replace(' ', '')+str(self.colradpy_run.data['atomic']['charge_state'])+'+'
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


    def select_closest_lines(self,selection,temperature,density,mass ,beta=0.0 ) -> requested_lines:
        if len(self.scaled_lumo_ergs) == 0:
            print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
            self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)
        arguments = []
        for ii in range(0,len(selection)):
            check = np.argmin( np.abs(self.wl_vac_nm - selection[ii] ))
            arguments.append(check) 
        arguments = np.array(arguments)
        density_index,temp_index,mass_index = self.get_param_indices(temperature,density,mass)
        return self.return_lines(arguments,temp_index,density_index,mass_index,beta)

    def select_strongest_n_lines(self,n_select:int,temperature,density,mass,beta=0.0) -> requested_lines:
        density_index,temp_index,mass_index = self.get_param_indices(temperature,density,mass)
        if len(self.scaled_lumo_ergs) == 0:
            print('scale_lumo_by_ion_mass() hasnt been called - scaling to unit solar mass')
            self.scale_lumo_by_ion_mass(mass_of_ion_solar_units=1.0)
        array_to_be_sorted = self.scaled_lumo_ergs[:,temp_index,density_index,mass_index].flatten()
        #array_to_be_sorted = self.wl_air_nm[:].flatten()
        arguments = np.argsort(array_to_be_sorted)[::-1] [0:n_select]
                    
        return self.return_lines(arguments,temp_index,density_index,mass_index,beta)



    def optical_depth(self,temperature,density,time_exp_days,velocity_c,mass_solarunits,ion_string='',printing=False,fraction_override=0.0):
        
        #select the densities closest to the user selectrion.
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))

        #for getting the populations
        populations = np.zeros(len(self.pops_normed[:,temp_index,density_index])+1)

        #ColRadPy's pops_no_norm as ground=1 and not included in the arry, include this and normalise.
        populations[1:] = self.colradpy_run.data['processed']['pops_no_norm'][:,0,temp_index,density_index]#self.pops_normed[:,temp_index,density_index]
        populations[0] = 1.0 #/ self.sum_pops[temp_index,density_index]
        populations/= np.sum(populations)
        
        #Setting the Sob density array
        num_lines = len(self.colradpy_run.data['rates']['excit']['col_transitions'])
        #print('optical depth nl:',num_lines)
        optical_depth = np.ones(num_lines)* SOBOLEV_CONST#Sobolev const = e^2 / (m_e c)
        optical_depth_test = np.zeros(num_lines)
        
        time_exp_cgs = time_exp_days * 24.0 * 3600.0
        volume = 4.0 * np.pi * (C_CGS* time_exp_cgs*velocity_c)**3 / 3.0
        
        num_ions = mass_solarunits * self.num_ions_in_a_solar_mass
        eps = 1.0 / (4.0 * np.pi)

        header = '       wlvac(nm),  transition,     E_j (cm-1),         Level j,           Pop(j),     E_i (cm-1),        Level i,           Pop(j), A_ij(s^-1),   SobDepth,     SobEsc'
        idd = 0 
        if volume != 0:
            idd = num_ions/volume 
        if fraction_override != 0.0 : 
            idd = fraction_override * self.density[density_index]
        if printing:
            print('Requested density: {:10.2e}'.format(density))
            print('Using density:     {:10.2e}'.format(self.density[density_index] ))
            print('Requested temp:    {:10.2f}'.format(temperature))
            print('Using temp:        {:10.2f}'.format(self.temp[temp_index] ))
            print('Volume:            {:10.2e}'.format(volume))
            print('Num ions:          {:10.2e}'.format(num_ions))
            print('Ion Density:       {:10.2e}'.format(idd))
        fvalue_array = np.zeros(num_lines)
        strings = []
        string_format = ' {:8.2f},     {:2} - {:2}, {:14.3f},  {:14}, {:10.2E}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        #print(header)
        self.escape_prob = np.ones_like(optical_depth)
        #ttt = np.argwhere(optical_depth > 1e-7)
        #self.escape_prob[ttt] = (1.0 - np.exp(-optical_depth[ttt])) / optical_depth[ttt]
        
        for ii in range(0, num_lines):
            
            upper = self.colradpy_run.data['rates']['excit']['col_transitions'][ii][0]
            lower = self.colradpy_run.data['rates']['excit']['col_transitions'][ii][1]
            correction = 1.0 - self.statistical_weights[lower-1]*populations[upper-1]\
                                /self.statistical_weights[upper-1]*populations[lower-1]
            
            #Convert A -value to f value
            f_value = self.avalues[upper-1,lower-1] * self.statistical_weights[upper-1]\
                        * self.wl_vac_ang_matrix[upper-1,lower-1]**2 / self.statistical_weights[lower-1]
            f_value /= OSCILLATOR_CONST
            fvalue_array[ii] = f_value

           # print(f_value)
            #calculate optical depth - put wavelength in MKS as so is the Sob factor... time_exp in secs?

            
            if fraction_override != 0.0:
                ion_density = populations[lower-1] * fraction_override * self.density[density_index]
            else:
                ion_density = populations[lower-1] * num_ions / volume 
            factor = f_value * self.wl_vac_ang_matrix[upper-1,lower-1] * 1e-8 * time_exp_cgs * ion_density * correction 
            optical_depth[ii] *= factor
            
            wlcm = self.wl_vac_ang_matrix[upper-1,lower-1] * 1e-8 
            eightpi = 8.0 * np.pi 
            test = self.avalues[upper-1,lower-1] * (wlcm**3) * time_exp_cgs * ion_density * self.statistical_weights[upper-1] * correction/ (eightpi*self.statistical_weights[lower-1])
            #print(test,optical_depth[ii])
            #optical_depth_test[ii] = (self.wl_air_nm[ii]*1e-7)**3 * eps*self.avalues[upper-1,lower-1] * 0.5 * populations[lower-1]*(self.statistical_weights[upper-1]/self.statistical_weights[lower-1] - populations[upper-1]/populations[lower-1])*time_exp_cgs *num_ions/volume
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
            wavel = self.wl_vac_ang_matrix[upper-1,lower-1] * 0.1 
            
            if optical_depth[ii] > 1e-7:
                self.escape_prob[ii] = (1.0 - np.exp(-optical_depth[ii])) / optical_depth[ii]

            

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
                                                               self.escape_prob[ii]
                                                        ))
            

        self.sobolov = optical_depth

        if printing:
            print(header)
            for ii in range(0,1000):
                    upper = self.colradpy_run.data['rates']['excit']['col_transitions'][ii][0]
                    lower = self.colradpy_run.data['rates']['excit']['col_transitions'][ii][1]
                    print('{:6}'.format(ion_string),strings[ii])

    def line_broadening_lumo_density(self,velocity_frac_speed_light,wavelength_array,temperature,density,mass):
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        mass_index = np.argmin(np.abs(self.mass - mass))
        #print(density_index,temp_index)
        broadened_spec = np.zeros_like(wavelength_array) 
        for (index,wavelength) in enumerate(self.wl_air_nm.flatten()):
            #print(index)
            if (wavelength != 0.0) and (wavelength<1000*wavelength_array[-1]):
                broadened_spec += self.scaled_lumo_ergs[index,temp_index,density_index,mass_index] * gaussian_kernel_lumo_density(wavelength,velocity_frac_speed_light,wavelength_array)
                #print(self.scaled_lumo_ergs[0],wavelength)

        return broadened_spec
    
    def line_broadening_lumo_density_fwhm(self,velocity_fwhm_cunits,wavelength_array,temperature,density):
        density_index = np.argmin(np.abs(self.density - density))
        temp_index = np.argmin(np.abs(self.temp - temperature))
        #print(density_index,temp_index)
        broadened_spec = np.zeros_like(wavelength_array) 
        #print(temp_index,density_index)

        for (index,wavelength) in enumerate(self.wl_air_nm.flatten()):
            if (wavelength != 0.0) and (wavelength<1000*wavelength_array[-1]):
                broadened_spec += self.scaled_lumo_ergs[index,temp_index,density_index] * gaussian_kernel_lumo_density_new(wavelength,velocity_fwhm_cunits,wavelength_array)
                #print(self.scaled_lumo_ergs[0],wavelength)

        return broadened_spec





