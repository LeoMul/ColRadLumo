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
        self. avalue_sum_check = np.inf
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

        self.header = '    wlvac(nm),    wlair(nm),  transition,     E_j (cm-1),          Level j,      E_i (cm-1),          Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        self.string_format = ' {:12.4f}, {:12.4f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        
        self.header_with_element_code = ' ElemCod, wlvac(nm), wlair(nm),  transition,     E_j (cm-1),         Level j,     E_i (cm-1),        Level i, A_ij(s^-1), pec cm^3/s,   L (ph/s),  L (erg/s)'
        self.string_format_with_element_code = '{:8}, {:12.4f}, {:12.4f},     {:2} - {:2}, {:14.3f},  {:14}, {:14.3f}, {:10}, {:10.2E}, {:10.2E}, {:10.2E}, {:10.2E}'
        
        
        if beta != 0:
            self.string_format = self.string_format + ', {:8.2f}'
        self.rule = 180*'-'
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