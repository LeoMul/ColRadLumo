from colradlumo_library import * 

#code to produce a spectrum, compared with at23vfi for now.

class spectrum:
    
    def __init__(self,
               central_temp:     float ,
               central_density:  float ,
               listOfadf04Files :list= [], 
               listOfVelocities: list = [],
               wlMax = 5000,
               wlMin = 500,
               wlNum = 1000,
               ionFraction = []
               ):
    
        #beta_22 = 0.04
        #beta_20 = 0.04
        #beta_45 = 0.07
        self.wlMin = wlMin 
        self.wlMax = wlMax 
     
        numSpec = len(listOfadf04Files)
        self.numSpec = numSpec
        self.ionFraction = ionFraction
        #assert (numSpec == len(listOfVelocities) )

        elementsChargeRoman = [] 
        elementsChargePlus = []
        nWavelengths=  wlNum
        self.nWavelengths = nWavelengths
        spectra = np.zeros([numSpec,nWavelengths])
        wavelengths = np.linspace(wlMin,wlMax,nWavelengths)
        colRadLumo_classes = []

        for ii,thisfile in enumerate(listOfadf04Files):
            thismodel = colradlumo_calc(thisfile,
                                        np.array([central_density]),
                                        np.array([central_temp]))
            thismodel.scale_lumo_by_ion_mass(np.array([1e0]))
            bs = thismodel.line_broadening_lumo_density_fwhm(
                listOfVelocities[ii],wavelengths,central_temp,central_density
            )
            spectra[ii,:] = bs
            elementsChargeRoman.append(thismodel.element + ' '+thismodel.charge_roman)
             #massLabel = fr"~M = {massStringVoodoo[:4]} $\times 10^{{{exponent}}} \mathrm{{ M}}_\odot$"

            elementChargePlusString = fr'{thismodel.element:>2}$^{{{thismodel.charge_plus}+}}$'
            
            elementsChargePlus.append(elementChargePlusString)

            colRadLumo_classes.append(thismodel)

        self.elementsChargeRoman = elementsChargeRoman 
        self.elementsChargePlus  = elementsChargePlus
        self.spectra = spectra 
        self.wavelengths = wavelengths 
        self.velocities = listOfVelocities 
        self.colRadLumo_classes = colRadLumo_classes
        self.continuum = np.zeros_like(wavelengths)
        
        #return elements,spectra,wavelengths,colRadLumo_classes





        
    def blackBodyAndAfterglow(self, bbT,agA,agB):
        wavelength_range = self.wavelengths
        T = bbT 
        a = agA
        b = agB 
        h_si = 6.63e-34 
        c_si = 3e8
        k_si = 1.38e-23
        blackbody = 2.0 * h_si * c_si * c_si* np.power(wavelength_range*1e-9,-5) * np.power(np.exp(h_si*c_si / (wavelength_range*1e-9*k_si*T))-1,-1)
        blackbody_cgs = 10.0 * blackbody
        R_photo_cgs = 6e15 
        blackbody_cgs_lumo = blackbody_cgs * R_photo_cgs*R_photo_cgs * 4.0 * np.pi * np.pi
        #units are erg cm-1 s-1
        #to convert to a-1
        blackbody_cgs_lumo = blackbody_cgs_lumo * 1e-8

        continuum = a*np.power(wavelength_range,b)*1e35 +blackbody_cgs_lumo   
        self.continuum = continuum
        #return continuum 
        

    def plot(self,
             dpi = 200,
             xlimits = [-1,-1],
             ylimits = [-1,-1],
             colormap = 'tab20',
             legend=True,
             cumulative=True,
             integralThresh = 0.005,
             massThresh = 1e-10,
             printInterestingLines=False,
             interestingThresh = 1e37,
             plotVFI=True,
             yscale='linear'):
        '''
        This routine was refactored by Google Gemini 3, 26.03.26.
        '''
        
        import matplotlib.pyplot as plt 
        import numpy as np
        from matplotlib import rc
        
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size': 12})

        xlimitsPlot = [self.wlMin, self.wlMax]
        ylimitsPlot = [0, 1.2]
        if xlimits != [-1, -1]: 
            xlimitsPlot = xlimits
        if ylimits != [-1, -1]: 
            ylimitsPlot = ylimits
        if yscale == 'log': 
            ylimitsPlot[0] = 1e-4


        ratios = self.integrals / self.sumIntegrals
        selected_indices = np.where(ratios > integralThresh)[0]
        
        if len(selected_indices) == 0:
            return None, None

        #Handle color to plot.
        cmap = plt.get_cmap(colormap)
        colors = [cmap(i) for i in np.linspace(0, 1, len(selected_indices))]
        color_mapping = {idx: colors[i] for i, idx in enumerate(selected_indices)}
        
        fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
        fig.subplots_adjust(right=0.7) 
        ax.set_yscale(yscale)

        if plotVFI:
            vfiwl, vfilum = getVFI()
            ax.step(vfiwl, vfilum * 1e-35, 'k', linewidth=0.2, where='mid', alpha=0.3)

        # 4. PLOTTING ORDER (The Fix)
        # We sort selected indices in REVERSE order (e.g., [50, 42, 10, 2])
        # This ensures the 'fullest' cumulative sum is the bottom layer.
        loop_order = sorted(selected_indices, reverse=True)

        handle_dict = {}

        # 5. Plotting Loop
        for ii in loop_order:
            color = color_mapping[ii]
            
            # Mass Label Logic
            massStringVoodoo = '{:8.2e}'.format(self.masses[ii])
            exphack = int(massStringVoodoo[-1])    
            massLabel = fr"~M = {massStringVoodoo[:4]} $\times 10^{{-{exphack}}} \mathrm{{ M}}_\odot$"
            full_label = self.elementsChargePlus[ii] + massLabel

            if cumulative:
                # Plotting from high index to low index prevents overshadowing
                ax.fill_between(self.wavelengths, self.cumulative[ii, :] * 1e-35, 
                                color=color, step='mid', zorder=-ii)
            else:
                thisSpec = self.spectra[ii, :] * self.masses[ii] * 1e-35
                ax.fill_between(self.wavelengths, thisSpec, color=color, step='mid', alpha=0.8)

            # Legend Line Handle
            line_handle, = ax.plot([], [], color=color, label=full_label, linewidth=1.5)
            handle_dict[ii] = line_handle

        # Plot total sum line
        ax.plot(self.wavelengths, self.total * 1e-35, color='k', linewidth=1, label='Total')

        # 6. Apply Limits & Labels
        ax.set_xlim(xlimitsPlot)
        ax.set_ylim(ylimitsPlot)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel(r'Luminosity (10$^{35}$ erg s$^{-1}$ Å$^{-1}$)')

        # 7. Legend: Restore order of the elements array
        if legend:
            sorted_keys = sorted(handle_dict.keys()) # Numerical order of indices
            final_handles = [handle_dict[k] for k in sorted_keys]
            # Add Total at the very top
            total_h, total_l = ax.get_legend_handles_labels()
            t_idx = total_l.index('Total')
            
            ax.legend([total_h[t_idx]] + final_handles, 
                      ['Total'] + [h.get_label() for h in final_handles],
                      loc="center left", bbox_to_anchor=(1.05, 0.67), frameon=False)

        return fig, ax
    
    def totalSpectrum(self,listOfMasses):
        self.masses = np.array(listOfMasses)
        ns = np.shape(self.spectra)[0]
        nw = np.shape(self.spectra)[1]
        self.total = np.zeros(nw)
        self.integrals = np.zeros(ns)
        self.cumulative = np.zeros_like(self.spectra)
        
        for ii in range(0,ns):
            contribution = self.spectra[ii,:] * listOfMasses[ii]
            self.integrals[ii] += 0.1*np.trapz(contribution,self.wavelengths) #0.1 to convert units back to nm^-1

        #self.order = np.argsort(self.integrals)[::-1]

        for ii in range(0,ns):
            contribution = self.spectra[ii,:] * listOfMasses[ii]
            self.total[:] += contribution
            self.cumulative[ii,:] = contribution 
            #Calculate cumulatives - this could probably be made more efficient with Numpy 
            for kk in range(0,ii):
                self.cumulative[ii,:] += self.spectra[kk,:] * listOfMasses[kk]
                
        self.total += self.continuum
        self.sumIntegrals = np.sum(self.integrals)
        
    def plotgemini(self,
             dpi = 200,
             xlimits = [-1,-1],
             ylimits = [-1,-1],
             colormap = 'tab20',
             legend=True,
             cumulative=True,
             integralThresh = 0.005,
             massThresh = 1e-10,
             printInterestingLines=False,
             interestingThresh = 1e37,
             plotVFI=True,
             yscale='linear'):
        
        import matplotlib.pyplot as plt 
        import numpy as np
        from matplotlib import rc
        from matplotlib.patches import Patch # For custom legend proxy
        
        # 1. Formatting
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size': 12})

        # 2. Axis Limit Logic
        xlimitsPlot = [self.wlMin, self.wlMax]
        ylimitsPlot = [0, 1.2]
        if xlimits != [-1, -1]: xlimitsPlot = xlimits
        if ylimits != [-1, -1]: ylimitsPlot = ylimits
        if yscale == 'log': ylimitsPlot[0] = 1e-4

        # 3. Filtering Logic
        ratios = self.integrals / self.sumIntegrals
        selected_indices = np.where(ratios > integralThresh)[0]
        excess_indices = np.where(ratios <= integralThresh)[0]
        
        cmap = plt.get_cmap(colormap)
        colors = [cmap(i) for i in np.linspace(0, 1, len(selected_indices))]
        color_mapping = {idx: colors[i] for i, idx in enumerate(selected_indices)}
        
        fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
        fig.subplots_adjust(right=0.7) 
        ax.set_yscale(yscale)

        if plotVFI:
            vfiwl, vfilum = getVFI()
            ax.step(vfiwl, vfilum * 1e-35, 'k', linewidth=0.2, where='mid', alpha=0.3)

        # 4. Calculate the "Excess" Spectrum
        excess_spectrum = np.zeros_like(self.total)
        for idx in excess_indices:
            excess_spectrum += self.spectra[idx, :] * self.masses[idx]
        
        # 5. Plotting the Striped Excess Layer
        # 'hatch' creates the stripes. 'edgecolor' sets the stripe color.
        excess_label = "Other Species"
        excess_color = 'whitesmoke' # Very light grey
        stripe_color = 'dimgrey' # Medium grey for stripes
        
        if cumulative:
            ax.fill_between(self.wavelengths, excess_spectrum * 1e-35, 
                            facecolor=excess_color, hatch='///', 
                            edgecolor=stripe_color, linewidth=0, step='mid', alpha=0.5)
        else:
            ax.fill_between(self.wavelengths, excess_spectrum * 1e-35, 
                            facecolor=excess_color, hatch='///', 
                            edgecolor=stripe_color, linewidth=0, step='mid', alpha=0.3)
        
        # Custom legend handle to show the stripes in the legend box
        excess_handle = Patch(facecolor=excess_color, hatch='///', 
                              edgecolor=stripe_color, label=excess_label)

        # 6. Plotting Selected Spectra (Reverse Order for Cumulative)
        loop_order = sorted(selected_indices, reverse=True)
        handle_dict = {}

        for ii in loop_order:
            color = color_mapping[ii]
            
            # Mass Label Logic
            massStringVoodoo = '{:8.2e}'.format(self.masses[ii])
            exphack = int(massStringVoodoo[-1])    
            massLabel = fr"~M = {massStringVoodoo[:4]} $\times 10^{{-{exphack}}} \mathrm{{ M}}_\odot$"
            full_label = self.elementsChargePlus[ii] + massLabel

            if cumulative:
                # Stack interesting species on top of the 'excess' baseline
                plot_data = (self.cumulative[ii, :] + excess_spectrum) * 1e-35
                ax.fill_between(self.wavelengths, plot_data, color=color, step='mid')
            else:
                thisSpec = self.spectra[ii, :] * self.masses[ii] * 1e-35
                ax.fill_between(self.wavelengths, thisSpec, color=color, step='mid', alpha=0.8)

            # Store handle as a line for the legend (as requested)
            line_handle, = ax.plot([], [], color=color, label=full_label, linewidth=1.5)
            handle_dict[ii] = line_handle

        # Plot total sum line
        total_handle, = ax.plot(self.wavelengths, self.total * 1e-35, color='k', linewidth=1, label='Total')

        # 7. Limits & Legend
        ax.set_xlim(xlimitsPlot)
        ax.set_ylim(ylimitsPlot)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel(r'Luminosity (10$^{35}$ erg s$^{-1}$ Å$^{-1}$)')

        if legend:
            sorted_ii = sorted(handle_dict.keys())
            # Final Order: Total -> Elements (Original Array Order) -> Striped Excess
            final_handles = [total_handle] + [handle_dict[k] for k in sorted_ii] + [excess_handle]
            final_labels = [h.get_label() for h in final_handles]

            ax.legend(final_handles, final_labels, loc="center left", 
                      bbox_to_anchor=(1.05, 0.67), frameon=False, fontsize=10)

        return fig, ax

def getVFI():
    #thanks to James Gillanders for the data reduction
    #and additionally the plotting. 
    #Much of the code in this function is his. 
    redshift = 0.0647

    spectrum_gillanders_29d = np.loadtxt('AT2023vfi_JWST_29d_fluxcal.txt')
    spectrum_gillanders_29d[:,0]/= (1+redshift) 
    
    #spectral reduction
    pixel_size = 80
    spectrum_gillanders_29d = spectrum_gillanders_29d[(spectrum_gillanders_29d[:,0] < 30026.9 - pixel_size) | (spectrum_gillanders_29d[:,0] > 30026.9 + pixel_size)]
    spectrum_gillanders_29d = spectrum_gillanders_29d[(spectrum_gillanders_29d[:,0] < 22242.3 - pixel_size) | (spectrum_gillanders_29d[:,0] > 22242.3 + pixel_size)]
    spectrum_gillanders_29d = spectrum_gillanders_29d[(spectrum_gillanders_29d[:,0] < 22907.8 - pixel_size) | (spectrum_gillanders_29d[:,0] > 22907.8 + pixel_size)]

    #nm
    spectrum_gillanders_29d[:,0] *= 0.1 #nm
    
    
    vfi_wl_nm = spectrum_gillanders_29d 
    
    d_l = 302 * 3.086e+24
    
    vfi_lum_ergsPerSec = spectrum_gillanders_29d[:,1] * d_l * d_l * 4.0 * np.pi #* 1e-35 
    return vfi_wl_nm,vfi_lum_ergsPerSec

