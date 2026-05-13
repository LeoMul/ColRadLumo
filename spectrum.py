'''

code to produce a spectrum, compared with at23vfi for now.

todo, implement the fortran code I am calling for the CRM into a proper module.

'''

import numpy as np 
class spectrum:
    
    def __init__(self,
               central_temp:     float ,
               central_density:  float ,
               listOfadf04Files :list= [], 
               listOfVelocities: list = [],
               listOfIonMasses: list = [],
               wlMax = 5000,
               wlMin = 500,
               wlNum = 1000,
               ionFraction = [],
               useSobolov = False,
               timeSinceExplosionDays =0.0 
               ):
    
        self.wlMin = wlMin 
        self.wlMax = wlMax 
        self.masses = listOfIonMasses
        numSpec = len(listOfadf04Files)
        self.numSpec = numSpec
        self.ionFraction = ionFraction
        #assert (numSpec == len(listOfVelocities) )
        self.central_temp    = central_temp
        self.central_density = central_density
        self.useSobolov      = useSobolov
        self.timeSinceExplosionDays = timeSinceExplosionDays
        elementsChargeRoman = [] 
        elementsChargePlus = []
        nWavelengths=  wlNum
        self.nWavelengths = nWavelengths
        spectra = np.zeros([numSpec,nWavelengths])
        wavelengths = np.linspace(wlMin,wlMax,nWavelengths)
        colRadLumo_classes = []

        for ii,thisfile in enumerate(listOfadf04Files):
            
            #thismodel = colradlumo_calc(thisfile,
            #                            np.array([central_density]),
            #                            np.array([central_temp]))
            
            rv = self.runLeoCRM(thisfile,listOfIonMasses[ii],listOfVelocities[ii],fractionOverride=-1.0)
            element = ''
            charge_plus = 99
            import os.path
            failure = os.path.isfile('ifail')
            if (not failure) : 

                bs = np.loadtxt('spectrum')[:,1]
                spectra[ii,:] = bs 
                gg = open('chargequick','r')
                line = gg.readline().split()
                gg.close()
                charge_plus = int(line[1])
                element = line[0]
                #print('debug leo',element,charge_plus)
                gg.close()
            else:
                os.system('rm ifail')
            #print(rv)
            
            charge_roman = int_to_roman(charge_plus+1)

            #if useSobolov:
            #    thismodel.convergeOpacity(central_temp,central_density,29,0.1,1e-3,False,False,0.5)
            
            #thismodel.scale_lumo_by_ion_mass(np.array([1e0]))
            #bs = thismodel.line_broadening_lumo_density_fwhm(
            #    listOfVelocities[ii],wavelengths,central_temp,central_density
            #)
            

            
            elementsChargeRoman.append(element + ' '+charge_roman)
             #massLabel = fr"~M = {massStringVoodoo[:4]} $\times 10^{{{exponent}}} \mathrm{{ M}}_\odot$"

            elementChargePlusString = fr'{element:>2}$^{{{charge_plus}+}}$'
            
            elementsChargePlus.append(elementChargePlusString)

            colRadLumo_classes.append(None)

        self.elementsChargeRoman = elementsChargeRoman 
        self.elementsChargePlus  = elementsChargePlus
        self.spectra = spectra 
        self.wavelengths = wavelengths 
        self.velocities = listOfVelocities 
        self.colRadLumo_classes = colRadLumo_classes
        self.continuum = np.zeros_like(wavelengths)
        
        #return elements,spectra,wavelengths,colRadLumo_classes


    def runLeoCRM(self,file,ionMassSolar,velocity,fractionOverride=-1.0):
        import os 
        crm_input = open('crm_input','w')
        crm_input.write(" &crm_input \n")
        crm_input.write(f" temperature            = {self.central_temp}\n")
        crm_input.write(f" density                = {self.central_density}\n")
        crm_input.write(f" sobolov                = .{self.useSobolov}.\n")
        crm_input.write(f" adf04path              ='{file}' \n")
        crm_input.write(f" timeSinceExplosionDays = {self.timeSinceExplosionDays} \n")
        crm_input.write(f" velocityExpansionC     = {velocity}\n")
        crm_input.write(f" massElementSolar       = {ionMassSolar}\n")
        crm_input.write(f" fractionOverride       = {fractionOverride}\n")
        crm_input.write(f" wlmin_nm               = {self.wlMin}\n")
        crm_input.write(f" wlmax_nm               = {self.wlMax}\n")
        crm_input.write(f" numwl                  = {self.nWavelengths}\n")
        crm_input.write(" &end \n")
        crm_input.close()
        
        return os.system('/Users/leomulholland/myCRM/bin/crm') 



        
    def blackBodyAndAfterglow(self, bbT=660,agA=3e3,agB=-1.31):
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

    
    def totalSpectrum(self):
        #self.masses = np.array(listOfMasses)
        ns = np.shape(self.spectra)[0]
        nw = np.shape(self.spectra)[1]
        self.total = np.zeros(nw)
        self.integrals = np.zeros(ns)
        self.cumulative = np.zeros_like(self.spectra)
        
        for ii in range(0,ns):
            contribution = self.spectra[ii,:] 
            self.integrals[ii] += 0.1*np.trapz(contribution,self.wavelengths) #0.1 to convert units back to nm^-1

        #self.order = np.argsort(self.integrals)[::-1]

        for ii in range(0,ns):
            contribution = self.spectra[ii,:] 
            self.total[:] += contribution
            self.cumulative[ii,:] = contribution 
            #Calculate cumulatives - this could probably be made more efficient with Numpy 
            for kk in range(0,ii):
                self.cumulative[ii,:] += self.spectra[kk,:] 
                
        self.total += self.continuum
        self.sumIntegrals = np.sum(self.integrals)
        
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
             yscale='linear',
             firstPeakColorOverride=False):
        
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
        if xlimits != [-1, -1]: 
            xlimitsPlot = xlimits
        if ylimits != [-1, -1]: 
            ylimitsPlot = ylimits
        if yscale == 'log': 
            ylimitsPlot[0] = 1e-4

        # 3. Filtering Logic
        ratios = self.integrals / self.sumIntegrals
        selected_indices = np.where(ratios > integralThresh)[0]
        excess_indices = np.where(ratios <= integralThresh)[0]
        
        cmap = plt.get_cmap(colormap)
        colors = [cmap(i) for i in np.linspace(0, 1, len(selected_indices))]
        color_mapping = {idx: colors[i] for i, idx in enumerate(selected_indices)}
        
        fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
        ax.set_yscale(yscale)

        if plotVFI:
            vfiwl, vfilum = getVFI()
            ax.step(vfiwl, vfilum * 1e-35, 'k', linewidth=0.2, where='mid', alpha=0.3)

        # 4. Calculate the "Excess" Spectrum
        excess_spectrum = np.zeros_like(self.total)
        for idx in excess_indices:
            excess_spectrum += self.spectra[idx, :] 
        
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
            if firstPeakColorOverride:
                try:
                    color = mycolordict[self.elementsChargeRoman[ii]]
                except:
                    print('Color match failed - using default',self.elementsChargeRoman[ii])
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
                thisSpec = self.spectra[ii, :] * 1e-35
                ax.fill_between(self.wavelengths, thisSpec, color=color, step='mid', alpha=0.8)

            # Store handle as a line for the legend (as requested)
            line_handle, = ax.plot([], [], color=color, label=full_label, linewidth=1.5)
            handle_dict[ii] = line_handle

        # Plot total sum line
        total_handle, = ax.plot(self.wavelengths, (self.total+excess_spectrum) * 1e-35, color='k', linewidth=1, label='Total')

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
                      bbox_to_anchor=(1.05, 0.63), frameon=False, fontsize=10)
            
            densStringVoodoo = '{:8.2e}'.format(self.central_density)
            exphack = int(densStringVoodoo[-1])
            densLabel = fr"$n_e = $ {densStringVoodoo[:4]} $\times 10^{{{exphack}}}$ cm$^{{-3}}$"
            tempLabel = fr"$T_e = $ {self.central_temp:5.0f} K"

            ax.text(1.08, 0.19, densLabel, 
                    transform=ax.transAxes, 
                    verticalalignment='top', 
                    horizontalalignment='left',
                    fontsize=10)
            ax.text(1.08, 0.15, tempLabel, 
                    transform=ax.transAxes, 
                    verticalalignment='top', 
                    horizontalalignment='left',
                    fontsize=10)
                        
        return fig, ax

def getVFI():
    #thanks to James Gillanders for the data reduction
    #and additionally the plotting. 
    #Much of the code in this function is his. 
    redshift = 0.0647

    spectrum_gillanders_29d = np.loadtxt('/Users/leomulholland/colrad_lumo/AT2023vfi_JWST_29d_fluxcal.txt')
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

mycolordict = {
    "Ga I"  :  (0.59607843, 0.87450980, 0.54117647, 1.00000000),
    "Ge I"  :  (0.78039216, 0.78039216, 0.78039216, 1.00000000),
    "Ge II" :  (0.76862745, 0.61176471, 0.58039216, 1.00000000),
    "As I"  :  (1.00000000, 0.49803922, 0.05490196, 1.00000000),
    "As II" :  (0.58039216, 0.40392157, 0.74117647, 1.00000000),
    "As III":  (1.00000000, 0.59607843, 0.58823529, 1.00000000),
    "Se I"  :  (0.85882353, 0.85882353, 0.55294118, 1.00000000),
    "Se II" :  (0.83921569, 0.15294118, 0.15686275, 1.00000000),
    "Se III":  (0.73725490, 0.74117647, 0.13333333, 1.00000000),
    "Se IV" :  (1.00000000, 0.73333333, 0.47058824, 1.00000000),
    "Br I"  :  (0.96862745, 0.71372549, 0.82352941, 1.00000000),
    "Br II" :  (0.89019608, 0.46666667, 0.76078431, 1.00000000),
    "Br III":  (0.54901961, 0.33725490, 0.29411765, 1.00000000),
    "Br IV" :  (0.49803922, 0.49803922, 0.49803922, 1.00000000),
    "Br V"  :  (0.09019608, 0.74509804, 0.81176471, 1.00000000),
    "Kr II" :  (0.17254902, 0.62745098, 0.17254902, 1.00000000),
    "Kr III":  (0.68235294, 0.78039216, 0.90980392, 1.00000000),
    "Kr IV" :  (0.61960784, 0.85490196, 0.89803922, 1.00000000),
    "Kr V"  :  (0.12156863, 0.46666667, 0.70588235, 1.00000000)
}

ROMAN = [
    (1000, "M"),
    ( 900, "CM"),
    ( 500, "D"),
    ( 400, "CD"),
    ( 100, "C"),
    (  90, "XC"),
    (  50, "L"),
    (  40, "XL"),
    (  10, "X"),
    (   9, "IX"),
    (   5, "V"),
    (   4, "IV"),
    (   1, "I"),
]

def int_to_roman(number):
    result = []
    for (arabic, roman) in ROMAN:
        (factor, number) = divmod(number, arabic)
        result.append(roman * factor)
        if number == 0:
            break
    return "".join(result)


def spectrumFromNonThermalOutput(nonthermaloutputpath: str,
                                 adf04filedir:str ='/Users/leomulholland/SeKrPaper/quickCalcs/firstPeak',
                                 wlmin = 100,
                                 wlmax = 5000,
                                 nwavelengths=2000,
                                 useSobolov = False):
    
    header = np.loadtxt(nonthermaloutputpath,max_rows=1)
    thermalElectronTemperature = header[0]
    imposedElectronDensity     = header[1]
    actualElectronDensity      = header[2]
    timeSinceExplosionDays     = header[3]
    velocityExpansionC         = header[4]
    
    density = imposedElectronDensity 
    #if the user let pynt calculate an electron density, use that instead.
    if density == 0: 
        density = actualElectronDensity
    
    nonthermaloutput = np.loadtxt(nonthermaloutputpath,usecols=(2,-1),dtype=str,delimiter=',',skiprows=1)

    keys = nonthermaloutput[:,0]
    masses = nonthermaloutput[:,1].astype(float)
    files = []
    velocities = []

    for ii,key in enumerate(keys):

        strkey = key.replace(' ','')
        #print(strkey)

        file = f'{adf04filedir}/adf04{strkey}'
        files.append(file)
        velocities.append(velocityExpansionC)
    
    thisSpectrum = spectrum(
        thermalElectronTemperature,
        density,
        files,
        velocities,
        masses,
        wlMax=wlmax,
        wlMin=wlmin,
        wlNum=nwavelengths,
        useSobolov = useSobolov,
        timeSinceExplosionDays = timeSinceExplosionDays
    )
    
    return thisSpectrum