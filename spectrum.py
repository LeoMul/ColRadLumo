from colradlumo_library import * 

#code to produce a spectrum, compared with at23vfi for now.

class spectrum:
    
    def __init__(self,
               central_temp:     float ,
               central_density:  float ,
               listOfadf04Files :list= [], 
               listOfVelocities: list = []
               ):
    
        #beta_22 = 0.04
        #beta_20 = 0.04
        #beta_45 = 0.07
        numSpec = len(listOfadf04Files)
        self.numSpec = numSpec
        #assert (numSpec == len(listOfVelocities) )

        elements = [] 

        nWavelengths=  1000
        self.nWavelengths = nWavelengths
        spectra = np.zeros([numSpec,nWavelengths])
        wavelengths = np.linspace(500,6000,nWavelengths)
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
            elements.append(thismodel.element + ' '+thismodel.charge_roman)
            colRadLumo_classes.append(thismodel)

        self.elements = elements 
        self.spectra = spectra 
        self.wavelengths = wavelengths 
        self.velocities = listOfVelocities 
        self.colRadLumo_classes = colRadLumo_classes
        self.continuum = np.zeros_like(wavelengths)
        
        #return elements,spectra,wavelengths,colRadLumo_classes




    def totalSpectrum(self,listOfMasses):
        self.masses = listOfMasses
        ns = np.shape(self.spectra)[0]
        nw = np.shape(self.spectra)[1]
        self.total = np.zeros(nw)
        for ii in range(0,ns):
            self.total[:] += self.spectra[ii,:] * listOfMasses[ii]
        self.total += self.continuum
        
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
        
    def plot(self,dpi = 200,xlimits = [500,5000],ylimits = [0,1.2],colormap = 'rainbow'):
        import matplotlib.pyplot as plt 
        from matplotlib.ticker import FormatStrFormatter
        from matplotlib import rc
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size': 14})

        cmap = plt.get_cmap(colormap)
        colors = [cmap(i) for i in np.linspace(0, 1, self.numSpec+1)]
        
        fig, ax = plt.subplots(figsize=(7,5),dpi=dpi)
        
        vfiwl,vfilum  = getVFI()
        
        plt.step(vfiwl,vfilum * 1e-35,'k',linewidth=0.2)
        
        
        for ii in range(0,self.numSpec):
            
            # FOrmatting scientific notation in matplotlib with latex is voodoo.
            massStringVoodoo = '{:8.2e}'.format(self.masses[ii])
            exphack = int(massStringVoodoo[-1])
#
            if (self.masses[ii] <= 1e-10):
                #exponent requires two digits...
                import sys 
                sys.exit()
                print('masses wont be formatted properly - debug @ leomulholland@hotmail.com')
#            #as such, the following three lines were greatly aided by Gemini.
#
            exponent = f"-{exphack}"
            massLabel = fr",~M = {massStringVoodoo[:4]} $\times 10^{{{exponent}}} \mathrm{{ M}}_\odot$"
            
            thisSpec = self.spectra[ii,:] * self.masses[ii] * 1e-35

            ax.fill_between(self.wavelengths,thisSpec,color = colors[ii+1])
            #print(massMantessa,exphack,massLabel1.format(massMantessa))
            print(massLabel)
            #massLabel = massLabel1.format(massMantessa) + massLabel2 + massLabel3.format(exphack) + massLabel4
            
            plt.plot(-1,-1,
                     color = colors[ii+1],
                     label = self.elements[ii] + massLabel )
        
        plt.plot(self.wavelengths,self.total* 1e-35,color = colors[0])
        
        plt.xlim(xlimits)
        plt.ylim(ylimits)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel(r'Luminosity (10$^3$$^5$ erg s$^{-1}$ Å$^{-1}$)')
        plt.legend(loc = 'upper left')
        return fig,ax



def getVFI():
    #thanks to James Gillanders for the data reduction
    
    vfi = np.loadtxt('AT2023vfi_JWST_29d_fluxcal.txt')
    d_l = 302 * 3.086e+24
    redshift = 0.0647
    vfi_wl_nm = 0.1*vfi[:,0]/(1+redshift) 
    vfi_lum_ergsPerSec = vfi[:,1] * d_l * d_l * 4.0 * np.pi #* 1e-35 

    
    
    return vfi_wl_nm,vfi_lum_ergsPerSec

