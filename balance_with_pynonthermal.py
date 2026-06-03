import pynonthermal
import numpy as np 
from astropy import units as u
from astropy import constants 
import periodictable
from scipy.interpolate import interp1d
import sys
import os 
import json 
import argparse
#get_eff_ionpot
''''
Charge state balance for a given set of elements.

Uses the pynonthermal code of Luke Shingles to obtain non-thermal ionization rates 
for a nebular gas.

Uses your recombination data, but I would prefer you used mine.

Self consistent iteration - optional.

Current todo's: 

Need to fix the electorn density malarky - e.g what is ued for ionization balancr and what is used for SF.

Possibly extend some index's by one - I think we can squeeze out an extra 
ion for the ionization balance. Its probably fine the way it is. 

Print out the degradation spectrum in each case.
Print out the total electron density from the ionization balance

Need to check that all of the units and input degradation is correct. 



'''
TOLERANCE = 5.E-2
MAXITER   = 20 
     
class input:
    '''
    Class for input, read from json.
    '''
    def __init__(self,
                 listOfAtomicNumbers=[],
                 pathsOfRecombinationData=[],
                 massesOfElements=[],
                 thermalElectronTemperature = None,
                 imposedElectronDensitySF    = None,
                 imposedElectronDensityRecombination    = None,
                 velocityExpansionC = None,
                 timeSinceExplosionDays = None, 
                 selfConsistent = False,
                 averageAtomicMass  = 140,
                 maxIonization      = 6,
                 depositionOverride = None,
                 velocityMaxForEfficiency   = 0.1,
                 depositionMode = "KasenBarnes"
                 ):        

        #Boring transfer of memory. 
        self.listOfAtomicNumbers        = listOfAtomicNumbers
        self.massesOfElements           = massesOfElements
        self.pathsOfRecombinationData   = pathsOfRecombinationData
        self.thermalElectronTemperature = thermalElectronTemperature
        self.imposedElectronDensitySF   = imposedElectronDensitySF
        self.imposedElectronDensityRecombination = imposedElectronDensityRecombination
        self.velocityExpansionC         = velocityExpansionC
        self.timeSinceExplosionDays     = timeSinceExplosionDays
        self.averageAtomicMass          = averageAtomicMass
        self.maxIonization              = maxIonization
        self.selfConsistent             = selfConsistent
        self.depositionOverride         = depositionOverride
        self.velocityMaxForEfficiency   = velocityMaxForEfficiency
        self.depositionMode             = depositionMode.lower() 
        
class nonThermalBalance:
    '''
    Non thermal ionization balance class.
    '''
    
    def __init__(self,input: input,outfile_suffix = ''):
        
        #Transfer memory we might need. 
        self.input                      = input
        self.listOfAtomicNumbers        = input.listOfAtomicNumbers
        self.pathsOfRecombinationData   = input.pathsOfRecombinationData
        self.thermalElectronTemperature = input.thermalElectronTemperature
        self.imposedElectronDensitySF     = input.imposedElectronDensitySF
        self.velocityExpansionC         = input.velocityExpansionC
        self.timeSinceExplosionDays     = input.timeSinceExplosionDays
        self.massesOfElements           = input.massesOfElements
        self.averageAtomicMass          = input.averageAtomicMass
        self.maxIonization              = input.maxIonization
        self.depositionOverride         = input.depositionOverride
        self.velocityMaxForEfficiency   = input.velocityMaxForEfficiency
        self.numberOfElements = len(input.listOfAtomicNumbers)
        self.depositionMode             = input.depositionMode
        self.imposedElectronDensityRecombination = input.imposedElectronDensityRecombination
        self.outfile = open(f'pynt-balance-{outfile_suffix}.out','w')
        
        #sanity checks on input data. 
        assert( len(input.massesOfElements)         == self.numberOfElements)
        assert( len(input.pathsOfRecombinationData) == self.numberOfElements * input.maxIonization)
        self.electronDensity = self.imposedElectronDensitySF
        self._setRecombinationRates()
        self._setInitialIonizationBalance()        
        self._setDepositionRateDensity()
        
    def _setRecombinationRates(self):
        #Store an recombination,ionization rate for each element
        self.recombinationRatesCoefficient = np.zeros([self.input.maxIonization, self.numberOfElements])
        self.ionizationRates = np.zeros_like(self.recombinationRatesCoefficient)
        
        counter = 0 
        for ii in range(0,self.numberOfElements):
            for jj in range(0,self.input.maxIonization):
                axelrodrate = 10 * (jj+1)**2 * RRAxelrod(self.thermalElectronTemperature, jj+1)
                if self.pathsOfRecombinationData[counter] == None:
                    rate = axelrodrate
                    self.outfile.write(f'no recomb data  found for {self.listOfAtomicNumbers[ii],jj+1} - using Axelrod rate of {rate:10.2e}\n')
                else:
                    data = np.loadtxt(self.pathsOfRecombinationData[counter])
                    rate = interp1d(data[:,0],data[:,1])(self.thermalElectronTemperature)
                    self.outfile.write(f'recomb data  found for {self.listOfAtomicNumbers[ii],jj+1} - using rate of {rate:10.2e} = {rate/axelrodrate:10.2e} Axelrod\n')
                
                
                
                self.recombinationRatesCoefficient[jj,ii] = rate 
                counter += 1 
        
        return None 
    
    def _setInitialIonizationBalance(self):
        #Initial Balance per element.  Will be discard later.
        #Assuming an equal ionization. The code will eventually get something better.
        balancePerElement = np.zeros(self.input.maxIonization)
        balancePerElement[:] = 1.0 
        balancePerElement[2] = 2.0
        balancePerElement[3] = 2.0
        balancePerElement[:] /= balancePerElement.sum()
        
        self.initBalance  = np.zeros(self.input.maxIonization * self.numberOfElements)
        self.initFraction = np.zeros_like(self.initBalance) 
        #Need total matter density for injection later.
        self.elementNumberDensities = np.zeros(self.numberOfElements)
        self.elementMassDensities   = np.zeros(self.numberOfElements)

        self.nparticles = np.zeros(self.numberOfElements)

        self.expansion_volume = 0.0
        for ii in range(0,self.numberOfElements):
            self.elementNumberDensities[ii],self.elementMassDensities[ii],self.nparticles [ii],self.expansion_volume = elementDensity(
                self.listOfAtomicNumbers[ii],
                self.massesOfElements[ii],
                self.velocityExpansionC,
                self.timeSinceExplosionDays
            )
            
        self.elementNumberDensityTotal = self.elementNumberDensities.sum()
        self.elementMassDensityTotal   = self.elementMassDensities  .sum()

        stride = self.input.maxIonization
        for ii in range(0,self.numberOfElements):
            self.initBalance [ii * stride : (ii+1) * stride ] = balancePerElement * self.elementNumberDensities[ii]
            self.initFraction[ii * stride : (ii+1) * stride ] = balancePerElement
                
        self.outfile.write(f'Total number densities: {self.elementNumberDensityTotal:10.2e}\n',)
        self.outfile.write('Element densities:\n')
        for aa,zz in enumerate(self.listOfAtomicNumbers):
            self.outfile.write(f'{zz:3} {self.elementNumberDensities[aa]:10.2e}\n')
        
        self.balance = self.initBalance.copy()
        self.ionFraction = self.initFraction.copy()
        return None 
            
    def _setDepositionRateDensity(self):
        
        '''
        Set deposition rate - eV /s / cm^3,
            =  Qdot * total matter density 
        '''
        match self.depositionMode.lower():
            
            case "kasenbarnes":
                self._kasenBarnesDeposition()
                
            case "artisdata":
                self._artisDataDeposition()
                
            case "override":
                if self.depositionOverride is None:
                    import sys
                    print('Override requested - but no deposition override. ') 
                    sys.exit()
                else:
                    self.depositionratedensity_ev = self.depositionOverride 
                    self.outfile.write('Using an user-override deposition density of {:10.3}\n'.format(self.depositionOverride))
            case _:
                import sys
                print('Deposition Mode not set. Valid options: KasenBarnes, ArtisData, Override.') 
                sys.exit()
        return None 
    
    def _kasenBarnesDeposition(self):
        self.heatingRate = heatingKasenBarnes(self.timeSinceExplosionDays, self.averageAtomicMass)
        self.thermalizationEfficiency = thermalizationKasenBarnes(self.timeSinceExplosionDays,velocity_max_c=self.velocityMaxForEfficiency)
        self.outfile.write(f'efficiency   of  {self.thermalizationEfficiency}\n')
        self.outfile.write(f'using a heating rate of {self.heatingRate} \n')
        self.depositionratedensity_ev  = self.heatingRate * self.elementNumberDensityTotal * self.thermalizationEfficiency
        try:
            self.depositionratedensity_ev = self.depositionratedensity_ev.value 
        except:
            pass
        return None 
    def _artisDataDeposition(self):
        from pathlib import Path
        ROOT_DIR = Path(__file__).parent
        TEXT_FILE = ROOT_DIR /'artisdata.dat'
        artisdata = np.loadtxt(TEXT_FILE)
        logdays      = np.log10( artisdata[:, 0] )
        logdep_per_g = np.log10( artisdata[:,-1] )
        from scipy.interpolate import interp1d
        interp = interp1d(logdays,logdep_per_g)
        self.depPerGram = 10 ** interp(np.log10(self.timeSinceExplosionDays))
        print(self.depPerGram,self.elementMassDensityTotal)
        self.depositionratedensity_ev = self.depPerGram * self.elementMassDensityTotal
        return None 
    
    def setNewBalanceDamped(self):
        self.outfile.write('damping new ion balance...\n')
        self.damp = 0.2
        stride = self.maxIonization
        self.ionFraction = (1.0 - self.damp) * self.ionFraction + self.damp*self.ionFractionOld 
        #renormalize
        self.actualElectronDensity = 0.0 
        for aa,z in enumerate(self.listOfAtomicNumbers):
            thisfrac = self.ionFraction[aa * stride : (aa+1) * stride ] 
            thisnorm = thisfrac.sum()
            thisfrac/= thisnorm
            self.ionFraction[aa * stride : (aa+1) * stride ] = thisfrac
            self.balance[aa * stride : (aa+1) * stride ]  = thisfrac * self.elementNumberDensities[aa]
            self.actualElectronDensity += np.sum ( np.arange(0,self.maxIonization,1,dtype=int) * thisfrac * self.elementNumberDensities[aa])
        return None
    
    def ionIter(self):
        stride = self.maxIonization
        
        self.actualElectronDensity = 0.0 
    
        for aa,Z in enumerate(self.listOfAtomicNumbers):
            
            thisBalance = ionizationBalance(
                self.ionizationRates[:,aa], 
                self.electronDensityForIonization * self.recombinationRatesCoefficient[:,aa],Z
                )[0:self.maxIonization]
            
            self.balance[aa * stride : (aa+1) * stride ]     = thisBalance * self.elementNumberDensities[aa]
            
            self.actualElectronDensity += np.sum ( np.arange(0,self.maxIonization,1,dtype=int) * thisBalance * self.elementNumberDensities[aa])
            
            self.ionFraction[aa * stride : (aa+1) * stride ] = thisBalance 
        
        return None 
    
    def ionizationBalance(self):
        #Calculate a new Ionization balance.
        self.balanceOld = self.balance.copy()
        self.ionFractionOld = self.ionFraction.copy() 
        if self.imposedElectronDensityRecombination is not None:
            self.electronDensityForIonization = self.imposedElectronDensityRecombination
        else:
            self.electronDensityForIonization = self.electronDensity
        self.outfile.write(f'Calculating Ionization Balance, using initial electron density {self.electronDensityForIonization:10.2e}\n')

        self.ionIter()
        
        #if self.imposedElectronDensitySF is None:
        #    counter = 0
        #    for _ in range(0,MAXITER):
        #        conv = (self.electronDensityForIonization - self.actualElectronDensity) / self.electronDensityForIonization
        #        self.outfile.write(f'    {self.electronDensityForIonization:10.2e} {self.actualElectronDensity:10.2e} {conv:10.2e}\n')
#
        #        self.electronDensityForIonization = self.actualElectronDensity 
        #        if (abs(conv) < TOLERANCE):
        #            break
        #        self.ionIter()
        #        counter += 1
        #    self.outfile.write(f'    Self consistency on electron density converged to {self.electronDensityForIonization:10.2e} in {counter:3} iterations.\n')
        
        #self.electronDensity = self.electronDensityForIonization
        self.outfile.write('Ionization iteration: \n')
        self.outfile.write(f' Electrons contributed by the part of the gas is: mycalc (afterThisIter): {self.actualElectronDensity:10.2e} pynt (before): {self.electronDensityPYNT:10.2e}\n')
        self.outfile.write(' The electron fraction x_e = {:10.2e}\n'.format(self.sf.get_n_e() / self.elementNumberDensityTotal))
        
        
        counter = 0 
        self.outfile.write('AN Charge  FracBefore FracAfter  %Change\n')
        for aa,atomicnumber in enumerate(self.listOfAtomicNumbers):
            for ii in range(0,self.maxIonization):
                cc = (self.ionFractionOld[counter]-self.ionFraction[counter])/self.ionFractionOld[counter]
                self.outfile.write(f'{atomicnumber:3} {ii:2}+  {self.ionFractionOld[counter]:10.2e} {self.ionFraction[counter]:10.2e} {cc:10.2e} {self.sf.get_eff_ionpot(atomicnumber,ii+1):10.2e}\n')
                counter += 1
        return None
    
    def checkConvergence(self):
        self.converged = True 
        
        if np.any(np.abs(self.ionFraction - self.ionFractionOld)/self.ionFractionOld > TOLERANCE):
            self.converged = False
        
        return self.converged
    
    def runSpencerFano(self):
        #Runs Luke Shingles' Spencer-Fano solver.
        #
        #First add all the ions in the gas to an array.
        self.ionizationRatesOld = self.ionizationRates
        ions = []
        counter = 0 
        for atomicNumber in self.listOfAtomicNumbers:
            for ii in range(0,self.input.maxIonization):
                ions.append(
                    (atomicNumber,ii+1,self.balance[counter])
                )
                counter += 1
        
        #Pass this array to initialize the Spencer-Fano solver.
        self.sf = pynonthermal.SpencerFanoSolver(emin_ev=1, emax_ev=3000, npts=400, verbose=True)
        for Z, ion_stage, n_ion in ions:
            #print('debug:',Z,ion_stage,n_ion)
            self.sf.add_ionisation(Z, ion_stage, n_ion)
            
        self.outfile.write(f"Calling SpencerFano with: Edep = {self.depositionratedensity_ev:10.2e} eV/s/cm3.\n")
        
        if self.depositionMode.lower() == 'kasenbarnes':
            self.outfile.write(f"                               = {self.thermalizationEfficiency:10.2e} * {self.heatingRate.value:10.2e} eV/s * {self.elementNumberDensityTotal:10.2e} /cm3 \n")
        
        #If the user wants their own density
        if self.imposedElectronDensitySF is not None:
            self.sf.solve(depositionratedensity_ev = self.depositionratedensity_ev ,override_n_e=self.imposedElectronDensitySF)
        else:
            self.sf.solve(depositionratedensity_ev = self.depositionratedensity_ev)
            self.electronDensity = self.sf.calculate_free_electron_density()

        #Call to analysis - I don't know what this does but it seems necessary. 
        self.sf.analyse_ntspectrum()
        self.electronDensityPYNT = self.sf.calculate_free_electron_density()
        
        self.engrid = self.sf.engrid
        self.yvec   = self.sf.yvec
        
        #Pass the calculated Ionization Rates back to the self class. 
        for aa in range(0,self.numberOfElements):
            for ii in range(0,self.input.maxIonization):
                self.ionizationRates[ii,aa] = self.sf.get_ionisation_ratecoeff(self.listOfAtomicNumbers[aa], ii+1)
        
        
        
        return None 


    def writeOutBalanceAndRates(self,fileSuffix=''):
        '''
        Writes out the non-thermal ionization rates calculated by this code.
        '''
        file = open(f'pynt-balance-{fileSuffix}.dat','w')
        if self.imposedElectronDensitySF is None:
            self.imposedElectronDensitySF = 0.0
        if self.imposedElectronDensityRecombination == 0.0:
            self.imposedElectronDensityRecombination = 0.0
        header = '{:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e}\n'.format(
                                    self.thermalElectronTemperature,
                                    self.imposedElectronDensitySF,
                                    self.imposedElectronDensityRecombination,
                                    self.actualElectronDensity,
                                    self.timeSinceExplosionDays,
                                    self.velocityExpansionC,
                                    self.depositionratedensity_ev,
                                    self.elementMassDensityTotal,
                                    sum(self.massesOfElements)
                                    )
        file.write(header)
        
        writeFormat = '{:>3}, {:4}, {:2}{:2}, {:>2}+->{:>2}+, {:12.6e}, {:12.6e}, {:12.6e}, {:12.6e}, {:12.6e}\n'
        counter = 0
        for aa,atomicNumber in enumerate(self.listOfAtomicNumbers):
            
            for ii in range(0,self.maxIonization):
                file.write(writeFormat.format(
                  str(periodictable.elements[atomicNumber]),
                  atomicNumber,
                  atomicNumber,atomicNumber-ii,ii,ii+1,self.sf.get_eff_ionpot(atomicNumber,ii+1),self.ionizationRates[ii,aa],
                  self.recombinationRatesCoefficient[ii,aa]*self.electronDensityForIonization,self.ionFraction[counter],self.ionFraction[counter] * self.massesOfElements[aa]
                ))
                counter += 1 
        
        file.close()
        
        file = open(f'pynt-deg-{fileSuffix}.dat','w')
        
        for ii in range(0,len(self.yvec)):
            file.write('{:14.6e} {:14.6e}\n'.format(
                self.engrid[ii],
                self.yvec  [ii]
            ))
        
        file.close()
        
        return None 
        
        

def elementDensity(elementNumber,ion_mass_solar,velocity_c,explosion_time_days):
    #Calculates rough density assuming a Homologous expansion of 0 to velocity_c
    fourthirdspi = 4.0 * np.pi  / 3.0 

    atomic_mass = periodictable.elements[elementNumber].mass * u.u
    ion_mass = ion_mass_solar * u.solMass
    
    v = velocity_c * constants.c.cgs
    t = explosion_time_days*u.day
    nparticles = (ion_mass / atomic_mass).to('') # convert Msun/dalton to a number
    
    vt = ((v * t).to('cm'))
    expansion_volume = fourthirdspi * vt ** 3
    ionNumberDensity = nparticles / expansion_volume
    ionMassDensity   = (ion_mass  / expansion_volume).to('g/cm^3')
    
    #print(nparticles,v,t,expansion_volume,ionDensity)
    
    return ionNumberDensity.value,ionMassDensity.value, nparticles.value,expansion_volume.value

def heatingKasenBarnes(time_exp_days,averageAtomMass = 140,powerLaw = -1.3):
    '''
    Radioactive Power of r-process material.
    Power law of:
    Q = 1e10 t_d^-1.3  erg/s/g,
    where td is the time since merger measured in days from:
        https://iopscience.iop.org/article/10.3847/1538-4357/ab06c2/pdf  
    and references within.
    This function converts this to the work per ion, assuming an average atomic mass of  
    averageAtomMass nuclear units (140 by default).
    
    Using <A> = 140 with t = 1d, , this returns a work per ion of ~ 1.45 eV/s/ion,
    which is 10 * the number that Luke Shingles gave me.. - but need to multiply
    by a thermalization efficiency.
    '''

    const = 1e10 * u.erg / (u.g * u.s) #from 
    
    averageAtomicMass = (averageAtomMass * u.u).to('g')
    
    heating  = (averageAtomicMass * const).to('eV/s')
    
    return heating * (time_exp_days ** powerLaw)

def thermalizationKasenBarnes(time_exp_days,
                              ejectaMassSolar = 0.05,
                              eta = 1.0,
                              velocity_max_c=0.1,
                              temporalindex = 1.5,
                              fastelectron_parition = 0.2):
    
    ejectaMassSolar = ejectaMassSolar
    
    t_e = 15.0 * (eta * ejectaMassSolar / 0.01)**0.667 * (velocity_max_c / 0.2) **-2
    
    return fastelectron_parition * (1.0 + time_exp_days/t_e) ** (-1. * temporalindex)


def ionizationBalance(ionization_rates, recombination_rates,atomicNumber):
    ''' 
    Coronal ionization balance. 
    I.e - only consider:
        - ionization    out of ground into all adjacent states
        - recombination out of ground into all adjacent states
    
    ionization_rates   [i] = rate of ionization    from i   to i+1 
    recombination_rates[i] = rate of recombination from i+1 to i    
    
    pads out the balance, by assuming the ionization decreases 
    geometrically from the last explicitly calculated and that the 
    recombination increases geometrically. 
    
    assuming a factor of 2 drop off/increase - the balance drops off with factor 4.
    '''    
    populations = np.zeros(len(ionization_rates)+1)
    populations[0] = 1.0
    for z in range(0,len(ionization_rates)):
        populations[z+1] = populations[z] * ionization_rates[z] / recombination_rates[z]
    
    pp = populations[-1]
    norm = populations.sum()
    current = pp 
    totalarray = [*populations]
    for z in range(len(ionization_rates),atomicNumber):
        current = current * 0.25 
        norm += current 
        totalarray.append(current)
    
    totalarray = np.array(totalarray)
    totalarray /= norm
    populations = populations / norm
    
    #self.outfile.write(f'Ion population: Z = {atomicNumber} {totalarray[:]}\n')
    
    #if populations[-1] > 1e-4:
    #    import sys
    #    print(f'Populations not converged for {atomicNumber} ', totalarray[:]) 
    #    sys.exit()    
        
    return populations

def RRAxelrod(temp_kelvin,i):
    #i = ionization
    #todo: check if this formula needs to be replaced by the other one from Axelrod's thesis.
    AXELROD_CONST = 1E-13 
    
    temp_1e4K = temp_kelvin / 1e4 
    
    return AXELROD_CONST * i*i * np.power(temp_1e4K,-0.5 )
'''
def getSFSolution(ionizationBalance):
    totalIonDensity = ionDensity(atomicNumber,massElement,velocityExpansionLight,timeSinceExplosionDays)
    ionBalanceDensity = totalIonDensity * ionizationBalance
    ions = []
    for ii,dens in enumerate(ionBalanceDensity):
        ions.append( (atomicNumber,ii+1,dens) )
    sf = pynonthermal.SpencerFanoSolver(emin_ev=1, emax_ev=3000, npts=2000, verbose=False)
    for Z, ion_stage, n_ion in ions:
        sf.add_ionisation(Z, ion_stage, n_ion.value)
    depositionratedensity_ev = heatingKasenBarnes(timeSinceExplosionDays).value * totalIonDensity.value
    sf.solve(depositionratedensity_ev = depositionratedensity_ev ,override_n_e=densityElectronsImposed)
    sf.analyse_ntspectrum()
    
    numIncludedInSF = len(ionBalanceInitial)
    ionizationRates = np.zeros(atomicNumber)
    for ii in range(0,numIncludedInSF):
        ionizationRates[ii] = sf.get_ionisation_ratecoeff(58, ii+1)
    for ii in range(numIncludedInSF,atomicNumber):
        ionizationRates[ii] = 0.5 * ionizationRates[ii-1]
    return sf, ionizationRates
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--json', help='Path of JSON for fitting.')
    args = parser.parse_args()
    if (not args.json):
        input_default = input()
        default = json.dumps(input_default.__dict__,indent=1)
        print(default)
    else:
        with open(args.json, 'r') as j:
            contents = json.loads(j.read())
        thisInput = input(**contents)
        ntb = nonThermalBalance(thisInput,outfile_suffix=args.json)
        
        ntb.outfile.write('Iter 1\n')
        ntb.runSpencerFano()
        ntb.ionizationBalance()
        
        if ntb.input.selfConsistent:
            nehistory = []
            for ii in range(0,MAXITER):
            
                ntb.outfile.write(f'Iter {ii+2}\n')
                
                ntb.runSpencerFano()

                ntb.ionizationBalance()
                nehistory.append(ntb.electronDensity)
                #ntb.setNewBalanceDamped()
                
                if len(nehistory) == 3: 
                    #Aitken Jump: 
                    n0,n1,n2 = nehistory[-3:]
                    
                    dd = n2 - 2.0 * n1 + n0 
                    if dd == 0: 
                        print(n2,n1,n0)
                        import sys 
                        sys.exit()
                    neNew = n2 - (n2 -n1)**2 / dd 
                    ntb.electronDensity = neNew
                    ntb.outfile.write(f'Aitken jump {n0:10.2e}{n1:10.2e}{n2:10.2e}{neNew:10.2e}\n')
                    ntb.ionizationBalance()
                    nehistory = []
                ntb.checkConvergence()

                if ntb.converged:
                    break
        
        
        
        ntb.writeOutBalanceAndRates(fileSuffix=args.json)
        
    return 0 

if __name__ == '__main__':
    
    main()