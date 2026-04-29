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


''''
Charge state balance for a given set of elements.

Uses the pynonthermal code of Luke Shingles to obtain non-thermal ionization rates 
for a nebular gas.

Uses your recombination data, but I would prefer you used mine.

Self consistent iteration - optional.

Current todo's: 

Have the code output the ionization fractions, for easier calculation of 
spectra later. 
Possibly extend some index's by one - I think we can squeeze out an extra 
ion for the ionization balance. Its probably fine the way it is. 
Print out the degradation spectrum in each case.
Print out the total electron density from the ionization balance, and 
compare with the imposed electron density. 

Need to check that all of the units and input degradation is correct. 

'''
TOLERANCE = 1.E-3 
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
                 imposedElectronDensity    = None,
                 velocityExpansionC = None,
                 timeSinceExplosionDays = None, 
                 selfConsistent = False,
                 averageAtomicMass = 140,
                 maxIonization     = 6 
                 ):        

        #Boring transfer of memory. 
        self.listOfAtomicNumbers        = listOfAtomicNumbers
        self.massesOfElements           = massesOfElements
        self.pathsOfRecombinationData   = pathsOfRecombinationData
        self.thermalElectronTemperature = thermalElectronTemperature
        self.imposedElectronDensity     = imposedElectronDensity
        self.velocityExpansionC         = velocityExpansionC
        self.timeSinceExplosionDays     = timeSinceExplosionDays
        self.averageAtomicMass          = averageAtomicMass
        self.maxIonization              = maxIonization
        self.selfConsistent             = selfConsistent
class nonThermalBalance:
    '''
    Non thermal ionization balance class.
    '''
    
    def __init__(self,input: input):
        
        #Transfer memory we might need. 
        self.input                      = input
        self.listOfAtomicNumbers        = input.listOfAtomicNumbers
        self.pathsOfRecombinationData   = input.pathsOfRecombinationData
        self.thermalElectronTemperature = input.thermalElectronTemperature
        self.imposedElectronDensity     = input.imposedElectronDensity
        self.velocityExpansionC         = input.velocityExpansionC
        self.timeSinceExplosionDays     = input.timeSinceExplosionDays
        self.massesOfElements           = input.massesOfElements
        self.averageAtomicMass          = input.averageAtomicMass
        self.maxIonization              = input.maxIonization
        
        self.numberOfElements = len(input.listOfAtomicNumbers)
        
        #sanity checks on input data. 
        assert( len(input.massesOfElements)         == self.numberOfElements)
        assert( len(input.pathsOfRecombinationData) == self.numberOfElements * input.maxIonization)

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
                
                if self.pathsOfRecombinationData[counter] == None:
                    rate = RRAxelrod(self.thermalElectronTemperature, jj+1)
                else:
                    data = np.loadtxt(self.pathsOfRecombinationData[counter])
                    rate = interp1d(data[:,0],data[:,1])(self.thermalElectronTemperature)
                
                self.recombinationRatesCoefficient[jj,ii] = rate 
                counter += 1 
        
        return None 
    
    def _setInitialIonizationBalance(self):
        #Initial Balance per element.  Will be discard later.
        #Assuming an equal ionization. The code will eventually get something better.
        balancePerElement = np.zeros(self.input.maxIonization)
        balancePerElement[:] = 1.0 
        balancePerElement[:] /= balancePerElement.sum()
        
        self.initBalance  = np.zeros(self.input.maxIonization * self.numberOfElements)
        self.initFraction = np.zeros_like(self.initBalance) 
        #Need total matter density for injection later.
        self.elementDensities = np.zeros(self.numberOfElements)
        for ii in range(0,self.numberOfElements):
            self.elementDensities[ii] = elementDensity(
                self.listOfAtomicNumbers[ii],
                self.massesOfElements[ii],
                self.velocityExpansionC,
                self.timeSinceExplosionDays
            ).value
            
        self.elementDensityTotal = self.elementDensities.sum()
        
        stride = self.input.maxIonization
        for ii in range(0,self.numberOfElements):
            self.initBalance [ii * stride : (ii+1) * stride ] = balancePerElement * self.elementDensities[ii]
            self.initFraction[ii * stride : (ii+1) * stride ] = balancePerElement
        print(self.initBalance)
        
        print(self.elementDensities,self.elementDensities.sum())
        
        self.balance = self.initBalance.copy()
        self.ionFraction = self.initFraction.copy()
        return None 
            
    def _setDepositionRateDensity(self):
        '''
        Set deposition rate - eV /s / cm^3,
            =  Qdot * total matter density 
        '''
        
        self.heatingRate = heatingKasenBarnes(
            self.timeSinceExplosionDays,
            self.averageAtomicMass
        )
        self.depositionratedensity_ev  = self.heatingRate * self.elementDensityTotal
        return None 
    
    def ionizationBalance(self):
        #Calculate a new Ionization balance
        self.balanceOld = self.balance.copy()
        self.ionFractionOld = self.ionFraction.copy() 
        
        stride = self.maxIonization
        for aa,Z in enumerate(self.listOfAtomicNumbers):
            
            thisBalance = ionizationBalance(
                self.ionizationRates[:,aa], 
                self.imposedElectronDensity * self.recombinationRatesCoefficient[:,aa]
                )[0:self.maxIonization]
            self.balance[aa * stride : (aa+1) * stride ]     = thisBalance * self.elementDensities[aa]
            self.ionFraction[aa * stride : (aa+1) * stride ] = thisBalance 

        print(self.ionFractionOld)
        print(self.ionFraction)
    
        return None
    
    def checkConvergence(self):
        self.converged = True 
        
        if np.any(np.abs(self.ionFraction - self.ionFractionOld) > TOLERANCE):
            self.converged = False
        
        return self.converged
    
    def runSpencerFano(self):
        #Runs Luke Shingles' Spencer-Fano solver.
        #
        #First add all the ions in the gas to an array.
        ions = []
        counter = 0 
        for atomicNumber in self.listOfAtomicNumbers:
            for ii in range(0,self.input.maxIonization):
                ions.append(
                    (atomicNumber,ii+1,self.balance[counter])
                )
                counter += 1
        
        #Pass this array to initialize the Spencer-Fano solver.
        sf = pynonthermal.SpencerFanoSolver(emin_ev=1, emax_ev=3000, npts=2000, verbose=False)
        for Z, ion_stage, n_ion in ions:
            #print('debug:',Z,ion_stage,n_ion)
            sf.add_ionisation(Z, ion_stage, n_ion)
        
        #If the user wants their own density
        if self.imposedElectronDensity is not None:
            sf.solve(depositionratedensity_ev = self.depositionratedensity_ev.value ,override_n_e=self.imposedElectronDensity)
        else:
            sf.solve(depositionratedensity_ev = self.depositionratedensity_ev.value)

        #Call to analysis - I don't know what this does but it seems necessary. 
        sf.analyse_ntspectrum()
        
        #Pass the calculated Ionization Rates back to the self class. 
        for aa in range(0,self.numberOfElements):
            for ii in range(0,self.input.maxIonization):
                self.ionizationRates[ii,aa] = sf.get_ionisation_ratecoeff(self.listOfAtomicNumbers[aa], ii+1)
        
        
        
        return None 


    def writeOutIonizationRates(self):
        '''
        Writes out the non-thermal ionization rates calculated by this code.
        '''
        file = open('nonThermalRates.dat','w')
        
        writeFormat = '{:>3}, {:4}, {:>2}+->{:>2}+, {:12.6e}\n'
        
        for aa,atomicNumber in enumerate(self.listOfAtomicNumbers):
            
            for ii in range(0,self.maxIonization):
                file.write(writeFormat.format(
                  str(periodictable.elements[atomicNumber]),atomicNumber,ii,ii+1,self.ionizationRates[ii,aa]  
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
    ionDensity = nparticles / expansion_volume
    
    
    #print(nparticles,v,t,expansion_volume,ionDensity)
    
    return ionDensity

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
    
    Using <A> = 140 with t = 1d, , this returns a work per ion of ~ 0.145 eV/s/ion,
    which is close to the number that Luke Shingles gave me. 
    '''

    const = 1e10 * u.erg / (u.g * u.s) #from 
    
    averageAtomicMass = (averageAtomMass * u.u).to('g')
    
    heating  = (averageAtomicMass * const).to('eV/s')
        
    return heating * time_exp_days ** powerLaw

def ionizationBalance(ionization_rates, recombination_rates):
    ''' 
    Coronal ionization balance. 
    I.e - only consider:
        - ionization    out of ground into all adjacent states
        - recombination out of ground into all adjacent states
    
    ionization_rates   [i] = rate of ionization    from i   to i+1 
    recombination_rates[i] = rate of recombination from i+1 to i    
    
    '''    
    populations = np.zeros(len(ionization_rates)+1)
    populations[0] = 1.0
    for z in range(0,len(ionization_rates)):
        populations[z+1] = populations[z] * ionization_rates[z] / recombination_rates[z]
    return populations / populations.sum()

def RRAxelrod(temp_kelvin,i):
    #i = ionization
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
        ntb = nonThermalBalance(thisInput)
        
        print('Iter 1')
        ntb.runSpencerFano()
        ntb.ionizationBalance()
        if ntb.input.selfConsistent:
            for ii in range(0,MAXITER):
            
                print(f'Iter {ii+2}')
                ntb.runSpencerFano()
                ntb.ionizationBalance()
                ntb.checkConvergence()

                if ntb.converged:
                    break
                
        ntb.writeOutIonizationRates()
        
    return 0 

if __name__ == '__main__':
    
    main()