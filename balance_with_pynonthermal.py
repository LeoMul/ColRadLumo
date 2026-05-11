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
TOLERANCE = 1.E-2
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
                 averageAtomicMass  = 140,
                 maxIonization      = 6,
                 depositionOverride = None,
                 velocityMaxForEfficiency   = 0.1 
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
        self.depositionOverride         = depositionOverride
        self.velocityMaxForEfficiency   = velocityMaxForEfficiency
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
        self.imposedElectronDensity     = input.imposedElectronDensity
        self.velocityExpansionC         = input.velocityExpansionC
        self.timeSinceExplosionDays     = input.timeSinceExplosionDays
        self.massesOfElements           = input.massesOfElements
        self.averageAtomicMass          = input.averageAtomicMass
        self.maxIonization              = input.maxIonization
        self.depositionOverride         = input.depositionOverride
        self.velocityMaxForEfficiency   = input.velocityMaxForEfficiency
        self.numberOfElements = len(input.listOfAtomicNumbers)
        
        self.outfile = open(f'pynt-balance-{outfile_suffix}.out','w')
        
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
        balancePerElement[:] /= balancePerElement.sum()
        
        self.initBalance  = np.zeros(self.input.maxIonization * self.numberOfElements)
        self.initFraction = np.zeros_like(self.initBalance) 
        #Need total matter density for injection later.
        self.elementDensities = np.zeros(self.numberOfElements)
        self.nparticles = np.zeros(self.numberOfElements)

        self.expansion_volume =0.0
        for ii in range(0,self.numberOfElements):
            self.elementDensities[ii],self.nparticles [ii],self.expansion_volume = elementDensity(
                self.listOfAtomicNumbers[ii],
                self.massesOfElements[ii],
                self.velocityExpansionC,
                self.timeSinceExplosionDays
            )
            
        self.elementDensityTotal = self.elementDensities.sum()
        
        stride = self.input.maxIonization
        for ii in range(0,self.numberOfElements):
            self.initBalance [ii * stride : (ii+1) * stride ] = balancePerElement * self.elementDensities[ii]
            self.initFraction[ii * stride : (ii+1) * stride ] = balancePerElement
                
        self.outfile.write(f'Total number densities: {self.elementDensityTotal:10.2e}\n',)
        self.outfile.write('Element densities:\n')
        for aa,zz in enumerate(self.listOfAtomicNumbers):
            self.outfile.write(f'{zz:3} {self.elementDensities[aa]:10.2e}\n')
        
        self.balance = self.initBalance.copy()
        self.ionFraction = self.initFraction.copy()
        return None 
            
    def _setDepositionRateDensity(self):
        '''
        Set deposition rate - eV /s / cm^3,
            =  Qdot * total matter density 
        '''
        if self.depositionOverride is None:

            self.heatingRate = heatingKasenBarnes(
                self.timeSinceExplosionDays,
                self.averageAtomicMass
            )
            self.thermalizationEfficiency = thermalizationKasenBarnes(self.timeSinceExplosionDays,velocity_max_c=self.velocityMaxForEfficiency)
            self.outfile.write(f'efficiency   of  {self.thermalizationEfficiency}\n')
            self.outfile.write(f'using a heating rate of {self.heatingRate} \n')
            self.depositionratedensity_ev  = self.heatingRate * self.elementDensityTotal * self.thermalizationEfficiency
            try:
                self.depositionratedensity_ev = self.depositionratedensity_ev.value 
            except:
                pass
        else:
            self.depositionratedensity_ev = self.depositionOverride 
            self.outfile.write('Using an user-override deposition density of {:10.3}\n'.format(self.depositionOverride))
        
        return None 
    
    def ionizationBalance(self):
        #Calculate a new Ionization balance
        self.balanceOld = self.balance.copy()
        self.ionFractionOld = self.ionFraction.copy() 
        
        stride = self.maxIonization
        
        self.actualElectronDensity = 0.0 
        
        for aa,Z in enumerate(self.listOfAtomicNumbers):
            
            thisBalance = ionizationBalance(
                self.ionizationRates[:,aa], 
                self.imposedElectronDensity * self.recombinationRatesCoefficient[:,aa],Z
                )[0:self.maxIonization]
            
            self.balance[aa * stride : (aa+1) * stride ]     = thisBalance * self.elementDensities[aa]
            
            self.actualElectronDensity += np.sum ( np.arange(0,self.maxIonization,1,dtype=int) * thisBalance * self.elementDensities[aa])
            
            self.ionFraction[aa * stride : (aa+1) * stride ] = thisBalance 

        #print(self.ionFractionOld)
        #print(self.ionFraction)
        
        self.outfile.write('Ionization iteration: \n')
        self.outfile.write(f' Electrons contributed by the part of the gas is: mycalc: {self.actualElectronDensity:10.2e} pynt: {self.electronDensityPYNT:10.2e}\n')
        
        
        
        counter = 0 
        for aa,atomicnumber in enumerate(self.listOfAtomicNumbers):
            for ii in range(0,self.maxIonization):
                self.outfile.write(f'{atomicnumber:3} {ii:2}+ {self.ionFraction[counter]:10.2e} {self.ionFractionOld[counter]:10.2e} \n')
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
        ions = []
        counter = 0 
        for atomicNumber in self.listOfAtomicNumbers:
            for ii in range(0,self.input.maxIonization):
                ions.append(
                    (atomicNumber,ii+1,self.balance[counter])
                )
                counter += 1
        
        #Pass this array to initialize the Spencer-Fano solver.
        sf = pynonthermal.SpencerFanoSolver(emin_ev=1, emax_ev=3000, npts=200, verbose=False)
        for Z, ion_stage, n_ion in ions:
            #print('debug:',Z,ion_stage,n_ion)
            sf.add_ionisation(Z, ion_stage, n_ion)
            
        self.outfile.write(f"Calling SpencerFano with: Edep = {self.depositionratedensity_ev:10.2e} eV/s/cm3.\n")
        if self.depositionOverride is None:
            self.outfile.write(f"                               = {self.thermalizationEfficiency:10.2e} * {self.heatingRate.value:10.2e} eV/s * {self.elementDensityTotal:10.2e} /cm3 \n")
        else:
            self.outfile.write(f"                               = user override.\n")
        
        #If the user wants their own density
        if self.imposedElectronDensity is not None:
            sf.solve(depositionratedensity_ev = self.depositionratedensity_ev ,override_n_e=self.imposedElectronDensity)
        else:
            sf.solve(depositionratedensity_ev = self.depositionratedensity_ev)

        #Call to analysis - I don't know what this does but it seems necessary. 
        sf.analyse_ntspectrum()
        self.electronDensityPYNT = sf.calculate_free_electron_density()
        #Pass the calculated Ionization Rates back to the self class. 
        for aa in range(0,self.numberOfElements):
            for ii in range(0,self.input.maxIonization):
                self.ionizationRates[ii,aa] = sf.get_ionisation_ratecoeff(self.listOfAtomicNumbers[aa], ii+1)
        
        
        
        return None 


    def writeOutBalanceAndRates(self,fileSuffix=''):
        '''
        Writes out the non-thermal ionization rates calculated by this code.
        '''
        file = open(f'pynt-balance-{fileSuffix}.dat','w')
        
        header = '{:13.7e} {:13.7e} {:13.7e} {:13.7e} {:13.7e}\n'.format(
                                    self.thermalElectronTemperature,
                                    self.imposedElectronDensity,
                                    self.actualElectronDensity,
                                    self.timeSinceExplosionDays,
                                    self.velocityExpansionC
                                    )
        file.write(header)
        
        writeFormat = '{:>3}, {:4}, {:2}{:2}, {:>2}+->{:>2}+, {:12.6e}, {:12.6e}, {:12.6e}, {:12.6e}\n'
        counter = 0
        for aa,atomicNumber in enumerate(self.listOfAtomicNumbers):
            
            for ii in range(0,self.maxIonization):
                file.write(writeFormat.format(
                  str(periodictable.elements[atomicNumber]),
                  atomicNumber,
                  atomicNumber,atomicNumber-ii,ii,ii+1,self.ionizationRates[ii,aa],
                  self.recombinationRatesCoefficient[ii,aa]*self.imposedElectronDensity,self.ionFraction[counter],self.ionFraction[counter] * self.massesOfElements[aa]
                ))
                counter += 1 
        
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
    
    return ionDensity.value, nparticles.value,expansion_volume.value

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
            for ii in range(0,MAXITER):
            
                ntb.outfile.write(f'Iter {ii+2}\n')
                ntb.runSpencerFano()
                ntb.ionizationBalance()
                ntb.checkConvergence()

                if ntb.converged:
                    break
        
        
        
        ntb.writeOutBalanceAndRates(fileSuffix=args.json)
        
    return 0 

if __name__ == '__main__':
    
    main()