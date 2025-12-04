import numpy as np 
from colradlumo_library import * 
def produce_lumo_for_jwst(adf04_path,density,temperature_kelvin,
                          amu=0.0,
                          opacity_override = 0.0,
                          time_exp=0.0,
                          velocity=0.0,
                          opacity_accelerator=False,
                          writepath=''):
    
    temperature_ev = np.array([temperature_kelvin]) / 11600 
    
    density_ar = np.array([density])
    
    
    colrad_lumo_run = colradlumo_calc(
        adf04_path,
        density_ar,
        temperature_ev,
        atomic_mass_number=amu
    )
    mass = np.array([1e0])

    if time_exp != 0:
        colrad_lumo_run.convergeOpacity(temperature_ev[0],
                                        density_ar[0],
                                        time_exp_days=time_exp,
                                        velocity_c=velocity,
                                        mass_solar=1.0,
                                        fraction_override=opacity_override,
                                        oscillator_breaker=opacity_accelerator)
        
    colrad_lumo_run.scale_lumo_by_ion_mass(mass)

    file_name = str(colrad_lumo_run.element) + str(colrad_lumo_run.charge_roman) + 'Dens={:8.2e}.Temp={:4.0f}.dat'.format(
        density,
        temperature_kelvin
    )
    
    if writepath != '':
        if not ('/' in writepath):
            writepath += '/'
        file_name = writepath + file_name
             
    g = open(file_name,'w')
    #colrad_lumo_run.colradpy_run.processed.
    format_string = '{:2} {:5} {:11.4f} {:11.4E}\n'
    
    lumo = colrad_lumo_run.scaled_lumo_ergs[:,0,0,0]
    
    args = np.argsort(lumo)[::-1]
    
    for ii in range(0 ,len(colrad_lumo_run.wl_air_nm)):
        lumo = colrad_lumo_run.scaled_lumo_ergs[args[ii],0,0,0]
        wl   = colrad_lumo_run.wl_vac_nm[args[ii]]
        g.write(format_string.format(
            colrad_lumo_run.element,
            colrad_lumo_run.charge_roman,
            wl,
            lumo
        )) 
        
    g.close()    
    
    
    return 0 