from colradlumo_library import * 

adf04_path = '/Users/leomulholland/TeIII_collection/TeIII_gamma/adf04/TeIII_LPM_QUB.adf04'
dens_cm_3 = 1.e6 
mass_number_te_iii = 127.6 #nuclear units

t = 10.5
Hotokezaka_density = 1e7 * (9.5/t)**3

dens_gfo = Hotokezaka_density
dens_grb = 3e5

temp_gfo = 2000/11600
temp_grb = 3000/11600

mass_te_iii = 0.001
mass_te_ii  = 0.001
mass_te_i   = 0.001
print('Masses: ')
print('M Te I   = {:6.4E} M_o '.format(mass_te_i  ))
print('M Te II  = {:6.4E} M_o '.format(mass_te_ii ))
print('M Te III = {:6.4E} M_o '.format(mass_te_iii))




te_i_lumo_gfo = colradlumo_calc('/Users/leomulholland/TeI_kappa/omegafiles/adf04/adf04',dens_gfo,temp_gfo,atomic_mass_number=mass_number_te_iii)
te_ii_lumo_gfo = colradlumo_calc('/Users/leomulholland/TeII_collection/TeII_epsilon/adas/adf04',dens_gfo,temp_gfo,atomic_mass_number=mass_number_te_iii)
te_iii_lumo_gfo = colradlumo_calc(adf04_path,dens_gfo,temp_gfo,atomic_mass_number=mass_number_te_iii)

te_i_lumo_grb = colradlumo_calc('/Users/leomulholland/TeI_kappa/omegafiles/adf04/adf04',dens_grb,temp_grb,atomic_mass_number=mass_number_te_iii)
te_ii_lumo_grb= colradlumo_calc('/Users/leomulholland/TeII_collection/TeII_epsilon/adas/adf04',dens_grb,temp_grb,atomic_mass_number=mass_number_te_iii)
te_iii_lumo_grb = colradlumo_calc(adf04_path,dens_grb,temp_grb,atomic_mass_number=mass_number_te_iii)


te_iii_lumo_gfo.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_i)
te_ii_lumo_gfo.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_ii)
te_i_lumo_gfo.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_iii)


te_iii_lumo_grb.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_i)
te_ii_lumo_grb.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_ii)
te_i_lumo_grb.scale_lumo_by_ion_mass(mass_of_ion_solar_units=mass_te_iii)


print('n_e = {:6.4E} cm^-3'.format(dens_gfo))
print('T_e = {:6.4f} eV \sim {:6.4f} K'.format(temp_gfo,temp_gfo*11600))

requested_lines = te_iii_lumo_gfo.select_strongest_n_lines(10)
requested_lines.display()
requested_lines = te_ii_lumo_gfo.select_strongest_n_lines(10)
requested_lines.display()
requested_lines = te_i_lumo_gfo.select_strongest_n_lines(10)
requested_lines.display()


print('n_e = {:6.4E} cm^-3'.format(dens_grb))
print('T_e = {:6.4f} eV \sim {:6.4f} K'.format(temp_grb,temp_grb*11600))

requested_lines = te_iii_lumo_grb.select_strongest_n_lines(10)
requested_lines.display()
requested_lines = te_ii_lumo_grb.select_strongest_n_lines(10)
requested_lines.display()
requested_lines = te_i_lumo_grb.select_strongest_n_lines(10)
requested_lines.display()