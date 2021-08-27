import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


delta = 0.33
log_energy = np.arange(2.15, 4.15, delta) # energy in GeV
energy = pow(10.,log_energy)
eff_area_const = 10000.*(12.6*100.*100.) #cm^2
eff_area = np.piecewise(energy, [energy>=100,energy>=214,energy>=457,energy>=1000], [67.0798*10000.,19522.2*10000.,95972.3*10000.,126646*10000.])
angular_size = 1.0
time = 100. # exposure in hours
bkg_rate = 1e-12*19.4*pow(energy/1000.,-3.5)*3600. # per hour per 2 degree FoV
crab_rate = 1e-12*37.5*pow(energy/1000.,-2.5)*3600.
geminga_rate = 1e-12*3.5*pow(energy/1000.,-2.2)*3600.

N_CR = eff_area*bkg_rate*time*(angular_size*angular_size)/(2.*2.)
Err_CR_Stat = pow(N_CR,0.5)

crab_percent = 0.1
N_gamma_Crab10 = eff_area*crab_rate*crab_percent*time
N_gamma_Geminga = eff_area*geminga_rate*time
N_normalization = eff_area*geminga_rate*time
#N_normalization = eff_area*crab_rate*crab_percent*time
Err_CR_Stat = pow(N_CR,0.5)

scale_index = 0.

Relative_Syst = [0.022,0.015,0.025,0.029]
#Err_CR_Syst = N_CR*0.01
Err_CR_Syst = N_CR*np.piecewise(energy, [energy>=100,energy>=214,energy>=457,energy>=1000], Relative_Syst)
Limit = 5.0*pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat*Err_CR_Stat,0.5)*pow(energy,scale_index)
Limit_Smooth = interp1d(energy, Limit, kind='cubic')
Limit_Err = 1.0*pow(Err_CR_Stat*Err_CR_Stat,0.5)*pow(energy,scale_index)
Limit_Err_Smooth = interp1d(energy, Limit_Err, kind='cubic')

Relative_Syst = [0.043,0.030,0.039,0.049]
#Err_CR_Syst = N_CR*0.03
Err_CR_Syst = N_CR*np.piecewise(energy, [energy>=100,energy>=214,energy>=457,energy>=1000], Relative_Syst)
Limit_Init = 5.0*pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat*Err_CR_Stat,0.5)*pow(energy,scale_index)
Limit_Init_Smooth = interp1d(energy, Limit_Init, kind='cubic')

Limit_Stat = 5.0*pow(Err_CR_Stat*Err_CR_Stat,0.5)*pow(energy,scale_index)
Limit_Stat_Smooth = interp1d(energy, Limit_Stat, kind='cubic')

fig, ax = plt.subplots()
#plt.plot(energy, N_gamma_Crab10*pow(energy,scale_index), color='r', linestyle='dashed', label='10% Crab in 1-deg FoV')
plt.plot(energy, N_gamma_Geminga*pow(energy,scale_index)/N_normalization, color='r', label='Geminga in 1-deg FoV')
plt.plot(energy, Limit/N_normalization, color='b', label='MIBE')
plt.plot(energy, Limit_Init/N_normalization, color='b',linestyle='dashed', label='ON/OFF')
plt.plot(energy, Limit_Stat/N_normalization, color='g', label='only stat. err.')
#plt.xlim(214., 10000.)
ax.axis('on')
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel('Unit in Geminga flux')
ax.legend(loc='best')
ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig("output_plots/Exclusion.png")

