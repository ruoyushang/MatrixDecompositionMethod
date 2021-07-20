import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


eff_area = 10000.*(100.*100.) #cm^2
angular_size = 1.0
time = 100. # exposure in hours
delta = 1.0
energy = np.arange(214, 10000, delta) # energy in GeV
bkg_rate = 1e-12*33.7*pow(energy/1000.,-3.4)*3600. # per hour per 2 degree FoV
crab_rate = 1e-12*37.5*pow(energy/1000.,-2.5)*3600.

N_CR = eff_area*bkg_rate*time*(angular_size*angular_size)/(2.*2.)
Err_CR_Stat = pow(N_CR,0.5)

crab_percent = 0.1
N_gamma_Crab10 = eff_area*crab_rate*crab_percent*time
Err_CR_Stat_Crab10 = pow(N_CR+N_gamma_Crab10,0.5)

#Err_CR_Syst = N_CR*0.01
Err_CR_Syst = N_CR*np.piecewise(energy, [energy>=214,energy>=457,energy>=1000], [0.014, 0.018, 0.037])
Limit = 5.0*pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab10*Err_CR_Stat_Crab10,0.5)
Limit_Err = 1.0*pow(Err_CR_Stat_Crab10*Err_CR_Stat_Crab10,0.5)

#Err_CR_Syst = N_CR*0.03
Err_CR_Syst = N_CR*np.piecewise(energy, [energy>=214,energy>=457,energy>=1000], [0.040, 0.048, 0.056])
Limit_Init = 5.0*pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab10*Err_CR_Stat_Crab10,0.5)

fig, ax = plt.subplots()
plt.plot(energy, N_gamma_Crab10*pow(energy,2), color='r', label='10% Crab in 1-deg FoV')
plt.plot(energy, Limit*pow(energy,2), color='b', label='MIBE')
plt.fill_between(energy, (Limit-Limit_Err)*pow(energy,2), (Limit+Limit_Err)*pow(energy,2), alpha=0.2, color='b')
plt.plot(energy, Limit_Init*pow(energy,2), color='b',linestyle='dashed', label='ON/OFF')
#plt.ylim(0, 10)
ax.axis('on')
ax.set_xlabel('Energy [GeV]')
ax.set_ylabel('$E^{2}$ Flux')
ax.legend(loc='best')
ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig("output_plots/Exclusion.png")

