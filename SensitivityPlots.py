import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

eff_area = 10000.*(100.*100.) #cm^2
angular_size = 1.0
energy_threshold = 214
bkg_rate = 1e-12*33.7*pow(energy_threshold/1000.,-3.4+1)*(1./(3.4-1))*3600. # per hour per 2 degree FoV
crab_rate = 1e-12*37.5*pow(energy_threshold/1000.,-2.5+1)*(1./(2.5-1))*3600.

delta = 0.01
time = np.arange(0.5, 300., delta) # exposure
N_CR = eff_area*bkg_rate*time*(angular_size*angular_size)/(2.*2.)
Err_CR_Stat = pow(N_CR,0.5)
crab_percent = 0.02
N_gamma_Crab02 = eff_area*crab_rate*crab_percent*time
Err_CR_Stat_Crab02 = pow(N_CR+N_gamma_Crab02,0.5)
crab_percent = 0.05
N_gamma_Crab05 = eff_area*crab_rate*crab_percent*time
Err_CR_Stat_Crab05 = pow(N_CR+N_gamma_Crab05,0.5)
crab_percent = 0.1
N_gamma_Crab10 = eff_area*crab_rate*crab_percent*time
Err_CR_Stat_Crab10 = pow(N_CR+N_gamma_Crab10,0.5)

Err_CR_Syst = N_CR*0.01
Sigma_Crab02 = N_gamma_Crab02/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab02*Err_CR_Stat_Crab02,0.5)
Sigma_Crab05 = N_gamma_Crab05/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab05*Err_CR_Stat_Crab05,0.5)
Sigma_Crab10 = N_gamma_Crab10/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab10*Err_CR_Stat_Crab10,0.5)

Err_CR_Syst = N_CR*0.03
Sigma_Crab02_Init = N_gamma_Crab02/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab02*Err_CR_Stat_Crab02,0.5)
Sigma_Crab05_Init = N_gamma_Crab05/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab05*Err_CR_Stat_Crab05,0.5)
Sigma_Crab10_Init = N_gamma_Crab10/pow(Err_CR_Syst*Err_CR_Syst+Err_CR_Stat_Crab10*Err_CR_Stat_Crab10,0.5)

time_data = [39.8,75.6,166.4,284.5]
data_sigma_crab02 = [1.4,2.3,4.0,4.2]
data_sigma_crab05 = [3.2,5.0,7.0,7.0]
data_sigma_crab10 = [6.0,9.5,12.0,12.0]

fig, ax = plt.subplots()
plt.plot(time, Sigma_Crab02, color='r', label='2% Crab')
plt.fill_between(time, Sigma_Crab02-1, Sigma_Crab02+1, alpha=0.2, color='r')
plt.plot(time, Sigma_Crab05, color='g', label='5% Crab')
plt.fill_between(time, Sigma_Crab05-1, Sigma_Crab05+1, alpha=0.2, color='g')
plt.plot(time, Sigma_Crab10, color='b', label='10% Crab')
plt.fill_between(time, Sigma_Crab10-1, Sigma_Crab10+1, alpha=0.2, color='b')
plt.plot(time, Sigma_Crab02_Init, color='r',linestyle='dashed')
plt.plot(time, Sigma_Crab05_Init, color='g',linestyle='dashed')
plt.plot(time, Sigma_Crab10_Init, color='b',linestyle='dashed')
plt.errorbar(time_data, data_sigma_crab02, yerr=1., color='r', marker='s', ls='none')
plt.errorbar(time_data, data_sigma_crab05, yerr=1., color='g', marker='s', ls='none')
plt.errorbar(time_data, data_sigma_crab10, yerr=1., color='b', marker='s', ls='none')
#plt.ylim(0, 10)
ax.axis('on')
ax.set_xlabel('Exposure [hour]')
ax.set_ylabel('Detection Significance')
ax.legend(loc='best')
#ax.set_yscale('log')
plt.savefig("output_plots/Significance_%sGeV.png"%(energy_threshold))

plt.clf()
fig, ax = plt.subplots()
Err_CR_Syst = N_CR*0.01
Relative_Stat = Err_CR_Stat/N_CR
Relative_Syst = Err_CR_Syst/N_CR
plt.plot(time, Relative_Stat, color='r', label='Stat err')
plt.plot(time, Relative_Syst, color='b', label='Syst err')
plt.ylim(0, 0.1)
ax.axis('on')
ax.set_xlabel('Exposure [hour]')
ax.set_ylabel('Relative Error')
ax.legend(loc='best')
plt.savefig("output_plots/RelativeErr_%sGeV.png"%(energy_threshold))

