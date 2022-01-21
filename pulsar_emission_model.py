
import array
import math
from array import *
from ROOT import *
from scipy import special
from scipy import interpolate
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

def IC_electron_energy_loss_function(E_e_GeV):
    E_e = E_e_GeV*1e9
    gamma_factor = E_e/m_e
    dEdt = -4./3.*sigma_thomson*speed_light*gamma_factor*gamma_factor*(U_cmb) # eV/s
    return dEdt

def Sync_electron_energy_loss_function(E_e_GeV):
    E_e = E_e_GeV*1e9
    gamma_factor = E_e/m_e
    dEdt = -4./3.*sigma_thomson*speed_light*gamma_factor*gamma_factor*(U_B) # eV/s
    return dEdt

def total_electron_energy_loss_function(E_e_GeV):
    return IC_electron_energy_loss_function(E_e_GeV)+Sync_electron_energy_loss_function(E_e_GeV)

def initial_electron_energy(dt_year,E_e_final_GeV):
    dt_sec = dt_year*365.*24.*60.*60.
    delta_inv_E = 4./3.*sigma_thomson*speed_light*pow(m_e,-2)*(U_cmb+U_B)*dt_sec # 1/eV
    inv_E_init = - delta_inv_E + 1./(E_e_final_GeV*1e9)
    E_init = 1./inv_E_init # eV
    return E_init/1e9 # GeV

def diffusion_coefficient(dt_year,E_e_final_GeV):
    E_e_GeV = initial_electron_energy(dt_year,E_e_final_GeV)
    return D0*pow(E_e_GeV,alpha)*year_to_sec # cm2/year

def diffusion_length(dt_year,E_e_final_GeV):
    pc_to_cm = 3.086e+18
    dt_year_limit = 1./(4./3.*sigma_thomson*speed_light*pow(m_e,-2)*(U_cmb+U_B)*E_e_final_GeV*1e9)/year_to_sec
    if dt_year>dt_year_limit:
        dt_year = dt_year_limit
    length = 2.*pow(quad(diffusion_coefficient, 0., dt_year, args=(E_e_final_GeV))[0],0.5)/pc_to_cm
    return length # pc

def electron_density_burst(t0_year,photon_location,source_location_initial,source_location_final,E_e_final_GeV):
    # electrons released at t0, radiate photons at t_now = pulsar age
    dt_year = pulsar_age_year-t0_year
    dt_year_limit = 1./(4./3.*sigma_thomson*speed_light*pow(m_e,-2)*(U_cmb+U_B)*E_e_final_GeV*1e9)/year_to_sec-1.
    if dt_year>dt_year_limit:
        return 0.
    source_location_at_t_0 = (source_location_initial-source_location_final)/pulsar_age_year*dt_year+source_location_final
    distance_to_source = pow(pow(photon_location[1]-source_location_at_t_0,2)+pow(photon_location[0],2),0.5)
    diffusion_length_at_t_now = diffusion_length(dt_year,E_e_final_GeV)
    energy_loss_at_t_now = total_electron_energy_loss_function(E_e_final_GeV)
    E_e_initial_GeV = initial_electron_energy(dt_year,E_e_final_GeV)
    energy_loss_at_t0 = total_electron_energy_loss_function(E_e_initial_GeV)
    Qt = injection_spectrum_time_dep(E_e_initial_GeV,t0_year)
    electron_density = Qt*exp(-pow(distance_to_source,2)/pow(diffusion_length_at_t_now,2))*pow(2.*3.14,-3/2)*pow(diffusion_length_at_t_now,-3)*energy_loss_at_t0/energy_loss_at_t_now
    #electron_density = pow(E_e_initial_GeV,-2)*exp(-pow(distance_to_source,2)/pow(diffusion_length_at_t_now,2))*pow(2.*3.14,-3/2)*pow(diffusion_length_at_t_now,-3)*energy_loss_at_t0/energy_loss_at_t_now
    return electron_density

def electron_density_continuum(photon_location_z,photon_location_x,t0_year,E_e_final_GeV):
    source_location_initial = proper_velocity_pc_per_year*pulsar_age_year
    source_location_final = 0.
    photon_location = [photon_location_z,photon_location_x]
    electron_density = quad(electron_density_burst, 0., t0_year, args=(photon_location,source_location_initial,source_location_final,E_e_final_GeV))[0]
    return electron_density

def electron_column_density_continuum(photon_location_x,t0_year,E_e_final_GeV):
    x_low = 0.
    x_up = 2.*diffusion_length(t0_year,E_e_final_GeV)
    electron_column_density = quad(electron_density_continuum, x_low, x_up, args=(photon_location_x,t0_year,E_e_final_GeV))[0]
    return 2.*electron_column_density

def electron_column_density_continuum_distance(photon_location,t0_year,E_e_final_GeV):
    return electron_column_density_continuum(photon_location,t0_year,E_e_final_GeV)

def injection_spectrum(E_e_GeV):
    E_cutoff_GeV = 400.*1e3
    return E_e_GeV*pow(E_e_GeV,-gamma_index)*exp(-E_e_GeV/E_cutoff_GeV)

def injection_normalization():
    E1_GeV = 0.1
    E_cutoff_GeV = 400.*1e3
    E_max_GeV = 10.*E_cutoff_GeV
    integral = quad(injection_spectrum, E1_GeV, E_max_GeV)[0]
    erg_to_GeV = 624.151
    spin_down_luminosity = 1.7*1e37*erg_to_GeV*year_to_sec # GeV/year
    Q0 = spin_down_luminosity/integral
    return Q0

def injection_spectrum_time_dep(E_e_GeV,t_year):
    tau0 = 15.*1e3 # year
    spectrum = Q0_normalization*pow(1+t_year/tau0,-2)*injection_spectrum(E_e_GeV)/E_e_GeV
    return spectrum

mag_field = 3. # muG
E_cmb = 6.6*1e-4 # eV
m_e = 0.511*1e6 # eV
sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
speed_light = 3.*1e8 # m/s
U_cmb = 2.6*1e5 # eV/m3
U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3

#D0 = 1e26 # cm2/s
D0 = 1e24 # cm2/s
alpha = 1/3
gamma_index = 2.
#gamma_index = 1.3
#gamma_index = 1.

#pulsar_age_year = 327*1e3 # year
#proper_velocity = 400 #km/s
pulsar_age_year = 20*1e3 # year
proper_velocity = 700 #km/s

km_to_pc = 3.24078e-14
year_to_sec = 365.*24.*60.*60.
proper_velocity_pc_per_year = proper_velocity*km_to_pc*year_to_sec
Q0_normalization = injection_normalization()
print ('Q0_normalization = %0.2e'%(Q0_normalization))
source_location_start = proper_velocity_pc_per_year*pulsar_age_year
print ('source_location_initial = %0.3e'%(source_location_start))

fig, ax = plt.subplots()

E_e_min = m_e*pow(100.*1e9/E_cmb,0.5)/1e9 # GeV
E_e_max = m_e*pow(10000.*1e9/E_cmb,0.5)/1e9 # GeV
log_energy = np.linspace(log10(E_e_min),log10(E_e_max),50)
xdata = pow(10.,log_energy)
tdata = np.linspace(0.,pulsar_age_year,50)
rdata = np.linspace(-0.25*source_location_start,source_location_start*1.25,40)

fig.clf()
axbig = fig.add_subplot()
vectorize_IC_loss = np.vectorize(IC_electron_energy_loss_function)
vectorize_Sync_loss = np.vectorize(Sync_electron_energy_loss_function)
vectorize_total_loss = np.vectorize(total_electron_energy_loss_function)
ydata_IC = vectorize_IC_loss(xdata)
ydata_Sync = vectorize_Sync_loss(xdata)
ydata_total = vectorize_total_loss(xdata)
axbig.plot(xdata, ydata_total,'k-',label='total')
axbig.plot(xdata, ydata_IC,'r-',label='IC')
axbig.plot(xdata, ydata_Sync,'g-',label='Sync.')
axbig.legend(loc='best')
axbig.set_xlabel('electron energy [GeV]')
axbig.set_ylabel('energy loss [eV/s]')
axbig.set_xscale('log')
plotname = 'EnergyLoss'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

E_e_1TeV = m_e*pow(10.*1e9/E_cmb,0.5)/1e9 # GeV
E_e_10TeV = m_e*pow(600.*1e9/E_cmb,0.5)/1e9 # GeV
E_e_100TeV = m_e*pow(4000.*1e9/E_cmb,0.5)/1e9 # GeV
fig.clf()
axbig = fig.add_subplot()
vectorize_Ee = np.vectorize(initial_electron_energy)
ydata_1TeV = vectorize_Ee(tdata,E_e_1TeV)
axbig.plot(tdata, ydata_1TeV,'k-',label='E_{ph} = 10 GeV')
ydata_10TeV = vectorize_Ee(tdata,E_e_10TeV)
axbig.plot(tdata, ydata_10TeV,'g-',label='E_{ph} = 0.6 TeV')
ydata_100TeV = vectorize_Ee(tdata,E_e_100TeV)
axbig.plot(tdata, ydata_100TeV,'r-',label='E_{ph} = 4 TeV')
axbig.legend(loc='best')
axbig.set_xlabel('\delta t [year]')
axbig.set_ylabel('electron energy [GeV]')
plt.ylim(1e1, 1e7)
axbig.set_yscale('log')
plotname = 'ElectronEnergy'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

t_cooling = -E_e_1TeV/(total_electron_energy_loss_function(E_e_1TeV)/1e9)/(365.*24.*60.*60.)
effective_length = proper_velocity_pc_per_year*t_cooling
print ('t cooling = %0.5e year at 0.6 TeV'%(t_cooling))
print ('effective_length = %0.3e pc'%(effective_length))
t_cooling = -E_e_100TeV/(total_electron_energy_loss_function(E_e_100TeV)/1e9)/(365.*24.*60.*60.)
effective_length = proper_velocity_pc_per_year*t_cooling
print ('t cooling = %0.5e year at 4 TeV'%(t_cooling))
print ('effective_length = %0.3e pc'%(effective_length))

fig.clf()
axbig = fig.add_subplot()
vectorize_L_ism = np.vectorize(diffusion_length)
ydata_1TeV = vectorize_L_ism(tdata,E_e_1TeV)
axbig.plot(tdata, ydata_1TeV,'k-',label='E_{ph} = 10 GeV')
ydata_10TeV = vectorize_L_ism(tdata,E_e_10TeV)
axbig.plot(tdata, ydata_10TeV,'g-',label='E_{ph} = 0.6 TeV')
ydata_100TeV = vectorize_L_ism(tdata,E_e_100TeV)
axbig.plot(tdata, ydata_100TeV,'r-',label='E_{ph} = 4 TeV')
axbig.legend(loc='best')
axbig.set_xlabel('\delta t [year]')
axbig.set_ylabel('diffusion length [pc]')
plotname = 'DiffusionLength'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
vectorize_density = np.vectorize(electron_column_density_continuum_distance)
ydata_1TeV = vectorize_density(rdata,pulsar_age_year,E_e_1TeV)
axbig.plot(rdata, 0.001*ydata_1TeV,'k-',label='E_{ph} = 10 GeV')
ydata_10TeV = vectorize_density(rdata,pulsar_age_year,E_e_10TeV)
axbig.plot(rdata, ydata_10TeV,'g-',label='E_{ph} = 0.6 TeV')
ydata_100TeV = 10.*vectorize_density(rdata,pulsar_age_year,E_e_100TeV)
axbig.plot(rdata, ydata_100TeV,'r-',label='E_{ph} = 4 TeV')
axbig.legend(loc='best')
axbig.set_xlabel('distance [pc]')
axbig.set_ylabel('electron density')
#axbig.set_yscale('log')
plotname = 'ElectronDensityContinuum'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
