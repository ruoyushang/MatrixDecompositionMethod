
import sys,ROOT
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
    return D0*pow(E_e_GeV/1000.,alpha)*year_to_sec # cm2/year

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
    distance_to_source = pow(pow(photon_location[0]-source_location_at_t_0[0],2)+pow(photon_location[1]-source_location_at_t_0[1],2)+pow(photon_location[2],2),0.5)
    distance_to_source_physical = distance_to_source*3.14/180.*pulsar_distance
    diffusion_length_at_t_now = diffusion_length(dt_year,E_e_final_GeV)
    energy_loss_at_t_now = total_electron_energy_loss_function(E_e_final_GeV)
    E_e_initial_GeV = initial_electron_energy(dt_year,E_e_final_GeV)
    energy_loss_at_t0 = total_electron_energy_loss_function(E_e_initial_GeV)
    Qt = injection_spectrum_time_dep(E_e_initial_GeV,t0_year)
    electron_density = Qt*np.exp(-pow(distance_to_source_physical,2)/pow(diffusion_length_at_t_now,2))*pow(2.*3.14,-3/2)*pow(diffusion_length_at_t_now,-3)*energy_loss_at_t0/energy_loss_at_t_now
    return electron_density

def electron_density_continuum(photon_location_z,photon_location_x,photon_location_y,t0_year,E_e_final_GeV):
    source_travel_distance = proper_velocity_pc_per_year*pulsar_age_year/pulsar_distance*180./3.14
    source_location_ref = np.array([PSR_tail_x,PSR_tail_y])
    source_location_final = np.array([PSR_head_x,PSR_head_y])
    source_ref_dx = source_location_ref[0]-source_location_final[0]
    source_ref_dy = source_location_ref[1]-source_location_final[1]
    source_ref_distance = pow(source_ref_dx*source_ref_dx+source_ref_dy*source_ref_dy,0.5)
    source_location_initial = np.array([source_location_final[0]+source_ref_dx*source_travel_distance/source_ref_distance,source_location_final[1]+source_ref_dy*source_travel_distance/source_ref_distance])
    photon_location = np.array([photon_location_x,photon_location_y,photon_location_z])
    electron_density = quad(electron_density_burst, 0., t0_year, args=(photon_location,source_location_initial,source_location_final,E_e_final_GeV))[0]
    return electron_density

def electron_column_density_continuum(photon_location_x,photon_location_y,t0_year,E_e_final_GeV):
    x_low = 0.
    x_up = 2.*diffusion_length(t0_year,E_e_final_GeV)/pulsar_distance*180./3.14
    electron_column_density = quad(electron_density_continuum, x_low, x_up, args=(photon_location_x,photon_location_y,t0_year,E_e_final_GeV))[0]
    return 2.*electron_column_density

def injection_spectrum(E_e_GeV):
    E_cutoff_GeV = 400.*1e3
    return E_e_GeV*pow(E_e_GeV,-gamma_index)*np.exp(-E_e_GeV/E_cutoff_GeV)

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

def smooth_map(x_ax,y_ax,z_map):

    z_map_smooth = []
    for ybin in range(0,len(y_ax)):
        z_along_x_axis = []
        for xbin in range(0,len(x_ax)):
            x = x_ax[xbin]
            y = y_ax[ybin]
            z_sum = 0.
            for ybin2 in range(0,len(y_ax)):
                for xbin2 in range(0,len(x_ax)):
                    x2 = x_ax[xbin2]
                    y2 = y_ax[ybin2]
                    z_sum += np.exp(-0.5*(pow(x-x2,2)+pow(y-y2,2))/pow(0.1,2))*z_map[ybin2][xbin2]
            z_along_x_axis += [z_sum]
        z_map_smooth += [z_along_x_axis]
    return z_map_smooth

def GetFluxMapData():

    InputFile = ROOT.TFile('output_pulsar/FluxSkymap.root')
    for ebin in range(0,len(energy_bin)-1):
        HistName = "hist_energy_flux_PWN_skymap_%s"%(ebin)
        hist_energy_flux_PWN_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_energy_flux_skymap_%s"%(ebin)
        hist_energy_flux_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_data_skymap_%s"%(ebin)
        hist_data_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_bkgd_skymap_%s"%(ebin)
        hist_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_syst_skymap_%s"%(ebin)
        hist_syst_skymap[ebin].Add(InputFile.Get(HistName))

def PlotDataMaps(name_tag):

    Skymap_nbins = 20
    zoomin_scale = 1
    MapCenter_x = PSR_head_x
    MapCenter_y = PSR_head_y
    Skymap_size = 2.
    hist_skymap = ROOT.TH2D("hist_skymap","",Skymap_nbins,MapCenter_x-Skymap_size,MapCenter_x+Skymap_size,Skymap_nbins,MapCenter_y-Skymap_size,MapCenter_y+Skymap_size)

    MapEdge_left  = hist_skymap.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_skymap.GetXaxis().GetBinLowEdge(Skymap_nbins+1)
    MapEdge_lower = hist_skymap.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_skymap.GetYaxis().GetBinLowEdge(Skymap_nbins+1)

    x_axis = np.linspace(MapEdge_left,MapEdge_right,Skymap_nbins)
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,Skymap_nbins)

    source_travel_distance = proper_velocity_pc_per_year*pulsar_age_year/pulsar_distance*180./3.14
    source_location_ref = np.array([PSR_tail_x,PSR_tail_y])
    source_location_final = np.array([PSR_head_x,PSR_head_y])
    source_ref_dx = source_location_ref[0]-source_location_final[0]
    source_ref_dy = source_location_ref[1]-source_location_final[1]
    source_ref_distance = pow(source_ref_dx*source_ref_dx+source_ref_dy*source_ref_dy,0.5)
    source_location_initial = np.array([source_location_final[0]+source_ref_dx*source_travel_distance/source_ref_distance,source_location_final[1]+source_ref_dy*source_travel_distance/source_ref_distance])

    E_ph = 0.6*1e3
    E_e = m_e*pow(E_ph*1e9/E_cmb,0.5)/1e9 # GeV
    grid_z = []
    for ybin in range(0,len(y_axis)):
        z_along_x_axis = []
        for xbin in range(0,len(x_axis)):
            x = x_axis[xbin]
            y = y_axis[ybin]
            z = electron_column_density_continuum(x,y,pulsar_age_year,E_e)
            z_along_x_axis += [z]
        grid_z += [z_along_x_axis]
    grid_z_smooth = smooth_map(x_axis,y_axis,grid_z)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.set_xlabel('X')
    axbig.set_ylabel('Y')
    axbig.imshow(grid_z_smooth, origin='lower', cmap='gray', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()))
    axbig.scatter(source_location_final[0], source_location_final[1], marker='^', c='r')
    axbig.scatter(source_location_initial[0], source_location_initial[1], marker='^', c='r')
    plotname = 'ElectronColumnDensity_0p6TeV_%s'%(name_tag)
    fig.savefig("output_pulsar/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

    E_ph = 4.0*1e3
    E_e = m_e*pow(E_ph*1e9/E_cmb,0.5)/1e9 # GeV
    grid_z = []
    for ybin in range(0,len(y_axis)):
        z_along_x_axis = []
        for xbin in range(0,len(x_axis)):
            x = x_axis[xbin]
            y = y_axis[ybin]
            z = electron_column_density_continuum(x,y,pulsar_age_year,E_e)
            z_along_x_axis += [z]
        grid_z += [z_along_x_axis]
    grid_z_smooth = smooth_map(x_axis,y_axis,grid_z)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.set_xlabel('X')
    axbig.set_ylabel('Y')
    axbig.imshow(grid_z_smooth, origin='lower', cmap='gray', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()))
    axbig.scatter(source_location_final[0], source_location_final[1], marker='^', c='r')
    axbig.scatter(source_location_initial[0], source_location_initial[1], marker='^', c='r')
    plotname = 'ElectronColumnDensity_4TeV_%s'%(name_tag)
    fig.savefig("output_pulsar/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()


E_cmb = 6.6*1e-4 # eV
m_e = 0.511*1e6 # eV
sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
speed_light = 3.*1e8 # m/s
U_cmb = 2.6*1e5 # eV/m3

#mag_field = 1. # muG
mag_field = 3. # muG
#mag_field = 6. # muG
U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3

D0 = 8.2*1e26 # cm2/s
alpha = 1/3

gamma_index = 2.
#gamma_index = 1.3
#gamma_index = 1.

#pulsar_age_year = 327*1e3 # year
#proper_velocity = 400 #km/s
pulsar_age_year = 20*1e3 # year
pulsar_distance = 3.2*1000. # pc

km_to_pc = 3.24078e-14
year_to_sec = 365.*24.*60.*60.
Q0_normalization = injection_normalization()
print ('Q0_normalization = %0.2e'%(Q0_normalization))

fig, ax = plt.subplots()

PSR_head_x = 286.98
PSR_head_y = 6.04
PSR_tail_x = 287.09
PSR_tail_y = 6.86
#PSR_head_x = 0.
#PSR_head_y = 0.
#PSR_tail_x = 0.
#PSR_tail_y = 1.


D0 = 8.2*1e26 # cm2/s
alpha = 1/3
proper_velocity = 2000 #km/s
proper_velocity_pc_per_year = proper_velocity*km_to_pc*year_to_sec
PlotDataMaps('high_diff')

D0 = 8.2*1e27 # cm2/s
alpha = 1/3
proper_velocity = 2000 #km/s
proper_velocity_pc_per_year = proper_velocity*km_to_pc*year_to_sec
PlotDataMaps('low_diff')


E_e_min = m_e*pow(100.*1e9/E_cmb,0.5)/1e9 # GeV
E_e_max = m_e*pow(10000.*1e9/E_cmb,0.5)/1e9 # GeV
log_energy = np.linspace(log10(E_e_min),log10(E_e_max),50)

E_e_10TeV = m_e*pow(600.*1e9/E_cmb,0.5)/1e9 # GeV
E_e_100TeV = m_e*pow(4000.*1e9/E_cmb,0.5)/1e9 # GeV

e_axis = pow(10.,log_energy)
t_axis = np.linspace(0.,pulsar_age_year,50)

fig.clf()
axbig = fig.add_subplot()
vectorize_IC_loss = np.vectorize(IC_electron_energy_loss_function)
vectorize_Sync_loss = np.vectorize(Sync_electron_energy_loss_function)
vectorize_total_loss = np.vectorize(total_electron_energy_loss_function)
ydata_IC = vectorize_IC_loss(e_axis)
ydata_Sync = vectorize_Sync_loss(e_axis)
ydata_total = vectorize_total_loss(e_axis)
axbig.plot(e_axis, ydata_total,'k-',label='total')
axbig.plot(e_axis, ydata_IC,'r-',label='IC')
axbig.plot(e_axis, ydata_Sync,'g-',label='Sync.')
axbig.legend(loc='best')
axbig.set_xlabel('electron energy [GeV]')
axbig.set_ylabel('energy loss [eV/s]')
axbig.set_xscale('log')
plotname = 'EnergyLoss'
fig.savefig("output_pulsar/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
vectorize_Ee = np.vectorize(initial_electron_energy)
ydata_10TeV = vectorize_Ee(t_axis,E_e_10TeV)
axbig.plot(t_axis, ydata_10TeV,'g-',label='E_{ph} = 0.6 TeV')
ydata_100TeV = vectorize_Ee(t_axis,E_e_100TeV)
axbig.plot(t_axis, ydata_100TeV,'r-',label='E_{ph} = 4 TeV')
axbig.legend(loc='best')
axbig.set_xlabel('\delta t [year]')
axbig.set_ylabel('electron energy [GeV]')
plt.ylim(1e1, 1e7)
axbig.set_yscale('log')
plotname = 'ElectronEnergy'
fig.savefig("output_pulsar/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()



