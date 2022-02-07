
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

def injection_normalization():
    E1_GeV = 0.1
    E_cutoff_GeV = 400.*1e3
    E_max_GeV = 10.*E_cutoff_GeV
    integral = quad(injection_spectrum, E1_GeV, E_max_GeV)[0] # GeV
    erg_to_GeV = 624.151
    year_to_sec = 365.*24.*60.*60.
    spin_down_luminosity = 1.7*1e37*erg_to_GeV*year_to_sec # GeV/year
    Q0 = spin_down_luminosity/integral # 1./GeV/year
    return Q0

def injection_spectrum(E_e_GeV):
    E_cutoff_GeV = 400.*1e3
    return E_e_GeV*pow(E_e_GeV,-gamma_index)*np.exp(-E_e_GeV/E_cutoff_GeV)

def injection_spectrum_time_dep(E_e_GeV,t_year=0.):
    tau0 = 15.*1e3 # year
    spectrum = Q0_normalization*pow(1+t_year/tau0,-2)*injection_spectrum(E_e_GeV)/E_e_GeV
    return spectrum # in unit of Q0, 1./GeV/year

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
    year_to_sec = 365.*24.*60.*60.
    E_e_GeV = initial_electron_energy(dt_year,E_e_final_GeV)
    return D0*pow(E_e_GeV/1000.,alpha)*year_to_sec # cm2/year

def diffusion_length(dt_year,E_e_final_GeV):
    pc_to_cm = 3.086e+18
    year_to_sec = 365.*24.*60.*60.
    dt_year_limit = 1./(4./3.*sigma_thomson*speed_light*pow(m_e,-2)*(U_cmb+U_B)*E_e_final_GeV*1e9)/year_to_sec
    if dt_year>dt_year_limit:
        dt_year = dt_year_limit
    length = 2.*pow(quad(diffusion_coefficient, 0., dt_year, args=(E_e_final_GeV))[0],0.5)/pc_to_cm
    return length # pc

def electron_density_burst(t0_year,photon_location,source_location_initial,source_location_final,E_e_final_GeV):
    # electrons released at t0, radiate photons at t_now = pulsar age
    year_to_sec = 365.*24.*60.*60.
    dt_year = pulsar_age_year-t0_year
    dt_year_limit = 1./(4./3.*sigma_thomson*speed_light*pow(m_e,-2)*(U_cmb+U_B)*E_e_final_GeV*1e9)/year_to_sec-1.
    if dt_year>dt_year_limit:
        return 0.
    source_location_at_t_0 = (source_location_initial-source_location_final)/pulsar_age_year*dt_year+source_location_final
    distance_to_source = pow(pow(photon_location[0]-source_location_at_t_0[0],2)+pow(photon_location[1]-source_location_at_t_0[1],2)+pow(photon_location[2],2),0.5)
    distance_to_source_physical = distance_to_source*3.14/180.*pulsar_distance
    diffusion_length_at_t_now = diffusion_length(dt_year,E_e_final_GeV)
    #km_to_pc = 3.24078e-14
    #larmor_radius = 33.36*E_e_final_GeV*(1./(mag_field*1e-6))*km_to_pc
    #diffusion_length_at_t_now = max(diffusion_length_at_t_now,larmor_radius)
    PSF_length = 0.1*3.14/180.*pulsar_distance
    diffusion_length_at_t_now = pow(pow(diffusion_length_at_t_now,2)+pow(PSF_length,2),0.5)
    energy_loss_at_t_now = total_electron_energy_loss_function(E_e_final_GeV)
    E_e_initial_GeV = initial_electron_energy(dt_year,E_e_final_GeV)
    energy_loss_at_t0 = total_electron_energy_loss_function(E_e_initial_GeV)
    Qt = injection_spectrum_time_dep(E_e_initial_GeV,t0_year) # 1./GeV/year
    electron_density_rate = Qt*np.exp(-pow(distance_to_source_physical,2)/pow(diffusion_length_at_t_now,2))*pow(2.*3.14,-3/2)*pow(diffusion_length_at_t_now,-3)*energy_loss_at_t0/energy_loss_at_t_now # 1./pc^3/GeV/year
    return electron_density_rate

def electron_density_continuum(photon_location_z,photon_location_x,photon_location_y,t_age,E_e_final_GeV):
    source_travel_distance = proper_velocity_pc_per_year*pulsar_age_year/pulsar_distance*180./3.14
    source_location_ref = np.array([PSR_tail_x,PSR_tail_y])
    source_location_final = np.array([PSR_head_x,PSR_head_y])
    source_ref_dx = source_location_ref[0]-source_location_final[0]
    source_ref_dy = source_location_ref[1]-source_location_final[1]
    source_ref_distance = pow(source_ref_dx*source_ref_dx+source_ref_dy*source_ref_dy,0.5)
    source_location_initial = np.array([source_location_final[0]+source_ref_dx*source_travel_distance/source_ref_distance,source_location_final[1]+source_ref_dy*source_travel_distance/source_ref_distance])
    photon_location = np.array([photon_location_x,photon_location_y,photon_location_z])
    electron_density = quad(electron_density_burst, 0., 1.0*t_age, args=(photon_location,source_location_initial,source_location_final,E_e_final_GeV))[0]
    return electron_density # 1./pc^3/GeV

def electron_column_density_continuum(E_e_final_GeV,photon_location_x,photon_location_y,t_age):
    x_low = 0.
    x_up = 2.*diffusion_length(t_age,E_e_final_GeV)/pulsar_distance*180./3.14 # degree
    electron_column_density = quad(electron_density_continuum, x_low, x_up, args=(photon_location_x,photon_location_y,t_age,E_e_final_GeV))[0]
    return 2.*electron_column_density # 1./pc^2/GeV

def electron_column_density_integrated(photon_location_x,photon_location_y,t_age,E_ph_GeV_low,E_ph_GeV_up):

    E_ph_GeV_avg = (E_ph_GeV_low+E_ph_GeV_up)/2.
    E_e_GeV_avg = m_e*pow(E_ph_GeV_avg*1e9/E_cmb,0.5)/1e9 # GeV
    E_e_GeV_up = m_e*pow(E_ph_GeV_up*1e9/E_cmb,0.5)/1e9 # GeV
    E_e_GeV_low = m_e*pow(E_ph_GeV_low*1e9/E_cmb,0.5)/1e9 # GeV
    electron_column_density = quad(electron_column_density_continuum, E_e_GeV_low, E_e_GeV_up, args=(photon_location_x,photon_location_y,t_age))[0]
    return electron_column_density # el/pc^2

def IC_photon_rate_column_density_integrated(pixel_size_deg2,E_ph_GeV_low,E_ph_GeV_up):

    E_ph_GeV_avg = (E_ph_GeV_low+E_ph_GeV_up)/2.
    E_e_GeV_avg = m_e*pow(E_ph_GeV_avg*1e9/E_cmb,0.5)/1e9 # GeV
    E_e_GeV_up = m_e*pow(E_ph_GeV_up*1e9/E_cmb,0.5)/1e9 # GeV
    E_e_GeV_low = m_e*pow(E_ph_GeV_low*1e9/E_cmb,0.5)/1e9 # GeV

    E1_factor_hat_up = E_ph_GeV_to_E1_factor_hat(E_ph_GeV_up,E_e_GeV_avg)
    E1_factor_hat_low = E_ph_GeV_to_E1_factor_hat(E_ph_GeV_low,E_e_GeV_avg)
    IC_photon_rate = IC_photon_rate_per_electron_integrated(E1_factor_hat_low,E1_factor_hat_up,E_e_GeV_avg) # 1/sec/el

    area_per_pixel = pixel_size_deg2*pow(3.14/180.*pulsar_distance,2) # pc^2
    pc_to_cm = 3.086e+18
    solid_angle_per_cm2_from_source = 1./(pow(pulsar_distance*pc_to_cm,2))  # 1/cm^2

    return IC_photon_rate*area_per_pixel*solid_angle_per_cm2_from_source # pc^2/sec/el/cm^2

def E1_factor_hat_to_E_ph_GeV(E1_factor_hat,E_e_GeV=1e5):

    gamma_factor_electron = E_e_GeV*1e9/m_e
    gamma_factor_cmb = 4.*E_cmb*gamma_factor_electron/m_e
    E1_factor = E1_factor_hat*gamma_factor_cmb/(1+gamma_factor_cmb)
    E_ph_GeV = E_e_GeV*E1_factor
    return E_ph_GeV

def E_ph_GeV_to_E1_factor_hat(E_ph_GeV,E_e_GeV=1e5):

    E1_factor = E_ph_GeV/E_e_GeV
    gamma_factor_electron = E_e_GeV*1e9/m_e
    gamma_factor_cmb = 4.*E_cmb*gamma_factor_electron/m_e
    E1_factor_hat = E1_factor*(1+gamma_factor_cmb)/gamma_factor_cmb
    return E1_factor_hat

def IC_photon_rate_per_electron(E1_factor_hat,E_e_GeV=1e5):

    gamma_factor_electron = E_e_GeV*1e9/m_e
    gamma_factor_cmb = 4.*E_cmb*gamma_factor_electron/m_e
    E1_factor = E1_factor_hat*gamma_factor_cmb/(1+gamma_factor_cmb)

    q_factor = E1_factor/(gamma_factor_cmb*(1.-E1_factor))
    F_factor = 2.*q_factor*np.log(q_factor)+(1.+2.*q_factor)*(1.-q_factor)+0.5*pow(q_factor*gamma_factor_cmb,2)/(1.+q_factor*gamma_factor_cmb)*(1.-q_factor)
    IC_photon_rate = 2*3.14*pow(radius_thomson,2)*m_e*speed_light*cmb_ph_density/(gamma_factor_electron*E_cmb)*F_factor #1/sec
    #return F_factor
    return IC_photon_rate # 1/sec/dE1

def IC_photon_rate_per_electron_integrated(E1_factor_hat_low,E1_factor_hat_up,E_e_GeV=1e5):
    IC_photon_rate = quad(IC_photon_rate_per_electron, E1_factor_hat_low, E1_factor_hat_up, args=(E_e_GeV))[0]
    return IC_photon_rate # 1/sec

def PlotAFunction(function,x_low,x_up,plot_name,logx=True,logy=True):
    x_axis = np.linspace(x_low,x_up,100)
    if logx:
        log_x_axis = np.linspace(log10(x_low),log10(x_up),100)
        x_axis = pow(10.,log_x_axis)
    fig.clf()
    axbig = fig.add_subplot()
    vectorized_func = np.vectorize(function)
    y_axis = vectorized_func(x_axis)
    axbig.plot(x_axis, y_axis,'k-')
    if logx: axbig.set_xscale('log')
    if logy: axbig.set_yscale('log')
    fig.savefig("output_plots/%s.png"%(plot_name),bbox_inches='tight')
    axbig.remove()

def PlotElectronColummnDensity(E_e_GeV_low,E_e_GeV_up,plot_tag):

    source_travel_distance = proper_velocity_pc_per_year*pulsar_age_year/pulsar_distance*180./3.14
    source_location_ref = np.array([PSR_tail_x,PSR_tail_y])
    source_location_final = np.array([PSR_head_x,PSR_head_y])
    source_ref_dx = source_location_ref[0]-source_location_final[0]
    source_ref_dy = source_location_ref[1]-source_location_final[1]
    source_ref_distance = pow(source_ref_dx*source_ref_dx+source_ref_dy*source_ref_dy,0.5)
    source_location_initial = np.array([source_location_final[0]+source_ref_dx*source_travel_distance/source_ref_distance,source_location_final[1]+source_ref_dy*source_travel_distance/source_ref_distance])

    Source_RA = 286.975
    Source_Dec = 6.269
    Skymap_nbins = 60
    MapEdge_left = Source_RA-3.
    MapEdge_right = Source_RA+3.
    MapEdge_lower = Source_Dec-3.
    MapEdge_upper = Source_Dec+3.
    pixel_size = ((MapEdge_right-MapEdge_left)/Skymap_nbins)*((MapEdge_upper-MapEdge_lower)/Skymap_nbins)

    x_axis = np.linspace(MapEdge_left,MapEdge_right,Skymap_nbins)
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,Skymap_nbins)
    grid_z = []
    IC_photon_rate = IC_photon_rate_column_density_integrated(pixel_size,E_e_GeV_low,E_e_GeV_up)
    for ybin in range(0,len(y_axis)):
        z_along_x_axis = []
        for xbin in range(0,len(x_axis)):
            print ('xbin = %s, ybin = %s'%(xbin,ybin))
            x = x_axis[xbin]
            y = y_axis[ybin]
            electron_column_density = electron_column_density_integrated(x,y,pulsar_age_year,E_e_GeV_low,E_e_GeV_up)
            z_along_x_axis += [electron_column_density*IC_photon_rate]
        grid_z += [z_along_x_axis]
    fig.clf()
    axbig = fig.add_subplot()
    axbig.set_xlabel('X')
    axbig.set_ylabel('Y')
    im = axbig.imshow(grid_z, origin='lower', cmap='gray', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()))
    cbar = fig.colorbar(im)
    cbar.set_label('electron density [$pc^{-2}$]')
    axbig.scatter(source_location_final[0], source_location_final[1], marker='^', c='r')
    axbig.scatter(source_location_initial[0], source_location_initial[1], marker='^', c='r')
    plotname = 'ElectronColumnDensity_%s'%(plot_tag)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

E_cmb = 6.6*1e-4 # eV
m_e = 0.511*1e6 # eV
sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
radius_thomson = pow(sigma_thomson*3./(8.*3.14),0.5) # m
speed_light = 3.*1e8 # m/s
U_cmb = 2.6*1e5 # eV/m3
cmb_ph_density = 440.*pow(100.,3) # 1/m3

#mag_field = 1. # muG
mag_field = 3. # muG
#mag_field = 6. # muG
U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3

D0 = 8.2*1e26 # cm2/s
alpha = 1/3
proper_velocity = 2000 #km/s
km_to_pc = 3.24078e-14
year_to_sec = 365.*24.*60.*60.
proper_velocity_pc_per_year = proper_velocity*km_to_pc*year_to_sec


gamma_index = 2.
pulsar_age_year = 20*1e3 # year
pulsar_distance = 3.2*1000. # pc
PSR_head_x = 286.98
PSR_head_y = 6.04
PSR_tail_x = 287.09
PSR_tail_y = 6.86

Q0_normalization = injection_normalization()
print ('Q0_normalization = %0.2e / year'%(Q0_normalization))

fig, ax = plt.subplots()
PlotAFunction(injection_spectrum_time_dep,1e3,1e5,'injection_spectrum')
PlotAFunction(IC_photon_rate_per_electron,0.01,0.99,'IC_photon_rate_per_electron',logx=False,logy=False)
PlotAFunction(E1_factor_hat_to_E_ph_GeV,0.,1.,'E1_factor_hat_to_E_ph_GeV',logx=False,logy=False)
PlotAFunction(E_ph_GeV_to_E1_factor_hat,1e2,1e4,'E_ph_GeV_to_E1_factor_hat',logy=False)
PlotElectronColummnDensity(600,1500,'1TeV')

