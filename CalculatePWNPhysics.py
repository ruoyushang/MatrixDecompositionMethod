import matplotlib.pyplot as plt
import math
import numpy as np

def FluxModel(input_t,input_d,input_Edot):

    #t_PSR = 342.0 # kyr
    #d_PSR = 0.25 # kpc # fit diffusion radius of 1.8 deg, using D = (8.2-5.9)*1e26
    #Edot = 27.8e+36 # erg/s
    #t_PSR = 20. # kyr
    #d_PSR = 3.0 # kpc # fit diffusion radius of 1.8 deg, using D = (8.2-5.9)*1e26
    #Edot = 2.8e+36 # erg/s
    t_PSR = input_t
    d_PSR = input_d
    Edot = input_Edot
    
    B_field = 5.0 # muG
    #B_field = 3.0 # muG
    #B_field = 1.2 # muG
    
    E_e_cutoff = 1000. # TeV
    el_spec_index = 2.2
    CMB_energy = 6.6e-4 # eV
    
    print ('math.log(Edot) = %0.2f'%(math.log10(Edot)))
    IC_luminosity = pow(10.,19.5+0.41*math.log10(Edot))
    print ('IC_luminosity = %0.2e erg/s'%(IC_luminosity))
    
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    kyr_to_sec = 1000.*365.*24.*60.*60.
    E_cmb = 6.6*1e-4 # eV
    m_e = 0.511*1e6 # eV
    sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
    speed_light = 3.*1e8 # m/s
    U_cmb = 2.6*1e5 # eV/m3
    U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(B_field*1e-6/1e4,2) # eV/m3
    gamma_factor = 1./(t_PSR*kyr_to_sec/m_e*speed_light*(4./3.*sigma_thomson*(U_cmb+U_B)))
    
    E_e_break = gamma_factor*m_e/1e12 # TeV
    E_sync_break = 0.2*(B_field/10.)*pow(E_e_break/1.,2) # eV
    E_IC_break = 5.*(CMB_energy/1e-3)*pow(E_e_break/1.,2)*1e9 # eV
    print ('E_e_break = %0.2f TeV'%(E_e_break))
    print ('E_sync_break = %0.2f eV'%(E_sync_break))
    print ('E_IC_break = %0.2f GeV'%(E_IC_break/1e9))
    
    E_sync_cutoff = 0.2*(B_field/10.)*pow(E_e_cutoff/1.,2) # eV
    E_IC_cutoff = 5.*(CMB_energy/1e-3)*pow(E_e_cutoff/1.,2)*1e9 # eV
    print ('E_sync_cutoff = %0.2f keV'%(E_sync_cutoff/1e3))
    print ('E_IC_cutoff = %0.2f TeV'%(E_IC_cutoff/1e12))
    
    diffusion_index = 1./3.
    diffusion_coefficient_break = 1e28*pow((E_e_break*1e12/1e9)/10.,diffusion_index)*pow(B_field/3.,-0.5) # cm2/s
    diffusion_radius_break = pow(t_PSR*kyr_to_sec*4.*diffusion_coefficient_break,0.5)/deg_to_cm
    print ('leptonic diffusion_radius_break = %0.2f deg'%(diffusion_radius_break))
    
    gamma_factor_cutoff = E_e_cutoff*1e12/m_e
    t_cooling_cutoff = m_e/speed_light/(4./3.*sigma_thomson*gamma_factor_cutoff*(U_cmb+U_B)) # sec
    print ('t_cooling_cutoff = %0.2f kyr'%(t_cooling_cutoff/kyr_to_sec))
    diffusion_coefficient_cutoff = 1e28*pow((E_e_cutoff*1e12/1e9)/10.,diffusion_index)*pow(B_field/3.,-0.5) # cm2/s
    diffusion_radius_cutoff = pow(t_cooling_cutoff*4.*diffusion_coefficient_cutoff,0.5)/deg_to_cm
    print ('leptonic diffusion_radius_cutoff = %0.2f deg'%(diffusion_radius_cutoff))
    
    lumi_ratio = 0.1*pow(B_field/10.,-2)
    print ('L_{gamma}/L_{X} = %0.2f'%(lumi_ratio))
    
    E_IC_start = E_IC_break*1e-4
    E_sync_start = E_IC_start/(5.*(CMB_energy/1e-3))*(0.2*(B_field/10.))*1e-9
    E_IC_end = E_IC_cutoff*1e1
    E_sync_end = E_IC_end/(5.*(CMB_energy/1e-3))*(0.2*(B_field/10.))*1e-9

    IC_flux = IC_luminosity/(4.*3.14*pow(d_PSR*1000.*3.086e18,2))/pow(E_IC_break/1e12,2)
    print ('IC_flux = %0.2e /erg/cm2/s'%(IC_flux))
    sync_luminosity = IC_luminosity/lumi_ratio
    sync_flux = sync_luminosity/(4.*3.14*pow(d_PSR*1000.*3.086e18,2))/pow(E_sync_break/1e12,2)
    print ('sync_flux = %0.2e /erg/cm2/s'%(sync_flux))

    flux_IC_break = IC_flux
    flux_sync_break = sync_flux
    flux_IC_start = flux_IC_break*pow(E_IC_start/E_IC_break,-(el_spec_index+1.)/2.)
    flux_sync_start = flux_sync_break*pow(E_sync_start/E_sync_break,-(el_spec_index+1.)/2.)
    flux_IC_cutoff = flux_IC_break*pow(E_IC_cutoff/E_IC_break,-(el_spec_index+2.)/2.)
    flux_sync_cutoff = flux_sync_break*pow(E_sync_cutoff/E_sync_break,-(el_spec_index+2.)/2.)
    flux_IC_end = flux_IC_cutoff*math.exp(-E_IC_end/E_IC_cutoff)
    flux_sync_end = flux_sync_cutoff*math.exp(-E_sync_end/E_sync_cutoff)
    
    IC_energy_axis = [E_IC_start,E_IC_break,E_IC_cutoff,E_IC_end]
    IC_flux_axis = [flux_IC_start,flux_IC_break,flux_IC_cutoff,flux_IC_end]
    Sync_energy_axis = [E_sync_start,E_sync_break,E_sync_cutoff,E_sync_end]
    Sync_flux_axis = [flux_sync_start,flux_sync_break,flux_sync_cutoff,flux_sync_end]
    IC_energy_axis = np.array(IC_energy_axis)
    IC_flux_axis = np.array(IC_flux_axis)
    Sync_energy_axis = np.array(Sync_energy_axis)
    Sync_flux_axis = np.array(Sync_flux_axis)

    return IC_energy_axis, IC_flux_axis*(IC_energy_axis/1e12)*(IC_energy_axis/1e12), Sync_energy_axis, Sync_flux_axis*(Sync_energy_axis/1e12)*(Sync_energy_axis/1e12)

IC_energy_axis_1kyr, IC_flux_axis_1kyr, Sync_energy_axis_1kyr, Sync_flux_axis_1kyr = FluxModel(1.,3.,1e36)
IC_energy_axis_10kyr, IC_flux_axis_10kyr, Sync_energy_axis_10kyr, Sync_flux_axis_10kyr = FluxModel(10.,3.,1e36)
IC_energy_axis_100kyr, IC_flux_axis_100kyr, Sync_energy_axis_100kyr, Sync_flux_axis_100kyr = FluxModel(100.,3.,1e36)

fig, ax = plt.subplots()
fig.clf()
axbig = fig.add_subplot()
axbig.plot(IC_energy_axis_1kyr, IC_flux_axis_1kyr,'b-',label='1 kyr')
axbig.plot(Sync_energy_axis_1kyr, Sync_flux_axis_1kyr,'b-')
axbig.plot(IC_energy_axis_10kyr, IC_flux_axis_10kyr,'g-',label='10 kyr')
axbig.plot(Sync_energy_axis_10kyr, Sync_flux_axis_10kyr,'g-')
axbig.plot(IC_energy_axis_100kyr, IC_flux_axis_100kyr,'r-',label='100 kyr')
axbig.plot(Sync_energy_axis_100kyr, Sync_flux_axis_100kyr,'r-')
axbig.axvspan(1e2, 1e3, alpha=0.2, color='orange')
axbig.axvspan(1e3, 1e5, alpha=0.2, color='red')
axbig.axvspan(0.5*1e7, 1e11, alpha=0.2, color='green')
axbig.axvspan(1e11, 1e13, alpha=0.2, color='blue')
axbig.axvspan(1e13, 1e15, alpha=0.2, color='purple')
axbig.legend(loc='best')
axbig.set_xlabel('photon energy [eV]')
axbig.set_ylabel('E^{2} Flux')
plotname = 'SimpleSpectrum'
axbig.set_xscale('log')
axbig.set_yscale('log')
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
