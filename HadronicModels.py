
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
from scipy import special
from scipy import interpolate
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import cycle
from scipy.integrate import quad

import CommonPlotFunctions

energy_index = CommonPlotFunctions.energy_index_scale

def SNR_radius_Sedov(t_age_year):

    radius = 1.54*1e19*pow(E_SN,1./5.)*pow(n_0,-1./5.)*pow(t_age_year/1000.,2./5.) # cm
    return radius

def SNR_age_Sedov(angular_radius,d_pc):

    radius_cm = angular_radius*(3.14/180.*d_pc*3.086e18)
    t_age_year = 1000.*pow(radius_cm/(1.54*1e19*pow(E_SN,1./5.)*pow(n_0,-1./5.)),5./2.)
    return t_age_year

def CR_energy_density(photon_z_cm,photon_x_deg,photon_y_deg,photon_E_low_GeV,photon_E_up_GeV,SN_data,D0_had):

    RA_SN = SN_data[0]
    Dec_SN = SN_data[1]
    d_SN = SN_data[2] # pc
    t_SN = SN_data[3] # year
    year_to_sec = 365.*24.*60.*60.
    deg_to_cm = 3.14/180.*d_SN*3.086e18

    photon_to_SN_cm = pow(pow((photon_x_deg-RA_SN)*deg_to_cm,2)+pow((photon_y_deg-Dec_SN)*deg_to_cm,2)+pow(photon_z_cm,2),0.5)

    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    photon_E_GeV = (photon_E_low_GeV+photon_E_up_GeV)/2.
    proton_E_GeV = (1./0.17)*2.*photon_E_GeV
    proton_E_up_GeV = (1./0.17)*2.*photon_E_up_GeV
    proton_E_low_GeV = (1./0.17)*2.*photon_E_low_GeV
    D0_at_E = D0_had*pow(proton_E_GeV/1000.,0.5) # cm2/s
    t_esc = t_Sedov*pow(proton_E_GeV/E_p_max,-1./2.48)
    t_diffusion = max(0.,t_SN-t_esc)
    diffusion_radius = pow(4.*D0_at_E*t_diffusion*year_to_sec,0.5) # cm
    SNR_radius = SNR_radius_Sedov(t_SN) # cm
    diffusion_radius = pow(pow(diffusion_radius,2)+pow(SNR_radius,2),0.5)
    vol_diff = 4.*3.14/3.*pow(diffusion_radius,3) # cm3
    attenuation = exp(-pow(photon_to_SN_cm,2)/pow(diffusion_radius,2))

    energy_fraction = math.log(proton_E_up_GeV/proton_E_low_GeV)/math.log(E_p_max/E_p_min)
    ECR_density = CR_efficiency*E_SN*energy_fraction/vol_diff*attenuation # 1/cm3
    return ECR_density

def CR_energy_column_density(photon_x_deg,photon_y_deg,photon_E_low_GeV,photon_E_up_GeV,SN_data,D0_had):

    RA_SN = SN_data[0]
    Dec_SN = SN_data[1]
    d_SN = SN_data[2] # pc
    t_SN = SN_data[3] # year
    year_to_sec = 365.*24.*60.*60.
    deg_to_cm = 3.14/180.*d_SN*3.086e18

    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    photon_E_GeV = (photon_E_low_GeV+photon_E_up_GeV)/2.
    proton_E_GeV = (1./0.17)*2.*photon_E_GeV
    D0_at_E = D0_had*pow(proton_E_GeV/1000.,0.5) # cm2/s
    t_esc = t_Sedov*pow(proton_E_GeV/E_p_max,-1./2.48)
    t_diffusion = max(0.,t_SN-t_esc)
    diffusion_radius = pow(4.*D0_at_E*t_diffusion*year_to_sec,0.5) # cm
    SNR_radius = SNR_radius_Sedov(t_SN) # cm
    diffusion_radius = pow(pow(diffusion_radius,2)+pow(SNR_radius,2),0.5)
    photon_z_low = 0.
    photon_z_up = 2.*diffusion_radius

    ECR_column_density = 2.*quad(CR_energy_density,photon_z_low,photon_z_up,args=(photon_x_deg,photon_y_deg,photon_E_low_GeV,photon_E_up_GeV,SN_data,D0_had))[0]
    return ECR_column_density # /cm2

def CalculateHadronicEmission(photon_x_bin,photon_y_bin,photon_E_low_GeV,photon_E_up_GeV,SN_data,hist_clouds,D0_had):

    photon_x_deg = hist_clouds.GetXaxis().GetBinCenter(photon_x_bin)
    photon_y_deg = hist_clouds.GetYaxis().GetBinCenter(photon_y_bin)
    bin_x_size = hist_clouds.GetXaxis().GetBinLowEdge(2)-hist_clouds.GetXaxis().GetBinLowEdge(1)
    bin_y_size = hist_clouds.GetYaxis().GetBinLowEdge(2)-hist_clouds.GetYaxis().GetBinLowEdge(1)

    RA_SN = SN_data[0]
    Dec_SN = SN_data[1]
    d_SN = SN_data[2] # pc
    t_SN = SN_data[3] # year
    year_to_sec = 365.*24.*60.*60.
    deg_to_cm = 3.14/180.*d_SN*3.086e18
    erg_to_GeV = 624.151

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19
    skymap_bin_area = bin_x_size*bin_y_size*deg_to_cm*deg_to_cm # cm2

    cloud_column_density = hist_clouds.GetBinContent(photon_x_bin,photon_y_bin)
    cloud_molecule_number = cloud_column_density*skymap_bin_area
    ECR_max_density = CR_energy_density(0.,photon_x_deg,photon_y_deg,photon_E_low_GeV,photon_E_up_GeV,SN_data,D0_had)
    A_factor = cloud_molecule_number*ECR_max_density*pow(d_SN/1000.,-2)
    photon_E_GeV = (photon_E_low_GeV+photon_E_up_GeV)/2.
    hadronic_emission = A_factor*f_Gamma*1e-10*pow(photon_E_GeV/1000.,1-Gamma_CR) # /cm2/s

    return hadronic_emission


fig, ax = plt.subplots()
energy_bin = CommonPlotFunctions.energy_bin

energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])
print ('energy_bin_cut_low = %s'%(energy_bin_cut_low))
print ('energy_bin_cut_up = %s'%(energy_bin_cut_up))

InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_data_skymap_sum"
nbins = InputFile.Get(HistName).GetNbinsX()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()

other_stars, other_star_coord = CommonPlotFunctions.GetGammaSourceInfo()
other_star_labels = []
other_star_markers = []
other_star_names = []
source_ra = (MapEdge_left+MapEdge_right)/2.
source_dec = (MapEdge_lower+MapEdge_upper)/2.
star_range = MapEdge_right-source_ra
for star in range(0,len(other_stars)):
    if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
    other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
    other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
    other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
    other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
    other_star_names += [other_stars[star]]

hist_mask_skymap_sum = ROOT.TH2D("hist_mask_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_excess_skymap_sum = ROOT.TH2D("hist_excess_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_data_skymap_sum = ROOT.TH2D("hist_data_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_bkgd_skymap_sum = ROOT.TH2D("hist_bkgd_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_syst_skymap_sum = ROOT.TH2D("hist_syst_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum_v2","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap = []
hist_flux_syst_skymap = []
hist_flux_normsyst_skymap = []
hist_cali_skymap = []
hist_data_skymap = []
hist_bkgd_skymap = []
hist_expo_skymap = []
hist_syst_skymap = []
hist_normsyst_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_flux_skymap += [ROOT.TH2D("hist_flux_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_flux_syst_skymap += [ROOT.TH2D("hist_flux_syst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_flux_normsyst_skymap += [ROOT.TH2D("hist_flux_normsyst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_cali_skymap += [ROOT.TH2D("hist_cali_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_data_skymap += [ROOT.TH2D("hist_data_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_bkgd_skymap += [ROOT.TH2D("hist_bkgd_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_expo_skymap += [ROOT.TH2D("hist_expo_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_syst_skymap += [ROOT.TH2D("hist_syst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_normsyst_skymap += [ROOT.TH2D("hist_normsyst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_mask_skymap_sum"
hist_mask_skymap_sum.Add(InputFile.Get(HistName))
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    HistName = "hist_energy_flux_skymap_%s"%(ebin)
    hist_flux_skymap[ebin].Add(InputFile.Get(HistName))
    hist_flux_skymap_sum.Add(InputFile.Get(HistName))
    HistName = "hist_energy_flux_syst_skymap_%s"%(ebin)
    hist_flux_syst_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_energy_flux_normsyst_skymap_%s"%(ebin)
    hist_flux_normsyst_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_cali_skymap_%s"%(ebin)
    hist_cali_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_data_skymap_%s"%(ebin)
    hist_data_skymap[ebin].Add(InputFile.Get(HistName))
    hist_data_skymap_sum.Add(InputFile.Get(HistName))
    hist_excess_skymap_sum.Add(InputFile.Get(HistName))
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    hist_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
    hist_bkgd_skymap_sum.Add(InputFile.Get(HistName))
    hist_excess_skymap_sum.Add(InputFile.Get(HistName),-1.)
    HistName = "hist_expo_skymap_%s"%(ebin)
    hist_expo_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_syst_skymap_%s"%(ebin)
    hist_syst_skymap[ebin].Add(InputFile.Get(HistName))
    hist_syst_skymap_sum.Add(InputFile.Get(HistName))
    HistName = "hist_normsyst_skymap_%s"%(ebin)
    hist_normsyst_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_mc_density = ROOT.TH2D("Hist_mc_density","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

pc_to_cm = 3.086e+18
CO_intensity_to_H_column_density = 2.*1e20
# Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_50_80_0th_moment.txt' # CO intensity (K km s^{-1} deg)
Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
Hist_mc_intensity.Scale(FITS_correction)
Hist_mc_column.Add(Hist_mc_intensity)
Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
print ('Hist_mc_column.Integral() = %s'%(Hist_mc_column.Integral()))
mc_depth = 1000.*(3.824-3.116)*pc_to_cm
Hist_mc_density.Add(Hist_mc_column)
Hist_mc_density.Scale(1./mc_depth)

RA_G41p1 = 286.9 # G41.1-0.3
Dec_G41p1= 7.1
d_G41p1 = 8.*1000. #pc
t_G41p1 = 1350. # years 
radius_G41p1 = 0.03 # deg
SNR_G41p1_data = [RA_G41p1,Dec_G41p1,d_G41p1,t_G41p1]

RA_G40p5 = 286.786
Dec_G40p5 = 6.498
#d_G40p5 = 8.7*1000. #pc, PSR J1907+0631
#t_G40p5 = 11.*1000. # years
d_G40p5 = 3.4*1000. #pc, Yang 2006
t_G40p5 = 11.6*1000. # years
#d_G40p5 = 5.5*1000. #pc, Downes 1980
#t_G40p5 = 38.59*1000. # years
radius_G40p5 = 0.18 # deg
SNR_G40p5_data = [RA_G40p5,Dec_G40p5,d_G40p5,t_G40p5]

proper_velocity = 1300 #km/s
km_to_pc = 3.24078e-14
year_to_sec = 365.*24.*60.*60.
proper_velocity_pc_per_year = proper_velocity*km_to_pc*year_to_sec
pulsar_age_year = 19.5*1e3 # year
pulsar_distance = 3.2*1000. # pc
PSR_head_x = 286.98
PSR_head_y = 6.04
travel_angle = 71.94
tail_length = proper_velocity_pc_per_year*pulsar_age_year/pulsar_distance*180./3.14
PSR_tail_x = tail_length*np.cos(travel_angle*3.14/180.)+PSR_head_x
PSR_tail_y = tail_length*np.sin(travel_angle*3.14/180.)+PSR_head_y

RA_J1908 = PSR_tail_x
Dec_J1908 = PSR_tail_y
#d_J1908 = 3.2*1000. #pc
d_J1908 = 2.0*1000. #pc
t_J1908 = pulsar_age_year # years
SNR_J1908_data = [RA_J1908,Dec_J1908,d_J1908,t_J1908]


E_SN = 1. # 10^{51} erg
B_SNR = 100. # micro G
B_ISM = 6. # micro G
M_ej = 1. # ejecta mass in solar mass unit
n_0 = 3.0 # cm^{-3} ambient density
t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
E_p_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit
E_p_min = 1. # GeV
CR_efficiency = 0.5
print ('t_Sedov = %0.2e years'%(t_Sedov))
print ('vel_init = %0.2e cm/s'%(vel_init))
print ('E_p_max = %0.2e GeV'%(E_p_max))

D0 = 1.*8.2*1e26
#radius_G41p1 = SNR_radius_Sedov(t_G41p1)/(3.14/180.*d_G41p1*3.086e18)
#radius_G40p5 = SNR_radius_Sedov(t_G40p5)/(3.14/180.*d_G40p5*3.086e18)
t_G41p1 = SNR_age_Sedov(radius_G41p1,d_G41p1)
t_G40p5 = SNR_age_Sedov(radius_G40p5,d_G40p5)
radius_J1908 = SNR_radius_Sedov(t_J1908)/(3.14/180.*d_J1908*3.086e18)
proton_E = (1./0.17)*2.*600.
D_of_E = D0*pow(proton_E/1000.,0.5) # cm2/s
t_esc = t_Sedov*pow(proton_E/E_p_max,-1./2.48)
t_diffusion_J1908 = max(0.,t_J1908-t_esc)
diffusion_radius_J1908 = pow(4.*D_of_E*t_diffusion_J1908*year_to_sec,0.5)/(3.14/180.*d_J1908*3.086e18) # deg
CR_distribution_radius_J1908 = pow(pow(diffusion_radius_J1908,2)+pow(radius_J1908,2),0.5)
print ('SNR G41.1-0.3 age = %0.2f kyr'%(t_G41p1/1000.))
print ('SNR G40.5-0.5 age = %0.2f kyr'%(t_G40p5/1000.))
print ('SNR J1908+06 radius = %0.2f deg'%(radius_J1908))
print ('CR_distribution_radius_J1908 = %0.2f deg'%(CR_distribution_radius_J1908))

#CR_source_data = SNR_G41p1_data
CR_source_data = SNR_G40p5_data
#CR_source_data = SNR_J1908_data

hist_hadronic_skymap_sum = ROOT.TH2D("hist_hadronic_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_hadronic_skymap = ROOT.TH2D("hist_hadronic_skymap","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_hadronic_skymap.Reset()
    for binx in range(0,hist_hadronic_skymap_sum.GetNbinsX()):
        for biny in range(0,hist_hadronic_skymap_sum.GetNbinsY()):
            J1908_emission = 0.
            J1908_emission = CalculateHadronicEmission(binx+1,biny+1,energy_bin[ebin],energy_bin[ebin+1],CR_source_data,Hist_mc_column,D0)
            hist_hadronic_skymap.SetBinContent(binx+1,biny+1,J1908_emission)
    hist_hadronic_skymap.Scale(pow(energy_bin[ebin]/1e3,energy_index-1))
    hist_hadronic_skymap_sum.Add(hist_hadronic_skymap)


hist_zscore = CommonPlotFunctions.GetSignificanceMap(hist_data_skymap_sum,hist_bkgd_skymap_sum,hist_syst_skymap_sum,False)
hist_zscore_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore)
hist_zscore_reflect.SetContour(3)
hist_zscore_reflect.SetContourLevel(0,3)
hist_zscore_reflect.SetContourLevel(1,4)
hist_zscore_reflect.SetContourLevel(2,5)
hist_zscore_reflect.SetLineColor(0)

Hist_mc_density_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_density)
CommonPlotFunctions.MatplotlibMap2D(Hist_mc_density_reflect,hist_zscore_reflect,fig,'RA','Dec','density [$1/cm^{3}$]','SkymapMolecularDensity.png')

Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,hist_zscore_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapMolecularColumn.png')

Hist_mc_intensity_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_intensity)
CommonPlotFunctions.MatplotlibMap2D(Hist_mc_intensity_reflect,hist_zscore_reflect,fig,'RA','Dec','Intensity','SkymapMolecularIntensity.png')

hist_hadronic_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_hadronic_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_hadronic_skymap_sum_reflect,hist_zscore_reflect,fig,'RA','Dec','hadronic emission','SkymapHadronicEmission.png')

hist_flux_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_flux_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_flux_skymap_sum_reflect,hist_zscore_reflect,fig,'RA','Dec','flux','SkymapFlux.png')

hist_UL_skymap = ROOT.TH2D("hist_UL_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_UL_skymap.Add(hist_flux_skymap_sum)
hist_UL_skymap.Add(hist_hadronic_skymap_sum,-1)
hist_UL_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_UL_skymap)
CommonPlotFunctions.MatplotlibMap2D(hist_UL_skymap_reflect,hist_zscore_reflect,fig,'RA','Dec','flux','SkymapUpperLimit.png')

Hist_fermi5 = ROOT.TH2D("Hist_fermi5","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_fermi5 = CommonPlotFunctions.GetSkyViewMap("MWL_maps/skv1826930706371_j1908_fermi5.txt", Hist_fermi5, True)
# 3-300 GeV Band 5, Atwood et al. 2009
Hist_fermi5_reflect = CommonPlotFunctions.reflectXaxis(Hist_fermi5)
CommonPlotFunctions.MatplotlibMap2D(Hist_fermi5_reflect,hist_zscore_reflect,fig,'RA','Dec','cnts$/s/cm^{2}/sr$','SkymapFermi5.png')

