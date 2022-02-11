
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

def SNR_radius_Sedov(t_age_year):

    radius = 1.54*1e19*pow(E_SN,1./5.)*pow(n_0,-1./5.)*pow(t_age_year/1000.,2./5.) # cm
    return radius

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
    D0_at_E = D0_had*pow(proton_E_GeV/1000.,0.5) # cm2/s
    t_esc = t_Sedov*pow(proton_E_GeV/E_p_max,-1./2.48)
    t_diffusion = max(0.,t_SN-t_esc)
    diffusion_radius = pow(4.*D0_at_E*t_diffusion*year_to_sec,0.5) # cm
    SNR_radius = SNR_radius_Sedov(t_SN) # cm
    diffusion_radius = pow(pow(diffusion_radius,2)+pow(SNR_radius,2),0.5)
    vol_diff = 4.*3.14/3.*pow(diffusion_radius,3) # cm3
    attenuation = exp(-pow(photon_to_SN_cm,2)/pow(diffusion_radius,2))

    energy_fraction = math.log(photon_E_up_GeV/photon_E_low_GeV)/math.log(E_p_max/E_p_min)
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

    #Gamma_CR = 2.1
    #f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    Gamma_CR = 2.3
    f_Gamma = 0.19
    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    skymap_bin_area = bin_x_size*bin_y_size*deg_to_cm*deg_to_cm

    cloud_density = hist_clouds.GetBinContent(binx,biny)
    ECR_column_density = CR_energy_column_density(photon_x_deg,photon_y_deg,photon_E_low_GeV,photon_E_up_GeV,SN_data,D0_had)
    ECR_in_a_pixel = skymap_bin_area*ECR_column_density 
    A_factor = cloud_density*ECR_in_a_pixel*pow(d_SN/1000.,-2)
    photon_E_GeV = (photon_E_low_GeV+photon_E_up_GeV)/2.
    hadronic_emission = A_factor*f_Gamma*1e-10*pow(photon_E_GeV/1000.,1-Gamma_CR) # /cm2/s

    return hadronic_emission

def hadronic_model(x,par):

    binx = hist_lep_model_global.GetXaxis().FindBin(x[0])
    biny = hist_lep_model_global.GetYaxis().FindBin(x[1])
    model_lep = hist_lep_model_global.GetBinContent(binx,biny)
    model_had = hist_had_model_global.GetBinContent(binx,biny)

    norm_lep = par[0]
    norm_had = par[1]
    return norm_lep*model_lep + norm_had*model_had

def FitHadronicModel(hist_flux,hist_flux_syst,hist_pwn,hist_psr,hist_lep_model,hist_had_model):

    global hist_lep_model_global
    hist_lep_model_global.Reset()
    hist_lep_model_global.Add(hist_lep_model)
    global hist_had_model_global
    hist_had_model_global.Reset()
    hist_had_model_global.Add(hist_had_model)

    hist_data = ROOT.TH2D("hist_data","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    hist_data.Add(hist_flux)
    for bx in range(0,hist_data.GetNbinsX()):
        for by in range(0,hist_data.GetNbinsY()):
            syst = hist_flux_syst.GetBinContent(bx+1,by+1)
            hist_data.SetBinError(bx+1,by+1,syst)

    scaling_factor_lep = hist_psr.Integral()/hist_lep_model.Integral()
    scaling_factor_had = 1.

    npar = 2
    simple_model_2d = ROOT.TF2('simple_model_2d',hadronic_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    simple_model_2d.SetParameter(0,scaling_factor_lep)
    simple_model_2d.SetParLimits(0,0.,2.0*scaling_factor_lep)
    simple_model_2d.SetParameter(1,scaling_factor_had)
    simple_model_2d.SetParLimits(1,0.,10.0*scaling_factor_had)
    hist_data.Fit('simple_model_2d')
    scaling_factor_lep = simple_model_2d.GetParameter(0)
    scaling_factor_had = simple_model_2d.GetParameter(1)

    return scaling_factor_lep, scaling_factor_had

def GetHadronicModelChi2(hist_flux,hist_flux_syst,hist_pwn,hist_psr):

    hist_lep_model = []
    hist_had_model = []
    for ebin in range(0,len(energy_bin_big)-1):
        hist_lep_model += [ROOT.TH2D("hist_lep_model_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_had_model += [ROOT.TH2D("hist_had_model_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

    #D_lep_par = [0,1,2,3,4,5,6,7,8,9]
    #D_had_par = [0,1,2,3,4,5,6,7,8,9]
    D_lep_par = [0,1]
    D_had_par = [0,1]
    chi2_array = np.full((len(energy_bin_big)-1,len(D_lep_par),len(D_had_par)),0.)
    D_lep_axis = [ [0 for x in range(len(D_lep_par))] for y in range(0,2)]
    D_had_axis = [ [0 for x in range(len(D_had_par))] for y in range(0,2)]

    for idx_l in range(0,len(D_lep_par)):
        InputFile = ROOT.TFile("output_pulsar_models/pulsar_skymap_D%s_V%s.root"%(D_lep_par[idx_l],0))
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        D0_lep = InfoTree.D0/1e26
        D_lep_axis[0][idx_l] = '%0.2f'%(D0_lep)
        D_lep_axis[1][idx_l] = D_lep_par[idx_l]
        D_had_axis[0][idx_l] = '%0.2f'%(D0_lep)
        D_had_axis[1][idx_l] = D_lep_par[idx_l]
        InputFile.Close()

    for idx_l in range(0,len(D_lep_par)):
        for idx_h in range(0,len(D_had_par)):

            energy_index = 2
            InputFile = ROOT.TFile("output_pulsar_models/pulsar_skymap_D%s_V%s.root"%(D_lep_par[idx_l],0))
            for ebin_big in range(energy_bin_cut_low,energy_bin_cut_up):
                hist_lep_model[ebin_big].Reset()
                for ebin in range(0,len(energy_bin)-1):
                    if energy_bin[ebin]<energy_bin_big[ebin_big]: continue
                    if energy_bin[ebin]>=energy_bin_big[ebin_big+1]: continue
                    HistName = "hist_pulsar_skymap_631_1585"
                    if ebin==2:
                        HistName = "hist_pulsar_skymap_631_1585"
                    elif ebin==3:
                        HistName = "hist_pulsar_skymap_1585_3981"
                    elif ebin==4:
                        HistName = "hist_pulsar_skymap_3981_10000"
                    elif ebin==5:
                        HistName = "hist_pulsar_skymap_10000_25118"
                    hist_lep_model[ebin_big].Add(InputFile.Get(HistName),pow(energy_bin[ebin]/1e3,energy_index-1))
            InputFile.Close()

            for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
                for binx in range(0,hist_had_model[ebin].GetNbinsX()):
                    for biny in range(0,hist_had_model[ebin].GetNbinsY()):
                        D0_had = float(D_had_axis[0][idx_h])*1e26
                        G40p5_emission = 0.
                        G41p1_emission = 0.
                        G40p5_emission = CalculateHadronicEmission(binx+1,biny+1,energy_bin_big[ebin],energy_bin_big[ebin+1],SNR_G40p5_data,Hist_mc_density_G40p5,D0_had)
                        #G41p1_emission = CalculateHadronicEmission(binx+1,biny+1,energy_bin_big[ebin],energy_bin_big[ebin+1],SNR_G41p1_data,Hist_mc_density_G41p1,D0_had)
                        hist_had_model[ebin].SetBinContent(binx+1,biny+1,G40p5_emission+G41p1_emission)
                hist_had_model[ebin].Scale(pow(energy_bin[ebin]/1e3,energy_index-1))
            
            for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
                hist_lep_model[ebin] = CommonPlotFunctions.Smooth2DMap(hist_lep_model[ebin],CommonPlotFunctions.smooth_size_spectroscopy,False,True)
                hist_had_model[ebin] = CommonPlotFunctions.Smooth2DMap(hist_had_model[ebin],CommonPlotFunctions.smooth_size_spectroscopy,False,True)
                scaling_factor_lep, scaling_factor_had = FitHadronicModel(hist_flux[ebin],hist_flux_syst[ebin],hist_pwn[ebin],hist_psr[ebin],hist_lep_model[ebin],hist_had_model[ebin])
                hist_lep_model[ebin].Scale(scaling_factor_lep)
                hist_had_model[ebin].Scale(scaling_factor_had)
                total_chi2 = 0.
                total_bins = 0.
                for binx in range(0,hist_flux[ebin].GetNbinsX()):
                    for biny in range(0,hist_flux[ebin].GetNbinsY()):
                        data_flux = hist_flux[ebin].GetBinContent(binx+1,biny+1)
                        data_flux_err = hist_flux[ebin].GetBinError(binx+1,biny+1)
                        data_flux_syst = hist_flux_syst[ebin].GetBinContent(binx+1,biny+1)
                        pwn_flux = hist_pwn[ebin].GetBinContent(binx+1,biny+1)
                        model_flux_lep = hist_lep_model[ebin].GetBinContent(binx+1,biny+1)
                        model_flux_had = hist_had_model[ebin].GetBinContent(binx+1,biny+1)
                        error_sq = (data_flux_err*data_flux_err+data_flux_syst*data_flux_syst)
                        if error_sq>0.:
                            flux_zscore = data_flux/pow(error_sq,0.5)
                            if flux_zscore>2.:
                                chi2 = pow(data_flux-pwn_flux-model_flux_lep-model_flux_had,2)/error_sq
                                total_chi2 += chi2
                                total_bins += 1.
                chi2_array[ebin,idx_l,idx_h] = total_chi2/total_bins

    return chi2_array, D_lep_axis, D_had_axis


fig, ax = plt.subplots()
energy_bin = CommonPlotFunctions.energy_bin
energy_bin_big = CommonPlotFunctions.energy_bin_big

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

hist_lep_model_global = ROOT.TH2D("hist_lep_model_global","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_had_model_global = ROOT.TH2D("hist_had_model_global","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_syst_skymap_sum = ROOT.TH2D("hist_flux_syst_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_pwn_skymap_sum = ROOT.TH2D("hist_pwn_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_psr_skymap_sum = ROOT.TH2D("hist_psr_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap = []
hist_flux_syst_skymap = []
hist_pwn_skymap = []
hist_psr_skymap = []
for ebin in range(0,len(energy_bin_big)-1):
    hist_flux_skymap += [ROOT.TH2D("hist_flux_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_flux_syst_skymap += [ROOT.TH2D("hist_flux_syst_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_pwn_skymap += [ROOT.TH2D("hist_pwn_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_psr_skymap += [ROOT.TH2D("hist_psr_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])

InputFile = ROOT.TFile("output_fitting/J1908_fit_skymap.root")
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    HistName = "hist_flux_skymap_big_E%s"%(ebin)
    hist_flux_skymap_sum.Add(InputFile.Get(HistName))
    hist_flux_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_flux_syst_skymap_big_E%s"%(ebin)
    hist_flux_syst_skymap_sum.Add(InputFile.Get(HistName))
    hist_flux_syst_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_fit_PWN_skymap_big_E%s"%(ebin)
    hist_pwn_skymap_sum.Add(InputFile.Get(HistName))
    hist_pwn_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_fit_PSR_ellipse_skymap_big_E%s"%(ebin)
    hist_psr_skymap_sum.Add(InputFile.Get(HistName))
    hist_psr_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()


Hist_mc_column_G40p5 = ROOT.TH2D("Hist_mc_column_G40p5","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_mc_column_G41p1 = ROOT.TH2D("Hist_mc_column_G41p1","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_mc_density_G40p5 = ROOT.TH2D("Hist_mc_density_G40p5","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
Hist_mc_density_G41p1 = ROOT.TH2D("Hist_mc_density_G41p1","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

pc_to_cm = 3.086e+18
MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_50_60_0th_moment.txt' # CO intensity (K km s^{-1} deg)
Hist_mc_column_G40p5 = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_column_G40p5, True)
Hist_mc_column_G40p5.Scale(2.*1e20) # H2 column density in unit of 1/cm2
print ('Hist_mc_column_G40p5.Integral() = %s'%(Hist_mc_column_G40p5.Integral()))
mc_depth_G40p5 = 1000.*(3.824-3.116)*pc_to_cm
Hist_mc_density_G40p5.Add(Hist_mc_column_G40p5)
Hist_mc_density_G40p5.Scale(1./mc_depth_G40p5)
MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_50_60_0th_moment.txt' # CO intensity (K km s^{-1} deg)
Hist_mc_column_G41p1 = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_column_G41p1, True)
Hist_mc_column_G41p1.Scale(2.*1e20) # H2 column density in unit of 1/cm2
mc_depth_G41p1 = 1000.*(8.294-7.586)*pc_to_cm
Hist_mc_density_G41p1.Add(Hist_mc_column_G41p1)
Hist_mc_density_G41p1.Scale(1./mc_depth_G41p1)

hist_hadronic_skymap = []
for ebin in range(0,len(energy_bin_big)-1):
    hist_hadronic_skymap += [ROOT.TH2D("hist_hadronic_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

RA_G41p1 = 286.9 # G41.1-0.3
Dec_G41p1= 7.1
d_G41p1 = 8.*1000. #pc
t_G41p1 = 1350. # years 
SNR_G41p1_data = [RA_G41p1,Dec_G41p1,d_G41p1,t_G41p1]

RA_G40p5 = 286.786
Dec_G40p5 = 6.498
d_G40p5 = 3.4*1000. #pc, Yang et al. 2006
t_G40p5 = 11.*1000. # years
SNR_G40p5_data = [RA_G40p5,Dec_G40p5,d_G40p5,t_G40p5]

E_SN = 1. # 10^{51} erg
B_SNR = 100. # micro G
B_ISM = 6. # micro G
M_ej = 1. # ejecta mass in solar mass unit
n_0 = 3.0 # cm^{-3} ambient density
t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
E_p_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit
E_p_min = 1. # GeV
CR_efficiency = 0.3
print ('t_Sedov = %0.2e years'%(t_Sedov))
print ('vel_init = %0.2e cm/s'%(vel_init))
print ('E_p_max = %0.2e GeV'%(E_p_max))

D0 = 8.2*1e26
hist_hadronic_skymap_sum = ROOT.TH2D("hist_hadronic_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
energy_index = 2
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    for binx in range(0,hist_hadronic_skymap[ebin].GetNbinsX()):
        for biny in range(0,hist_hadronic_skymap[ebin].GetNbinsY()):
            G40p5_emission = 0.
            G41p1_emission = 0.
            G40p5_emission = CalculateHadronicEmission(binx+1,biny+1,energy_bin_big[ebin],energy_bin_big[ebin+1],SNR_G40p5_data,Hist_mc_density_G40p5,D0)
            #G41p1_emission = CalculateHadronicEmission(binx+1,biny+1,energy_bin_big[ebin],energy_bin_big[ebin+1],SNR_G41p1_data,Hist_mc_density_G41p1,D0)
            hist_hadronic_skymap[ebin].SetBinContent(binx+1,biny+1,G40p5_emission+G41p1_emission)
    hist_hadronic_skymap[ebin].Scale(pow(energy_bin_big[ebin]/1e3,energy_index-1))
    hist_hadronic_skymap_sum.Add(hist_hadronic_skymap[ebin])

chi2_array_pub, D_lep_axis, D_had_axis = GetHadronicModelChi2(hist_flux_skymap,hist_flux_syst_skymap,hist_pwn_skymap,hist_psr_skymap)
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    plt.imshow(chi2_array_pub[ebin,:,:])
    plt.xticks(D_had_axis[1], D_had_axis[0])
    plt.yticks(D_lep_axis[1], D_lep_axis[0])
    plt.xlabel("SNR diffusion coefficient 1e26 $\mathrm{cm}^{2}\mathrm{s}^{-1}$")
    plt.ylabel("pulsar diffusion coefficient 1e26 $\mathrm{cm}^{2}\mathrm{s}^{-1}$")
    plt.colorbar()
    plotname = 'HadronicModelChi2_E%s'%(ebin)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    plt.clf()


canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 650, 600)
pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
pad1.SetBottomMargin(0.10)
pad1.SetRightMargin(0.20)
pad1.SetLeftMargin(0.10)
pad1.SetTopMargin(0.10)
pad1.SetBorderMode(0)
pad1.Draw()
pad1.cd()

Hist_mc_density_G40p5_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_density_G40p5)
x1 = Hist_mc_density_G40p5_reflect.GetXaxis().GetXmin()
x2 = Hist_mc_density_G40p5_reflect.GetXaxis().GetXmax()
y1 = Hist_mc_density_G40p5_reflect.GetYaxis().GetXmin()
y2 = Hist_mc_density_G40p5_reflect.GetYaxis().GetXmax()
IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 505, "+")
raLowerAxis.SetLabelSize(Hist_mc_density_G40p5_reflect.GetXaxis().GetLabelSize())
raLowerAxis.Draw()
Hist_mc_density_G40p5_reflect.GetXaxis().SetLabelOffset(999)
Hist_mc_density_G40p5_reflect.GetXaxis().SetTickLength(0)
Hist_mc_density_G40p5_reflect.Draw("COL4Z")
raLowerAxis.Draw()
canvas.SaveAs('output_plots/SkymapMolecularDensity_G40p5.png')

hist_hadronic_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_hadronic_skymap_sum)
hist_hadronic_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_hadronic_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_hadronic_skymap_sum_reflect.Draw("COL4Z")
raLowerAxis.Draw()
canvas.SaveAs('output_plots/SkymapHadronicEmission.png')

