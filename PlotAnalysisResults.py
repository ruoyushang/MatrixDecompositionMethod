
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

import CommonPlotFunctions

fig, ax = plt.subplots()

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")
np.set_printoptions(precision=2)

skymap_zoomin_scale = CommonPlotFunctions.skymap_zoomin_scale
smooth_size_spectroscopy = CommonPlotFunctions.smooth_size_spectroscopy
Skymap_size_x = CommonPlotFunctions.Skymap_size_x
Skymap_size_y = CommonPlotFunctions.Skymap_size_y
Skymap_nbins_x = CommonPlotFunctions.Skymap_nbins_x
Skymap_nbins_y = CommonPlotFunctions.Skymap_nbins_y
energy_bin = CommonPlotFunctions.energy_bin
calibration_radius = CommonPlotFunctions.calibration_radius
elev_range = CommonPlotFunctions.elev_range
energy_index_scale = CommonPlotFunctions.energy_index_scale
Smoothing = CommonPlotFunctions.Smoothing
doGalacticCoord = CommonPlotFunctions.doGalacticCoord
N_bins_for_deconv = CommonPlotFunctions.N_bins_for_deconv
energy_fine_bin = CommonPlotFunctions.energy_fine_bin
gamma_hadron_dim_ratio_w = CommonPlotFunctions.gamma_hadron_dim_ratio_w
gamma_hadron_dim_ratio_l = CommonPlotFunctions.gamma_hadron_dim_ratio_l
gamma_hadron_low_end = CommonPlotFunctions.gamma_hadron_low_end
MSCW_blind_cut = CommonPlotFunctions.MSCW_blind_cut
MSCL_blind_cut = CommonPlotFunctions.MSCL_blind_cut
Skymap_nzones_x = CommonPlotFunctions.Skymap_nzones_x
Skymap_nzones_y = CommonPlotFunctions.Skymap_nzones_y

method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

folder_path = CommonPlotFunctions.folder_path

energy_bin_cut_low = 0
energy_bin_cut_up = 6

theta2_bins = [0,3]
#theta2_bins = [0,4]

#distance_sq_limit = 3.24
distance_sq_limit = 2.56
#distance_sq_limit = 100.

# all time
mjd_tag = []
mjd_tag += ['']
#mjd_tag += ['_MJD53613to54343'] #V4 2005-2007 Aug 31
#mjd_tag += ['_MJD54343to55074'] #V4 2007-2009 Aug 31
#mjd_tag += ['_MJD55074to55804'] #V5 2009-2011 Aug 31
#mjd_tag += ['_MJD55804to56535'] #V5+V6 2011-2013 Aug 31
#mjd_tag += ['_MJD56535to57265'] #V6 2013-2015 Aug 31
#mjd_tag += ['_MJD57265to57996'] #V6 2015-2017 Aug 31
#mjd_tag += ['_MJD57996to58726'] #V6 2017-2019 Aug 31
#mjd_tag += ['_MJD58726to59457'] #V6 2019-2021 Aug 31

ONOFF_tag = 'ON'
sample_list = []


def PrepreSample():
    global ONOFF_tag
    global sample_list
    for dt in range(0,len(data_type)):
        if sys.argv[1]=='%s_%s'%(observation_name,data_type[dt]):
            if data_type[dt]=='OFF':
                ONOFF_tag = 'OFF'
            else:
                ONOFF_tag = 'ON'
            ONOFF_tag += '_Model0'
            sample_list = []
            for de in range(0,len(data_epoch)):
                sample_list += ['%s_%s'%(data_epoch[de],data_type[dt])]
    

data_type = ['OFF','ON','Imposter1','Imposter2','Imposter3','Imposter4','Imposter5','Imposter6']

observation_name = 'MGRO_J1908'
data_epoch = ['MGRO_J1908_V5','MGRO_J1908_V6']
if 'Crab' in sys.argv[1]:
    observation_name = 'Crab'
    data_epoch = ['CrabV5','CrabV6']
if 'Crab_Elev70' in sys.argv[1]:
    observation_name = 'Crab_Elev70'
    data_epoch = ['Crab_Elev70_V5','Crab_Elev70_V6']
    elev_range = [70,90]
if 'Crab_Elev60' in sys.argv[1]:
    observation_name = 'Crab_Elev60'
    data_epoch = ['Crab_Elev60_V5','Crab_Elev60_V6']
    elev_range = [60,80]
if 'Crab_Elev50' in sys.argv[1]:
    observation_name = 'Crab_Elev50'
    data_epoch = ['Crab_Elev50_V5','Crab_Elev50_V6']
    elev_range = [50,70]
if 'Crab_Elev30' in sys.argv[1]:
    observation_name = 'Crab_Elev30'
    data_epoch = ['Crab_Elev30_V5','Crab_Elev30_V6']
    elev_range = [30,50]
if 'Crab_Offset_1p0' in sys.argv[1]:
    observation_name = 'Crab_Offset_1p0'
    data_epoch = ['Crab_Offset_1p0_V6']
if 'Crab_Offset_1p5' in sys.argv[1]:
    observation_name = 'Crab_Offset_1p5'
    data_epoch = ['Crab_Offset_1p5_V6']
if 'CrabRHV' in sys.argv[1]:
    observation_name = 'CrabRHV'
    data_epoch = ['CrabRHV_V6']
if 'NGC1275' in sys.argv[1]:
    observation_name = 'NGC1275'
    data_epoch = ['NGC1275V6']
if '1ES0229' in sys.argv[1]:
    observation_name = '1ES0229'
    data_epoch = ['1ES0229V5','1ES0229V6']
if 'LHAASO_J2032' in sys.argv[1]:
    observation_name = 'LHAASO_J2032'
    data_epoch = ['LHAASO_J2032_V6','LHAASO_J2032_V5']
if 'LHAASO_J2032_Fall2017' in sys.argv[1]:
    observation_name = 'LHAASO_J2032_Fall2017'
    data_epoch = ['LHAASO_J2032_Fall2017_V6','LHAASO_J2032_Fall2017_V5']
if 'LHAASO_J2032_Baseline' in sys.argv[1]:
    observation_name = 'LHAASO_J2032_Baseline'
    data_epoch = ['LHAASO_J2032_Baseline_V6','LHAASO_J2032_Baseline_V5']
if 'LHAASO_J1843' in sys.argv[1]:
    observation_name = 'LHAASO_J1843'
    data_epoch = ['LHAASO_J1843_V6']
if 'LHAASO_J1929' in sys.argv[1]:
    observation_name = 'LHAASO_J1929'
    data_epoch = ['LHAASO_J1929_V6']
if 'LHAASO_J2108' in sys.argv[1]:
    observation_name = 'LHAASO_J2108'
    data_epoch = ['LHAASO_J2108_V6']
if 'PSR_J1841_m0345' in sys.argv[1]:
    observation_name = 'PSR_J1841_m0345'
    data_epoch = ['PSR_J1841_m0345_V6','PSR_J1841_m0345_V5']
if 'PSR_J1856_p0245' in sys.argv[1]:
    observation_name = 'PSR_J1856_p0245'
    data_epoch = ['PSR_J1856_p0245_V6','PSR_J1856_p0245_V5']
if 'PSR_J1849_m0003' in sys.argv[1]:
    observation_name = 'PSR_J1849_m0003'
    data_epoch = ['PSR_J1849_m0003_V6','PSR_J1849_m0003_V5']
if 'PSR_J1938_p2213' in sys.argv[1]:
    observation_name = 'PSR_J1938_p2213'
    data_epoch = ['PSR_J1938_p2213_V6','PSR_J1938_p2213_V5']
if 'PSR_J2021_p3651' in sys.argv[1]:
    observation_name = 'PSR_J2021_p3651'
    data_epoch = ['PSR_J2021_p3651_V6','PSR_J2021_p3651_V5']
if 'PSR_J2021_p4026' in sys.argv[1]:
    observation_name = 'PSR_J2021_p4026'
    data_epoch = ['PSR_J2021_p4026_V6','PSR_J2021_p4026_V5']
if 'PSR_J2032_p4127_Baseline' in sys.argv[1]:
    observation_name = 'PSR_J2032_p4127_Baseline'
    data_epoch = ['PSR_J2032_p4127_Baseline_V6','PSR_J2032_p4127_Baseline_V5']
if 'PSR_J2032_p4127_Fall2017' in sys.argv[1]:
    observation_name = 'PSR_J2032_p4127_Fall2017'
    data_epoch = ['PSR_J2032_p4127_Fall2017_V6']
if 'PSR_J0633_p0632' in sys.argv[1]:
    observation_name = 'PSR_J0633_p0632'
    data_epoch = ['PSR_J0633_p0632_V6','PSR_J0633_p0632_V5']
if 'PSR_B0611_p22' in sys.argv[1]:
    observation_name = 'PSR_B0611_p22'
    data_epoch = ['PSR_B0611_p22_V6','PSR_B0611_p22_V5']
if 'IC443HotSpot' in sys.argv[1]:
    observation_name = 'IC443HotSpot'
    data_epoch = ['IC443HotSpotV5','IC443HotSpotV6']
if 'WComae' in sys.argv[1]:
    observation_name = 'WComae'
    data_epoch = ['WComaeV5','WComaeV6']
if 'Geminga' in sys.argv[1]:
    observation_name = 'Geminga'
    data_epoch = ['GemingaV5','GemingaV6']
if 'CasA' in sys.argv[1]:
    observation_name = 'CasA'
    data_epoch = ['CasAV6']
if 'M82' in sys.argv[1]:
    observation_name = 'M82'
    data_epoch = ['M82V5','M82V6']
if 'UrsaMajorII' in sys.argv[1]:
    observation_name = 'UrsaMajorII'
    data_epoch = ['UrsaMajorIIV6']
if 'UrsaMinor' in sys.argv[1]:
    observation_name = 'UrsaMinor'
    data_epoch = ['UrsaMinorV5','UrsaMinorV6']
if 'RGB_J0710_p591' in sys.argv[1]:
    observation_name = 'RGB_J0710_p591'
    data_epoch = ['RGB_J0710_p591_V5','RGB_J0710_p591_V6']
if 'GammaCygni' in sys.argv[1]:
    observation_name = 'GammaCygni'
    data_epoch = ['GammaCygniV5','GammaCygniV6']
if '1ES0502' in sys.argv[1]:
    observation_name = '1ES0502'
    data_epoch = ['1ES0502V5','1ES0502V6']
if '1ES0414' in sys.argv[1]:
    observation_name = '1ES0414'
    data_epoch = ['1ES0414V5']
if '3C273' in sys.argv[1]:
    observation_name = '3C273'
    data_epoch = ['3C273V5','3C273V6']
if 'BLLac' in sys.argv[1]:
    observation_name = 'BLLac'
    data_epoch = ['BLLacV5','BLLacV6']
if 'Draco' in sys.argv[1]:
    observation_name = 'Draco'
    data_epoch = ['DracoV5','DracoV6']
if '3C264' in sys.argv[1]:
    observation_name = '3C264'
    data_epoch = ['3C264V6']
if 'PG1553' in sys.argv[1]:
    observation_name = 'PG1553'
    data_epoch = ['PG1553V5','PG1553V6']
if 'OJ287' in sys.argv[1]:
    observation_name = 'OJ287'
    data_epoch = ['OJ287V6']
if 'H1426' in sys.argv[1]:
    observation_name = 'H1426'
    data_epoch = ['H1426V6']
if 'PKS1424' in sys.argv[1]:
    observation_name = 'PKS1424'
    data_epoch = ['PKS1424V5','PKS1424V6']
if 'Segue1' in sys.argv[1]:
    observation_name = 'Segue1'
    data_epoch = ['Segue1V5','Segue1V6']
if 'Boomerang' in sys.argv[1]:
    observation_name = 'Boomerang'
    data_epoch = ['BoomerangV5','BoomerangV6']
if 'Tycho' in sys.argv[1]:
    observation_name = 'Tycho'
    data_epoch = ['TychoV5','TychoV6']
if 'MGRO_J2019' in sys.argv[1]:
    observation_name = 'MGRO_J2019'
    data_epoch = ['MGRO_J2019_V5','MGRO_J2019_V6']
if 'HESS_J1825' in sys.argv[1]:
    observation_name = 'HESS_J1825'
    data_epoch = ['HESS_J1825_V6']
if '1ES1011' in sys.argv[1]:
    observation_name = '1ES1011'
    data_epoch = ['1ES1011V6']
if '1ES0647' in sys.argv[1]:
    observation_name = '1ES0647'
    data_epoch = ['1ES0647V6']
if 'SNR_G150p3Plus04p5' in sys.argv[1]:
    observation_name = 'SNR_G150p3Plus04p5'
    data_epoch = ['SNR_G150p3Plus04p5_V6']
if 'V_V725_Tau' in sys.argv[1]:
    observation_name = 'V_V725_Tau'
    data_epoch = ['V_V725_Tau_V6']

Gal_l_center = [30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230]
for glc in range(0,len(Gal_l_center)):
    if 'GalacticPlane_All_l%s'%(Gal_l_center[glc]) in sys.argv[1]:
        observation_name = 'GalacticPlane_All_l%s'%(Gal_l_center[glc])
        data_epoch = ['GalacticPlane_All_l%s_V6'%(Gal_l_center[glc]),'GalacticPlane_All_l%s_V5'%(Gal_l_center[glc])]

PrepreSample()

root_file_tags = []

for elev in range(0,len(elev_range)-1):
    elev_tag = '_TelElev%sto%s'%(elev_range[elev],elev_range[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        #theta2_tag = '_Y%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

print ('Get %s'%(root_file_tags[0]))

#for elev in range(0,len(elev_range)-1):
#    # all time
#    root_file_tags += [method_tag+'_TelElev%sto%s_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # 1ES 1215 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56372to56464_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # WComae and 1ES 1218 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD54400to54700_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # Geminga
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56000to56500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58000to58700_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58700to59000_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # MGRO J2031
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to57997_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57997to59000_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # 3 slices of 1ES 0229
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to56500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56501to57500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57501to59500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]

selection_tag = root_file_tags[0]

PercentCrab = ''


MSCW_chi2_upper = -0.5
MSCL_chi2_upper = -0.5
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))+MSCL_blind_cut
MSCW_plot_lower = -gamma_hadron_low_end*gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))-MSCW_blind_cut
MSCL_plot_lower = -gamma_hadron_low_end*gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))-MSCL_blind_cut
#MSCW_plot_lower = -0.5
#MSCL_plot_lower = -0.5
ErecS_lower_cut = 0
ErecS_upper_cut = 0

n_good_matches = 0
exposure_hours = 0.
exposure_hours_ref = 0.
NSB_mean_data = 0.
Zenith_mean_data = 0.
NSB_RMS_data = 0.
Zenith_RMS_data = 0.
source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
n_control_samples = 1
if sys.argv[1]=='SgrA':
    n_control_samples = 1
MJD_Start = 2147483647
MJD_End = 0
Data_runlist_NSB = ROOT.std.vector("double")(10)
Data_runlist_L3Rate = ROOT.std.vector("double")(10)
roi_name = ROOT.std.vector("string")(10)
roi_ra = ROOT.std.vector("double")(10)
roi_dec = ROOT.std.vector("double")(10)
roi_radius_inner = ROOT.std.vector("double")(10)
roi_radius_outer = ROOT.std.vector("double")(10)
max_chi2_diff2_position = ROOT.std.vector("double")(10)
max_chi2_diff2_position_this_energy = 0.

Syst_MDM = 0.02
Syst_Init = 0.02
Syst_Redu = 0.02
Syst_Corr = 0.02
Syst_Clos = 0.02

integration_radii = [0.1,0.5,1.0,1.5,2.0]

energy_syst = []
energy_syst += [0.121]
energy_syst += [0.017]
energy_syst += [0.019]
energy_syst += [0.036]
energy_syst += [0.080]
energy_syst += [0.095]
energy_syst += [0.095]


elev_bins = [25,35,45,55,65,75,85]
MJD_bins = [53613,55074,56535,57996,59457]

other_stars = []
other_star_coord = []

bright_star_ra = []
bright_star_dec = []
bright_star_brightness = []
faint_star_ra = []
faint_star_dec = []
faint_star_brightness = []

#NRGBs = 5
#NCont = 104
#stops = [0.00,0.40,0.50,0.60,1.00]
#red =   [0.00,1.00,1.00,1.00,1.00]
#green = [0.00,1.00,1.00,1.00,0.00]
#blue =  [1.00,1.00,1.00,1.00,0.00]
#ROOT.TColor.CreateGradientColorTable(NRGBs,array('d',stops),array('d',red),array('d',green),array('d',blue),NCont)
#ROOT.gStyle.SetNumberContours(NCont)

def set_palette(name, ncontours=999):
    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    ROOT.TColor.CreateGradientColorTable(len(stops),array('d',stops),array('d',red),array('d',green),array('d',blue),ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def EstimateHadronicAngularDiffusionLength(photon_energy,mag_field,d_SNR,t_SNR,esc_param):

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    M_ej = 1. # ejecta mass in solar mass unit
    B_ISM = mag_field # micro G
    deg_to_cm = 3.14/180.*d_SNR*1000.*3.086e18

    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    E_vis = (1./0.17)*2.*photon_energy  # visible CR energy threshold

    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
    E_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit

    t_esc = t_Sedov*pow(E_vis/E_max,-1./esc_param)
    t_diff = t_SNR-t_esc
    if t_diff<=0.: return 0.1
    diffusion_ism = 1e28*pow(E_vis/10.,0.5)*pow(B_ISM/3.,-0.5) # cm2/s
    diffusion_radius = pow(t_diff*(4.*diffusion_ism*365.*24.*60.*60.),0.5)/deg_to_cm

    return diffusion_radius

def Smooth2DMap(Hist_Old,smooth_size,addLinearly,normalized):
    return CommonPlotFunctions.Smooth2DMap(Hist_Old,smooth_size,addLinearly,normalized)

def FindExtension(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range):
    return CommonPlotFunctions.FindExtension(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range)

def EstimateLeptonicDiffusionCoefficient(photon_energy,mag_field,d_PSR,t_PSR,diff_coeff,diff_coeff_index):
    # photon_energy in GeV
    # mag_field in muG
    # d_PSR in kpc
    # t_PSR in yr
    # diff_coeff in D = D0 (E/4TeV)^{diff_coeff_index}
    E_cmb = 6.6*1e-4 # eV
    m_e = 0.511*1e6 # eV
    E_e = m_e*pow(photon_energy*1e9/E_cmb,0.5) # eV
    gamma_factor = E_e/m_e
    sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
    speed_light = 3.*1e8 # m/s
    U_cmb = 2.6*1e5 # eV/m3
    U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3
    t_cooling = m_e/speed_light/(4./3.*sigma_thomson*gamma_factor*(U_cmb+U_B)) # sec
    t_cooling = min(t_cooling,t_PSR*365.*24.*60.*60.)
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    D_ism = diff_coeff*pow(E_e/(1e12),diff_coeff_index)

    return D_ism

def EstimateDiffusionCoefficient(photon_energy,angular_extent,mag_field,d_PSR,t_PSR,diff_coeff_index):
    # photon_energy in GeV
    # angular_extent in deg
    # mag_field in muG
    # d_PSR in kpc
    E_cmb = 6.6*1e-4 # eV
    m_e = 0.511*1e6 # eV
    E_e = m_e*pow(photon_energy*1e9/E_cmb,0.5) # eV
    gamma_factor = E_e/m_e
    sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
    speed_light = 3.*1e8 # m/s
    U_cmb = 2.6*1e5 # eV/m3
    U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3
    t_cooling = m_e/speed_light/(4./3.*sigma_thomson*gamma_factor*(U_cmb+U_B)) # sec
    t_cooling = min(t_cooling,t_PSR*365.*24.*60.*60.)
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    source_radius = angular_extent*deg_to_cm
    D_ism = 1./(4.*t_cooling)*pow(source_radius,2) # cm2/s
    print ('E_ph = %0.2f GeV'%(photon_energy))
    print ('E_e = %0.2f GeV'%(E_e/1e9))
    print ('t_cooling = %0.2f year'%(t_cooling/(365.*24.*60.*60.)))
    print ('D_ism = %0.3f x10^{28} cm2/s at E_e = %0.2f TeV'%(D_ism/1e28,E_e/1e12))
    D_ism_1TeV = D_ism*pow(1e12/E_e,diff_coeff_index)
    print ('D_ism = %0.3f x10^{28} cm2/s at E_e = 1 TeV'%(D_ism_1TeV/1e28))
    D_ism_100TeV = D_ism*pow(1e14/E_e,diff_coeff_index)
    print ('D_ism = %0.3f x10^{28} cm2/s at E_e = 100 TeV'%(D_ism_100TeV/1e28))

def ConvertGalacticToRaDec(l, b):
    return CommonPlotFunctions.ConvertGalacticToRaDec(l, b)

def ConvertRaDecToGalactic(ra, dec):
    return CommonPlotFunctions.ConvertRaDecToGalactic(ra, dec)

def GetBrightStarInfo(file_path):

    global bright_star_ra
    global bright_star_dec
    global bright_star_brightness
    global faint_star_ra
    global faint_star_dec
    global faint_star_brightness

    print ('GetBrightStarInfo Read file: %s'%(file_path))
    if not os.path.isfile(file_path): 
        print ('No such a file!')
        return
    bright_star_ra = []
    bright_star_dec = []
    bright_star_brightness = []
    faint_star_ra = []
    faint_star_dec = []
    faint_star_brightness = []
    InputFile = ROOT.TFile(file_path)
    StarTree = InputFile.Get("StarTree")
    print ('StarTree.GetEntries() = %s'%(StarTree.GetEntries()))
    for entry in range(0,StarTree.GetEntries()):
        StarTree.GetEntry(entry)
        distance_to_center = pow(pow(source_ra-StarTree.star_ra,2)+pow(source_dec-StarTree.star_dec,2),0.5)
        bright_star_ra += [StarTree.star_ra]
        bright_star_dec += [StarTree.star_dec]
        bright_star_brightness += [StarTree.star_brightness]
    FaintStarTree = InputFile.Get("FaintStarTree")
    for entry in range(0,FaintStarTree.GetEntries()):
        FaintStarTree.GetEntry(entry)
        distance_to_center = pow(pow(source_ra-FaintStarTree.faint_star_ra,2)+pow(source_dec-FaintStarTree.faint_star_dec,2),0.5)
        faint_star_ra += [FaintStarTree.faint_star_ra]
        faint_star_dec += [FaintStarTree.faint_star_dec]
        faint_star_brightness += [FaintStarTree.faint_star_brightness]
    InputFile.Close()

def GetGammaSourceInfo():

    global other_stars
    global other_star_coord

    other_stars, other_star_coord = CommonPlotFunctions.GetGammaSourceInfo(None) 

def ResetStackedShowerHistograms():

    Hist_EffArea_Sum.Reset()

    Hist2D_OnData_Point_Sum.Reset()
    Hist2D_OnData_Ring_Sum.Reset()
    Hist2D_OnData_Sum.Reset()
    Hist2D_OnBkgd_Sum.Reset()
    Hist2D_OnGamma_Sum.Reset()
    Hist2D_OnDark_Sum.Reset()
    Hist_OnData_MSCL_Sum.Reset()
    Hist_OnBkgd_MSCL_Sum.Reset()
    Hist_OnDark_MSCL_Sum.Reset()
    Hist_OnData_MSCW_Sum.Reset()
    Hist_OnBkgd_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_wGamma_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_woGamma_MSCW_Sum.Reset()
    Hist_OnGamma_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_wGamma_MSCL_Sum.Reset()
    Hist_OnBkgd_Unblind_woGamma_MSCL_Sum.Reset()
    Hist_OnGamma_MSCL_Sum.Reset()
    Hist_OnDark_MSCW_Sum.Reset()

    Hist2D_Rank0_Data_Sum.Reset()
    Hist2D_Rank1_Data_Sum.Reset()
    Hist2D_Rank2_Data_Sum.Reset()
    Hist2D_Rank3_Data_Sum.Reset()
    Hist2D_Rank4_Data_Sum.Reset()
    Hist2D_Rank0_Dark_Sum.Reset()
    Hist2D_Rank1_Dark_Sum.Reset()
    Hist2D_Rank2_Dark_Sum.Reset()
    Hist2D_Rank3_Dark_Sum.Reset()
    Hist2D_Rank4_Dark_Sum.Reset()

    Hist2D_U_Proj_Sum.Reset()
    Hist2D_V_Proj_Sum.Reset()
    Hist2D_Coeff_Data_Sum.Reset()
    Hist2D_Coeff_Bkgd_Sum.Reset()

    Hist_SystErr_MSCL.Reset()
    Hist_SystErr_MSCW.Reset()
    Hist_SystErr_Energy.Reset()
    Hist_SystErr_Theta2_Sum.Reset()

    Hist_OnData_SR_XYoff_Sum.Reset()
    Hist_OnDark_SR_XYoff_Sum.Reset()
    Hist_OnData_CR_XYoff_Sum.Reset()
    Hist_OnData_R2off_Sum.Reset()
    Hist_OnData_ISR_R2off_Sum.Reset()
    Hist_OnData_Yoff_Sum.Reset()
    Hist_OnBkgd_R2off_Sum.Reset()
    Hist_OnBkgd_Yoff_Sum.Reset()
    Hist_OnBkgd_Yoff_Raw_Sum.Reset()
    Hist_OnData_Theta2_Sum.Reset()
    Hist_OnBkgd_Theta2_Sum.Reset()
    Hist_OnRFoV_Theta2_Sum.Reset()
    Hist_OnDark_Theta2_Sum.Reset()
    Hist_OnData_Energy_Sum.Reset()
    Hist_OnBkgd_Energy_Sum.Reset()
    Hist_OnRFoV_Energy_Sum.Reset()
    Hist_OnDark_Energy_Sum.Reset()
    Hist_OnData_Zenith_Sum.Reset()
    Hist_OnBkgd_Zenith_Sum.Reset()
    Hist_OnData_Height_Sum.Reset()
    Hist_OnBkgd_Height_Sum.Reset()
    Hist_OnData_Depth_Sum.Reset()
    Hist_OnBkgd_Depth_Sum.Reset()
    Hist_OnData_Rcore_Sum.Reset()
    Hist_OnBkgd_Rcore_Sum.Reset()
    Hist_OnData_Energy_CamCenter_Sum.Reset()
    Hist_OnBkgd_Energy_CamCenter_Sum.Reset()

    Hist_OnData_Skymap_Sum.Reset()
    Hist_OnDark_Skymap_Sum.Reset()
    Hist_OnBkgd_Skymap_Sum.Reset()
    Hist_OnData_Skymap_Galactic_Sum.Reset()
    Hist_OnBkgd_Skymap_Galactic_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_Sum.Reset()
    Hist_OnBkgd_Skymap_Galactic_Syst_MDM.Reset()

    Hist_OnData_Skymap_ProjX_Sum.Reset()
    Hist_OnDark_Skymap_ProjX_Sum.Reset()
    Hist_OnBkgd_Skymap_ProjX_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_ProjX_Sum.Reset()
    Hist_OnData_Skymap_ProjY_Sum.Reset()
    Hist_OnDark_Skymap_ProjY_Sum.Reset()
    Hist_OnBkgd_Skymap_ProjY_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_ProjY_Sum.Reset()

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Reset() 
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Reset() 
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Reset()
        Hist_OnData_RoI_X_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_X_Sum[nth_roi].Reset()
        Hist_OnData_RoI_Y_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_Y_Sum[nth_roi].Reset()
        Hist_OnData_RoI_MJD_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Reset()

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Reset()
        Hist_OffData_MSCL_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCL_Sum[nth_sample].Reset()
        Hist_OffData_MSCW_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCW_Sum[nth_sample].Reset()
        Hist_OffData_Energy_Sum[nth_sample].Reset()
        Hist_OffBkgd_Energy_Sum[nth_sample].Reset()
        Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].Reset()
        Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].Reset()

def GetSourceInfo(file_list):

    global N_bins_for_deconv
    global MSCW_blind_cut
    global MSCL_blind_cut
    global MSCW_chi2_upper
    global MSCL_chi2_upper
    #global MSCW_plot_lower
    #global MSCL_plot_lower
    global n_good_matches
    global exposure_hours
    global exposure_hours_ref
    global NSB_mean_data
    global Zenith_mean_data
    global NSB_RMS_data
    global Zenith_RMS_data
    #global Skymap_size
    #global Skymap_nbins
    global source_ra
    global source_dec
    global source_l
    global source_b
    global n_control_samples
    global MJD_Start
    global MJD_End
    global Data_runlist_NSB
    global Data_runlist_L3Rate
    global roi_name
    global roi_ra
    global roi_dec
    global roi_radius_inner
    global roi_radius_outer
    global max_chi2_diff2_position

    n_good_matches = 0
    exposure_hours = 0.
    exposure_hours_ref = 0.
    NSB_mean_data = 0.
    Zenith_mean_data = 0.
    NSB_RMS_data = 0.
    Zenith_RMS_data = 0.
    for path in range(0,len(file_list)):
        print ('Read file: %s'%(file_list[path]))
        if not os.path.isfile(file_list[path]):
            print ('Cannot get source info from file: %s'%(file_list[path]))
            continue
        print ('Getting source info from file: %s'%(file_list[path]))
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.SetBranchAddress('Data_runlist_NSB',ROOT.AddressOf(Data_runlist_NSB))
        InfoTree.SetBranchAddress('Data_runlist_L3Rate',ROOT.AddressOf(Data_runlist_L3Rate))
        InfoTree.SetBranchAddress('roi_name',ROOT.AddressOf(roi_name))
        InfoTree.SetBranchAddress('roi_ra',ROOT.AddressOf(roi_ra))
        InfoTree.SetBranchAddress('roi_dec',ROOT.AddressOf(roi_dec))
        InfoTree.SetBranchAddress('roi_radius_inner',ROOT.AddressOf(roi_radius_inner))
        InfoTree.SetBranchAddress('roi_radius_outer',ROOT.AddressOf(roi_radius_outer))
        InfoTree.GetEntry(0)
        Hist_NSB.Reset()
        Hist_L3Rate.Reset()
        for entry in range(0,len(Data_runlist_NSB)):
            Hist_NSB.Fill(Data_runlist_NSB[entry])
        for entry in range(0,len(Data_runlist_L3Rate)):
            Hist_L3Rate.Fill(Data_runlist_L3Rate[entry])
        #N_bins_for_deconv = InfoTree.N_bins_for_deconv
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        MSCW_chi2_upper = InfoTree.MSCW_chi2_upper
        MSCL_chi2_upper = InfoTree.MSCL_chi2_upper
        #MSCW_plot_lower = InfoTree.MSCW_plot_lower
        #MSCL_plot_lower = InfoTree.MSCL_plot_lower
        #Skymap_size = InfoTree.Skymap_size
        #Skymap_nbins = InfoTree.Skymap_nbins
        source_ra = InfoTree.mean_tele_point_ra
        source_dec = InfoTree.mean_tele_point_dec
        if doGalacticCoord:
            source_ra = InfoTree.mean_tele_point_l
            source_dec = InfoTree.mean_tele_point_b
        source_l = InfoTree.mean_tele_point_l
        source_b = InfoTree.mean_tele_point_b
        exposure_hours += InfoTree.exposure_hours_usable
        MJD_Start = min(InfoTree.MJD_Start,MJD_Start)
        MJD_End = max(InfoTree.MJD_End,MJD_End)
        print ('len(roi_ra) = %s'%(len(roi_ra)))

        NewInfoTree = InputFile.Get("NewInfoTree")
        NewInfoTree.SetBranchAddress('max_chi2_diff2_position',ROOT.AddressOf(max_chi2_diff2_position))
        NewInfoTree.GetEntry(0)
        NSB_mean_data = NewInfoTree.NSB_mean_data
        Zenith_mean_data = NewInfoTree.Zenith_mean_data
        NSB_RMS_data = NewInfoTree.NSB_RMS_data
        Zenith_RMS_data = NewInfoTree.Zenith_RMS_data
        #if 'Theta20' in file_list[path] and 'Up' in file_list[path]:
        #if 'Y0' in file_list[path]:

        HistName = "Hist_EffArea"
        Hist_EffArea.Reset()
        Hist_EffArea.Add(InputFile.Get(HistName))
        Hist_EffArea.Scale(InfoTree.exposure_hours_usable*3600.)
        Hist_EffArea_Sum.Add(Hist_EffArea)

        InputFile.Close()

def MergeHistogram(hist_sum,hist_input):

    bin_size_x_1 = hist_sum.GetXaxis().GetBinCenter(1)-hist_sum.GetXaxis().GetBinCenter(0)
    bin_size_y_1 = hist_sum.GetYaxis().GetBinCenter(1)-hist_sum.GetYaxis().GetBinCenter(0)
    bin_size_x_2 = hist_input.GetXaxis().GetBinCenter(1)-hist_input.GetXaxis().GetBinCenter(0)
    bin_size_y_2 = hist_input.GetYaxis().GetBinCenter(1)-hist_input.GetYaxis().GetBinCenter(0)

    bin_size_x_ratio = bin_size_x_2/bin_size_x_1
    bin_size_y_ratio = bin_size_y_2/bin_size_y_1

    hist_output = hist_sum.Clone()
    for binx in range(0,hist_sum.GetNbinsX()):
        for biny in range(0,hist_sum.GetNbinsY()):
            old_content = hist_sum.GetBinContent(binx+1,biny+1)
            old_error = pow(hist_sum.GetBinError(binx+1,biny+1),2)
            x = hist_sum.GetXaxis().GetBinCenter(binx+1)
            y = hist_sum.GetYaxis().GetBinCenter(biny+1)
            new_binx = hist_input.GetXaxis().FindBin(x)-1
            new_biny = hist_input.GetYaxis().FindBin(y)-1
            new_content = hist_input.GetBinContent(new_binx+1,new_biny+1)/(bin_size_x_ratio*bin_size_y_ratio) 
            new_error = pow(hist_input.GetBinError(new_binx+1,new_biny+1),2)/(bin_size_x_ratio*bin_size_y_ratio) 
            hist_output.SetBinContent(binx+1,biny+1,old_content+new_content)
            hist_output.SetBinError(binx+1,biny+1,pow(old_error+new_error,0.5))

    return hist_output

def GetShowerHistogramsFromFile(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    global MSCW_chi2_upper
    global MSCL_chi2_upper
    global max_chi2_diff2_position_this_energy

    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    MSCW_chi2_upper = InfoTree.MSCW_chi2_upper
    MSCL_chi2_upper = InfoTree.MSCL_chi2_upper
    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(-1.*MSCL_blind_cut+0.01)-1
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut+0.01)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(-1.*MSCW_blind_cut+0.01)-1
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut+0.01)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnData.Reset()
    Hist2D_OnData.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnData,InputFile.Get(HistName))
    #Hist2D_OnData.Add(HistTemp)
    HistName = "Hist_OnData_Point_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnData_Point.Reset()
    Hist2D_OnData_Point.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_Ring_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnData_Ring.Reset()
    Hist2D_OnData_Ring.Add(InputFile.Get(HistName))

    HistName = "Hist_OnDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnDark.Reset()
    Hist2D_OnDark.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnDark,InputFile.Get(HistName))
    #Hist2D_OnDark.Add(HistTemp)

    HistName = "Hist_OnBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnBkgd.Reset()
    Hist2D_OnBkgd.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnBkgd,InputFile.Get(HistName))
    #Hist2D_OnBkgd.Add(HistTemp)

    HistName = "Hist_Gamma_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Getting histogram %s'%(HistName))
    Hist2D_OnGamma.Reset()
    Hist2D_OnGamma.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnGamma,InputFile.Get(HistName))
    #Hist2D_OnGamma.Add(HistTemp)

    Hist2D_Rank0_Data.Reset()
    Hist2D_Rank1_Data.Reset()
    Hist2D_Rank2_Data.Reset()
    Hist2D_Rank3_Data.Reset()
    Hist2D_Rank4_Data.Reset()
    Hist2D_Rank0_Dark.Reset()
    Hist2D_Rank1_Dark.Reset()
    Hist2D_Rank2_Dark.Reset()
    Hist2D_Rank3_Dark.Reset()
    Hist2D_Rank4_Dark.Reset()

    Hist2D_U_Proj.Reset()
    Hist2D_V_Proj.Reset()
    Hist2D_Coeff_Data.Reset()
    Hist2D_Coeff_Bkgd.Reset()

    HistName = "Hist_Rank0_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank0_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank1_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank1_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank2_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank2_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank3_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank3_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank4_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank4_Data.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank0_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank0_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank1_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank1_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank2_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank2_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank3_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank3_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank4_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank4_Dark.Add(InputFile.Get(HistName))

    HistName = "Hist_U_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_U_Proj.Add(InputFile.Get(HistName))

    HistName = "Hist_V_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_V_Proj.Add(InputFile.Get(HistName))

    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Coeff_Data.Add(InputFile.Get(HistName))

    HistName = "Hist_Coeff_Bkgd_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Coeff_Bkgd.Add(InputFile.Get(HistName))

    #if Hist1D_Data_Rank0_LeftVector.Integral()==0:
    Hist_VVV_Eigenvalues.Reset()
    Hist_Bkgd_Chi2.Reset()
    Hist_Bkgd_Chi2_Diff.Reset()
    Hist_Bkgd_Chi2_Diff2.Reset()
    Hist_Bkgd_Converge_Blind.Reset()
    Hist_Bkgd_Converge_Unblind.Reset()
    Hist1D_Data_Rank0_LeftVector.Reset()
    Hist1D_Data_Rank1_LeftVector.Reset()
    Hist1D_Data_Rank2_LeftVector.Reset()
    Hist1D_Data_Rank3_LeftVector.Reset()
    Hist1D_Data_Rank0_RightVector.Reset()
    Hist1D_Data_Rank1_RightVector.Reset()
    Hist1D_Data_Rank2_RightVector.Reset()
    Hist1D_Data_Rank3_RightVector.Reset()
    Hist1D_Bkgd_Rank0_LeftVector.Reset()
    Hist1D_Bkgd_Rank1_LeftVector.Reset()
    Hist1D_Bkgd_Rank2_LeftVector.Reset()
    Hist1D_Bkgd_Rank3_LeftVector.Reset()
    Hist1D_Bkgd_Rank0_RightVector.Reset()
    Hist1D_Bkgd_Rank1_RightVector.Reset()
    Hist1D_Bkgd_Rank2_RightVector.Reset()
    Hist1D_Bkgd_Rank3_RightVector.Reset()
    Hist1D_Dark_Rank0_LeftVector.Reset()
    Hist1D_Dark_Rank1_LeftVector.Reset()
    Hist1D_Dark_Rank2_LeftVector.Reset()
    Hist1D_Dark_Rank3_LeftVector.Reset()
    Hist1D_Dark_Rank0_RightVector.Reset()
    Hist1D_Dark_Rank1_RightVector.Reset()
    Hist1D_Dark_Rank2_RightVector.Reset()
    Hist1D_Dark_Rank3_RightVector.Reset()
    HistName = "Hist_VVV_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_VVV_Eigenvalues.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_Diff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2_Diff.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_Diff2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2_Diff2.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Converge_Blind_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Converge_Blind.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Converge_Unblind_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Converge_Unblind.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank3_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank3_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank3_RightVector.Add(InputFile.Get(HistName))

    Hist_VVV_Eigenvalues.SetMinimum(1e-2)
    #MakeOneHistPlot(Hist_VVV_Eigenvalues,'entry','eigenvalue','VVV_Eigenvalue_%s_%s'%(sample_list[0],ErecS_lower_cut_int),True)
    #MakeOneHistPlot(Hist_Bkgd_Chi2,'log10 #alpha','#chi^{2} in CR','Bkgd_Chi2_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)
    #MakeOneHistPlot(Hist_Bkgd_Chi2_Diff,'log10 #alpha','#chi^{2} 1st derivative in CR','Bkgd_Chi2_Diff_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)
    #MakeOneHistPlot(Hist_Bkgd_Chi2_Diff2,'log10 #alpha','#chi^{2} 2nd derivative in CR','Bkgd_Chi2_Diff2_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)

    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank0_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank0_LeftVector]
    #legends += ['inital vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank0_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{1}}','Rank0_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank1_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank1_LeftVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank1_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{2}}','Rank1_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank2_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank2_LeftVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank2_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{3}}','Rank2_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank0_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank0_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank0_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{1}}','Rank0_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank1_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank1_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank1_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{2}}','Rank1_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank2_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank2_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank2_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{3}}','Rank2_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)

    #if Hist2D_OnData.Integral()<1600. or Hist2D_OnDark.Integral()<1600.:
    #    Hist2D_OnData.Reset()
    #    Hist2D_OnDark.Reset()
    #    Hist2D_OnBkgd.Reset()
    #    Hist2D_OnGamma.Reset()
    #    Hist2D_Rank0_Data.Reset()
    #    Hist2D_Rank1_Data.Reset()
    #    Hist2D_Rank2_Data.Reset()
    #    Hist2D_Rank3_Data.Reset()
    #    Hist2D_Rank4_Data.Reset()
    #    Hist2D_Rank0_Dark.Reset()
    #    Hist2D_Rank1_Dark.Reset()
    #    Hist2D_Rank2_Dark.Reset()
    #    Hist2D_Rank3_Dark.Reset()
    #    Hist2D_Rank4_Dark.Reset()
    #    Hist2D_U_Proj.Reset()
    #    Hist2D_V_Proj.Reset()
    #    Hist2D_Coeff_Data.Reset()
    #    Hist2D_Coeff_Bkgd.Reset()

    Hist_OnData_MSCL.Reset()
    Hist_OnData_MSCL.Add(Hist2D_OnData.ProjectionX("Hist1D_OnData_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnData_MSCW.Reset()
    Hist_OnData_MSCW.Add(Hist2D_OnData.ProjectionY("Hist1D_OnData_MSCW",bin_lower_x,bin_upper_x))

    Hist_OnDark_MSCL.Reset()
    Hist_OnDark_MSCL.Add(Hist2D_OnDark.ProjectionX("Hist1D_OnDark_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnDark_MSCW.Reset()
    Hist_OnDark_MSCW.Add(Hist2D_OnDark.ProjectionY("Hist1D_OnDark_MSCW",bin_lower_x,bin_upper_x))

    Hist_OnBkgd_MSCL.Reset()
    Hist_OnBkgd_MSCL.Add(Hist2D_OnBkgd.ProjectionX("Hist1D_OnBkgd_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnBkgd_MSCW.Reset()
    Hist_OnBkgd_MSCW.Add(Hist2D_OnBkgd.ProjectionY("Hist1D_OnBkgd_MSCW",bin_lower_x,bin_upper_x))
    Hist_OnGamma_MSCW.Reset()
    Hist_OnGamma_MSCW.Add(Hist2D_OnGamma.ProjectionY("Hist1D_OnGamma_MSCW",bin_lower_x,bin_upper_x))
    Hist_OnGamma_MSCL.Reset()
    Hist_OnGamma_MSCL.Add(Hist2D_OnGamma.ProjectionX("Hist1D_OnGamma_MSCL",bin_lower_y,bin_upper_y))

    for nth_sample in range(0,n_control_samples):

        HistName = "Hist_OffData2_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print ('Getting histogram %s'%(HistName))
        Hist2D_OffData[nth_sample].Reset()
        Hist2D_OffData[nth_sample].Add(InputFile.Get(HistName))
        #HistTemp = MergeHistogram(Hist2D_OffData[nth_sample],InputFile.Get(HistName))
        #Hist2D_OffData[nth_sample].Add(HistTemp)

        #if Hist2D_OffData[nth_sample].Integral()<1600.:
        #    Hist2D_OffData[nth_sample].Reset()
        if math.isnan(Hist2D_OffData[nth_sample].Integral()):
            Hist2D_OffData[nth_sample].Reset()

        Hist_OffData_MSCL[nth_sample].Reset()
        Hist_OffData_MSCL[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionX("Hist1D_OffData_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffData_MSCW[nth_sample].Reset()
        Hist_OffData_MSCW[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionY("Hist1D_OffData_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

    InputFile.Close()


def StackShowerHistograms():

    Hist2D_OnData_Point_Sum.Add(Hist2D_OnData_Point)
    Hist2D_OnData_Ring_Sum.Add(Hist2D_OnData_Ring)
    Hist2D_OnData_Sum.Add(Hist2D_OnData)
    Hist2D_OnDark_Sum.Add(Hist2D_OnDark)
    Hist2D_OnBkgd_Sum.Add(Hist2D_OnBkgd)
    Hist2D_OnGamma_Sum.Add(Hist2D_OnGamma)

    if Hist2D_U_Proj_Sum.Integral()==0.:
        Hist2D_U_Proj_Sum.Add(Hist2D_U_Proj)
        Hist2D_V_Proj_Sum.Add(Hist2D_V_Proj)
        Hist2D_Coeff_Data_Sum.Add(Hist2D_Coeff_Data)
        Hist2D_Coeff_Bkgd_Sum.Add(Hist2D_Coeff_Bkgd)

    Hist2D_Rank0_Data_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank1_Data)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank1_Data)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank2_Data)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank1_Data)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank2_Data)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank3_Data)

    #Hist2D_Rank0_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    #Hist2D_Rank0_Data_Sum.Add(Hist2D_OnData)
    #Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    #Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    #Hist2D_Rank1_Data_Sum.Add(Hist2D_OnData)
    #Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    #Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    #Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    #Hist2D_Rank2_Data_Sum.Add(Hist2D_OnData)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank3_Data,-1.)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_OnData)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank3_Data,-1.)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank4_Data,-1.)
    #Hist2D_Rank4_Data_Sum.Add(Hist2D_OnData)

    Hist2D_Rank0_Dark_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank1_Data)
    Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank2_Data)
    Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank3_Data)
    #Hist2D_Rank0_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank0_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank3_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank3_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank4_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_OnDark)

    Hist_OnData_MSCL_Sum.Add(Hist_OnData_MSCL)
    Hist_OnData_MSCW_Sum.Add(Hist_OnData_MSCW)
    Hist_OnDark_MSCL_Sum.Add(Hist_OnDark_MSCL)
    Hist_OnDark_MSCW_Sum.Add(Hist_OnDark_MSCW)
    Hist_OnBkgd_MSCL_Sum.Add(Hist_OnBkgd_MSCL)
    Hist_OnBkgd_MSCW_Sum.Add(Hist_OnBkgd_MSCW)
    Hist_OnBkgd_Unblind_wGamma_MSCW_Sum.Add(Hist_OnBkgd_Unblind_wGamma_MSCW)
    Hist_OnBkgd_Unblind_woGamma_MSCW_Sum.Add(Hist_OnBkgd_Unblind_woGamma_MSCW)
    Hist_OnGamma_MSCW_Sum.Add(Hist_OnGamma_MSCW)
    Hist_OnBkgd_Unblind_wGamma_MSCL_Sum.Add(Hist_OnBkgd_Unblind_wGamma_MSCL)
    Hist_OnBkgd_Unblind_woGamma_MSCL_Sum.Add(Hist_OnBkgd_Unblind_woGamma_MSCL)
    Hist_OnGamma_MSCL_Sum.Add(Hist_OnGamma_MSCL)

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Add(Hist2D_OffData[nth_sample])

        Hist_OffData_MSCL_Sum[nth_sample].Add(Hist_OffData_MSCL[nth_sample])
        Hist_OffData_MSCW_Sum[nth_sample].Add(Hist_OffData_MSCW[nth_sample])
        Hist_OffBkgd_MSCL_Sum[nth_sample].Add(Hist_OffBkgd_MSCL[nth_sample])
        Hist_OffBkgd_MSCW_Sum[nth_sample].Add(Hist_OffBkgd_MSCW[nth_sample])

def CalculateSignificance(s,b,err):
    if (b*b+(s+b)*err*err)==0.: return 0.
    if (s+b)*(b+err*err)==0.: return 0.
    if ((s+b)*(b+err*err)/(b*b+(s+b)*err*err))<=0.: return 0.
    first_term = (s+b)*math.log((s+b)*(b+err*err)/(b*b+(s+b)*err*err))
    if err>0. and b>0:
        second_term = b*b/(err*err)*math.log(1.+err*err*s/(b*(b+err*err)))
    else: 
        second_term = 0.
    result = 0.
    if first_term>second_term: result = pow(2*(first_term-second_term),0.5)
    else: result = pow(2*(-first_term+second_term),0.5)
    if s>0: return result
    else: return -1.*result

def Hist2DIntegralAndError(Hist):

    return CommonPlotFunctions.Hist2DIntegralAndError(Hist)

def IntegralAndError(Hist,bin1,bin2):
    
    integral = 0
    error = 0
    for b in range(1,Hist.GetNbinsX()+1):
        if b<bin1: continue
        if b>bin2: continue
        integral += Hist.GetBinContent(b)
        error += pow(Hist.GetBinError(b),2)
    error = pow(error,0.5)
    error = max(error,1.)
    if math.isnan(integral) or math.isnan(error):
        integral = 0
        error = 1.
    return integral, error

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeMultiplePlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetLeftMargin(0.2)
    pad3.SetBorderMode(1)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.8)
    pad2.SetBottomMargin(0.2)
    pad2.SetLeftMargin(0.2)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad2.SetGrid()
    pad2.Draw()
    pad3.Draw()

    pad2.cd()
    if logy: pad2.SetLogy()

    min_heigh = 0
    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,len(Hists)):
        mean += [0]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].GetXaxis().SetTitleOffset(0.8)
            Hists[h].GetXaxis().SetTitleSize(0.06)
            Hists[h].GetXaxis().SetLabelSize(0.06)
            Hists[h].GetYaxis().SetLabelSize(0.06)
            Hists[h].GetYaxis().SetTitleOffset(1.2)
            Hists[h].GetYaxis().SetTitleSize(0.06)
            Hists[h].GetXaxis().SetTitle(title_x)
            Hists[h].GetYaxis().SetTitle(title_y)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
            if min_heigh > Hists[h].GetMinimum(): 
                min_heigh = Hists[h].GetMinimum()

    gap = 0.1*(max_heigh-min_heigh)
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(2)
        Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy,normalized):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.1)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

    min_heigh = 0
    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,len(Hists)):
        mean += [0]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].GetXaxis().SetTitleOffset(0.8)
            Hists[h].GetXaxis().SetTitleSize(0.06)
            Hists[h].GetXaxis().SetLabelSize(0.06)
            Hists[h].GetYaxis().SetLabelSize(0.06)
            Hists[h].GetYaxis().SetTitleOffset(1.2)
            Hists[h].GetYaxis().SetTitleSize(0.06)
            Hists[h].GetXaxis().SetTitle(title_x)
            Hists[h].GetYaxis().SetTitle(title_y)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
            if min_heigh > Hists[h].GetMinimum(): 
                min_heigh = Hists[h].GetMinimum()

    gap = 0.1*(max_heigh-min_heigh)
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        if normalized: Hists[0].DrawNormalized("E")
        else: Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        if normalized: Hists[max_hist].DrawNormalized("E")
        else: Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(2)
        if normalized: Hists[h].DrawNormalized("E same")
        else: Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

#def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
#    
#    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
#    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
#    pad3.SetBottomMargin(0.0)
#    pad3.SetTopMargin(0.03)
#    pad3.SetBorderMode(1)
#    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
#    pad1.SetBottomMargin(0.1)
#    pad1.SetTopMargin(0.0)
#    pad1.SetBorderMode(0)
#    pad1.SetGrid()
#    pad1.Draw()
#    pad3.Draw()
#
#    pad1.cd()
#    if logy: pad1.SetLogy()
#
#    min_heigh = 0
#    max_heigh = 0
#    max_hist = 0
#    mean = []
#    rms = []
#    amp = []
#    for h in range(0,len(Hists)):
#        mean += [0]
#        rms += [0]
#        amp += [0]
#        if Hists[h]!=0:
#            Hists[h].GetXaxis().SetTitleOffset(0.8)
#            Hists[h].GetXaxis().SetTitleSize(0.06)
#            Hists[h].GetXaxis().SetLabelSize(0.06)
#            Hists[h].GetYaxis().SetLabelSize(0.06)
#            Hists[h].GetYaxis().SetTitleOffset(1.2)
#            Hists[h].GetYaxis().SetTitleSize(0.06)
#            Hists[h].GetXaxis().SetTitle(title_x)
#            Hists[h].GetYaxis().SetTitle(title_y)
#            if max_heigh < Hists[h].GetMaximum(): 
#                max_heigh = Hists[h].GetMaximum()
#                max_hist = h
#            if min_heigh > Hists[h].GetMinimum(): 
#                min_heigh = Hists[h].GetMinimum()
#
#    gap = 0.1*(max_heigh-min_heigh)
#    if not y_max==0. and not y_min==0.:
#        Hists[0].SetMaximum(y_max)
#        Hists[0].SetMinimum(y_min)
#        Hists[0].DrawNormalized("E")
#    else:
#        if not logy:
#            Hists[max_hist].SetMaximum(max_heigh+gap)
#            Hists[max_hist].SetMinimum(min_heigh-gap)
#        Hists[max_hist].DrawNormalized("E")
#
#    for h in range(0,len(Hists)):
#        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
#        #if Hists[h]!=0:
#        Hists[h].SetLineColor(colors[h])
#        Hists[h].SetLineWidth(2)
#        Hists[h].DrawNormalized("E same")
#
#    pad3.cd()
#    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
#    legend.SetTextFont(42)
#    legend.SetBorderSize(0)
#    legend.SetTextSize(0.1)
#    legend.SetFillColor(0)
#    legend.SetFillStyle(0)
#    legend.SetLineColor(0)
#    legend.Clear()
#    for h in range(0,len(Hists)):
#        if Hists[h]!=0:
#            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
#    legend.Draw("SAME")
#
#    if logx: 
#        pad1.SetLogx()
#
#    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def GetCrabFlux():

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    func_crab.SetParameters(37.5,-2.467,-0.16)

    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1
    count_crab_total = 0.
    for binx in range(0,Hist_EffArea_Sum.GetNbinsX()):
        if binx+1<bin_lower: continue
        if binx+1>bin_upper: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        flux_crab_this_energy = func_crab.Eval(Hist_EffArea_Sum.GetBinCenter(binx+1))
        count_crab_total += flux_crab_this_energy*Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE
    return count_crab_total

def GetBackgroundFluxInterval(lower_cut,upper_cut):

    bin_lower = Hist_OnData_Energy.FindBin(lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(upper_cut)-1
    flux_bkgd_total = 0.
    for binx in range(0,Hist_EffArea_Sum.GetNbinsX()):
        if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
        if binx+1<bin_lower: continue
        if binx+1>bin_upper: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        bkg_this_energy, bkg_err_this_energy = IntegralAndError(Hist_OnBkgd_Energy_Sum,binx+1,binx+2)
        flux_bkgd_this_energy = bkg_this_energy
        flux_bkgd_total += flux_bkgd_this_energy
    return flux_bkgd_total

def MakeCrabUnitSpectrumPlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    n_roi = len(Hist_data)
    pads = []

    nth_pad = 0
    pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,1-0.7,1,1)]
    pads[nth_pad].SetBorderMode(1)
    pads[nth_pad].SetTopMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.0)
    pads[nth_pad].SetTopMargin(0.2)
    pads[nth_pad].Draw("same")
    nth_pad = 1
    pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,0,1,1-0.7)]
    pads[nth_pad].SetBorderMode(1)
    pads[nth_pad].SetTopMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.2)
    pads[nth_pad].Draw("same")

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    Func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    Func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_data_energy = []
    Hist_bkgd_energy = []
    Hist_ratio_energy = []
    for roi in range(0,n_roi):

        Hist_data_energy += [ROOT.TH1D("Hist_data_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            Hist_data_energy[roi].SetBinContent(binx+1,Hist_data[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_data_energy[roi].GetNbinsX()):
            Hist_data_energy[roi].SetBinError(binx+1,pow(Hist_data_energy[roi].GetBinContent(binx+1),0.5))

        Hist_bkgd_energy += [ROOT.TH1D("Hist_bkgd_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            Hist_bkgd_energy[roi].SetBinContent(binx+1,Hist_bkgd[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_energy[roi].GetNbinsX()):
            Hist_bkgd_energy[roi].SetBinError(binx+1,Hist_bkgd_energy[roi].GetBinContent(binx+1)*Syst_MDM)

        Hist_ratio_energy += [ROOT.TH1D("Hist_ratio_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        Hist_ratio_energy[roi].Add(Hist_data_energy[roi])
        Hist_ratio_energy[roi].Add(Hist_bkgd_energy[roi],-1.)
        Hist_ratio_energy[roi].Divide(Hist_bkgd_energy[roi])
        for binx in range(0,Hist_data_energy[roi].GetNbinsX()):
            lower_cut = Hist_data[roi].GetBinLowEdge(binx+1)
            upper_cut = Hist_data[roi].GetBinLowEdge(binx+2)
            deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
            bkgd_flux = 0.
            if Hist_EffArea_Sum.GetBinContent(binx+1)!=0.:
                bkgd_flux = Hist_bkgd_energy[roi].GetBinContent(binx+1)/(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
            crab_flux = func_crab.Eval(Hist_EffArea_Sum.GetBinCenter(binx+1))
            Hist_ratio_energy[roi].SetBinContent(binx+1,Hist_ratio_energy[roi].GetBinContent(binx+1)*bkgd_flux/crab_flux)
            Hist_ratio_energy[roi].SetBinError(binx+1,Hist_ratio_energy[roi].GetBinError(binx+1)*bkgd_flux/crab_flux)

        pads[0].cd()
        Hist_ratio_energy[roi].SetLineColor(1)
        Hist_ratio_energy[roi].SetLineWidth(3)
        Hist_ratio_energy[roi].SetMinimum(-0.1)
        Hist_ratio_energy[roi].Draw("E")

        pads[1].cd()
        Hist_data_energy[roi].GetXaxis().SetTitle("MJD")
        Hist_data_energy[roi].SetLineColor(1)
        Hist_data_energy[roi].SetLineWidth(3)
        Hist_data_energy[roi].Draw("E")
        fill_color = 38
        stack = ROOT.THStack("stack", "")
        set_histStyle( Hist_bkgd_energy[roi] , fill_color)
        stack.Add( Hist_bkgd_energy[roi] )
        stack.Draw("hist same")
        Hist_data_energy[roi].Draw("E same")

        pads[0].SetLogx()
        pads[1].SetLogx()
        pads[1].SetLogy()

        c_both.SaveAs('output_plots/%s_RoI%s_%s.png'%(name,roi,selection_tag))

def GetFermiFluxJ1908(energy_index):
    energies = [pow(10.,10.6),pow(10.,10.86),pow(10.,11.12),pow(10.,11.37)]
    fluxes = [pow(10.,-11.31),pow(10.,-11.29),pow(10.,-11.12),pow(10.,-10.88)]
    flux_errs = [pow(10.,-11.31),pow(10.,-11.29),pow(10.,-11.12),pow(10.,-10.88)]
    flux_errs_up = [pow(10.,-11.17),pow(10.,-11.12),pow(10.,-10.96),pow(10.,-10.76)]
    flux_errs_low = [pow(10.,-11.53),pow(10.,-11.57),pow(10.,-11.38),pow(10.,-11.05)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e9
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiUpperLimitFluxGeminga(energy_index):
    energies = [pow(10.,2.09),pow(10.,2.35),pow(10.,2.61),pow(10.,2.87)]
    fluxes = [pow(10.,-7.35),pow(10.,-7.23),pow(10.,-7.34),pow(10.,-7.18)]
    fluxes_err = [pow(10.,-7.35),pow(10.,-7.23),pow(10.,-7.34),pow(10.,-7.18)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        fluxes_err[entry] = fluxes[entry]*0.3

    return energies, fluxes, fluxes_err

def GetFermiHAWCFluxGeminga(energy_index):
    energies = [pow(10.,3.90),pow(10.,4.60)]
    fluxes = [pow(10.,-8.13),pow(10.,-8.36)]
    flux_errs = [0.,0.]
    flux_errs_up = [pow(10.,-8.06),pow(10.,-8.30)]
    flux_errs_low = [pow(10.,-8.24),pow(10.,-8.47)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxGeminga(energy_index):
    energies = [pow(10.,1.03),pow(10.,1.30),pow(10.,1.56),pow(10.,1.82)]
    fluxes = [pow(10.,-7.27),pow(10.,-7.29),pow(10.,-7.41),pow(10.,-7.35)]
    flux_errs = [pow(10.,-7.27),pow(10.,-7.29),pow(10.,-7.41),pow(10.,-7.35)]
    flux_errs_up = [pow(10.,-7.14),pow(10.,-7.16),pow(10.,-7.29),pow(10.,-7.23)]
    flux_errs_low = [pow(10.,-7.46),pow(10.,-7.49),pow(10.,-7.58),pow(10.,-7.50)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHAWCFluxJ1908(energy_index):
    energies = [pow(10.,12.00),pow(10.,12.20),pow(10.,12.47),pow(10.,12.74),pow(10.,13.01),pow(10.,13.29),pow(10.,13.55),pow(10.,13.78),pow(10.,14.00),pow(10.,14.21)]
    fluxes = [pow(10.,-10.60),pow(10.,-10.62),pow(10.,-10.64),pow(10.,-10.74),pow(10.,-10.86),pow(10.,-11.03),pow(10.,-11.18),pow(10.,-11.36),pow(10.,-11.37),pow(10.,-11.55)]
    flux_errs = [pow(10.,-10.60),pow(10.,-10.62),pow(10.,-10.64),pow(10.,-10.74),pow(10.,-10.86),pow(10.,-11.03),pow(10.,-11.18),pow(10.,-11.36),pow(10.,-11.37),pow(10.,-11.55)]
    flux_errs_up = [pow(10.,-10.57),pow(10.,-10.58),pow(10.,-10.61),pow(10.,-10.70),pow(10.,-10.82),pow(10.,-10.99),pow(10.,-11.14),pow(10.,-11.29),pow(10.,-11.28),pow(10.,-11.38)]
    flux_errs_low = [pow(10.,-10.64),pow(10.,-10.65),pow(10.,-10.68),pow(10.,-10.77),pow(10.,-10.89),pow(10.,-11.07),pow(10.,-11.22),pow(10.,-11.43),pow(10.,-11.48),pow(10.,-11.81)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e9
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs


def gaussian_func(x,a,b,c):
    return a*np.exp(-0.5*pow((x-b)/c,2))

def power_law_func(x,a,b):
    return a*pow(10,-12)*pow(x*1./1000.,b)

def flux_crab_func(x):
    # Crab https://arxiv.org/pdf/1508.06442.pdf
    return 37.5*pow(10,-12)*pow(x*1./1000.,-2.467-0.16*log(x/1000.))
def flux_veritas_boomerang_func(x):
    # Boomerang  TeV^{-1}cm^{-2}s^{-1}
    return 2.46*pow(10,-15)*pow(x*1./20000.,-2.29)
def flux_veritas_gamma_cygni_func(x):
    # VER J2019+407  TeV^{-1}cm^{-2}s^{-1}
    return 1.5*pow(10,-12)*pow(x*1./1000.,-2.37)
    #return 5.01*pow(10,-13)*pow(x*1./1500.,-2.79)
def flux_hawc_gamma_cygni_func(x):
    # 3HWC J2020+403  TeV^{-1}cm^{-2}s^{-1}
    return 11.4*pow(10,-15)*pow(x*1./7000.,-3.11)
def flux_2hwc_J1928_point_func(x):
    # 2HWC J1928+177 TeV^{-1}cm^{-2}s^{-1}
    return 10.0*pow(10,-15)*pow(x*1./7000.,-2.56)
def flux_2hwc_J1844_point_func(x):
    # 2HWC J1844-032 TeV^{-1}cm^{-2}s^{-1}
    return 46.8*pow(10,-15)*pow(x*1./7000.,-2.64)
def flux_2hwc_J1844_extend_func(x):
    # 2HWC J1844-032 TeV^{-1}cm^{-2}s^{-1}
    return 92.8*pow(10,-15)*pow(x*1./7000.,-2.51)
def flux_hess_J1843_func(x):
    # HESS J1843-033 TeV^{-1}cm^{-2}s^{-1}
    return 9.14*pow(10,-13)*pow(x*1./1870.,-2.15)
def flux_veritas_j1908_func(x):
    # MGRO J1908  TeV^{-1}cm^{-2}s^{-1}
    return 4.23*pow(10,-12)*pow(x*1./1000.,-2.2)
def flux_hess_j1908_func(x):
    # MGRO J1908  TeV^{-1}cm^{-2}s^{-1}
    return 4.14*pow(10,-12)*pow(x*1./1000.,-2.1)
def flux_hawc_j1908_func(x):
    # MGRO J1908  TeV^{-1}cm^{-2}s^{-1}
    return 0.95*pow(10,-13)*pow(x*1./10000.,-2.46-0.11*log(x/10000.))
def flux_j1857_func(x):
    # HESS J1857+026
    return 5.37*pow(10,-12)*pow(x*1./1000.,-2.16)
def flux_ic443_func(x):
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    #return 0.838*pow(10,-12)*pow(x*1./1000.,-2.99)
    # IC 443 https://arxiv.org/pdf/1512.01911.pdf
    return 9.92*pow(10,-12)*pow(x*1./550.,-2.8)
    #return 9.92*pow(10,-13)*pow(x*1./1000.,-2.8)
def flux_ic443_hawc_func(x):
    # 3HWC J0617+224 https://arxiv.org/pdf/2007.08582.pdf
    return 4.5*pow(10,-15)*pow(x*1./7000.,-3.05)
def flux_1es1218_func(x):
    # 1ES 1218 https://arxiv.org/pdf/0810.0301.pdf
    return 7.5*pow(10,-12)*pow(x*1./1000.,-3.08)
def flux_geminga_func(x):
    return 48.7*pow(10,-15)*pow(x*1./7000.,-2.23)  #https://arxiv.org/pdf/1702.02992.pdf table 3
    #return 5.9*pow(10,-15)*pow(x*1./7000.,-2.52)  #https://arxiv.org/pdf/2007.08582.pdf table 2

def MakeSpectrumInNonCrabUnit(ax_local,hist_data,hist_bkgd,hist_syst,radii,legends,title,doCalibrate,doBkgFlux,E_index):
    
    global calibration

    Hist_Flux = []
    for nth_roi in range(0,len(hist_data)):
        Hist_Flux += [ROOT.TH1D("Hist_Flux_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        if doBkgFlux:
            Hist_Flux[nth_roi].Add(hist_bkgd[nth_roi])
            Hist_Flux[nth_roi].Add(hist_syst[nth_roi])
            for binx in range(energy_bin_cut_low,energy_bin_cut_up):
                deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
                scale = 0.
                if Hist_EffArea_Sum.GetBinContent(binx+1)>0.:
                    scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
                Hist_Flux[nth_roi].SetBinContent(binx+1,Hist_Flux[nth_roi].GetBinContent(binx+1)*scale)
                Hist_Flux[nth_roi].SetBinError(binx+1,Hist_Flux[nth_roi].GetBinError(binx+1)*scale)
        else:
            if doCalibrate:
                Hist_Flux[nth_roi].Add(hist_data[nth_roi])
                Hist_Flux[nth_roi].Add(hist_bkgd[nth_roi],-1.)
                Hist_Flux[nth_roi].Add(hist_syst[nth_roi],-1.)
                Hist_Flux[nth_roi].Divide(hist_bkgd[nth_roi])
                for binx in range(energy_bin_cut_low,energy_bin_cut_up):
                    print ('=================================')
                    print ('energy %s'%(Hist_Flux[nth_roi].GetBinLowEdge(binx+1)))
                    print ('rate %s'%(Hist_Flux[nth_roi].GetBinContent(binx+1)))
                    print ('bkgd %s'%(hist_bkgd[nth_roi].GetBinContent(binx+1)))
                    print ('error %s'%(hist_syst[nth_roi].GetBinError(binx+1)))
            else:
                Hist_Flux[nth_roi].Add(hist_data[nth_roi])
                Hist_Flux[nth_roi].Add(hist_bkgd[nth_roi],-1.)
                Hist_Flux[nth_roi].Add(hist_syst[nth_roi],-1.)
                for binx in range(energy_bin_cut_low,energy_bin_cut_up):
                    deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
                    scale = 0.
                    if Hist_EffArea_Sum.GetBinContent(binx+1)>0.:
                        scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
                    Hist_Flux[nth_roi].SetBinContent(binx+1,Hist_Flux[nth_roi].GetBinContent(binx+1)*scale)
                    Hist_Flux[nth_roi].SetBinError(binx+1,Hist_Flux[nth_roi].GetBinError(binx+1)*scale)

    roi_xdata = []
    roi_ydata = []
    roi_error = []
    roi_zeros = []
    for nth_roi in range(0,len(hist_data)):
        xdata = []
        ydata = []
        error = []
        zeros = []
        for binx in range(energy_bin_cut_low,energy_bin_cut_up):
            xdata += [Hist_Flux[nth_roi].GetBinLowEdge(binx+1)]
            zeros += [0.]
            if doCalibrate:
                ydata += [Hist_Flux[nth_roi].GetBinContent(binx+1)*calibration[binx]*pow(radii[nth_roi],2)/(0.4*0.4)]
                error += [Hist_Flux[nth_roi].GetBinError(binx+1)*calibration[binx]*pow(radii[nth_roi],2)/(0.4*0.4)]
            else:
                ydata += [Hist_Flux[nth_roi].GetBinContent(binx+1)]
                error += [Hist_Flux[nth_roi].GetBinError(binx+1)]
        roi_xdata += [xdata]
        roi_ydata += [ydata]
        roi_error += [error]
        roi_zeros += [zeros]

    for nth_roi in range(0,len(roi_ydata)):
        for binx in range(0,len(roi_ydata[nth_roi])):
            roi_ydata[nth_roi][binx] = pow(roi_xdata[nth_roi][binx],E_index)*roi_ydata[nth_roi][binx]
            roi_error[nth_roi][binx] = pow(roi_xdata[nth_roi][binx],E_index)*roi_error[nth_roi][binx]

    cycol = cycle('krgbcmy')
    for nth_roi in range(0,len(hist_data)):
        n_1sigma = 0
        for binx in range(0,len(roi_ydata[nth_roi])):
            if roi_error[nth_roi][binx]==0.: continue
            if roi_ydata[nth_roi][binx]/roi_error[nth_roi][binx]>1.:
                n_1sigma += 1
        next_color = next(cycol)
        if doUpperLimit:
            ax_local.fill_between(roi_xdata[nth_roi], roi_zeros[nth_roi], roi_error[nth_roi], alpha=0.2, color='r')
        if 'Crab' in legends[nth_roi] or n_1sigma<4 or doUpperLimit:
            ax_local.errorbar(roi_xdata[nth_roi], roi_ydata[nth_roi], roi_error[nth_roi], color=next_color, marker='s', ls='none', label='%s'%(legends[nth_roi]))
        else:
            start = (roi_ydata[nth_roi][0]/pow(10,-12), -2.)
            popt, pcov = curve_fit(power_law_func, np.array(roi_xdata[nth_roi]), np.array(roi_ydata[nth_roi]), p0=start, sigma=np.array(error))
            ax_local.plot(np.array(roi_xdata[nth_roi]), power_law_func(np.array(roi_xdata[nth_roi]), *popt),color=next_color)
            flux_fit = power_law_func(np.array(roi_xdata[nth_roi]), *popt)
            residual = np.array(roi_ydata[nth_roi]) - flux_fit
            chisq = np.sum((residual/np.array(roi_error[nth_roi]))**2)
            dof = len(roi_xdata[nth_roi])-2
            ax_local.errorbar(roi_xdata[nth_roi], roi_ydata[nth_roi], roi_error[nth_roi], color=next_color, marker='s', ls='none', label='%s, $\Gamma$ = %0.1f, $\chi^{2}/dof = %0.1f$'%(legends[nth_roi],popt[1],chisq/dof))

    log_energy = np.linspace(log10(Hist_Flux[0].GetBinLowEdge(1)),log10(Hist_Flux[0].GetBinLowEdge(Hist_Flux[0].GetNbinsX())),50)
    xdata = pow(10.,log_energy)
    for nth_roi in range(0,len(hist_data)):
        if 'background flux' in legends[nth_roi]:
            print ('background flux = %s'%(np.array(roi_ydata[nth_roi])))
        if 'Crab' in legends[nth_roi] and nth_roi==0:
            vectorize_f = np.vectorize(flux_crab_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='1508.06442')
            if doCalibrate:
                xdata_array = []
                for binx in range(0,Hist_Flux[nth_roi].GetNbinsX()):
                    xdata_array += [Hist_Flux[0].GetBinLowEdge(binx+1)]
                #ydata = pow(np.array(xdata_array),E_index)*vectorize_f(xdata_array)
                #calibration_new = []
                #for binx in range(0,Hist_Flux[nth_roi].GetNbinsX()):
                #    if Hist_Flux[0].GetBinContent(binx+1)>0.:
                #        calibration_new += [ydata[binx]/Hist_Flux[0].GetBinContent(binx+1)]
                #    else:
                #        calibration_new += [0.]
                #print ('calibration_new = %s'%(calibration_new))
        if 'VHE region' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_veritas_j1908_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='1404.7185 (VERITAS)')
        if 'HAWC region' in legends[nth_roi]:
            vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
            ydata_hawc = pow(xdata,E_index)*vectorize_f_hawc(xdata)
            ax_local.plot(xdata, ydata_hawc,'r-',label='1909.08609 (HAWC)')
            vectorize_f = np.vectorize(flux_veritas_j1908_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'b-',label='1404.7185 (VERITAS)')
            vectorize_f = np.vectorize(flux_hess_j1908_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'g-',label='0904.3409 (HESS)')
        if 'IC 443' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_ic443_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='0905.3291')
            vectorize_f = np.vectorize(flux_ic443_hawc_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'g-',label='2007.08582 (HAWC)')
        if '1ES 1218+304' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_1es1218_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='0810.0301')
        if 'Geminga Pulsar' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_geminga_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='HAWC extrapolation')
        if 'J1856+0245' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_j1857_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='From Aleksic et al. (2014)')
        if 'J2019+407' in legends[nth_roi]:
            vectorize_f = np.vectorize(flux_veritas_gamma_cygni_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'r-',label='From Weinstein et al. (2012)')
            vectorize_f = np.vectorize(flux_hawc_gamma_cygni_func)
            ydata = pow(xdata,E_index)*vectorize_f(xdata)
            ax_local.plot(xdata, ydata,'g-',label='3HWC J2020+403')

    ax_local.legend(loc='best')
    ax_local.set_xlabel('Energy [GeV]')
    ax_local.set_ylabel('Flux [$\mathrm{TeV}^{-1}\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    ax_local.set_xscale('log')
    ax_local.set_yscale('log')
    return(ax_local)


def MakeStarEffectPlot(Hist_data,Hist_bkgd,legends,colors,name,range_lower,range_upper):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad = ROOT.TPad("pad","pad",0,0,1,1)
    pad.SetBottomMargin(0.1)
    pad.SetTopMargin(0.1)
    pad.SetLeftMargin(0.15)
    pad.SetBorderMode(1)
    pad.Draw()
    pad.cd()

    norm_bin_low_target = Hist_data[0].FindBin(range_lower)
    norm_bin_up_target = Hist_data[0].FindBin(range_upper)-1

    n_roi = 0
    star_brightness, star_effect, star_brightness_err, star_effect_err = array( 'd' ), array( 'd' ) , array( 'd' ) , array( 'd' )
    for roi in range(0,len(Hist_data)):
        if not 'b-mag' in legends[roi]: continue
        star_brightness.append(float(legends[roi].strip('b-mag ')))
        star_brightness_err.append(0.)
        data_count = Hist_data[roi].Integral(norm_bin_low_target,norm_bin_up_target)
        bkgd_count = Hist_bkgd[roi].Integral(norm_bin_low_target,norm_bin_up_target)
        if (bkgd_count>0.):
            star_effect.append( (data_count-bkgd_count)/bkgd_count )
            star_effect_err.append( pow(data_count,0.5)/bkgd_count )
            n_roi += 1

    gr = ROOT.TGraphErrors( n_roi, star_brightness, star_effect , star_brightness_err, star_effect_err )
    gr.SetLineColor( 1 )
    gr.SetLineWidth( 4 )
    gr.SetMarkerColor( 1 )
    gr.SetMarkerStyle( 21 )
    gr.GetXaxis().SetTitle( 'b-mag' )
    gr.GetYaxis().SetTitle( 'star syst.' )
    gr.Draw( 'AP' )

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeChi2Plot(Hists,legends,colors,stack_it,title,name,doSum,range_lower,range_upper,syst):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,0.8)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
    pad2.SetBottomMargin(0.39)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad1.SetGrid()
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    norm_bin_low_target = 1
    norm_bin_up_target = Hists[0].GetNbinsX()
    if 'MSCW' in name or 'MSCL' in name:
        norm_bin_low_target = Hists[0].FindBin(range_lower)
        norm_bin_up_target = Hists[0].FindBin(range_upper)-1

    max_heigh = 0
    max_hist = 0
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
    #Hists[max_hist].Draw("E")
    Hists[0].Draw("E")

    Hist_Bkgd = Hists[1].Clone()
    Hist_Bkgd.Reset()
    Hist_Syst = Hists[1].Clone()
    Hist_Syst.Reset()
    found_syst_hist = False
    for h in range(0,len(Hists)):
        if stack_it[h] or 'syst' in legends[h]:
            Hist_Syst.Add(Hists[h])
        if 'syst' in legends[h]:
            found_syst_hist = True
        if colors[h]!=0 and colors[h]!=1:
            Hist_Bkgd.Add(Hists[h].Clone())
    #if not found_syst_hist:
    #    for binx in range(0,Hist_Syst.GetNbinsX()):
    #        syst_err = Syst_MDM*Hist_Bkgd.GetBinContent(binx+1)
    #        stat_err = Hist_Bkgd.GetBinError(binx+1)
    #        Hist_Syst.SetBinError(binx+1,pow(syst_err*syst_err+stat_err*stat_err,0.5))

    fill_color = [0,0,46,0,38,30]
    if doSum:
        stack = ROOT.THStack("stack", "")
        for h in range(1,len(Hists)):
            if stack_it[h]:
                set_histStyle( Hists[h] , fill_color[colors[h]])
                stack.Add( Hists[h] )
        stack.Draw("hist same")
    if found_syst_hist:
        Hist_Syst.SetFillColor(1)
        Hist_Syst.SetFillStyle(3004)
        Hist_Syst.SetMarkerSize(0)
        Hist_Syst.Draw("e2 same")

    #for h in range(0,len(Hists)):
    #    if colors[h]==0: continue
    #    Hists[h].SetLineWidth(3)
    #    Hists[h].Draw("E same")
    Hists[0].SetLineWidth(3)
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.6,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    #legend.SetTextSize(0.2)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hists[0],legends[0],"pl")
    for h in range(1,len(Hists)):
        if 'syst' in legends[h]: continue
        if doSum:
            legend.AddEntry(Hists[h],legends[h],"f")
        else:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")
    #lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hists[0],norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkgd,norm_bin_low_target,norm_bin_up_target)
    syst_err = Syst_MDM
    if not 'MDM' in name: syst_err = 0.
    err_bkg = pow(err_bkg*err_bkg+(syst_err*predict_bkg)*(syst_err*predict_bkg),0.5)
    Sig = 1.*CalculateSignificance(data_SR-predict_bkg,predict_bkg,err_bkg)
    lumilab2 = ROOT.TLatex(0.5,0.60,'Excess = %0.1f#pm%0.1f'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5)) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab1 = ROOT.TLatex(0.5,0.40,'Data = %0.1f#pm%0.1f'%(data_SR,err_SR) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab3 = ROOT.TLatex(0.5,0.20,'Bkg = %0.1f#pm%0.1f'%(predict_bkg,err_bkg) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0: 
        sbratio = (data_SR-predict_bkg)/(predict_bkg)
    if not data_SR-predict_bkg==0 and not predict_bkg==0:
        sbratio_err = sbratio*pow(pow(pow(err_SR*err_SR+err_bkg*err_bkg,0.5)/(data_SR-predict_bkg),2)+pow(err_bkg/predict_bkg,2),0.5)
    else: sbratio_err = 0
    #lumilab4 = ROOT.TLatex(0.15,0.20,'S/B = %0.3f#pm%0.3f'%(sbratio,sbratio_err) )
    #lumilab4.SetNDC()
    #lumilab4.SetTextSize(0.15)
    #lumilab4.Draw()

    pad2.cd()
    Hist_Band = Hist_Bkgd.Clone()
    if found_syst_hist:
        Hist_Band.Reset()
        Hist_Band.Add(Hist_Syst)
    Hist_Band.Divide(Hist_Bkgd)
    Hist_Band.SetFillColor(1)
    if not found_syst_hist:
        Hist_Band.SetFillColor(0)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    Hist_Band.GetYaxis().SetTitleSize(0.10)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(1.2)
    Hist_Band.SetMinimum(0.8)
    if 'MJD' in title:
        Hist_Band.SetMaximum(3.0)
        Hist_Band.SetMinimum(0.)
    Hist_Band.Draw("e2")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    Hist_Ratio = Hists[0].Clone()
    Hist_Ratio.Divide(Hist_Bkgd)
    Hist_Ratio.SetLineWidth(2)
    Hist_Ratio.SetLineColor(1)
    Hist_Ratio.Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def PlotsStackedHistograms(tag):

    if energy_bin[energy_bin_cut_low]>=1000.:
        Hist_OnData_Skymap_ProjX_Sum.Rebin(2)
        Hist_OnDark_Skymap_ProjX_Sum.Rebin(2)
        Hist_OnBkgd_Skymap_ProjX_Sum.Rebin(2)
        Hist_OnData_Skymap_ProjY_Sum.Rebin(2)
        Hist_OnDark_Skymap_ProjY_Sum.Rebin(2)
        Hist_OnBkgd_Skymap_ProjY_Sum.Rebin(2)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjX_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Skymap_ProjX_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnBkgd_Skymap_Syst_ProjX_Sum]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Skymap_ProjX_%s'%(tag)
    title = 'RA'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjX_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Skymap_ProjX_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Dark_Skymap_ProjX_%s'%(tag)
    title = 'RA'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjY_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Skymap_ProjY_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnBkgd_Skymap_Syst_ProjY_Sum]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Skymap_ProjY_%s'%(tag)
    title = 'Dec'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjY_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Skymap_ProjY_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Dark_Skymap_ProjY_%s'%(tag)
    title = 'Dec'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_MSCW]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_wGamma_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnGamma_MSCW_Sum]
    legends += ['#gamma fixed at 10% of bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Unblind_wGamma_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_wGamma_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnGamma_MSCL_Sum]
    legends += ['#gamma fixed at 10% of bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Unblind_wGamma_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_woGamma_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Unblind_woGamma_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_woGamma_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Unblind_woGamma_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_MSCW_Sum]
    legends += ['init. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_MSCW_Dark_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_MSCL]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_MSCL_Sum]
    legends += ['init. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_MSCL_Dark_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    #if energy_bin[energy_bin_cut_low]>=1000.:
    Hist_OnData_Theta2_Sum.Rebin(2)
    Hist_OnBkgd_Theta2_Sum.Rebin(2)
    Hist_OnRFoV_Theta2_Sum.Rebin(2)
    Hist_OnDark_Theta2_Sum.Rebin(2)
    Hist_SystErr_Theta2_Sum.Rebin(2)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Yoff_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Yoff_Sum]
    legends += ['bkg. predict.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Yoff_MDM_%s'%(tag)
    title = 'Y off'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_R2off_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_R2off_Sum]
    legends += ['bkg. predict.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_R2off_MDM_%s'%(tag)
    title = 'R2 off'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    xdata = []
    ydata = []
    ydata_err = []
    ybkgd = []
    ybkgd_err = []
    for binx in range(0,Hist_OnData_Theta2_Sum.GetNbinsX()):
        if Hist_OnData_Theta2_Sum.GetBinCenter(binx+1)>4.: continue
        xdata += [Hist_OnData_Theta2_Sum.GetBinCenter(binx+1)]
        ydata += [Hist_OnData_Theta2_Sum.GetBinContent(binx+1)]
        ydata_err += [Hist_OnData_Theta2_Sum.GetBinError(binx+1)]
        ybkgd += [Hist_OnBkgd_Theta2_Sum.GetBinContent(binx+1)]
        bkgd_err = Hist_OnBkgd_Theta2_Sum.GetBinError(binx+1)
        syst_err = Hist_SystErr_Theta2_Sum.GetBinError(binx+1)
        ybkgd_err += [pow(bkgd_err*bkgd_err+syst_err*syst_err,0.5)]
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(xdata,ydata,ydata_err,color='k',marker='s',ls='none',label='data')
    axbig.errorbar(xdata,ybkgd,ybkgd_err,color='r',marker='s',ls='none',label='background')
    axbig.legend(loc='best')
    axbig.set_xlabel('squared angle from source location $\\theta^{2}$')
    axbig.set_ylabel('event count')
    plotname = 'Stack_Theta2_%s'%(tag)
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    axbig.remove()

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Theta2_Sum]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Theta2_MDM_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)
    plotname = 'SignificanceCurve_Theta2_MDM_%s'%(tag)
    title = 'analysis cut on #theta^{2}'
    MakeAccumulatedSignificancePlot(Hist_OnData_Theta2_Sum,Hist_OnBkgd_Theta2_Sum,Hist_SystErr_Theta2_Sum,0,title,'significance',plotname)
    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnRFoV_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Theta2_RFoV_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)
    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Theta2_Dark_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_X_Sum[nth_roi].Rebin(4)
        Hist_OnBkgd_RoI_X_Sum[nth_roi].Rebin(4)
        Hist_OnData_RoI_Y_Sum[nth_roi].Rebin(4)
        Hist_OnBkgd_RoI_Y_Sum[nth_roi].Rebin(4)
    for nth_roi in range(0,len(roi_ra)):
        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_X_Sum[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_X_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_X_MDM_%s'%(nth_roi,tag)
        title = 'rotated X axis'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Y_Sum[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Y_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_Y_MDM_%s'%(nth_roi,tag)
        title = 'rotated Y axis'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    #if energy_bin[energy_bin_cut_low]>=300.:
    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Rebin(2)
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Rebin(2)
    for nth_roi in range(1,len(roi_ra)):
        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Theta2_Sum[nth_roi]]
        legends += ['obs. data (%s)'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        Hists += [Hist_SumSyst_RoI_Theta2[nth_roi]]
        legends += ['syst. error']
        colors += [0]
        stack_it += [False]
        plotname = 'Stack_RoI%s_Theta2_MDM_%s'%(nth_roi,tag)
        title = 'squared angle from source location #theta^{2}'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        legends += ['obs. data (%s)'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_Energy_MDM_%s'%(nth_roi,tag)
        title = 'energy [GeV]'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

        mjd_start = 2147483647
        mjd_end = 0
        for binx in range(0,Hist_OnData_RoI_MJD_Sum[nth_roi].GetNbinsX()):
            bin_content = Hist_OnData_RoI_MJD_Sum[nth_roi].GetBinContent(binx+1)
            bin_lower_edge = Hist_OnData_RoI_MJD_Sum[nth_roi].GetBinLowEdge(binx+1)
            if bin_content>0:
                if mjd_start>bin_lower_edge:
                    mjd_start = bin_lower_edge
                if mjd_end<bin_lower_edge:
                    mjd_end = bin_lower_edge
        hist_mjd_data_plot = ROOT.TH1D("hist_mjd_data_plot","",10,mjd_start,mjd_end)
        hist_mjd_bkgd_plot = ROOT.TH1D("hist_mjd_bkgd_plot","",10,mjd_start,mjd_end)
        for binx in range (0,Hist_OnData_RoI_MJD_Sum[nth_roi].GetNbinsX()):
            mjd_bin_center = Hist_OnData_RoI_MJD_Sum[nth_roi].GetBinCenter(binx+1)
            mjd_data = Hist_OnData_RoI_MJD_Sum[nth_roi].GetBinContent(binx+1)
            mjd_bkgd = Hist_OnBkgd_RoI_MJD_Sum[nth_roi].GetBinContent(binx+1)
            mjd_bin = hist_mjd_data_plot.FindBin(mjd_bin_center)
            mjd_data_old_content = hist_mjd_data_plot.GetBinContent(mjd_bin)
            mjd_bkgd_old_content = hist_mjd_bkgd_plot.GetBinContent(mjd_bin)
            hist_mjd_data_plot.SetBinContent(mjd_bin,mjd_data_old_content+mjd_data)
            hist_mjd_bkgd_plot.SetBinContent(mjd_bin,mjd_bkgd_old_content+mjd_bkgd)
        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [hist_mjd_data_plot]
        legends += ['obs. data (%s)'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [hist_mjd_bkgd_plot]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_MJD_MDM_%s'%(nth_roi,tag)
        title = 'MJD'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hist_data = []
    Hist_bkgd = []
    Hist_syst = []
    radii = []
    legends = []
    for nth_roi in range(1,len(roi_ra)):
        if (roi_name[nth_roi] in exclude_roi): continue
        Hist_data += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        Hist_bkgd += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        Hist_syst += [Hist_SumSyst_RoI_Energy[nth_roi]]
        radii += [roi_radius_outer[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
    title = 'energy [GeV]'
    #fig.clf()
    #axbig = fig.add_subplot()
    #MakeSpectrumInNonCrabUnit(axbig,Hist_data,Hist_bkgd,Hist_syst,radii,legends,title,True,False,0.)
    #plotname = 'Flux_Calibrate_%s'%(tag)
    #fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    #axbig.remove()
    #fig.clf()
    #axbig = fig.add_subplot()
    #MakeSpectrumInNonCrabUnit(axbig,Hist_data,Hist_bkgd,Hist_syst,radii,legends,title,True,False,2.)
    #plotname = 'FluxE2_Calibrate_%s'%(tag)
    #fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    #axbig.remove()
    fig.clf()
    axbig = fig.add_subplot()
    MakeSpectrumInNonCrabUnit(axbig,Hist_data,Hist_bkgd,Hist_syst,radii,legends,title,False,False,0.)
    plotname = 'Flux_NoCalibrate_%s'%(tag)
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    axbig.remove()

    RoI_zscore = []
    for nth_roi in range(0,len(Hist_data)):
        data_total = 0.
        bkgd_total = 0.
        stat_total = 0.
        syst_total = 0.
        for binx in range(0,Hist_data[nth_roi].GetNbinsX()):
            data_total += Hist_data[nth_roi].GetBinContent(binx+1)
            bkgd_total += Hist_bkgd[nth_roi].GetBinContent(binx+1)
            stat_total += pow(Hist_bkgd[nth_roi].GetBinError(binx+1),2)
            syst_total += pow(Hist_syst[nth_roi].GetBinError(binx+1),2)
        syst_total = pow(syst_total+stat_total,0.5)
        zscore = CalculateSignificance(data_total-bkgd_total,bkgd_total,syst_total)
        RoI_zscore += [zscore]
    fig.clf()
    axbig = fig.add_subplot()
    ind = np.arange(len(Hist_data))
    width = 0.35
    rects1 = axbig.bar(ind, RoI_zscore, width, color='#FF9848')
    axbig.set_xlim(-width,len(ind)+width)
    axbig.set_ylabel("z score", fontsize=18)
    xTickMarks = legends
    axbig.set_xticks(ind+0.5*width)
    xtickNames = axbig.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    fig.subplots_adjust(bottom=0.15)
    plotname = 'RoIZscore_MDM_%s'%(tag)
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    axbig.remove()

    #Hist_data = []
    #Hist_bkgd = []
    #Hist_syst = []
    #legends = []
    #Hist_data += [Hist_OnData_RoI_Energy_Sum[0]]
    #Hist_bkgd += [Hist_OnBkgd_RoI_Energy_Sum[0]]
    #Hist_syst += [Hist_SumSyst_RoI_Energy[0]]
    #legends += ['background flux']
    #title = 'energy [GeV]'
    #fig.clf()
    #axbig = fig.add_subplot()
    #MakeSpectrumInNonCrabUnit(axbig,Hist_data,Hist_bkgd,Hist_syst,radii,legends,title,False,True,2.7)
    #plotname = 'FluxE27_MDM_%s'%(tag)
    #fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    #axbig.remove()

    #Hist_data_energy = []
    #Hist_bkgd_energy = []
    #legends = []
    #colors = []
    #for nth_roi in range(0,len(roi_ra)):
    #    legends += ['%s'%(roi_name[nth_roi])]
    #    colors += [nth_roi+1]
    #    Hist_data_energy += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
    #    Hist_bkgd_energy += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
    #plotname = 'Spectrum_RoI%s_CrabUnit_MDM_%s'%(nth_roi,tag)
    #title = 'GeV'
    #MakeCrabUnitSpectrumPlot(Hist_data_energy,Hist_bkgd_energy,legends,colors,title,plotname)

    #Hist_data_theta2 = []
    #Hist_bkgd_theta2 = []
    #legends = []
    #colors = []
    #for nth_roi in range(0,len(roi_ra)):
    #    legends += ['%s'%(roi_name[nth_roi])]
    #    colors += [nth_roi+1]
    #    Hist_data_theta2 += [Hist_OnData_RoI_Theta2_Sum[nth_roi]]
    #    Hist_bkgd_theta2 += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
    #plotname = 'StarEffect_RoI%s_MDM_%s'%(nth_roi,tag)
    #MakeStarEffectPlot(Hist_data_theta2,Hist_bkgd_theta2,legends,colors,plotname,0.,0.075)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Energy]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Energy_MDM_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnRFoV_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Energy]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Energy_RFoV_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_CamCenter_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Energy_CamCenter_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Energy_CamCenter_MDM_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Energy_Dark_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Zenith_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Zenith_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Zenith_MDM_%s'%(tag)
    title = 'Zenith [degree]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Height_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Height_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Height_MDM_%s'%(tag)
    title = 'Height [km]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Depth_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Depth_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Depth_MDM_%s'%(tag)
    title = 'Depth [km]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Rcore_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Rcore_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Rcore_MDM_%s'%(tag)
    title = 'R core [m]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)


    #for nth_sample in range(0,n_control_samples):

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_MSCW_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_MSCW_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_MSCW_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'MSCW'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCW_plot_lower,MSCW_blind_cut,-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_MSCL_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_MSCL_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_MSCL_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'MSCL'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCL_plot_lower,MSCL_blind_cut,-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_Energy_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_Energy_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_Energy_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'energy [GeV]'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,pow(10,4.0),-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    #Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
    #    #Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
    #    Hists += [Hist_OffData_CameraFoV_Theta2_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_CameraFoV_Theta2_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'squared angle from camera center #theta^{2}'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,1.,-1)

def CalculateSystErrorIndividual_v3():

    Syst_MDM_single = 0.

    binx_low_target = Hist_OnData_MSCL.FindBin(MSCL_plot_lower)
    binx_up_target = Hist_OnData_MSCL.FindBin(MSCL_plot_upper)-1
    binx_blind_target = Hist_OnData_MSCL.FindBin(MSCL_blind_cut)-1
    biny_low_target = Hist_OnData_MSCW.FindBin(MSCW_plot_lower)
    biny_up_target = Hist_OnData_MSCW.FindBin(MSCW_plot_upper)-1
    biny_blind_target = Hist_OnData_MSCW.FindBin(MSCW_blind_cut)-1

    #Total_OnBkgd = Hist2D_OnBkgd.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
    Total_OnBkgd = Hist2D_OnBkgd.Integral()
    if Total_OnBkgd==0.: return Syst_MDM_single

    Hist2D_Diff = Hist2D_OnBkgd.Clone()
    #Hist2D_Diff = Hist2D_OnDark.Clone()
    Hist2D_Diff.Add(Hist2D_OnData,-1.)
    for binx in range(0,Hist2D_Diff.GetNbinsX()):
        for biny in range(0,Hist2D_Diff.GetNbinsY()):
            if binx<=binx_blind_target and biny<=biny_blind_target: continue
            #if binx>binx_blind_target and biny>biny_blind_target: continue
            Syst_MDM_single += abs(Hist2D_Diff.GetBinContent(binx+1,biny+1))
    Syst_MDM_single = 2.*(pow(Syst_MDM_single,1.0))/Total_OnBkgd
    #Syst_MDM_single = 1.*pow(max(0.,Syst_MDM_single-Hist2D_OnData.Integral()+Total_OnBkgd),0.5)/Total_OnBkgd

    return Syst_MDM_single

def CalculateSystError_v3():

    global Syst_MDM
    Syst_MDM = 0.

    binx_low_target = Hist_OnData_MSCL.FindBin(MSCL_plot_lower)
    binx_up_target = Hist_OnData_MSCL.FindBin(MSCL_plot_upper)-1
    binx_blind_target = Hist_OnData_MSCL.FindBin(MSCL_blind_cut)-1
    biny_low_target = Hist_OnData_MSCW.FindBin(MSCW_plot_lower)
    biny_up_target = Hist_OnData_MSCW.FindBin(MSCW_plot_upper)-1
    biny_blind_target = Hist_OnData_MSCW.FindBin(MSCW_blind_cut)-1

    #Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
    Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral()
    if Total_OnBkgd==0.: return

    Hist2D_Diff = Hist2D_OnBkgd_Sum.Clone()
    #Hist2D_Diff = Hist2D_OnDark_Sum.Clone()
    Hist2D_Diff.Add(Hist2D_OnData_Sum,-1.)
    for binx in range(0,Hist2D_Diff.GetNbinsX()):
        for biny in range(0,Hist2D_Diff.GetNbinsY()):
            if binx<=binx_blind_target and biny<=biny_blind_target: continue
            #if binx>binx_blind_target and biny>biny_blind_target: continue
            Syst_MDM += abs(Hist2D_Diff.GetBinContent(binx+1,biny+1))
    Syst_MDM = 2.*(pow(Syst_MDM,1.0))/Total_OnBkgd
    #Syst_MDM = 1.*pow(max(0.,Syst_MDM-Hist2D_OnData_Sum.Integral()+Total_OnBkgd),0.5)/Total_OnBkgd

def Theta2HistScale(Hist,scale,scale_err):

    for b in range(1,Hist.GetNbinsX()+1):
        old_content = Hist.GetBinContent(b)
        old_error = Hist.GetBinError(b)
        new_content = old_content*scale
        new_error = 0
        if old_content>0 and scale>0:
            #new_error = new_content*(old_error/old_content)
            new_error = pow(pow(old_content*scale_err,2)+pow(old_error*scale,2),0.5)
        Hist.SetBinContent(b,new_content)

def NormalizeEnergyHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Energy.Reset()
    Hist_OnData_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Energy.Reset()
    Hist_OnBkgd_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnRFoV_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnRFoV_Energy.Reset()
    Hist_OnRFoV_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_SR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_Energy.Reset()
    Hist_OnDark_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Energy_CamCenter_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Energy_CamCenter.Reset()
    Hist_OnData_Energy_CamCenter.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Energy_CamCenter_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Energy_CamCenter.Reset()
    Hist_OnBkgd_Energy_CamCenter.Add(InputFile.Get(HistName))

    HistName = "Hist_OnData_SR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Zenith.Reset()
    Hist_OnData_Zenith.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Zenith.Reset()
    Hist_OnBkgd_Zenith.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Height_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Height.Reset()
    Hist_OnData_Height.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Height_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Height.Reset()
    Hist_OnBkgd_Height.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Depth_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Depth.Reset()
    Hist_OnData_Depth.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Depth_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Depth.Reset()
    Hist_OnBkgd_Depth.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Rcore_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Rcore.Reset()
    Hist_OnData_Rcore.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Rcore_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Rcore.Reset()
    Hist_OnBkgd_Rcore.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print ('read histogram: %s'%(HistName))
        Hist_OnData_RoI_Energy[nth_roi].Reset()
        Hist_OnData_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
        Hist_OnBkgd_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_NormSyst_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_NormSyst_RoI_Energy[nth_roi].Reset()
        Hist_NormSyst_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_ShapeSyst_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_ShapeSyst_RoI_Energy[nth_roi].Reset()
        Hist_ShapeSyst_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_SR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_MJD[nth_roi].Reset()
        Hist_OnData_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_MJD[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))

    InputFile.Close()

    #if Hist2D_OnData.Integral()<1600.:
    #    Hist_OnData_Energy.Reset()
    #    Hist_OnBkgd_Energy.Reset()
    #    Hist_OnRFoV_Energy.Reset()
    #    Hist_OnData_Energy_CamCenter.Reset()
    #    Hist_OnBkgd_Energy_CamCenter.Reset()
    #    Hist_OnDark_Energy.Reset()
    #    Hist_OnData_Zenith.Reset()
    #    Hist_OnBkgd_Zenith.Reset()
    #    for nth_roi in range(0,len(roi_ra)):
    #        Hist_OnData_RoI_Energy[nth_roi].Reset()
    #        Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
    #        Hist_OnData_RoI_MJD[nth_roi].Reset()
    #        Hist_OnBkgd_RoI_MJD[nth_roi].Reset()

def NormalizeTheta2Histograms(FilePath,ebin):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_SR_XYoff.Reset()
    Hist_OnData_SR_XYoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_SR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_SR_XYoff.Reset()
    Hist_OnDark_SR_XYoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_CR_XYoff.Reset()
    Hist_OnData_CR_XYoff.Add(InputFile.Get(HistName))
    Hist_Data_CR_XYoff[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Yoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Yoff.Reset()
    Hist_OnData_Yoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_R2off.Reset()
    Hist_OnData_R2off.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_ISR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_ISR_R2off.Reset()
    Hist_OnData_ISR_R2off.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_R2off.Reset()
    Hist_OnBkgd_R2off.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Yoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Yoff.Reset()
    Hist_OnBkgd_Yoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Yoff_Raw_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Yoff_Raw.Reset()
    Hist_OnBkgd_Yoff_Raw.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Theta2.Reset()
    Hist_OnData_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Theta2.Reset()
    Hist_OnBkgd_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnRFoV_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnRFoV_Theta2.Reset()
    Hist_OnRFoV_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_Theta2.Reset()
    Hist_OnDark_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_NormSyst_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_NormSyst_Theta2.Reset()
    Hist_NormSyst_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_ShapeSyst_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_ShapeSyst_Theta2.Reset()
    Hist_ShapeSyst_Theta2.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Theta2[nth_roi].Reset()
        Hist_OnData_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()
        Hist_OnBkgd_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_NormSyst_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_NormSyst_RoI_Theta2[nth_roi].Reset()
        Hist_NormSyst_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_ShapeSyst_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_ShapeSyst_RoI_Theta2[nth_roi].Reset()
        Hist_ShapeSyst_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_SR_Skymap_RoI_X_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_X[nth_roi].Reset()
        Hist_OnData_RoI_X[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_Skymap_RoI_X_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_X[nth_roi].Reset()
        Hist_OnBkgd_RoI_X[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_SR_Skymap_RoI_Y_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Y[nth_roi].Reset()
        Hist_OnData_RoI_Y[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_Skymap_RoI_Y_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Y[nth_roi].Reset()
        Hist_OnBkgd_RoI_Y[nth_roi].Add(InputFile.Get(HistName))

    InputFile.Close()

    #if Hist2D_OnData.Integral()<1600.:
    #    Hist_OnData_SR_XYoff.Reset()
    #    Hist_OnDark_SR_XYoff.Reset()
    #    Hist_OnData_CR_XYoff.Reset()
    #    Hist_OnData_Yoff.Reset()
    #    Hist_OnBkgd_Yoff.Reset()
    #    Hist_OnBkgd_Yoff_Raw.Reset()
    #    Hist_OnData_Theta2.Reset()
    #    Hist_OnBkgd_Theta2.Reset()
    #    Hist_OnRFoV_Theta2.Reset()
    #    Hist_OnDark_Theta2.Reset()
    #    for nth_roi in range(0,len(roi_ra)):
    #        Hist_OnData_RoI_Theta2[nth_roi].Reset()
    #        Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()

def StackEnergyHistograms():

    Hist_OnData_Energy_Sum.Add(Hist_OnData_Energy)
    Hist_OnBkgd_Energy_Sum.Add(Hist_OnBkgd_Energy)
    Hist_OnRFoV_Energy_Sum.Add(Hist_OnRFoV_Energy)
    Hist_OnData_Energy_CamCenter_Sum.Add(Hist_OnData_Energy_CamCenter)
    Hist_OnBkgd_Energy_CamCenter_Sum.Add(Hist_OnBkgd_Energy_CamCenter)
    Hist_OnDark_Energy_Sum.Add(Hist_OnDark_Energy)
    Hist_OnData_Zenith_Sum.Add(Hist_OnData_Zenith)
    Hist_OnBkgd_Zenith_Sum.Add(Hist_OnBkgd_Zenith)
    Hist_OnData_Height_Sum.Add(Hist_OnData_Height)
    Hist_OnBkgd_Height_Sum.Add(Hist_OnBkgd_Height)
    Hist_OnData_Depth_Sum.Add(Hist_OnData_Depth)
    Hist_OnBkgd_Depth_Sum.Add(Hist_OnBkgd_Depth)
    Hist_OnData_Rcore_Sum.Add(Hist_OnData_Rcore)
    Hist_OnBkgd_Rcore_Sum.Add(Hist_OnBkgd_Rcore)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Add(Hist_OnData_RoI_Energy[nth_roi])
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Energy[nth_roi])
        Hist_OnData_RoI_MJD_Sum[nth_roi].Add(Hist_OnData_RoI_MJD[nth_roi])
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Add(Hist_OnBkgd_RoI_MJD[nth_roi])

def StackTheta2Histograms():

    Hist_OnData_SR_XYoff_Sum.Add(Hist_OnData_SR_XYoff)
    Hist_OnDark_SR_XYoff_Sum.Add(Hist_OnDark_SR_XYoff)
    Hist_OnData_CR_XYoff_Sum.Add(Hist_OnData_CR_XYoff)
    Hist_OnData_R2off_Sum.Add(Hist_OnData_R2off)
    Hist_OnData_ISR_R2off_Sum.Add(Hist_OnData_ISR_R2off)
    Hist_OnData_Yoff_Sum.Add(Hist_OnData_Yoff)
    Hist_OnBkgd_R2off_Sum.Add(Hist_OnBkgd_R2off)
    Hist_OnBkgd_Yoff_Sum.Add(Hist_OnBkgd_Yoff)
    Hist_OnBkgd_Yoff_Raw_Sum.Add(Hist_OnBkgd_Yoff_Raw)
    Hist_OnData_Theta2_Sum.Add(Hist_OnData_Theta2)
    Hist_OnBkgd_Theta2_Sum.Add(Hist_OnBkgd_Theta2)
    Hist_OnRFoV_Theta2_Sum.Add(Hist_OnRFoV_Theta2)
    Hist_OnDark_Theta2_Sum.Add(Hist_OnDark_Theta2)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Add(Hist_OnData_RoI_Theta2[nth_roi])
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Theta2[nth_roi])
        Hist_OnData_RoI_X_Sum[nth_roi].Add(Hist_OnData_RoI_X[nth_roi])
        Hist_OnBkgd_RoI_X_Sum[nth_roi].Add(Hist_OnBkgd_RoI_X[nth_roi])
        Hist_OnData_RoI_Y_Sum[nth_roi].Add(Hist_OnData_RoI_Y[nth_roi])
        Hist_OnBkgd_RoI_Y_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Y[nth_roi])

def NormalizeCameraFoVHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    for nth_sample in range(0,n_control_samples):

        HistName = "Hist_OffData_SR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffData_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))
        HistName = "Hist_OffData_CR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))

        #if Hist2D_OffData[nth_sample].Integral()<1600.:
        #    Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
        #    Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()

        data_total, data_err = IntegralAndError(Hist_OffData_Energy[nth_sample],bin_lower,bin_upper)
        if data_total==0.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
            return

        old_integral = Hist_OffData_CameraFoV_Theta2[nth_sample].Integral()
        data_scale = 0
        data_scale_err = 0
        if not data_total==0 and not old_integral==0:
            data_scale = data_total/old_integral
            data_scale_err = data_scale*(data_err/data_total)
        else:
            data_scale = 0
            data_scale_err = 0
        if not data_total==0:
            Theta2HistScale(Hist_OffData_CameraFoV_Theta2[nth_sample],data_scale,data_scale_err)
        else:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Scale(0)

        bkg_total, bkg_err = IntegralAndError(Hist_OffBkgd_Energy[nth_sample],bin_lower,bin_upper)
        if bkg_total==0.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
            return

        old_integral = Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Integral()
        bkgd_scale = 0
        bkgd_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            bkgd_scale = bkg_total/old_integral
            bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        else:
            bkgd_scale = 0
            bkgd_scale_err = 0
        if not bkg_total==0:
            Theta2HistScale(Hist_OffBkgd_CameraFoV_Theta2[nth_sample],bkgd_scale,bkgd_scale_err)
        else:
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Scale(0)

    InputFile.Close()

def StackCameraFoVHistograms():

    for nth_sample in range(0,n_control_samples):
        Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].Add(Hist_OffData_CameraFoV_Theta2[nth_sample])
        Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].Add(Hist_OffBkgd_CameraFoV_Theta2[nth_sample])

def RaDecHistScale(Hist,scale,scale_err):

    for bx in range(1,Hist.GetNbinsX()+1):
        for by in range(1,Hist.GetNbinsY()+1):
            old_content = Hist.GetBinContent(bx,by)
            old_error = Hist.GetBinError(bx,by)
            new_content = old_content*scale
            new_error = 0
            if scale>0:
                new_error = old_error*scale
            Hist.SetBinContent(bx,by,new_content)
            Hist.SetBinError(bx,by,new_error)

def NormalizeSkyMapHistograms(FilePath,ebin):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_Data_Skymap"
    Hist_Data_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Elev_Skymap"
    Hist_Data_Elev_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Azim_Skymap"
    Hist_Data_Azim_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_NSB_Skymap"
    Hist_Data_NSB_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_MJD_Skymap"
    Hist_Data_MJD_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_ElevNSB"
    Hist_Data_ElevNSB.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_ElevAzim"
    Hist_Data_ElevAzim.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_ElevNSB"
    Hist_Dark_ElevNSB.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_ElevAzim"
    Hist_Dark_ElevAzim.Add(InputFile.Get(HistName))

    HistName = "Hist_OnData_EffArea_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_EffArea_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_ISR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Expo_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap.Reset()
    Hist_OnData_Skymap.Add(InputFile.Get(HistName))
    Hist_Data_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_Skymap.Reset()
    Hist_OnDark_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap_Galactic.Reset()
    Hist_OnData_Skymap_Galactic.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap.Reset()
    Hist_OnBkgd_Skymap.Add(InputFile.Get(HistName))
    Hist_Bkgd_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    Hist_Expo_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnRFoV_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Rfov_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_NormSyst_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Syst_Norm.Reset()
    Hist_OnBkgd_Skymap_Syst_Norm.Add(InputFile.Get(HistName))
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_Norm.GetNbinsX()):
        for biny in range(0,Hist_OnBkgd_Skymap_Syst_Norm.GetNbinsY()):
            old_content = Hist_NormSyst_Energy_Skymap[ebin].GetBinContent(binx+1,biny+1)
            new_content = Hist_OnBkgd_Skymap_Syst_Norm.GetBinContent(binx+1,biny+1)
            Hist_NormSyst_Energy_Skymap[ebin].SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content*new_content,0.5))
            #Hist_NormSyst_Energy_Skymap[ebin].SetBinContent(binx+1,biny+1,old_content+new_content)
    for xy_bin in range(0,len(integration_radii)):
        HistName = "Hist_ShapeSyst_Skymap_ErecS%sto%s_Bin%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int,xy_bin)
        Hist_OnBkgd_Skymap_Syst_Shape.Reset()
        Hist_OnBkgd_Skymap_Syst_Shape.Add(InputFile.Get(HistName))
        for binx in range(0,Hist_OnBkgd_Skymap_Syst_Shape.GetNbinsX()):
            for biny in range(0,Hist_OnBkgd_Skymap_Syst_Shape.GetNbinsY()):
                old_content = Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin].GetBinContent(binx+1,biny+1)
                new_content = Hist_OnBkgd_Skymap_Syst_Shape.GetBinContent(binx+1,biny+1)
                Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin].SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content*new_content,0.5))
    HistName = "Hist_OnData_CR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Galactic.Reset()
    Hist_OnBkgd_Skymap_Galactic.Add(InputFile.Get(HistName))

    InputFile.Close()

    #if Hist2D_OnData.Integral()<1600.:
    #    Hist_OnData_Skymap.Reset()
    #    Hist_OnDark_Skymap.Reset()
    #    Hist_OnData_Skymap_Galactic.Reset()
    #    Hist_OnBkgd_Skymap.Reset()
    #    Hist_OnBkgd_Skymap_Galactic.Reset()

def StackSkymapHistograms(ebin):

    Hist_OnData_Skymap_Sum.Add(Hist_OnData_Skymap)
    Hist_OnDark_Skymap_Sum.Add(Hist_OnDark_Skymap)
    Hist_OnData_Skymap_Galactic_Sum.Add(Hist_OnData_Skymap_Galactic)
    Hist_OnBkgd_Skymap_Sum.Add(Hist_OnBkgd_Skymap)
    Hist_OnBkgd_Skymap_Galactic_Sum.Add(Hist_OnBkgd_Skymap_Galactic)

    #Syst_MDM_single = CalculateSystErrorIndividual_v3()
    Syst_MDM_single = Syst_MDM
    print ('Syst_MDM_single = %s'%(Syst_MDM_single))

    for binx in range(0,Hist_OnBkgd_Skymap_Syst_Sum.GetNbinsX()):
        for biny in range(0,Hist_OnBkgd_Skymap_Syst_Sum.GetNbinsY()):
            old_content = Hist_OnBkgd_Skymap_Syst_Sum.GetBinContent(binx+1,biny+1)
            new_content_norm = Hist_OnBkgd_Skymap_Syst_Norm.GetBinContent(binx+1,biny+1)
            new_content_shape = Hist_OnBkgd_Skymap_Syst_Shape.GetBinContent(binx+1,biny+1)
            Hist_OnBkgd_Skymap_Syst_Sum.SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content_norm*new_content_norm+new_content_shape*new_content_shape,0.5))
    for binx in range(0,Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetNbinsX()):
        for biny in range(0,Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetNbinsY()):
            old_content = Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetBinContent(binx+1,biny+1)
            new_content = Hist_OnBkgd_Skymap_Galactic.GetBinContent(binx+1,biny+1)*Syst_MDM_single
            Hist_OnBkgd_Skymap_Galactic_Syst_MDM.SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content*new_content,0.5))

    for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
        old_err = Hist_SystErr_MSCL.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_MSCL.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_MSCL.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
        old_err = Hist_SystErr_MSCW.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_MSCW.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_MSCW.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_Energy.GetNbinsX()):
        old_err = Hist_SystErr_Energy.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_Energy.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_Energy.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_Theta2_Sum.GetNbinsX()):
        old_err = Hist_SystErr_Theta2_Sum.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_NormSyst_Theta2.GetBinContent(binx+1),2)+pow(Hist_ShapeSyst_Theta2.GetBinContent(binx+1),2)
        Hist_SystErr_Theta2_Sum.SetBinError(binx+1,pow(new_err,0.5))

    for nth_roi in range(0,len(roi_ra)):
        print ('Hist_NormSyst_RoI_Energy[nth_roi]:')
        Hist_NormSyst_RoI_Energy[nth_roi].Print("All")
        print ('Hist_ShapeSyst_RoI_Energy[nth_roi]:')
        Hist_ShapeSyst_RoI_Energy[nth_roi].Print("All")
        for binx in range(0,Hist_SumSyst_RoI_Energy[nth_roi].GetNbinsX()):
            old_err = Hist_SumSyst_RoI_Energy[nth_roi].GetBinError(binx+1)
            new_err = pow(old_err,2)+pow(Hist_NormSyst_RoI_Energy[nth_roi].GetBinContent(binx+1),2)+pow(Hist_ShapeSyst_RoI_Energy[nth_roi].GetBinContent(binx+1),2)
            Hist_SumSyst_RoI_Energy[nth_roi].SetBinError(binx+1,pow(new_err,0.5))

    for nth_roi in range(0,len(roi_ra)):
        for binx in range(0,Hist_SumSyst_RoI_Theta2[nth_roi].GetNbinsX()):
            old_err = Hist_SumSyst_RoI_Theta2[nth_roi].GetBinError(binx+1)
            new_err = pow(old_err,2)+pow(Hist_NormSyst_RoI_Theta2[nth_roi].GetBinContent(binx+1),2)+pow(Hist_ShapeSyst_RoI_Theta2[nth_roi].GetBinContent(binx+1),2)
            Hist_SumSyst_RoI_Theta2[nth_roi].SetBinError(binx+1,pow(new_err,0.5))

def CorrectLEE(z_score,n_tests,threshold_sigma):

    if abs(z_score)<threshold_sigma:
        return z_score
    new_z_score = 0.
    if z_score>0.:
        p_value = 1.-st.norm.cdf(z_score)
        p_value = min(p_value*n_tests,0.9999)
        new_z_score = st.norm.ppf(1.-p_value)
    else:
        p_value = st.norm.cdf(z_score)
        p_value = min(p_value*n_tests,0.9999)
        new_z_score = st.norm.ppf(p_value)
    return min(abs(new_z_score),z_score)

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn):

    return CommonPlotFunctions.GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn)

def reflectXaxis(hist):

    return CommonPlotFunctions.reflectXaxis(hist)

def FillSkymapHoles(hist_map, map_resolution):

    return CommonPlotFunctions.FillSkymapHoles(hist_map, map_resolution)

def ConvertRaDecToGalacticMap(hist_map_radec,hist_map_galactic):

    hist_map_galactic.Reset()
    for bx1 in range(1,hist_map_galactic.GetNbinsX()+1):
        for by1 in range(1,hist_map_galactic.GetNbinsY()+1):
            locationx1 = hist_map_galactic.GetXaxis().GetBinCenter(bx1)
            locationy1 = hist_map_galactic.GetYaxis().GetBinCenter(by1)
            locationx2, locationy2 = ConvertGalacticToRaDec(locationx1,locationy1)
            bx2 = hist_map_radec.GetXaxis().FindBin(locationx2)
            by2 = hist_map_radec.GetYaxis().FindBin(locationy2)
            content = hist_map_radec.GetBinContent(bx2,by2)
            hist_map_galactic.SetBinContent(bx1,by1,content)
    return hist_map_galactic

def GetSkyViewMap(map_file, hist_map, isRaDec):

    return CommonPlotFunctions.GetSkyViewMap(map_file, hist_map, isRaDec)

def GetGalacticCoordMap(map_file, hist_map, isRaDec):

    return CommonPlotFunctions.GetGalacticCoordMap(map_file, hist_map, isRaDec)

def GetCOSkymap(hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open('MWL_maps/CO_skymap_v2.txt')
    if isRaDec: 
        inputFile = open('MWL_maps/CO_skymap_RaDec_v2.txt')
    for line in inputFile:
        sig = float(line.split(' ')[0])
        l = float(line.split(' ')[1])
        b = float(line.split(' ')[2])
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,0.))

    map_resolution = 0.2
    #map_resolution = 0.4
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def OutputToTxtFile(hist_input,name):

    file1 = open(name,"w")
    for binx in range(0,hist_input.GetNbinsX()):
        for biny in range(0,hist_input.GetNbinsY()):
            file1.write('x: %s y: %s sigma: %s \n'%(hist_input.GetXaxis().GetBinCenter(binx+1),hist_input.GetYaxis().GetBinCenter(biny+1),hist_input.GetBinContent(binx+1,biny+1)))
    file1.close()

def Count5SigmaBins(Hist_sigma,Hist_excess,threshold_sigma):

    skymap_bin_size_x = Hist_sigma.GetXaxis().GetBinCenter(2)-Hist_sigma.GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = Hist_sigma.GetYaxis().GetBinCenter(2)-Hist_sigma.GetYaxis().GetBinCenter(1)
    unit_area = skymap_bin_size_x*skymap_bin_size_y
    n_tests = Hist_sigma.GetNbinsX()*Hist_sigma.GetNbinsY()
    global_p_value = 1.
    nbins = 0
    for binx in range(Hist_sigma.GetNbinsX()):
        for biny in range(Hist_sigma.GetNbinsY()):
            Sig = Hist_sigma.GetBinContent(binx+1,biny+1)
            if Sig>threshold_sigma:
                p_value = 1.-st.norm.cdf(Sig)
                p_value_correct = min(p_value*n_tests,1.0)
                global_p_value = global_p_value*p_value_correct
                nbins += 1
    if nbins==0:
        return Hist_sigma.GetMaximum()
    return st.norm.ppf(1.-global_p_value)

def VariableSkymapBins(Hist_Data_input,Hist_Bkgd_input,Hist_Syst_input):

    threshold_sigma = 2.5
    max_nbins = 60
    #list_nbins = [60,40,30,20,15,12,10,8,6,5,4,3,2,1]
    list_nbins = [30,20,15,12,10,6,5,4,3,2,1]

    hist_maxsig = ROOT.TH1D("hist_maxsig","",len(list_nbins),0,len(list_nbins))
    list_maxsig = []
    hist_maxsig_location = ROOT.TH1D("hist_maxsig_location","",len(list_nbins),0,len(list_nbins))
    list_maxsig_location = []
    list_binsize = []

    max_nbins_5sigma = 0
    best_nbins = 10
    for nbins in list_nbins:
        Hist_Data = reflectXaxis(Hist_Data_input)
        Hist_Bkgd = reflectXaxis(Hist_Bkgd_input)
        Hist_Syst = reflectXaxis(Hist_Syst_input)
        Hist_Data.Rebin2D(int(max_nbins/nbins),int(max_nbins/nbins))
        Hist_Bkgd.Rebin2D(int(max_nbins/nbins),int(max_nbins/nbins))
        Hist_Syst.Rebin2D(int(max_nbins/nbins),int(max_nbins/nbins))
        MapEdge_left = Hist_Data.GetXaxis().GetBinLowEdge(1)
        MapEdge_right = Hist_Data.GetXaxis().GetBinLowEdge(Hist_Data.GetNbinsX()+1)
        MapEdge_lower = Hist_Data.GetYaxis().GetBinLowEdge(1)
        MapEdge_upper = Hist_Data.GetYaxis().GetBinLowEdge(Hist_Data.GetNbinsY()+1)
        MapCenter_x = (MapEdge_right+MapEdge_left)/2.
        MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.

        Hist_Zscore = GetSignificanceMap(Hist_Data,Hist_Bkgd,Hist_Syst,False)
        max_sigma = max(abs(Hist_Zscore.GetMaximum()),abs(Hist_Zscore.GetMinimum()))
        max_sigma = CorrectLEE(max_sigma,Hist_Zscore.GetNbinsX()*Hist_Zscore.GetNbinsY(),1.0)
        list_maxsig += [max_sigma]
        max_sigma_x, max_sigma_y = FindHistAbsMaxBinXY(Hist_Zscore)
        max_sigma_x = max_sigma_x-MapCenter_x
        max_sigma_y = max_sigma_y-MapCenter_y
        list_maxsig_location += [pow(max_sigma_x*max_sigma_x+max_sigma_y*max_sigma_y,0.5)]
        skymap_bin_size_x = Hist_Data.GetXaxis().GetBinCenter(2)-Hist_Data.GetXaxis().GetBinCenter(1)
        skymap_bin_size_y = Hist_Data.GetYaxis().GetBinCenter(2)-Hist_Data.GetYaxis().GetBinCenter(1)
        list_binsize += [skymap_bin_size_x]

    return list_maxsig, list_binsize

def GetEffectiveAreaWeight(map_x,map_y,hist_effarea,hist_expo,exposure_hours,ebin):

    flux_calibration = [1.844566095034703, 2.026946094149479, 2.4147010122690085, 2.308880464973313, 2.0379949496452765, 1.6185467800446822]

    binx = hist_effarea.GetXaxis().FindBin(map_x)
    biny = hist_effarea.GetYaxis().FindBin(map_y)
    effarea_max = hist_effarea.GetMaximum()
    effarea = hist_effarea.GetBinContent(binx,biny)
    effarea = max(0.3*effarea_max,effarea)*10000.*exposure_hours*3600.
    #exposure = hist_expo.GetBinContent(binx,biny)
    exposure = hist_expo.GetMaximum()

    effarea_weight = exposure/effarea*flux_calibration[ebin]

    return effarea_weight

def GetFluxCalibration(map_x,map_y,energy):

    return 1.

    binx = Hist_Data_Skymap.GetXaxis().FindBin(map_x)
    biny = Hist_Data_Skymap.GetYaxis().FindBin(map_y)
    data_count = Hist_Data_Skymap.GetBinContent(binx,biny)
    data_elev = Hist_Data_Elev_Skymap.GetBinContent(binx,biny)
    if data_count==0.: return 0.
    avg_elev = data_elev/data_count
    
    flux_calibration = [3.715619485492126e-08, 3.4686074327294743e-09, 2.5491285892959776e-10, 1.6854055449950507e-11]
    #flux_calibration = []
    #if avg_elev<=85. and avg_elev>75.:
    #    flux_calibration = [3.352236564473935e-08, 2.658740853181333e-09, 2.0324882825007042e-10, 1.6339503640997207e-11]
    #if avg_elev<=75. and avg_elev>65.:
    #    flux_calibration = [3.749034537077448e-08, 3.161714293121568e-09, 2.1383421146594705e-10, 1.6093279456176878e-11]
    #if avg_elev<=65. and avg_elev>55.:
    #    flux_calibration = [3.914394112202656e-08, 3.8960057300817696e-09, 2.4197334576912557e-10, 1.5785941060327452e-11]
    #if avg_elev<=55. and avg_elev>45.:
    #    flux_calibration = [2.920184347950322e-08, 4.934389551239988e-09, 2.9907377075292634e-10, 1.8132248320894532e-11]
    #if avg_elev<=45. and avg_elev>35.:
    #    flux_calibration = [0.0, 5.58195550720075e-09, 4.624524584170672e-10, 2.808212850301435e-11]

    return flux_calibration[energy]


def MakeSpectrumIndexSkymap(exposure_in_hours,hist_data,hist_bkgd,hist_rfov,hist_normsyst,hist_syst,hist_expo,hist_effarea,title_x,title_y,name,zoomin_scale):

    global calibration
    global calibration_radius
    global Hist_MWL_global
    global current_gamma_energy
    global next_gamma_energy
    global doReferenceFlux

    isRaDec = False
    if 'RaDec' in name: isRaDec = True

    other_star_labels = []
    other_star_markers = []
    other_star_names = []
    bright_star_labels = []
    bright_star_markers = []
    faint_star_labels = []
    faint_star_markers = []
    other_star_labels_gal = []
    other_star_markers_gal = []
    other_star_names_gal = []
    bright_star_labels_gal = []
    bright_star_markers_gal = []
    faint_star_labels_gal = []
    faint_star_markers_gal = []
    star_range = 3.0/zoomin_scale
    for star in range(0,len(other_stars)):
        if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
        other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
        other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
        other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
        other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
        other_star_names += [other_stars[star]]
    for star in range(0,len(bright_star_ra)):
        bright_star_markers += [ROOT.TMarker(-bright_star_ra[star],bright_star_dec[star],30)]
        bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
        bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
        bright_star_labels += [ROOT.TLatex(-bright_star_ra[star]-0.15,bright_star_dec[star]+0.15,'b-mag %s'%(bright_star_brightness[star]))]
        bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
        bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
    for star in range(0,len(faint_star_ra)):
        faint_star_markers += [ROOT.TMarker(-faint_star_ra[star],faint_star_dec[star],30)]
        faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
        faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
        faint_star_labels += [ROOT.TLatex(-faint_star_ra[star]-0.15,faint_star_dec[star]+0.15,'b-mag %s'%(faint_star_brightness[star]))]
        faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
        faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)
    for star in range(0,len(other_stars)):
        gal_l, gal_b = ConvertRaDecToGalactic(other_star_coord[star][0],other_star_coord[star][1])
        if pow(source_l-gal_l,2)+pow(source_b-gal_b,2)>star_range*star_range: continue
        other_star_markers_gal += [ROOT.TMarker(-gal_l,gal_b,2)]
        other_star_labels_gal += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,other_stars[star])]
        other_star_markers_gal[len(other_star_markers_gal)-1].SetMarkerSize(1.5)
        other_star_labels_gal[len(other_star_labels_gal)-1].SetTextSize(0.03)
        other_star_names_gal += [other_stars[star]]
    for star in range(0,len(bright_star_ra)):
        gal_l, gal_b = ConvertRaDecToGalactic(bright_star_ra[star],bright_star_dec[star])
        bright_star_markers_gal += [ROOT.TMarker(-gal_l,gal_b,30)]
        bright_star_markers_gal[len(bright_star_markers_gal)-1].SetMarkerColor(2)
        bright_star_markers_gal[len(bright_star_markers_gal)-1].SetMarkerSize(1.5)
        bright_star_labels_gal += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(bright_star_brightness[star]))]
        bright_star_labels_gal[len(bright_star_labels_gal)-1].SetTextColor(2)
        bright_star_labels_gal[len(bright_star_labels_gal)-1].SetTextSize(0.02)
    for star in range(0,len(faint_star_ra)):
        gal_l, gal_b = ConvertRaDecToGalactic(faint_star_ra[star],faint_star_dec[star])
        faint_star_markers_gal += [ROOT.TMarker(-gal_l,gal_b,30)]
        faint_star_markers_gal[len(faint_star_markers_gal)-1].SetMarkerColor(3)
        faint_star_markers_gal[len(faint_star_markers_gal)-1].SetMarkerSize(1.5)
        faint_star_labels_gal += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(faint_star_brightness[star]))]
        faint_star_labels_gal[len(faint_star_labels_gal)-1].SetTextColor(3)
        faint_star_labels_gal[len(faint_star_labels_gal)-1].SetTextSize(0.02)

    for star in range(0,len(other_star_markers)):
        other_star_markers[star].SetMarkerSize(1.5*zoomin_scale)
        other_star_labels[star].SetTextSize(0.03)
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].SetMarkerSize(1.5*zoomin_scale)
        bright_star_labels[star].SetTextSize(0.03)
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].SetMarkerSize(1.5*zoomin_scale)
        faint_star_labels[star].SetTextSize(0.03)

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 650, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.10)
    pad1.SetRightMargin(0.20)
    pad1.SetLeftMargin(0.10)
    pad1.SetTopMargin(0.10)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()

    title_offset = 1.4
    max_nbins = 60

    MapEdge_left = hist_data[0].GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_data[0].GetXaxis().GetBinLowEdge(hist_data[0].GetNbinsX()+1)
    MapEdge_lower = hist_data[0].GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_data[0].GetYaxis().GetBinLowEdge(hist_data[0].GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.

    hist_zscore_skymap = []
    hist_excess_skymap = []
    hist_cali_skymap = []
    hist_data_skymap = []
    hist_flux_skymap = []
    hist_energy_flux_skymap = []
    hist_energy_flux_syst_skymap = []
    hist_energy_flux_normsyst_skymap = []
    hist_syst_skymap = []
    hist_normsyst_skymap = []
    hist_bkgd_skymap = []
    hist_rfov_skymap = []
    hist_effarea_skymap = []
    hist_expo_skymap = []
    hist_expo_hour_skymap = []
    hist_expo_rebin_skymap = []
    for ebin in range(0,len(energy_bin)-1):
        hist_flux_skymap += [ROOT.TH2D("hist_flux_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_energy_flux_skymap += [ROOT.TH2D("hist_energy_flux_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_energy_flux_syst_skymap += [ROOT.TH2D("hist_energy_flux_syst_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_energy_flux_normsyst_skymap += [ROOT.TH2D("hist_energy_flux_normsyst_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_zscore_skymap += [ROOT.TH2D("hist_zscore_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_excess_skymap += [ROOT.TH2D("hist_excess_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_cali_skymap += [ROOT.TH2D("hist_cali_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_data_skymap += [ROOT.TH2D("hist_data_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_syst_skymap += [ROOT.TH2D("hist_syst_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_normsyst_skymap += [ROOT.TH2D("hist_normsyst_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_bkgd_skymap += [ROOT.TH2D("hist_bkgd_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_rfov_skymap += [ROOT.TH2D("hist_rfov_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_effarea_skymap += [ROOT.TH2D("hist_effarea_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_expo_skymap += [ROOT.TH2D("hist_expo_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_expo_hour_skymap += [ROOT.TH2D("hist_expo_hour_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_expo_rebin_skymap += [ROOT.TH2D("hist_expo_rebin_skymap_%s"%(ebin),"",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
    for ebin in range(0,len(energy_bin)-1):
        for bx in range(0,hist_data[0].GetNbinsX()):
            for by in range(0,hist_data[0].GetNbinsY()):
                bx_center = hist_data[0].GetXaxis().GetBinCenter(bx+1)
                by_center = hist_data[0].GetYaxis().GetBinCenter(by+1)
                bx2 = hist_syst_skymap[ebin].GetXaxis().FindBin(bx_center)
                by2 = hist_syst_skymap[ebin].GetYaxis().FindBin(by_center)
                new_content = hist_syst[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_syst_skymap[ebin].GetBinContent(bx2,by2)
                hist_syst_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_normsyst[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_normsyst_skymap[ebin].GetBinContent(bx2,by2)
                hist_normsyst_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_data[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_data_skymap[ebin].GetBinContent(bx2,by2)
                hist_data_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_error = hist_data[ebin].GetBinError(bx+1,by+1)
                old_error = hist_data_skymap[ebin].GetBinError(bx2,by2)
                hist_data_skymap[ebin].SetBinError(bx2,by2,pow(new_error*new_error+old_error*old_error,0.5))
                new_content = hist_bkgd[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_bkgd_skymap[ebin].GetBinContent(bx2,by2)
                hist_bkgd_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_error = hist_bkgd[ebin].GetBinError(bx+1,by+1)
                old_error = hist_bkgd_skymap[ebin].GetBinError(bx2,by2)
                hist_bkgd_skymap[ebin].SetBinError(bx2,by2,pow(new_error*new_error+old_error*old_error,0.5))
                new_content = hist_rfov[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_rfov_skymap[ebin].GetBinContent(bx2,by2)
                hist_rfov_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_error = hist_rfov[ebin].GetBinError(bx+1,by+1)
                old_error = hist_rfov_skymap[ebin].GetBinError(bx2,by2)
                hist_rfov_skymap[ebin].SetBinError(bx2,by2,pow(new_error*new_error+old_error*old_error,0.5))
                new_content = hist_effarea[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_effarea_skymap[ebin].GetBinContent(bx2,by2)
                hist_effarea_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_expo[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_expo_skymap[ebin].GetBinContent(bx2,by2)
                hist_expo_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_error = hist_expo[ebin].GetBinError(bx+1,by+1)
                old_error = hist_expo_skymap[ebin].GetBinError(bx2,by2)
                hist_expo_skymap[ebin].SetBinError(bx2,by2,pow(new_error*new_error+old_error*old_error,0.5))
    for ebin in range(0,len(energy_bin)-1):
        if not Smoothing:
            continue
        hist_bkgd_skymap[ebin] = Smooth2DMap(hist_bkgd_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_rfov_skymap[ebin] = Smooth2DMap(hist_rfov_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_data_skymap[ebin] = Smooth2DMap(hist_data_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_effarea_skymap[ebin] = Smooth2DMap(hist_effarea_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_expo_skymap[ebin] = Smooth2DMap(hist_expo_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_syst_skymap[ebin] = Smooth2DMap(hist_syst_skymap[ebin],smooth_size_spectroscopy,False,True)
        hist_normsyst_skymap[ebin] = Smooth2DMap(hist_normsyst_skymap[ebin],smooth_size_spectroscopy,False,True)
    for ebin in range(0,len(energy_bin)-1):
        expo_rebin = 1
        #if ebin==3:
        #    expo_rebin = 2
        #if ebin==4:
        #    expo_rebin = 4
        #if ebin==5:
        #    expo_rebin = 5
        hist_expo_rebin_skymap[ebin].Add(hist_expo_skymap[ebin])
        hist_expo_rebin_skymap[ebin].Rebin2D(expo_rebin,expo_rebin)
        hist_expo_rebin_skymap[ebin].Scale(1./(expo_rebin*expo_rebin))
        for bx in range(0,hist_expo_skymap[ebin].GetNbinsX()):
            for by in range(0,hist_expo_skymap[ebin].GetNbinsY()):
                cell_x = hist_expo_skymap[ebin].GetXaxis().GetBinCenter(bx+1)
                cell_y = hist_expo_skymap[ebin].GetYaxis().GetBinCenter(by+1)
                rebin_x = hist_expo_rebin_skymap[ebin].GetXaxis().FindBin(cell_x)
                rebin_y = hist_expo_rebin_skymap[ebin].GetYaxis().FindBin(cell_y)
                expo_content = hist_expo_rebin_skymap[ebin].GetBinContent(rebin_x,rebin_y)
                hist_expo_skymap[ebin].SetBinContent(bx+1,by+1,expo_content)

    flux_index = 2.75
    energy_index = energy_index_scale
    skymap_bin_size_x = hist_data_skymap[0].GetXaxis().GetBinCenter(2)-hist_data_skymap[0].GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = hist_data_skymap[0].GetYaxis().GetBinCenter(2)-hist_data_skymap[0].GetYaxis().GetBinCenter(1)
    for ebin in range(0,len(energy_bin)-1):
        hist_zscore_skymap[ebin] = GetSignificanceMap(hist_data_skymap[ebin],hist_bkgd_skymap[ebin],hist_syst_skymap[ebin],False)
        hist_excess_skymap[ebin].Add(hist_data_skymap[ebin])
        hist_excess_skymap[ebin].Add(hist_bkgd_skymap[ebin],-1)
        data_sum = hist_data_skymap[ebin].Integral()
        for bx in range(0,hist_data_skymap[0].GetNbinsX()):
            for by in range(0,hist_data_skymap[0].GetNbinsY()):
                data_content = hist_data_skymap[ebin].GetBinContent(bx+1,by+1)
                bkgd_content = hist_bkgd_skymap[ebin].GetBinContent(bx+1,by+1)
                if data_content==0.: continue
                if bkgd_content==0.: continue
                cell_x = hist_data_skymap[ebin].GetXaxis().GetBinCenter(bx+1)
                cell_y = hist_data_skymap[ebin].GetYaxis().GetBinCenter(by+1)
                rebin_x = hist_expo_skymap[ebin].GetXaxis().FindBin(cell_x)
                rebin_y = hist_expo_skymap[ebin].GetYaxis().FindBin(cell_y)
                expo_content_max = hist_expo_skymap[ebin].GetMaximum()
                expo_content = hist_expo_skymap[ebin].GetBinContent(rebin_x,rebin_y)
                expo_content = max(expo_content,0.5*expo_content_max)

                #map_x = source_ra
                #map_y = source_dec
                map_x = hist_data_skymap[ebin].GetXaxis().GetBinCenter(bx+1)
                map_y = hist_data_skymap[ebin].GetYaxis().GetBinCenter(by+1)
                flux_calibration = GetFluxCalibration(map_x,map_y,ebin)
                hist_cali_skymap[ebin].SetBinContent(bx+1,by+1,flux_calibration)
                correction = flux_calibration*(skymap_bin_size_x*skymap_bin_size_y)/(3.14*calibration_radius*calibration_radius)
                stat_data_err = hist_data_skymap[ebin].GetBinError(bx+1,by+1)/expo_content*correction
                stat_bkgd_err = hist_bkgd_skymap[ebin].GetBinError(bx+1,by+1)/expo_content*correction
                stat_err = pow(stat_data_err*stat_data_err+stat_bkgd_err*stat_bkgd_err,0.5)
                syst_err = hist_syst_skymap[ebin].GetBinContent(bx+1,by+1)/expo_content*correction
                syst_err = max(syst_err,1e-4*correction)
                normsyst_err = hist_normsyst_skymap[ebin].GetBinContent(bx+1,by+1)/expo_content*correction
                normsyst_err = max(normsyst_err,1e-4*correction)
                flux_content = (data_content-bkgd_content)/expo_content*correction

                #map_x = hist_data_skymap[ebin].GetXaxis().GetBinCenter(bx+1)
                #map_y = hist_data_skymap[ebin].GetYaxis().GetBinCenter(by+1)
                #correction = GetEffectiveAreaWeight(map_x,map_y,hist_effarea_skymap[ebin],hist_expo_skymap[ebin],exposure_hours,ebin)
                #stat_data_err = hist_data_skymap[ebin].GetBinError(bx+1,by+1)*correction
                #stat_bkgd_err = hist_bkgd_skymap[ebin].GetBinError(bx+1,by+1)*correction
                #stat_err = pow(stat_data_err*stat_data_err+stat_bkgd_err*stat_bkgd_err,0.5)
                #syst_err = hist_syst_skymap[ebin].GetBinContent(bx+1,by+1)*correction
                #syst_err = max(syst_err,1e-4*correction)
                #normsyst_err = hist_normsyst_skymap[ebin].GetBinContent(bx+1,by+1)*correction
                #normsyst_err = max(normsyst_err,1e-4*correction)
                #flux_content = (data_content-bkgd_content)*correction

                hist_flux_skymap[ebin].SetBinContent(bx+1,by+1,flux_content)
                hist_flux_skymap[ebin].SetBinError(bx+1,by+1,stat_err)
                hist_energy_flux_skymap[ebin].SetBinContent(bx+1,by+1,flux_content*pow(energy_bin[ebin]/1e3,energy_index))
                hist_energy_flux_skymap[ebin].SetBinError(bx+1,by+1,stat_err*pow(energy_bin[ebin]/1e3,energy_index))
                hist_energy_flux_syst_skymap[ebin].SetBinContent(bx+1,by+1,syst_err*pow(energy_bin[ebin]/1e3,energy_index))
                hist_energy_flux_normsyst_skymap[ebin].SetBinContent(bx+1,by+1,normsyst_err*pow(energy_bin[ebin]/1e3,energy_index))

    
    list_eaxis = []
    list_fdata = []
    list_error = []
    list_zeros = []
    list_uplim = []
    legends = []
    for nth_roi in range(1,len(roi_ra)):
        if (roi_name[nth_roi] in exclude_roi): continue
        legends += ['%s'%(roi_name[nth_roi])]
    list_eaxis, list_fdata, list_stat_error = GetRegionIntegral(hist_energy_flux_skymap)
    list_eaxis, list_syst_error, list_syst_empty = GetRegionIntegral(hist_energy_flux_syst_skymap)
    for nth_roi in range(0,len(list_eaxis)):
        zeros = []
        uplim = []
        total_error = []
        for ebin in range(0,len(list_eaxis[nth_roi])):
            #total_error += [pow(pow(list_stat_error[nth_roi][ebin],2)+pow(list_syst_error[nth_roi][ebin],2),0.5)]
            total_error += [list_stat_error[nth_roi][ebin]]
            zeros += [0.]
            uplim += [0]
        list_uplim += [uplim]
        list_zeros += [zeros]
        list_error += [total_error]
    print ('list_fdata = %s'%(list_fdata))
    print ('list_error = %s'%(list_error))
    print ('list_uplim = %s'%(list_uplim))
    print ('list_eaxis = %s'%(list_eaxis))

    cycol = cycle('krgbcmy')
    fig.clf()
    axbig = fig.add_subplot()
    for nth_roi in range(0,len(list_fdata)):
        n_1sigma = 0
        eaxis_fit = []
        fdata_fit = []
        error_fit = []
        for binx in range(0,len(list_fdata[nth_roi])):
            if list_error[nth_roi][binx]==0.: continue
            if list_fdata[nth_roi][binx]/list_error[nth_roi][binx]>1.:
                n_1sigma += 1
                eaxis_fit += [list_eaxis[nth_roi][binx]]
                fdata_fit += [list_fdata[nth_roi][binx]]
                error_fit += [list_error[nth_roi][binx]]
            else:
                list_uplim[nth_roi][binx] = 1
                list_fdata[nth_roi][binx] = 1.*list_error[nth_roi][binx]
                list_error[nth_roi][binx] = 0.3*list_fdata[nth_roi][binx]
        uplims = np.array(list_uplim[nth_roi], dtype=bool)
        next_color = next(cycol)
        doFit = True
        if 'Crab' in legends[nth_roi]:
            doFit = False
        if 'Geminga' in legends[nth_roi]:
            doFit = False
        if 'HAWC' in legends[nth_roi]:
            doFit = False
        if n_1sigma<3 or doUpperLimit:
            doFit = False
        if not doFit:
            axbig.errorbar(list_eaxis[nth_roi],list_fdata[nth_roi],list_error[nth_roi],color=next_color,marker='_',uplims=uplims,ls='none',label='%s'%(legends[nth_roi]))
            print ('=================================================')
            print ('Region %s'%(legends[nth_roi]))
            print ('list_fdata[nth_roi] = %s'%(list_fdata[nth_roi]))
            print ('list_error[nth_roi] = %s'%(list_error[nth_roi]))
        else:
            start = (fdata_fit[0]/pow(10,-12), -2.)
            popt, pcov = curve_fit(power_law_func,np.array(eaxis_fit),np.array(fdata_fit),p0=start,sigma=np.array(error_fit))
            axbig.plot(np.array(list_eaxis[nth_roi]),power_law_func(np.array(list_eaxis[nth_roi]),*popt),color=next_color)
            flux_fit = power_law_func(np.array(eaxis_fit), *popt)
            residual = np.array(fdata_fit) - flux_fit
            chisq = np.sum((residual/np.array(error_fit))**2)
            dof = len(eaxis_fit)-2
            axbig.errorbar(list_eaxis[nth_roi],list_fdata[nth_roi],list_error[nth_roi],color=next_color,marker='_',uplims=uplims,ls='none',label='%s, $\Gamma$ = %0.1f, $\chi^{2}/dof = %0.1f$'%(legends[nth_roi],popt[1]-energy_index,chisq/dof))
        if doUpperLimit:
            axbig.fill_between(list_eaxis[nth_roi], list_zeros[nth_roi], list_error[nth_roi], alpha=0.2, color='r')

    #if doReferenceFlux and len(list_eaxis[0])>=1:
    #    log_energy = np.linspace(log10(list_eaxis[0][0]),log10(list_eaxis[0][len(list_eaxis[0])-1]),50)
    #    xdata = pow(10.,log_energy)
    #    for nth_roi in range(0,len(list_fdata)):
    #        if 'Crab' in legends[nth_roi]:
    #            if not doReferenceFlux: 
    #                continue
    #            vectorize_f = np.vectorize(flux_crab_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='1508.06442')
    #            xdata_array = []
    #            for binx in range(0,len(list_fdata[nth_roi])):
    #                xdata_array += [list_eaxis[nth_roi][binx]]
    #            ydata = pow(np.array(xdata_array)/1e3,energy_index)*vectorize_f(xdata_array)
    #            calibration_new = []
    #            for binx in range(0,len(list_fdata[nth_roi])):
    #                if list_fdata[nth_roi][binx]>0.:
    #                    calibration_new += [ydata[binx]/list_fdata[nth_roi][binx]]
    #                else:
    #                    calibration_new += [0.]
    #            print ('=======================================================================')
    #            print ('new flux_calibration = %s'%(calibration_new))
    #            print ('=======================================================================')
    #            doReferenceFlux = False
    #        if 'VHE region' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_veritas_j1908_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='1404.7185 (VERITAS)')
    #        if 'HAWC region' in legends[nth_roi]:
    #            log_energy = np.linspace(log10(1e2),log10(1e5),50)
    #            xdata_ref = pow(10.,log_energy)
    #            vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
    #            ydata_hawc = pow(xdata_ref/1e3,energy_index)*vectorize_f_hawc(xdata_ref)
    #            axbig.plot(xdata_ref, ydata_hawc,'r-',label='1909.08609 (HAWC)')
    #            axbig.fill_between(xdata_ref, ydata_hawc-0.15*ydata_hawc, ydata_hawc+0.15*ydata_hawc, alpha=0.2, color='r')
    #            # HAWC systematic uncertainty, The Astrophysical Journal 881, 134. Fig 13
    #            #vectorize_f = np.vectorize(flux_veritas_j1908_func)
    #            #ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            #axbig.plot(xdata, ydata,'b-',label='1404.7185 (VERITAS)')
    #            #vectorize_f = np.vectorize(flux_hess_j1908_func)
    #            #ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            #axbig.plot(xdata, ydata,'g-',label='0904.3409 (HESS)')
    #            HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1908(energy_index)
    #            Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1908(energy_index)
    #            axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='g',marker='s',ls='none',label='Fermi')
    #            axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC')
    #        if 'IC 443' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_ic443_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='0905.3291')
    #            vectorize_f = np.vectorize(flux_ic443_hawc_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'g-',label='2007.08582 (HAWC)')
    #        if '1ES 1218+304' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_1es1218_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='0810.0301')
    #        if 'Geminga Pulsar' in legends[nth_roi]:
    #            #log_energy = np.linspace(log10(3.2e3),log10(93e3),50)
    #            #xdata_ref = pow(10.,log_energy)
    #            #vectorize_f = np.vectorize(flux_geminga_func)
    #            #ydata = pow(xdata_ref/1e3,energy_index)*vectorize_f(xdata_ref)
    #            #axbig.plot(xdata_ref, ydata,'r-',label='HAWC')
    #            #axbig.fill_between(xdata_ref, ydata-0.17*ydata, ydata+0.34*ydata, alpha=0.2, color='r')
    #            HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetFermiHAWCFluxGeminga(energy_index)
    #            axbig.plot(HAWC_energies, HAWC_fluxes,'g-',label='HAWC')
    #            axbig.fill_between(HAWC_energies, np.array(HAWC_fluxes)-np.array(HAWC_flux_errs), np.array(HAWC_fluxes)+np.array(HAWC_flux_errs), alpha=0.2, color='g')
    #            Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxGeminga(energy_index)
    #            axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='r',marker='_',ls='none',label='Fermi')
    #            Fermi_UL_energies, Fermi_UL_fluxes, Fermi_UL_err = GetFermiUpperLimitFluxGeminga(energy_index)
    #            uplims = np.array([1,1,1,1], dtype=bool)
    #            axbig.errorbar(Fermi_UL_energies,Fermi_UL_fluxes,Fermi_UL_err,color='r',marker='_',ls='none',uplims=uplims)
    #        if 'J1856+0245' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_j1857_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='From Aleksic et al. (2014)')
    #        if 'J2019+407' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_veritas_gamma_cygni_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='From Weinstein et al. (2012)')
    #        if 'Boomerang' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_veritas_boomerang_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'r-',label='2005.13699')
    #        if 'J1929+1745' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_2hwc_J1928_point_func)
    #            hawc_log_energy = np.linspace(log10(2000.),log10(56000.),50)
    #            hawc_xdata = pow(10.,hawc_log_energy)
    #            ydata = pow(hawc_xdata/1e3,energy_index)*vectorize_f(hawc_xdata)
    #            axbig.plot(hawc_xdata, ydata,'r-',label='2HWC J1928+177 (point-like)')
    #        if 'J1843-0338' in legends[nth_roi]:
    #            vectorize_f = np.vectorize(flux_2hwc_J1844_point_func)
    #            hawc_log_energy = np.linspace(log10(2000.),log10(56000.),50)
    #            hawc_xdata = pow(10.,hawc_log_energy)
    #            ydata = pow(hawc_xdata/1e3,energy_index)*vectorize_f(hawc_xdata)
    #            axbig.plot(hawc_xdata, ydata,'r-',label='2HWC J1844-032 (point-like)')
    #            vectorize_f = np.vectorize(flux_2hwc_J1844_extend_func)
    #            hawc_log_energy = np.linspace(log10(2000.),log10(56000.),50)
    #            hawc_xdata = pow(10.,hawc_log_energy)
    #            ydata = pow(hawc_xdata/1e3,energy_index)*vectorize_f(hawc_xdata)
    #            axbig.plot(hawc_xdata, ydata,'m-',label='2HWC J1844-032 (0.6-deg ext.)')
    #            vectorize_f = np.vectorize(flux_hess_J1843_func)
    #            ydata = pow(xdata/1e3,energy_index)*vectorize_f(xdata)
    #            axbig.plot(xdata, ydata,'g-',label='HESS J1843-033')

    #axbig.legend(loc='best')
    #axbig.set_xlabel('Energy [GeV]')
    #axbig.set_ylabel('$E^{%s}$ Flux [$\mathrm{TeV}^{%s}\mathrm{cm}^{-2}\mathrm{s}^{-1}$]'%(energy_index,-1+energy_index))
    #axbig.set_xscale('log')
    #axbig.set_yscale('log')
    #plotname = 'FluxFromMap'
    #fig.savefig("output_plots/%s_%s.png"%(plotname,sys.argv[1]),bbox_inches='tight')
    #axbig.remove()

    # energy inclusive histograms
    hist_data_skymap_sum = ROOT.TH2D("hist_data_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_bkgd_skymap_sum = ROOT.TH2D("hist_bkgd_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_effarea_skymap_sum = ROOT.TH2D("hist_effarea_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_expo_skymap_sum = ROOT.TH2D("hist_expo_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_syst_skymap_sum = ROOT.TH2D("hist_syst_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_rate_skymap_sum = ROOT.TH2D("hist_rate_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_mask_skymap_sum = ROOT.TH2D("hist_mask_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_flux_syst_skymap_sum = ROOT.TH2D("hist_flux_syst_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_expected_all_gamma_ray_skymap = ROOT.TH2D("hist_expected_all_gamma_ray_skymap","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    for ebin in range(0,len(energy_bin)-1):
        for bx in range(0,hist_data_skymap[0].GetNbinsX()):
            for by in range(0,hist_data_skymap[0].GetNbinsY()):
                data_new_content = hist_data_skymap[ebin].GetBinContent(bx+1,by+1)
                data_old_content = hist_data_skymap_sum.GetBinContent(bx+1,by+1)
                hist_data_skymap_sum.SetBinContent(bx+1,by+1,data_new_content+data_old_content)
                data_new_error = hist_data_skymap[ebin].GetBinError(bx+1,by+1)
                data_old_error = hist_data_skymap_sum.GetBinError(bx+1,by+1)
                hist_data_skymap_sum.SetBinError(bx+1,by+1,pow(data_new_error*data_new_error+data_old_error*data_old_error,0.5))
                bkgd_new_content = hist_bkgd_skymap[ebin].GetBinContent(bx+1,by+1)
                bkgd_old_content = hist_bkgd_skymap_sum.GetBinContent(bx+1,by+1)
                hist_bkgd_skymap_sum.SetBinContent(bx+1,by+1,bkgd_new_content+bkgd_old_content)
                bkgd_new_error = hist_bkgd_skymap[ebin].GetBinError(bx+1,by+1)
                bkgd_old_error = hist_bkgd_skymap_sum.GetBinError(bx+1,by+1)
                hist_bkgd_skymap_sum.SetBinError(bx+1,by+1,pow(bkgd_new_error*bkgd_new_error+bkgd_old_error*bkgd_old_error,0.5))
                effarea_new_content = hist_effarea_skymap[ebin].GetBinContent(bx+1,by+1)
                effarea_old_content = hist_effarea_skymap_sum.GetBinContent(bx+1,by+1)
                hist_effarea_skymap_sum.SetBinContent(bx+1,by+1,effarea_new_content+effarea_old_content)
                expo_new_content = hist_expo_skymap[ebin].GetBinContent(bx+1,by+1)
                expo_old_content = hist_expo_skymap_sum.GetBinContent(bx+1,by+1)
                hist_expo_skymap_sum.SetBinContent(bx+1,by+1,expo_new_content+expo_old_content)
                expo_new_error = hist_expo_skymap[ebin].GetBinError(bx+1,by+1)
                expo_old_error = hist_expo_skymap_sum.GetBinError(bx+1,by+1)
                hist_expo_skymap_sum.SetBinError(bx+1,by+1,pow(expo_new_error*expo_new_error+expo_old_error*expo_old_error,0.5))
                syst_new_content = hist_syst_skymap[ebin].GetBinContent(bx+1,by+1)
                syst_old_content = hist_syst_skymap_sum.GetBinContent(bx+1,by+1)
                hist_syst_skymap_sum.SetBinContent(bx+1,by+1,pow(syst_new_content*syst_new_content+syst_old_content*syst_old_content,0.5))
    hist_zscore_skymap_sum = GetSignificanceMap(hist_data_skymap_sum,hist_bkgd_skymap_sum,hist_syst_skymap_sum,False)
    hist_rate_skymap_sum.Add(hist_data_skymap_sum)
    hist_rate_skymap_sum.Add(hist_bkgd_skymap_sum,-1.)
    hist_rate_skymap_sum.Scale(1./exposure_in_hours)
    bin_count = 0.
    avg_bkgd = 0.
    for bx in range(0,hist_data_skymap[0].GetNbinsX()):
        for by in range(0,hist_data_skymap[0].GetNbinsY()):
            data_content = hist_data_skymap_sum.GetBinContent(bx+1,by+1)
            bkgd_content = hist_bkgd_skymap_sum.GetBinContent(bx+1,by+1)
            if data_content!=0.:
                avg_bkgd += bkgd_content
                bin_count += 1.
    if bin_count>0.: avg_bkgd = avg_bkgd/bin_count
    for bx in range(0,hist_data_skymap[0].GetNbinsX()):
        for by in range(0,hist_data_skymap[0].GetNbinsY()):
            data_content = hist_data_skymap_sum.GetBinContent(bx+1,by+1)
            bkgd_content = hist_bkgd_skymap_sum.GetBinContent(bx+1,by+1)
            zscore = hist_zscore_skymap_sum.GetBinContent(bx+1,by+1)
            if bkgd_content<0.2*avg_bkgd:
                hist_rate_skymap_sum.SetBinContent(bx+1,by+1,0.)
            if zscore<3.0:
                hist_rate_skymap_sum.SetBinContent(bx+1,by+1,0.)

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for bx in range(0,hist_data_skymap[0].GetNbinsX()):
            for by in range(0,hist_data_skymap[0].GetNbinsY()):
                flux_content = hist_energy_flux_skymap[ebin].GetBinContent(bx+1,by+1)
                flux_error = hist_energy_flux_skymap[ebin].GetBinError(bx+1,by+1)
                old_content = hist_flux_skymap_sum.GetBinContent(bx+1,by+1)
                old_error = hist_flux_skymap_sum.GetBinError(bx+1,by+1)
                hist_flux_skymap_sum.SetBinContent(bx+1,by+1,old_content+flux_content)
                hist_flux_skymap_sum.SetBinError(bx+1,by+1,pow(old_error*old_error+flux_error*flux_error,0.5))
                syst_content = hist_energy_flux_syst_skymap[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_flux_syst_skymap_sum.GetBinContent(bx+1,by+1)
                #hist_flux_syst_skymap_sum.SetBinContent(bx+1,by+1,pow(old_content*old_content+syst_content*syst_content,0.5))
                hist_flux_syst_skymap_sum.SetBinContent(bx+1,by+1,old_content+syst_content)
    for bx in range(0,hist_flux_skymap_sum.GetNbinsX()):
        for by in range(0,hist_flux_skymap_sum.GetNbinsY()):
            cell_x = hist_flux_skymap_sum.GetXaxis().GetBinCenter(bx+1)
            cell_y = hist_flux_skymap_sum.GetYaxis().GetBinCenter(by+1)
            distance_sq = pow(cell_x-MapCenter_x,2) + pow(cell_y-MapCenter_y,2)
            #if distance_sq>distance_sq_limit: 
            #    hist_flux_skymap_sum.SetBinContent(bx+1,by+1,0.)
            #    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            #        hist_energy_flux_skymap[ebin].SetBinContent(bx+1,by+1,0.)

    hist_contour = hist_zscore_skymap_sum.Clone()
    hist_zscore_skymap_reflect = reflectXaxis(hist_zscore_skymap_sum)
    hist_contour_reflect = reflectXaxis(hist_contour)
    hist_contour_reflect.SetContour(3)
    hist_contour_reflect.SetContourLevel(0,3)
    hist_contour_reflect.SetContourLevel(1,4)
    hist_contour_reflect.SetContourLevel(2,5)
    hist_zscore_skymap_reflect.GetYaxis().SetTitle(title_y)
    hist_zscore_skymap_reflect.GetXaxis().SetTitle(title_x)
    hist_zscore_skymap_reflect.GetZaxis().SetTitle('Significance')
    hist_zscore_skymap_reflect.GetZaxis().SetTitleOffset(title_offset)
    hist_zscore_skymap_reflect.SetMaximum(5)
    hist_zscore_skymap_reflect.SetMinimum(-5)
    hist_zscore_skymap_reflect.Draw("COL4Z")
    hist_contour_reflect.Draw("CONT3 same")
    hist_zscore_skymap_reflect.GetXaxis().SetLabelOffset(999)
    hist_zscore_skymap_reflect.GetXaxis().SetTickLength(0)
    x1 = hist_zscore_skymap_reflect.GetXaxis().GetXmin()
    x2 = hist_zscore_skymap_reflect.GetXaxis().GetXmax()
    y1 = hist_zscore_skymap_reflect.GetYaxis().GetXmin()
    y2 = hist_zscore_skymap_reflect.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 505, "+")
    raLowerAxis.SetLabelSize(hist_zscore_skymap_reflect.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        mycircles[nth_roi].SetLineWidth(2)
        if nth_roi==0: continue
        if (roi_name[nth_roi] in exclude_roi): continue
        mycircles[nth_roi].Draw("same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    #for star in range(0,len(faint_star_markers)):
    #    faint_star_markers[star].Draw("same")
    #    faint_star_labels[star].Draw("same")
    print ('faint_star_ra = %s'%(faint_star_ra))
    print ('faint_star_dec = %s'%(faint_star_dec))
    canvas.SaveAs('output_plots/SkymapZscore_%s_E%sto%s.png'%(sys.argv[1],energy_bin_cut_low,energy_bin_cut_up))

    hist_rel_syst_skymap_sum = hist_syst_skymap_sum.Clone()
    hist_rel_syst_skymap_sum.Divide(hist_bkgd_skymap_sum)
    hist_rel_syst_skymap_sum_reflect = reflectXaxis(hist_rel_syst_skymap_sum)
    hist_rel_syst_skymap_sum_reflect.Draw("COL4Z")
    hist_rel_syst_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
    hist_rel_syst_skymap_sum_reflect.GetXaxis().SetTickLength(0)
    raLowerAxis.Draw()
    Hist_Data_Elev_Skymap.Divide(Hist_Data_Skymap)
    Hist_Data_Azim_Skymap.Divide(Hist_Data_Skymap)
    Hist_Data_NSB_Skymap.Divide(Hist_Data_Skymap)
    Hist_Data_MJD_Skymap.Divide(Hist_Data_Skymap)
    canvas.SaveAs('output_plots/SkymapRelSyst_%s_%s.png'%(name,selection_tag))
    hist_elev_skymap = Hist_Data_Elev_Skymap.Clone()
    hist_elev_skymap_reflect = reflectXaxis(hist_elev_skymap)
    hist_elev_skymap_reflect.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapElev_%s_%s.png'%(name,selection_tag))
    hist_azim_skymap = Hist_Data_Azim_Skymap.Clone()
    hist_azim_skymap_reflect = reflectXaxis(hist_azim_skymap)
    hist_azim_skymap_reflect.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapAzim_%s_%s.png'%(name,selection_tag))
    hist_nsb_skymap = Hist_Data_NSB_Skymap.Clone()
    hist_nsb_skymap_reflect = reflectXaxis(hist_nsb_skymap)
    hist_nsb_skymap_reflect.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapNSB_%s_%s.png'%(name,selection_tag))
    hist_mjd_skymap = Hist_Data_MJD_Skymap.Clone()
    hist_mjd_skymap_reflect = reflectXaxis(hist_mjd_skymap)
    hist_mjd_skymap_reflect.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapMJD_%s_%s.png'%(name,selection_tag))
    Hist_Data_ElevNSB.Draw("COL4Z")
    canvas.SaveAs('output_plots/ElevNSB_Data_%s_%s.png'%(name,selection_tag))
    Hist_Data_ElevAzim.Draw("COL4Z")
    canvas.SaveAs('output_plots/ElevAzim_Data_%s_%s.png'%(name,selection_tag))
    Hist_Dark_ElevNSB.Draw("COL4Z")
    canvas.SaveAs('output_plots/ElevNSB_Dark_%s_%s.png'%(name,selection_tag))
    Hist_Dark_ElevAzim.Draw("COL4Z")
    canvas.SaveAs('output_plots/ElevAzim_Dark_%s_%s.png'%(name,selection_tag))

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_expo_hour_skymap[ebin].Reset()
        hist_expo_hour_skymap[ebin].Add(hist_bkgd_skymap[ebin])
        bin_size_0 = hist_bkgd_skymap[ebin].GetXaxis().GetBinLowEdge(2)-hist_bkgd_skymap[ebin].GetXaxis().GetBinLowEdge(1)
        bin_size_1 = Hist_Data_CR_XYoff[ebin].GetXaxis().GetBinLowEdge(2)-Hist_Data_CR_XYoff[ebin].GetXaxis().GetBinLowEdge(1)
        bkgd_rate = Hist_Data_CR_XYoff[ebin].GetMaximum()*(bin_size_0*bin_size_0)/(exposure_hours*bin_size_1*bin_size_1)
        hist_expo_hour_skymap[ebin].Scale(1./bkgd_rate)

    hist_expo_hour_skymap_sum = ROOT.TH2D("hist_expo_hour_skymap_sum","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_expo_hour_skymap_sum.Add(hist_bkgd_skymap_sum)
    bin_size_0 = hist_bkgd_skymap_sum.GetXaxis().GetBinLowEdge(2)-hist_bkgd_skymap_sum.GetXaxis().GetBinLowEdge(1)
    bin_size_1 = Hist_OnData_CR_XYoff_Sum.GetXaxis().GetBinLowEdge(2)-Hist_OnData_CR_XYoff_Sum.GetXaxis().GetBinLowEdge(1)
    bkgd_rate = Hist_OnData_CR_XYoff_Sum.GetMaximum()*(bin_size_0*bin_size_0)/(exposure_hours*bin_size_1*bin_size_1)
    hist_expo_hour_skymap_sum.Scale(1./bkgd_rate)
    hist_expo_hour_skymap_sum_reflect = reflectXaxis(hist_expo_hour_skymap_sum)
    hist_expo_hour_skymap_sum_reflect.GetYaxis().SetTitle(title_y)
    hist_expo_hour_skymap_sum_reflect.GetXaxis().SetTitle(title_x)
    hist_expo_hour_skymap_sum_reflect.GetZaxis().SetTitle('exposure [hour]')
    hist_expo_hour_skymap_sum_reflect.GetZaxis().SetTitleOffset(title_offset)
    hist_expo_hour_skymap_sum_reflect.Draw("COL4Z")
    #hist_contour_reflect.Draw("CONT3 same")
    hist_expo_hour_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
    hist_expo_hour_skymap_sum_reflect.GetXaxis().SetTickLength(0)
    raLowerAxis.Draw()
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    canvas.SaveAs('output_plots/SkymapExpoHour_%s_%s.png'%(name,selection_tag))

    hist_rate_skymap_sum_reflect = reflectXaxis(hist_rate_skymap_sum)
    hist_rate_skymap_sum_reflect.GetYaxis().SetTitle(title_y)
    hist_rate_skymap_sum_reflect.GetXaxis().SetTitle(title_x)
    hist_rate_skymap_sum_reflect.GetZaxis().SetTitle('event rate [ evt / hour / %0.1f-deg ]'%(smooth_size_spectroscopy))
    hist_rate_skymap_sum_reflect.GetZaxis().SetTitleOffset(title_offset)
    hist_rate_skymap_sum_reflect.Draw("COL4Z")
    hist_contour_reflect.Draw("CONT3 same")
    hist_rate_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
    hist_rate_skymap_sum_reflect.GetXaxis().SetTickLength(0)
    raLowerAxis.Draw()
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        mycircles[nth_roi].SetLineWidth(2)
        if nth_roi==0: continue
        if (roi_name[nth_roi] in exclude_roi): continue
        mycircles[nth_roi].Draw("same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    #for star in range(0,len(faint_star_markers)):
    #    faint_star_markers[star].Draw("same")
    #    faint_star_labels[star].Draw("same")
    canvas.SaveAs('output_plots/SkymapRate_%s_%s.png'%(name,selection_tag))
    hist_flux_skymap_sum_reflect = reflectXaxis(hist_flux_skymap_sum)
    hist_flux_skymap_sum_reflect.GetYaxis().SetTitle(title_y)
    hist_flux_skymap_sum_reflect.GetXaxis().SetTitle(title_x)
    hist_flux_skymap_sum_reflect.GetZaxis().SetTitle('Flux (cm^{-2} s^{-1})')
    hist_flux_skymap_sum_reflect.GetZaxis().SetTitleOffset(title_offset)
    hist_flux_skymap_sum_reflect.Draw("COL4Z")
    hist_contour_reflect.Draw("CONT3 same")
    hist_flux_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
    hist_flux_skymap_sum_reflect.GetXaxis().SetTickLength(0)
    raLowerAxis.Draw()
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        mycircles[nth_roi].SetLineWidth(2)
        if nth_roi==0: continue
        if (roi_name[nth_roi] in exclude_roi): continue
        mycircles[nth_roi].Draw("same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    #for star in range(0,len(faint_star_markers)):
    #    faint_star_markers[star].Draw("same")
    #    faint_star_labels[star].Draw("same")
    canvas.SaveAs('output_plots/SkymapFlux_%s_%s.png'%(name,selection_tag))
    CommonPlotFunctions.MatplotlibMap2D(hist_flux_skymap_sum_reflect,None,fig,'RA','Dec','flux','SkymapFlux2_%s_%s.png'%(name,selection_tag))

    hist_zscore_skymap_galactic = ROOT.TH2D("hist_zscore_skymap_galactic","",int(Skymap_nbins_x/zoomin_scale),source_l-MapSize_x/zoomin_scale,source_l+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),source_b-MapSize_y/zoomin_scale,source_b+MapSize_y/zoomin_scale)
    #ConvertRaDecToGalacticMap(hist_zscore_skymap_sum,hist_zscore_skymap_galactic)
    ConvertRaDecToGalacticMap(hist_flux_skymap_sum,hist_zscore_skymap_galactic)
    hist_zscore_skymap_galactic_reflect = reflectXaxis(hist_zscore_skymap_galactic)
    hist_zscore_skymap_galactic_reflect.GetYaxis().SetTitle('gal. latitude')
    hist_zscore_skymap_galactic_reflect.GetXaxis().SetTitle('gal. longitude')
    hist_zscore_skymap_galactic_reflect.GetZaxis().SetTitle('Flux [$\mathrm{TeV}^{1}\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    hist_zscore_skymap_galactic_reflect.GetZaxis().SetTitleOffset(title_offset)
    #hist_zscore_skymap_galactic_reflect.SetMaximum(5)
    #hist_zscore_skymap_galactic_reflect.SetMinimum(-5)
    hist_zscore_skymap_galactic_reflect.Draw("COL4Z")
    hist_zscore_skymap_galactic_reflect.GetXaxis().SetLabelOffset(999)
    hist_zscore_skymap_galactic_reflect.GetXaxis().SetTickLength(0)
    x1 = hist_zscore_skymap_galactic_reflect.GetXaxis().GetXmin()
    x2 = hist_zscore_skymap_galactic_reflect.GetXaxis().GetXmax()
    y1 = hist_zscore_skymap_galactic_reflect.GetYaxis().GetXmin()
    y2 = hist_zscore_skymap_galactic_reflect.GetYaxis().GetXmax()
    IncValues_galactic = ROOT.TF1( "IncValues_galactic", "-x", -x2, -x1 )
    galLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues_galactic", 505, "+")
    galLowerAxis.SetLabelSize(hist_zscore_skymap_galactic_reflect.GetXaxis().GetLabelSize())
    galLowerAxis.Draw()
    for star in range(0,len(other_star_markers_gal)):
        other_star_markers_gal[star].Draw("same")
        other_star_labels_gal[star].SetTextAngle(45)
        other_star_labels_gal[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].SetMarkerSize(1.5*zoomin_scale)
        bright_star_labels[star].SetTextSize(0.03)
    canvas.SaveAs('output_plots/SkymapZscoreGal_%s_E%sto%s.png'%(sys.argv[1],energy_bin_cut_low,energy_bin_cut_up))

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_effarea_skymap[ebin].Divide(hist_expo_skymap[ebin])
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_effarea_skymap_reflect = reflectXaxis(hist_effarea_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_effarea_skymap_reflect,None,fig,'RA','Dec','effective area','SkymapEffArea_%s_E%s.png'%(sys.argv[1],ebin))

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        #hist_zscore_skymap_reflect = reflectXaxis(hist_zscore_skymap[ebin])
        hist_zscore_skymap_reflect = reflectXaxis(hist_excess_skymap[ebin])
        hist_zscore_skymap_reflect.GetYaxis().SetTitle(title_y)
        hist_zscore_skymap_reflect.GetXaxis().SetTitle(title_x)
        hist_zscore_skymap_reflect.GetZaxis().SetTitleOffset(title_offset)
        hist_zscore_skymap_reflect.GetXaxis().SetLabelOffset(999)
        hist_zscore_skymap_reflect.GetXaxis().SetTickLength(0)
        hist_zscore_skymap_reflect.Draw("COL4Z")
        raLowerAxis.Draw()
        mycircles = []
        for nth_roi in range(0,len(roi_ra)):
            mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
            mycircles[nth_roi].SetFillStyle(0)
            mycircles[nth_roi].SetLineColor(2)
            mycircles[nth_roi].SetLineWidth(2)
            if nth_roi==0: continue
            if (roi_name[nth_roi] in exclude_roi): continue
            mycircles[nth_roi].Draw("same")
        for star in range(0,len(other_star_markers)):
            other_star_markers[star].Draw("same")
            other_star_labels[star].Draw("same")
        for star in range(0,len(bright_star_markers)):
            bright_star_markers[star].Draw("same")
            bright_star_labels[star].Draw("same")
        #for star in range(0,len(faint_star_markers)):
        #    faint_star_markers[star].Draw("same")
        #    faint_star_labels[star].Draw("same")
        canvas.SaveAs('output_plots/SkymapExcess_%s_%s_E%s.png'%(name,selection_tag,ebin))

        hist_bkgd_skymap_reflect = reflectXaxis(hist_bkgd_skymap[ebin])
        hist_bkgd_skymap_reflect.GetYaxis().SetTitle(title_y)
        hist_bkgd_skymap_reflect.GetXaxis().SetTitle(title_x)
        hist_bkgd_skymap_reflect.GetZaxis().SetTitleOffset(title_offset)
        hist_bkgd_skymap_reflect.GetXaxis().SetLabelOffset(999)
        hist_bkgd_skymap_reflect.GetXaxis().SetTickLength(0)
        hist_bkgd_skymap_reflect.Draw("COL4Z")
        raLowerAxis.Draw()
        mycircles = []
        for nth_roi in range(0,len(roi_ra)):
            mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
            mycircles[nth_roi].SetFillStyle(0)
            mycircles[nth_roi].SetLineColor(2)
            mycircles[nth_roi].SetLineWidth(2)
            if nth_roi==0: continue
            if (roi_name[nth_roi] in exclude_roi): continue
            mycircles[nth_roi].Draw("same")
        for star in range(0,len(other_star_markers)):
            other_star_markers[star].Draw("same")
            other_star_labels[star].Draw("same")
        for star in range(0,len(bright_star_markers)):
            bright_star_markers[star].Draw("same")
            bright_star_labels[star].Draw("same")
        canvas.SaveAs('output_plots/SkymapBkgd_%s_%s_E%s.png'%(name,selection_tag,ebin))

    if doMWLMap:
        #Hist_MWL_global = GetCOSkymap(Hist_MWL_global, isRaDec)
        #Hist_MWL_global = GetHawcSkymap(Hist_MWL_global, isRaDec)
        #MWL_label = 'HAWC significance map'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv25081611580573_NineYear_INTEGRAL_IBIS_17_35keV_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = '9 year Integral IBIS 17-35 keV'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv25082422748501_SwiftBAT_70_Month_14_195keV_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = 'Swift BAT 70 months 14-195 keV'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv25082676066561_ROSAT_Xray_Band7_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = 'ROSAT X-ray Band 7'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv23774123153055_fermi5_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = 'Fermi 5'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv23774832380255_CO_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = 'CO Galactic Plane Survey'
        #Hist_MWL_global = GetSkyViewMap("MWL_maps/skv25103829773371_Effelsberg_Bonn_HI_J1908.txt", Hist_MWL_global, isRaDec)
        #MWL_label = 'Effelsberg-Bonn HI Survey'
        #MWL_map_file += ['MWL_maps/skv23774832380255_CO_J1908.txt']
        #MWL_map_title += ['CO intensity (K km s^{-1} deg)']
        #MWL_map_file += ['MWL_maps/skv790080564098_fermi_band4.txt']
        #MWL_map_title += ['cnts/s/cm^{2}/sr']
        #MWL_map_file += ['MWL_maps/skv790080564098_fermi_band5.txt']
        #MWL_map_title += ['cnts/s/cm^{2}/sr']

        RA_SN = 286.786
        Dec_SN= 6.498
        d_SN = 7.9 #kpc
        #t_SN = 20.*1000. # lower limit by Downes et al. (1980)
        t_SN = 11.*1000. # age of PSR J1907+0631 in year
        r_SN = 0.43 # Yang et al. 2006, ChJAA, 6, 210.

        M_ej = 1. # ejecta mass in solar mass unit
        n_0 = 3.0 # cm^{-3} ambient density

        if sys.argv[1]=='MGRO_J1908_ON':

            #MWL_map_file = 'MWL_maps/FGN_04000+0000_2x2_12CO_v1.00_cube_53_78_0th_moment.txt' # CO intensity (K km s^{-1} deg)
            ##MWL_map_file = 'MWL_maps/FGN_04000+0000_2x2_12CO_v1.00_cube_72_75_0th_moment.txt'
            #MWL_map_title = 'H_{2} column density (1/cm^{2})'
            #Hist_MWL_global = GetGalacticCoordMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(2.*1e20) # H2 column density
            ##MWL_label = 'FUGIN project'

            #MWL_map_file = 'MWL_maps/skv23774832380255_CO_J1908.txt' # CO intensity (K km s^{-1} deg)
            #MWL_map_title = 'H_{2} column density (1/cm^{2})'
            #Hist_MWL_global = GetSkyViewMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(2.*1e20) # H2 column density

            # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
            # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
            #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_incl_0th_moment.txt' # CO intensity (K km s^{-1} deg)
            if J1908_model=='plerion':
                MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_20_40_0th_moment.txt' # CO intensity (K km s^{-1} deg)
                #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_45_65_0th_moment.txt' # CO intensity (K km s^{-1} deg)
            if J1908_model=='SNR_3p4kpc_PSR_moving' or J1908_model=='SNR_3p4kpc_PSR_stationary':
                MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_45_65_0th_moment.txt' # CO intensity (K km s^{-1} deg)
            if J1908_model=='SNR_7p9kpc_PSR_moving' or J1908_model=='SNR_7p9kpc_PSR_stationary':
                MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_50_75_0th_moment.txt' # CO intensity (K km s^{-1} deg)
            MWL_map_title = 'H_{2} column density (1/cm^{2})'
            Hist_MWL_global = GetGalacticCoordMap(MWL_map_file, Hist_MWL_global, isRaDec)
            Hist_MWL_global.Scale(2.*1e20*1000.) # H2 column density, the source FITS file has a mistake in velocity km/s -> m/s

            #MWL_map_file = 'MWL_maps/VGPS_cont_MOS041_ds9.txt' # CO intensity (K km s^{-1} deg)
            #MWL_map_title = 'H_{2} column density (1/cm^{2})'
            #Hist_MWL_global = GetGalacticCoordMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(1.823*1e18*(25.)) # optical thin

        elif sys.argv[1]=='GammaCygni_ON':

            #MWL_map_file = 'MWL_maps/skv4052116888244_gammaCygni.txt' # CO intensity (K km s^{-1} deg)
            #MWL_map_title = 'H_{2} column density (1/cm^{2})'
            #Hist_MWL_global = GetSkyViewMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(2.*1e20) # H2 column density

            MWL_map_file = 'MWL_maps/CGPS_1420MHz.txt' # CO intensity (K km s^{-1} deg)
            MWL_map_title = 'H_{2} column density (1/cm^{2})'
            Hist_MWL_global = GetGalacticCoordMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(1.823*1e18*(9254.-8059.)/1000.) # optical thin
            Hist_MWL_global.Scale(1.823*1e18) # optical thin

        elif sys.argv[1]=='IC443HotSpot_ON':

            MWL_map_file = 'MWL_maps/CGPS_1420MHz_IC443.txt' # CO intensity (K km s^{-1} deg)
            MWL_map_title = 'H_{2} column density (1/cm^{2})'
            Hist_MWL_global = GetGalacticCoordMap(MWL_map_file, Hist_MWL_global, isRaDec)
            #Hist_MWL_global.Scale(1.823*1e18*(9254.-8059.)/1000.) # optical thin
            Hist_MWL_global.Scale(1.823*1e18) # optical thin

        Hist_MWL_global_reflect = reflectXaxis(Hist_MWL_global)
        Hist_MWL_global_reflect.GetYaxis().SetTitle(title_y)
        Hist_MWL_global_reflect.GetXaxis().SetTitle(title_x)
        Hist_MWL_global_reflect.GetZaxis().SetTitle(MWL_map_title)
        Hist_MWL_global_reflect.GetZaxis().SetTitleOffset(title_offset+0.1)
        Hist_MWL_global_reflect.Draw("COL4Z")
        hist_contour_reflect.Draw("CONT3 same")
        Hist_MWL_global_reflect.GetXaxis().SetLabelOffset(999)
        Hist_MWL_global_reflect.GetXaxis().SetTickLength(0)
        raLowerAxis.Draw()
        mycircles = []
        for nth_roi in range(0,len(roi_ra)):
            mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
            mycircles[nth_roi].SetFillStyle(0)
            mycircles[nth_roi].SetLineColor(2)
            mycircles[nth_roi].SetLineWidth(2)
            if nth_roi==0: continue
            if (roi_name[nth_roi] in exclude_roi): continue
            mycircles[nth_roi].Draw("same")
        for star in range(0,len(other_star_markers)):
            other_star_markers[star].Draw("same")
            other_star_labels[star].Draw("same")
        for star in range(0,len(bright_star_markers)):
            bright_star_markers[star].Draw("same")
            bright_star_labels[star].Draw("same")
        #for star in range(0,len(faint_star_markers)):
        #    faint_star_markers[star].Draw("same")
        #    faint_star_labels[star].Draw("same")
        canvas.SaveAs('output_plots/SkymapMWL_%s_%s.png'%(name,selection_tag))



    # spectral-index map
    hist_index_skymap = ROOT.TH2D("hist_index_skymap","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_chisq_skymap = ROOT.TH2D("hist_chisq_skymap","",int(Skymap_nbins_x/zoomin_scale),MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,int(Skymap_nbins_y/zoomin_scale),MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    #for binx in range(0,hist_index_skymap.GetNbinsX()):
    #    for biny in range(0,hist_index_skymap.GetNbinsY()):
    #        hist_index_skymap.SetBinContent(binx+1,biny+1,-99.)
    #        hist_chisq_skymap.SetBinContent(binx+1,biny+1,-99.)
    #        eaxis = []
    #        fdata = []
    #        error = []
    #        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    #            eaxis += [energy_bin[ebin]]
    #            fdata += [hist_flux_skymap[ebin].GetBinContent(binx+1,biny+1)]
    #            error += [hist_flux_skymap[ebin].GetBinError(binx+1,biny+1)]
    #        n_2sigma = 0.
    #        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    #            data_content = hist_data_skymap[ebin].GetBinContent(binx+1,biny+1)
    #            bkgd_content = hist_bkgd_skymap[ebin].GetBinContent(binx+1,biny+1)
    #            data_error = hist_data_skymap[ebin].GetBinError(binx+1,biny+1)
    #            bkgd_error = hist_bkgd_skymap[ebin].GetBinError(binx+1,biny+1)
    #            total_error = pow(data_error*data_error+bkgd_error*bkgd_error,0.5)
    #            if data_content==0.: continue
    #            if (data_content-bkgd_content)/pow(data_content,0.5)>2.0: n_2sigma += 1
    #        if n_2sigma<3: continue
    #        start = (fdata[0]/pow(10,-12), -2.)
    #        popt, pcov = curve_fit(power_law_func, np.array(eaxis), np.array(fdata), p0=start, sigma=np.array(error))
    #        flux_fit = power_law_func(np.array(eaxis), *popt)
    #        residual = np.array(fdata) - flux_fit
    #        chisq = np.sum((residual/np.array(error))**2)
    #        dof = len(eaxis)-2
    #        hist_index_skymap.SetBinContent(binx+1,biny+1,popt[1])
    #        hist_chisq_skymap.SetBinContent(binx+1,biny+1,chisq/dof)
    #        if chisq/dof<4.: continue
    #        ax.cla()
    #        ax.plot(np.array(eaxis), power_law_func(np.array(eaxis), *popt),color='b')
    #        ax.errorbar(eaxis, fdata, error, color='b', marker='s', ls='none')
    #        ax.set_xscale('log')
    #        ax.set_yscale('log')
    #        plotname = 'Skymap_flux_fit_cell_%s_%s'%(binx,biny)
    #        fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    #        ax.cla()
    #hist_index_skymap_reflect = reflectXaxis(hist_index_skymap)
    #hist_index_skymap_reflect.GetYaxis().SetTitle(title_y)
    #hist_index_skymap_reflect.GetXaxis().SetTitle(title_x)
    #hist_index_skymap_reflect.GetZaxis().SetTitle('spectral index')
    #hist_index_skymap_reflect.GetZaxis().SetTitleOffset(title_offset)
    #hist_index_skymap_reflect.SetMaximum(-2.)
    #hist_index_skymap_reflect.SetMinimum(-3.)
    #hist_index_skymap_reflect.Draw("COL4Z")
    #hist_contour_reflect.Draw("CONT3 same")
    #hist_index_skymap_reflect.GetXaxis().SetLabelOffset(999)
    #hist_index_skymap_reflect.GetXaxis().SetTickLength(0)
    #raLowerAxis.Draw()
    #mycircles = []
    #for nth_roi in range(0,len(roi_ra)):
    #    mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
    #    mycircles[nth_roi].SetFillStyle(0)
    #    mycircles[nth_roi].SetLineColor(2)
    #    mycircles[nth_roi].SetLineWidth(2)
    #    if nth_roi==0: continue
    #    if (roi_name[nth_roi] in exclude_roi): continue
    #    mycircles[nth_roi].Draw("same")
    #for star in range(0,len(other_star_markers)):
    #    other_star_markers[star].Draw("same")
    #    other_star_labels[star].Draw("same")
    #for star in range(0,len(bright_star_markers)):
    #    bright_star_markers[star].Draw("same")
    #    bright_star_labels[star].Draw("same")
    ##for star in range(0,len(faint_star_markers)):
    ##    faint_star_markers[star].Draw("same")
    ##    faint_star_labels[star].Draw("same")
    #canvas.SaveAs('output_plots/SkymapIndex_%s_%s.png'%(name,selection_tag))
    #hist_chisq_skymap_reflect = reflectXaxis(hist_chisq_skymap)
    #hist_chisq_skymap_reflect.GetYaxis().SetTitle(title_y)
    #hist_chisq_skymap_reflect.GetXaxis().SetTitle(title_x)
    #hist_chisq_skymap_reflect.GetZaxis().SetTitle('#chi^{2}/DoF')
    #hist_chisq_skymap_reflect.GetZaxis().SetTitleOffset(title_offset)
    #hist_chisq_skymap_reflect.SetMaximum(4.)
    #hist_chisq_skymap_reflect.SetMinimum(0.)
    #hist_chisq_skymap_reflect.Draw("COL4Z")
    #hist_contour_reflect.Draw("CONT3 same")
    #hist_chisq_skymap_reflect.GetXaxis().SetLabelOffset(999)
    #hist_chisq_skymap_reflect.GetXaxis().SetTickLength(0)
    #raLowerAxis.Draw()
    #mycircles = []
    #for nth_roi in range(0,len(roi_ra)):
    #    mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
    #    mycircles[nth_roi].SetFillStyle(0)
    #    mycircles[nth_roi].SetLineColor(2)
    #    mycircles[nth_roi].SetLineWidth(2)
    #    if nth_roi==0: continue
    #    if (roi_name[nth_roi] in exclude_roi): continue
    #    mycircles[nth_roi].Draw("same")
    #for star in range(0,len(other_star_markers)):
    #    other_star_markers[star].Draw("same")
    #    other_star_labels[star].Draw("same")
    #for star in range(0,len(bright_star_markers)):
    #    bright_star_markers[star].Draw("same")
    #    bright_star_labels[star].Draw("same")
    ##for star in range(0,len(faint_star_markers)):
    ##    faint_star_markers[star].Draw("same")
    ##    faint_star_labels[star].Draw("same")
    #canvas.SaveAs('output_plots/SkymapChisq_%s_%s.png'%(name,selection_tag))


    if doMWLMap:
        total_obs_flux_integral, total_model_flux_integral = FitCloudRegionFlux(hist_data_skymap,hist_bkgd_skymap,hist_expo_skymap,hist_syst_skymap,Hist_MWL_global)

    if not doExtentFit:

        RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()
        profile_center_x, profile_center_y, profile_center_z = RA_PSR, Dec_PSR, d_PSR
        profile, profile_err, theta2 = FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,MapSize_y/zoomin_scale)
        fig.clf()
        axbig = fig.add_subplot()
        axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
        axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
        axbig.set_xlabel('angular distance from source [degree]')
        new_tick_locations = axbig.get_xticks()
        axbig.legend(loc='best')
        plotname = 'ProfileVsTheta2'
        fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
        axbig.remove()


        #if sys.argv[1]=='MGRO_J1908_ON':
        #    for bx in range(0,hist_flux_skymap_sum.GetNbinsX()):
        #        for by in range(0,hist_flux_skymap_sum.GetNbinsY()):
        #            cell_x = hist_flux_skymap_sum.GetXaxis().GetBinCenter(bx+1)
        #            cell_y = hist_flux_skymap_sum.GetYaxis().GetBinCenter(by+1)
        #            distance_sq = pow(cell_x-MapCenter_x,2) + pow(cell_y-MapCenter_y,2)
        #    center_radec = [287.05, 6.39] #3HWC J1908+063
        #    for bx in range(0,hist_flux_skymap_sum.GetNbinsX()):
        #        for by in range(0,hist_flux_skymap_sum.GetNbinsY()):
        #            hist_mask_skymap_sum.SetBinContent(bx+1,by+1,0.)
        #            cell_x = hist_flux_skymap_sum.GetXaxis().GetBinCenter(bx+1)
        #            cell_y = hist_flux_skymap_sum.GetYaxis().GetBinCenter(by+1)
        #            distance_sq = pow(cell_x-center_radec[0],2) + pow(cell_y-center_radec[1],2)
        #            if distance_sq>distance_sq_limit: 
        #                hist_mask_skymap_sum.SetBinContent(bx+1,by+1,1.)
        #    exclusion_radec = []
        #    exclusion_radec += [[286.975, 6.03777777778]] #PSR J1907+0602
        #    exclusion_radec += [[286.63, 7.3821666666666665]] #PSR J1906+0722
        #    exclusion_radec += [[286.786, 6.498]] #G40.5-0.5
        #    exclusion_radec += [[286.89166666666665, 7.133333333333334]] #G41.1-0.3
        #    exclusion_radec += [[287.05, 6.39]] #3HWC J1908+063
        #    exclusion_radec += [[288.404,4.930]] #SS 433 e1
        #    exclusion_radec += [[287.654,5.037]] #SS 433 w1
        #    for region in range(0,len(exclusion_radec)):
        #        for bx in range(0,hist_flux_skymap_sum.GetNbinsX()):
        #            for by in range(0,hist_flux_skymap_sum.GetNbinsY()):
        #                cell_x = hist_flux_skymap_sum.GetXaxis().GetBinCenter(bx+1)
        #                cell_y = hist_flux_skymap_sum.GetYaxis().GetBinCenter(by+1)
        #                distance_sq = pow(cell_x-exclusion_radec[region][0],2) + pow(cell_y-exclusion_radec[region][1],2)
        #                if distance_sq<pow(0.5,2): 
        #                    hist_mask_skymap_sum.SetBinContent(bx+1,by+1,1.)
        #    output_file = ROOT.TFile("output_fitting/J1908_skymap.root","recreate")
        #    hist_mask_skymap_sum.Write()
        #    hist_flux_skymap_sum.Write()
        #    hist_flux_syst_skymap_sum.Write()
        #    hist_data_skymap_sum.Write()
        #    hist_bkgd_skymap_sum.Write()
        #    hist_expo_skymap_sum.Write()
        #    hist_syst_skymap_sum.Write()
        #    for ebin in range(0,len(energy_bin)-1):
        #        hist_energy_flux_skymap[ebin].Write()
        #        hist_energy_flux_syst_skymap[ebin].Write()
        #        hist_energy_flux_normsyst_skymap[ebin].Write()
        #        hist_cali_skymap[ebin].Write()
        #        hist_data_skymap[ebin].Write()
        #        hist_bkgd_skymap[ebin].Write()
        #        hist_expo_skymap[ebin].Write()
        #        hist_syst_skymap[ebin].Write()
        #        hist_normsyst_skymap[ebin].Write()
        #    output_file.Close();

        if int(sys.argv[2])==0 and int(sys.argv[3])==len(energy_bin)-1:
            output_file = ROOT.TFile("output_fitting/%s_skymap_%s.root"%(sys.argv[1],folder_path),"recreate")
            Hist_Data_Elev_Skymap.Write()
            Hist_Data_Azim_Skymap.Write()
            Hist_Data_NSB_Skymap.Write()
            Hist_Data_MJD_Skymap.Write()
            hist_expo_hour_skymap_sum.Write()
            for ebin in range(0,len(energy_bin)-1):
                hist_energy_flux_skymap[ebin].Write()
                hist_energy_flux_syst_skymap[ebin].Write()
                hist_data_skymap[ebin].Write()
                hist_bkgd_skymap[ebin].Write()
                hist_rfov_skymap[ebin].Write()
                hist_expo_skymap[ebin].Write()
                hist_expo_hour_skymap[ebin].Write()
                hist_effarea_skymap[ebin].Write()
            output_file.Close();

    else:

        center_x, center_y = source_ra, source_dec
        d_PSR = 1.
        t_PSR = 1.
        if sys.argv[1]=='MGRO_J1908_ON':
            center_x, center_y = 286.975, 6.03777777778 #PSR J1907+0602
            #center_x, center_y = 286.786, 6.498 #SNR G40.5-0.5
            d_PSR = 3.2 #kpc
            t_PSR = 19.5*1000. #yr
        if sys.argv[1]=='Geminga_ON':
            center_x, center_y = 98.476, 17.770 #Geminga Pulsar
        if sys.argv[1]=='GammaCygni_ON':
            center_x, center_y = 305.208, 40.43 # SNR G78.2+2.1
            #center_x, center_y = 305.020, 40.757 # VER J2019+407
            #center_x, center_y = 305.377, 40.448 # PSR J2021+4026

        current_gamma_energy = energy_bin[energy_bin_cut_low]
        next_gamma_energy = energy_bin[energy_bin_cut_up]

        skymap_fitting_model = 'hybrid'
        if sys.argv[1]=='MGRO_J1908_ON':
            skymap_fitting_model = 'hybrid'
            #skymap_fitting_model = 'lepton'
            #skymap_fitting_model = 'hadron'
        elif sys.argv[1]=='GammaCygni_ON':
            skymap_fitting_model = 'hadron'
        elif sys.argv[1]=='IC443HotSpot_ON':
            skymap_fitting_model = 'hadron'
        else:
            skymap_fitting_model = 'lepton'

        if sys.argv[1]=='MGRO_J1908_ON':

            print ('=============================================')

            RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()
            RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()

            if not skymap_fitting_model=='lepton':
                #RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()
                deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
                hydrogen_to_solar_mass = 8.4144035*1e-58
                skymap_bin_area = 3.14*pow(smooth_size_spectroscopy*deg_to_cm,2)
                Hist_MWL_global_mass_reflect = Hist_MWL_global_reflect.Clone()
                Hist_MWL_global_mass_reflect.Scale(skymap_bin_area/3.14*2.8*hydrogen_to_solar_mass) # in solar mass
                MWL_map_title = 'H_{2} mass'
                Hist_MWL_global_mass_reflect.GetZaxis().SetTitle(MWL_map_title)
                Hist_MWL_global_mass_reflect.Draw("COL4Z")
                hist_contour_reflect.Draw("CONT3 same")
                Hist_MWL_global_mass_reflect.GetXaxis().SetLabelOffset(999)
                Hist_MWL_global_mass_reflect.GetXaxis().SetTickLength(0)
                raLowerAxis.Draw()
                mycircles = []
                for nth_roi in range(0,len(roi_ra)):
                    mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius_outer[nth_roi])]
                    mycircles[nth_roi].SetFillStyle(0)
                    mycircles[nth_roi].SetLineColor(2)
                    mycircles[nth_roi].SetLineWidth(2)
                    if nth_roi==0: continue
                    if (roi_name[nth_roi] in exclude_roi): continue
                    mycircles[nth_roi].Draw("same")
                for star in range(0,len(other_star_markers)):
                    other_star_markers[star].Draw("same")
                    other_star_labels[star].Draw("same")
                canvas.SaveAs('output_plots/SkymapMWL_mass_%s_%s.png'%(name,selection_tag))


    return hist_index_skymap

def GetRegionIntegral(hist_skymap):

    list_eaxis = []
    list_fdata = []
    list_error = []
    for nth_roi in range(1,len(roi_ra)):
        if (roi_name[nth_roi] in exclude_roi): continue
        eaxis = []
        fdata = []
        error = []
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            flux_sum = 0.
            flux_stat_err = 0.
            flux_syst_err = 0.
            for bx in range(0,hist_skymap[0].GetNbinsX()):
                for by in range(0,hist_skymap[0].GetNbinsY()):
                    radius_inner = roi_radius_inner[nth_roi]
                    radius_outer = roi_radius_outer[nth_roi]
                    bin_ra = hist_skymap[0].GetXaxis().GetBinCenter(bx+1)
                    bin_dec = hist_skymap[0].GetYaxis().GetBinCenter(by+1)
                    distance = pow(pow(bin_ra-roi_ra[nth_roi],2) + pow(bin_dec-roi_dec[nth_roi],2),0.5)
                    if sys.argv[1]=='Crab_ON':
                        radius_inner = 0.
                        radius_outer = calibration_radius
                    if distance<radius_inner: continue
                    if distance>radius_outer: continue
                    flux_sum += hist_skymap[ebin].GetBinContent(bx+1,by+1)
                    flux_stat_err += pow(hist_skymap[ebin].GetBinError(bx+1,by+1),2)
            radius_inner = roi_radius_inner[nth_roi]
            radius_outer = roi_radius_outer[nth_roi]
            eaxis += [energy_bin[ebin]]
            fdata += [flux_sum]
            error += [pow(flux_stat_err,0.5)]
        list_eaxis += [eaxis]
        list_fdata += [fdata]
        list_error += [error]

    return list_eaxis, list_fdata, list_error

def FindLocalMaximum(Hist_sig, init_x, init_y):

    max_sig = 0.
    for bx in range(0,Hist_sig.GetNbinsX()):
        for by in range(0,Hist_sig.GetNbinsY()):
            bin_x = Hist_sig.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_sig.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_x-init_x,2)+pow(bin_y-init_y,2),0.5)
            if distance < 0.1:
                if max_sig < Hist_sig.GetBinContent(bx+1,by+1):
                    max_sig = Hist_sig.GetBinContent(bx+1,by+1)
    return max_sig

def FindHistAbsMaxBinXY(hist):

    max_content = 0.
    bin_x = 0.
    bin_y = 0.
    for bx in range(0,hist.GetNbinsX()):
        for by in range(0,hist.GetNbinsY()):
            if abs(hist.GetBinContent(bx+1,by+1))>max_content:
                max_content = abs(hist.GetBinContent(bx+1,by+1))
                bin_x = hist.GetXaxis().GetBinCenter(bx+1)
                bin_y = hist.GetYaxis().GetBinCenter(by+1)
    return bin_x, bin_y

def FindHistMaxBinXY(hist):

    max_content = 0.
    bin_x = 0.
    bin_y = 0.
    for bx in range(0,hist.GetNbinsX()):
        for by in range(0,hist.GetNbinsY()):
            if hist.GetBinContent(bx+1,by+1)>max_content:
                max_content = hist.GetBinContent(bx+1,by+1)
                bin_x = hist.GetXaxis().GetBinCenter(bx+1)
                bin_y = hist.GetYaxis().GetBinCenter(by+1)
    return bin_x, bin_y

def GetExtentionRMS(Hist_data, Hist_bkgd, Hist_exposure, roi_x, roi_y, roi_size):

    Hist_Excess = Hist_data.Clone()
    Hist_Excess.Add(Hist_bkgd,-1.)
    Hist_Excess.Divide(Hist_exposure)

    max_exposure = Hist_exposure.GetMaximum()
    if max_exposure==0.:
        return 0.

    weighted_distance_sq = 0.
    total_weight = 0.
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_x = Hist_Excess.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_Excess.GetYaxis().GetBinCenter(by+1)
            distance_sq = pow(bin_x-roi_x,2)+pow(bin_y-roi_y,2)
            if distance_sq>pow(roi_size,2): continue
            weight = pow(Hist_Excess.GetBinContent(bx+1,by+1),2)
            total_weight += weight
            weighted_distance_sq += weight*distance_sq

    if total_weight==0.: return 0.
    return pow(weighted_distance_sq/total_weight,0.5)

def hybrid_model_2d(x, par):

    flux_integral = 0.

    lep_par = [par[0],par[2],par[3],par[4]]
    flux_integral += leptonic_model_2d(x, lep_par)
    had_par = [par[1],par[2],par[3],par[4]]
    flux_integral += hadronic_model_2d(x, had_par)

    return flux_integral

def leptonic_model_2d(x, par):

    RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    n_segments = 10
    pulsar_initial_x = RA_PSR_initial
    pulsar_initial_y = Dec_PSR_initial
    pulsar_final_x = RA_PSR
    pulsar_final_y = Dec_PSR
    skymap_bin_area = 3.14*pow(smooth_size_spectroscopy*deg_to_cm,2)

    spectral_index = 2.2
    el_brightness = par[0]
    diffusion_coefficient_1TeV = par[1]
    RA_PSR_initial = par[2]
    Dec_PSR_initial = par[3]
 
    distance_sq = pow(x[0]-source_ra,2) + pow(x[1]-source_dec,2)
    if distance_sq>1.5*1.5: return 0.

    flux_integral = 0.
    for segment in range(0,n_segments):
        # segment=0 is the final pulsar position
        diff_time = t_PSR/n_segments*(segment+1)
        mag_field = 6.
        E_cmb = 6.6*1e-4 # eV
        m_e = 0.511*1e6 # eV
        sigma_thomson = 6.65*1e-29 # Thomson cross section in m2
        speed_light = 3.*1e8 # m/s
        U_cmb = 2.6*1e5 # eV/m3
        U_B = 6.24*1e18/(8.*3.14*1e-7)*pow(mag_field*1e-6/1e4,2) # eV/m3

        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            E_e = m_e*pow(energy_bin[ebin]*1e9/E_cmb,0.5) # eV
            gamma_factor = E_e/m_e
            t_cooling = m_e/speed_light/(4./3.*sigma_thomson*gamma_factor*(U_cmb+U_B)) # sec
            energy_loss = exp(-diff_time*365.*24.*60.*60./t_cooling)
            t_cooling = min(diff_time*365.*24.*60.*60.,t_cooling)

            diffusion_coefficient = EstimateLeptonicDiffusionCoefficient(energy_bin[ebin],mag_field,d_PSR,t_PSR,diffusion_coefficient_1TeV,0.5)
            diffusion_radius = pow(4.*diffusion_coefficient*t_cooling,0.5)/deg_to_cm
            diffusion_radius = pow(smooth_size_spectroscopy*smooth_size_spectroscopy+diffusion_radius*diffusion_radius,0.5)
            pulsar_current_x = (pulsar_initial_x-pulsar_final_x)/t_PSR*diff_time+pulsar_final_x
            pulsar_current_y = (pulsar_initial_y-pulsar_final_y)/t_PSR*diff_time+pulsar_final_y
            distance = pow(pow(x[0]-pulsar_current_x,2)+pow(x[1]-pulsar_current_y,2),0.5)
            vol_diff = 4.*3.14/3.*pow(diffusion_radius*deg_to_cm,3)
            attenuation = exp(-pow(distance/diffusion_radius,2))
            line_of_sight = diffusion_radius
            flux_integral += el_brightness*energy_loss*pow(energy_bin[ebin]/energy_bin[energy_bin_cut_low],-spectral_index+1)*attenuation*diffusion_radius/line_of_sight*pow(3.14,0.5)/2.*math.erf(line_of_sight/diffusion_radius)

    return flux_integral

def hadronic_model_2d(x, par):

    CR_efficiency = par[0]
    diffusion_coefficient_1TeV = par[1]
    RA_SN = par[2]
    Dec_SN = par[3]

    distance_sq = pow(x[0]-source_ra,2) + pow(x[1]-source_dec,2)
    if distance_sq>1.5*1.5: return 0.

    bin_ra = Hist_MWL_global.GetXaxis().FindBin(x[0])
    bin_dec = Hist_MWL_global.GetYaxis().FindBin(x[1])

    flux_integral = 0.
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):

        hadronic_emission = CalculateHadronicEmission(Hist_MWL_global,bin_ra,bin_dec,ebin,CR_efficiency,diffusion_coefficient_1TeV,RA_SN,Dec_SN)
        flux_integral += hadronic_emission

    return flux_integral

def EstimateCosmicRayBrightness(hist_flux,roi_x,roi_y,roi_radius,diffusion_radius):

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()
    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    skymap_bin_area = 3.14*pow(smooth_size_spectroscopy*deg_to_cm,2)

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19

    total_flux = 0.
    total_flux_err = 0.
    total_mass = 0.
    for binx in range(0,hist_flux.GetNbinsX()):
        for biny in range(0,hist_flux.GetNbinsY()):
            cell_x = hist_flux.GetXaxis().GetBinCenter(binx+1)
            cell_y = hist_flux.GetYaxis().GetBinCenter(biny+1)
            distance_to_roi = pow(pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2),0.5)
            if distance_to_roi>roi_radius: continue
            total_flux += hist_flux.GetBinContent(binx+1,biny+1)
            total_flux_err += hist_flux.GetBinError(binx+1,biny+1)
            total_mass += Hist_MWL_global.GetBinContent(binx+1,biny+1)*skymap_bin_area

    vol_diff = 4.*3.14/3.*pow(diffusion_radius*deg_to_cm,3)
    distance_to_SNR = pow(pow(RA_SN-roi_x,2)+pow(Dec_SN-roi_y,2),0.5)
    attenuation = exp(-pow(distance_to_SNR/diffusion_radius,2))
    line_of_sight = diffusion_radius
    A_par = total_flux/(f_Gamma*1e-10*(pow(current_gamma_energy/1000.,1-Gamma_CR)-pow(next_gamma_energy/1000.,1-Gamma_CR)))
    ECR_density_local_avg = A_par/(total_mass*pow(d_SN,-2))
    ECR_density_in_SNR = ECR_density_local_avg/(attenuation*diffusion_radius*pow(3.14,0.5)/2.*math.erf(line_of_sight/diffusion_radius))*line_of_sight
    CR_brightness = ECR_density_in_SNR*vol_diff

    CR_brightness = max(0.,CR_brightness)
    CR_brightness_err = CR_brightness*total_flux_err/total_flux

    return CR_brightness, CR_brightness_err

def disk_2d_to_1D(x, B):
    return B
def gaussian_2d_to_1D(x, rms, A):
    #return A * np.exp( -0.5*(x/rms)**2 )
    #return A * 1./(x+0.085*rms) * np.exp(-1.54*pow(x/rms,1.52))
    return A * 1./(x+0.06*rms) * np.exp(-1.*pow(x/rms,2))
def gaussian_disk_2d_to_1D(x, rms, A):
    #return A * 1./(x+0.085*rms) * np.exp(-1.54*pow(x/rms,1.52))
    return A * 1.22/(rms*pow(3.14,3/2)*(x+0.06*rms)) * np.exp(-1.*pow(x/rms,2))
def gaussian_2d(x, y, x0, y0, rms, A, B):
    return A * np.exp( -0.5*((x-x0)/rms)**2 -0.5*((y-y0)/rms)**2) + B
def _gaussian_2d(M, *args):
    # This is the callable that is passed to curve_fit. M is a (2,N) array
    # where N is the total number of data points in Z, which will be ravelled
    # to one dimension.
    # https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
    x, y = M
    arr = np.zeros(x.shape)
    n_param = 5
    for i in range(len(args)//n_param):
       arr += gaussian_2d(x, y, *args[i*n_param:i*n_param+n_param])
    return arr

def CalculateHadronicEmission(hist_column_density_map,bin_ra,bin_dec,e_bin,CR_efficiency,diffusion_coefficient_1TeV,RA_SN_input, Dec_SN_input):

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()
    RA_SN = RA_SN_input
    Dec_SN = Dec_SN_input

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit
    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
    E_p_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit

    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    E_p_min_vis = (1./0.17)*2.*energy_bin[e_bin]  # visible CR energy threshold
    E_p_max_vis = (1./0.17)*2.*energy_bin[e_bin+1]  # visible CR energy threshold
    f_vis = log(E_p_max_vis/E_p_min_vis)/log(E_p_max/1.)

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19
    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    bin_x_size = hist_column_density_map.GetXaxis().GetBinLowEdge(2)-hist_column_density_map.GetXaxis().GetBinLowEdge(1)
    bin_y_size = hist_column_density_map.GetYaxis().GetBinLowEdge(2)-hist_column_density_map.GetYaxis().GetBinLowEdge(1)
    skymap_bin_area = bin_x_size*bin_y_size*deg_to_cm*deg_to_cm

    #t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    #CR_efficiency = 3.14/5.*(1.-pow(t_Sedov/t_SN,0.4))
    #esc_param = 2.0

    binx = hist_column_density_map.GetXaxis().FindBin(bin_ra)
    biny = hist_column_density_map.GetYaxis().FindBin(bin_dec)
    distance = pow(pow(bin_ra-RA_SN,2)+pow(bin_dec-Dec_SN,2),0.5)
    MC_column_density = hist_column_density_map.GetBinContent(binx,biny)
    MC_mass = MC_column_density*skymap_bin_area

    #diffusion_radius = DeriveSupernovaDiffusionLength(esc_param,t_SN,d_SN,n_0,energy_bin[e_bin])
    #diffusion_radius = pow(smooth_size_spectroscopy*smooth_size_spectroscopy+diffusion_radius*diffusion_radius,0.5)
    diffusion_coefficient = diffusion_coefficient_1TeV*pow((1./0.17)*2.*energy_bin[e_bin]/1000.,0.5)
    diffusion_radius = pow(4.*diffusion_coefficient*t_SN*365.*24.*60.*60.,0.5)/deg_to_cm
    vol_diff = 4.*3.14/3.*pow(diffusion_radius*deg_to_cm,3)
    attenuation = exp(-(distance*distance)/pow(diffusion_radius,2))
    #vol_diff = 4.*3.14/3.*pow(diffusion_radius*deg_to_cm,3)
    #attenuation = exp(-(distance*distance+R_esc*R_esc)/pow(diffusion_radius,2))
    #if distance>diffusion_radius:
    #    A_esc = 4.*3.14*R_esc*R_esc
    #    vol_diff = diffusion_radius*A_esc*distance/R_esc*pow(deg_to_cm,3)
    #    attenuation = exp(-(distance*distance+R_esc*R_esc)/pow(diffusion_radius,2))*math.sinh(2.*distance*R_esc/diffusion_radius)

    ECR_density_in_SNR = CR_efficiency*E_SN/vol_diff
    ECR_density_local = ECR_density_in_SNR*attenuation
    line_of_sight = diffusion_radius
    ECR_density_local_avg = ECR_density_in_SNR*attenuation*diffusion_radius/line_of_sight*pow(3.14,0.5)/2.*math.erf(line_of_sight/diffusion_radius)
    A_par = MC_mass*ECR_density_local_avg*pow(d_SN,-2)
    hadronic_emission = A_par*(f_Gamma*1e-10*(pow(energy_bin[e_bin]/1000.,1-Gamma_CR)-pow(energy_bin[e_bin+1]/1000.,1-Gamma_CR)))

    return hadronic_emission

def FitCloudRegionFlux(hist_data_skymap,hist_bkgd_skymap,hist_expo_skymap,hist_syst_skymap,hist_column_density_map):

    # PSR J1907+0602
    PSR_ra = 286.975
    PSR_dec = 6.03777777778
    PSR_radius = 0.2

    # CO north
    #region_ra = 286.786
    #region_dec = 7.1
    #region_radius = 0.2
    # PSR tail 2
    region_ra = 287.2
    region_dec = 6.6
    region_radius = 0.2

    hist_hadronic_map = hist_column_density_map.Clone()
    hist_hadronic_map.Reset()

    energy_array = []
    cloud_mass = []
    psr_flux_integral = []
    psr_flux_integral_err = []
    obs_flux_integral = []
    obs_flux_integral_err = []
    had_model_flux_integral = []
    for ebin in range(0,len(energy_bin)-1):
        energy_array += [energy_bin[ebin]]
        cloud_mass += [0.]
        psr_flux_integral += [0.]
        psr_flux_integral_err += [0.]
        obs_flux_integral += [0.]
        obs_flux_integral_err += [0.]
        had_model_flux_integral += [0.]

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):

        data_sum = 0.
        bkgd_sum = 0.
        expo_sum = 0.
        syst_sum = 0.
        for bx in range(0,hist_column_density_map.GetNbinsX()):
            for by in range(0,hist_column_density_map.GetNbinsY()):

                bin_ra = hist_column_density_map.GetXaxis().GetBinCenter(bx+1)
                bin_dec = hist_column_density_map.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(bin_ra-PSR_ra,2) + pow(bin_dec-PSR_dec,2)
                if distance_sq>PSR_radius*PSR_radius: continue

                data_content = hist_data_skymap[ebin].GetBinContent(bx+1,by+1)
                bkgd_content = hist_bkgd_skymap[ebin].GetBinContent(bx+1,by+1)
                expo_content = hist_expo_skymap[ebin].GetBinContent(bx+1,by+1)
                syst_err = hist_syst_skymap[ebin].GetBinContent(bx+1,by+1)

                data_sum += data_content
                bkgd_sum += bkgd_content
                expo_sum += expo_content
                syst_sum += syst_err

        map_x = region_ra
        map_y = region_dec
        flux_calibration = GetFluxCalibration(map_x,map_y,ebin)
        correction = flux_calibration*pow(region_radius,2)/(calibration_radius*calibration_radius)
        flux_content = 0.
        flux_err = 0.
        if expo_sum>0.:
            flux_content = (data_sum-bkgd_sum)/expo_sum*correction
            flux_err = syst_sum/expo_sum*correction
        psr_flux_integral[ebin] = flux_content*pow(energy_bin[ebin]/1e3,1)
        psr_flux_integral_err[ebin] = flux_err*pow(energy_bin[ebin]/1e3,1)

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):

        flux_index = 2.75
        data_sum = 0.
        bkgd_sum = 0.
        expo_sum = 0.
        syst_sum = 0.
        for bx in range(0,hist_column_density_map.GetNbinsX()):
            for by in range(0,hist_column_density_map.GetNbinsY()):

                bin_ra = hist_column_density_map.GetXaxis().GetBinCenter(bx+1)
                bin_dec = hist_column_density_map.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(bin_ra-region_ra,2) + pow(bin_dec-region_dec,2)
                if distance_sq>region_radius*region_radius: continue

                data_content = hist_data_skymap[ebin].GetBinContent(bx+1,by+1)
                bkgd_content = hist_bkgd_skymap[ebin].GetBinContent(bx+1,by+1)
                expo_content = hist_expo_skymap[ebin].GetBinContent(bx+1,by+1)
                syst_err = hist_syst_skymap[ebin].GetBinContent(bx+1,by+1)

                data_sum += data_content
                bkgd_sum += bkgd_content
                expo_sum += expo_content
                syst_sum += syst_err

        map_x = region_ra
        map_y = region_dec
        flux_calibration = GetFluxCalibration(map_x,map_y,ebin)
        correction = flux_calibration*pow(region_radius,2)/(calibration_radius*calibration_radius)
        flux_content = 0.
        flux_err = 0.
        if expo_sum>0.:
            flux_content = (data_sum-bkgd_sum)/expo_sum*correction
            flux_err = syst_sum/expo_sum*correction
        #obs_flux_integral[ebin] = flux_content/(flux_index-1)*(energy_bin[ebin]/1000.-pow(energy_bin[ebin+1]/1000.,1-flux_index)/pow(energy_bin[ebin]/1000.,-flux_index))
        #obs_flux_integral_err[ebin] = flux_err/(flux_index-1)*(energy_bin[ebin]/1000.-pow(energy_bin[ebin+1]/1000.,1-flux_index)/pow(energy_bin[ebin]/1000.,-flux_index))
        obs_flux_integral[ebin] = flux_content*pow(energy_bin[ebin]/1e3,1)
        obs_flux_integral_err[ebin] = flux_err*pow(energy_bin[ebin]/1e3,1)

        RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()
        deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
        bin_x_size = hist_column_density_map.GetXaxis().GetBinLowEdge(2)-hist_column_density_map.GetXaxis().GetBinLowEdge(1)
        bin_y_size = hist_column_density_map.GetYaxis().GetBinLowEdge(2)-hist_column_density_map.GetYaxis().GetBinLowEdge(1)
        skymap_bin_area = bin_x_size*bin_y_size*deg_to_cm*deg_to_cm
        Gamma_CR = 2.1
        f_Gamma = 0.9
        #Gamma_CR = 2.2
        #f_Gamma = 0.43
        #Gamma_CR = 2.3
        #f_Gamma = 0.19
        E_SN = 1.
        M_ej = 1. # ejecta mass in solar mass unit
        t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
        CR_efficiency = 3.14/5.*(1.-pow(t_Sedov/t_SN,0.4))
        esc_param = 2.0
        diffusion_coefficient_1TeV = 8.2*1e26

        for bx in range(0,hist_column_density_map.GetNbinsX()):
            for by in range(0,hist_column_density_map.GetNbinsY()):

                bin_ra = hist_column_density_map.GetXaxis().GetBinCenter(bx+1)
                bin_dec = hist_column_density_map.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(bin_ra-region_ra,2) + pow(bin_dec-region_dec,2)
                if distance_sq>region_radius*region_radius: continue

                MC_column_density = hist_column_density_map.GetBinContent(bx+1,by+1)
                MC_mass = MC_column_density*skymap_bin_area
                cloud_mass[ebin] += MC_mass

                hadronic_emission = CalculateHadronicEmission(hist_column_density_map,bin_ra,bin_dec,ebin,CR_efficiency,diffusion_coefficient_1TeV,RA_SN,Dec_SN)

                had_model_flux_integral[ebin] += hadronic_emission

    total_psr_flux_integral = 0.
    total_obs_flux_integral = 0.
    total_model_flux_integral = 0.
    for ebin in range(0,len(energy_bin)-1):
        total_psr_flux_integral += psr_flux_integral[ebin]
        total_obs_flux_integral += obs_flux_integral[ebin]
        total_model_flux_integral += had_model_flux_integral[ebin]

    hydrogen_to_solar_mass = 8.4144035*1e-58
    print ('number of molecules = %0.2e (%0.2e solar mass)'%(cloud_mass[energy_bin_cut_low],cloud_mass[energy_bin_cut_low]*hydrogen_to_solar_mass))
    print ('psr_flux_integral = %0.2e /cm2/s'%(total_psr_flux_integral))
    print ('obs_flux_integral = %0.2e /cm2/s'%(total_obs_flux_integral))
    print ('had_model_flux_integral = %0.2e /cm2/s'%(total_model_flux_integral))

    #return energy_array,obs_flux_integral,obs_flux_integral_err,had_model_flux_integral
    return total_obs_flux_integral, total_model_flux_integral

def flat_disk_model(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    bkgd_scale = Hist_bkgd_global.GetBinContent(binx,biny)

    disk_center_x = par[0]
    disk_center_y = par[1]
    disk_radius = par[2]
    disk_norm = par[3]
    distance_sq = pow(x[0]-disk_center_x,2) + pow(x[1]-disk_center_y,2)
    if distance_sq<disk_radius:
        return disk_norm*bkgd_scale
    return 0.

def flat_ellipse_model(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    bkgd_scale = Hist_bkgd_global.GetBinContent(binx,biny)

    ellipse_center_x1 = par[0]
    ellipse_center_y1 = par[1]
    ellipse_center_x2 = par[2]
    ellipse_center_y2 = par[3]
    ellipse_semi_major_ax = par[4]
    ellipse_norm = par[5]
    if pow(pow(x[0]-ellipse_center_x1,2)+pow(x[1]-ellipse_center_y1,2),0.5) + pow(pow(x[0]-ellipse_center_x2,2)+pow(x[1]-ellipse_center_y2,2),0.5) < 2.*ellipse_semi_major_ax:
        return ellipse_norm*bkgd_scale
    return 0.


def gaussian_1d(x,a,b):

    RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()
    profile_center_x = RA_PSR
    profile_center_y = Dec_PSR
    profile_expo, profile_expo_err, theta2_expo = FindExtension(Hist_bkgd_global,None,profile_center_x,profile_center_y,2.5)
    expo = np.interp(x, theta2_expo, profile_expo)

    return a*np.exp(-0.5*x*x/b*b)*expo


def FitHybridModel2D(Hist_flux):

    MapEdge_left = Hist_flux.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_flux.GetXaxis().GetBinLowEdge(Hist_flux.GetNbinsX()+1)
    MapEdge_lower = Hist_flux.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_flux.GetYaxis().GetBinLowEdge(Hist_flux.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = Hist_flux.GetXaxis().GetBinLowEdge(2)-Hist_flux.GetXaxis().GetBinLowEdge(1)

    RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()

    lep_diffusion_radius = 0.2
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    vol_diff = 4.*3.14/3.*pow(lep_diffusion_radius*deg_to_cm,3)
    total_flux = 0.
    total_area = 0.
    for binx in range(0,Hist_flux.GetNbinsX()):
        for biny in range(0,Hist_flux.GetNbinsY()):
            cell_x = Hist_flux.GetXaxis().GetBinCenter(binx+1)
            cell_y = Hist_flux.GetYaxis().GetBinCenter(biny+1)
            distance = pow(pow(cell_x-RA_PSR,2)+pow(cell_y-Dec_PSR,2),0.5)
            if distance>lep_diffusion_radius: continue
            total_flux += Hist_flux.GetBinContent(binx+1,biny+1)
            total_area += 1.
    E_el_vis = total_flux/100.

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit
    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year

    fix_lep_r_diff = False
    D_ism = 0.082*1e28
    D_ism_upper = (0.082+0.209)*1e28
    D_ism_lower = (0.082-0.059)*1e28

    #diffusion_radius = 1.0
    #esc_param = DeriveSupernovaEscapeParameter(diffusion_radius,t_SN,d_SN,n_0)
    #esc_param = 2.0
    #diffusion_radius = DeriveSupernovaDiffusionLength(esc_param,t_SN,d_SN,n_0,current_gamma_energy)
    #CR_brightness, CR_brightness_err = EstimateCosmicRayBrightness(Hist_flux,286.786,7.1,0.5,diffusion_radius)
    #CR_efficiency = CR_brightness
    CR_efficiency = 3.14/5.*(1.-pow(t_Sedov/t_SN,0.4))

    print ('Initial E_el_vis = %0.3f e-16'%(E_el_vis*1e16))
    print ('Initial CR_efficiency = %0.3f'%(CR_efficiency))
    print ('Initial D_ism = %0.3f e28 cm2/s'%(D_ism/1e28))

    npar = 5
    hybrid_flux_2d = ROOT.TF2('hybrid_flux_2d',hybrid_model_2d,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    hybrid_flux_2d.SetParameter(0,E_el_vis)
    hybrid_flux_2d.SetParLimits(0,0.5*E_el_vis,1.5*E_el_vis)
    hybrid_flux_2d.SetParameter(1,CR_efficiency)
    hybrid_flux_2d.SetParLimits(1,0.5*CR_efficiency,1.5*CR_efficiency)
    hybrid_flux_2d.SetParameter(2,D_ism)
    hybrid_flux_2d.SetParLimits(2,D_ism_lower,D_ism_upper)
    hybrid_flux_2d.FixParameter(2,D_ism)
    hybrid_flux_2d.SetParameter(3,RA_PSR_initial)
    hybrid_flux_2d.FixParameter(3,RA_PSR_initial)
    hybrid_flux_2d.SetParameter(4,Dec_PSR_initial)
    hybrid_flux_2d.FixParameter(4,Dec_PSR_initial)

    Hist_flux_new = Hist_flux.Clone()

    Hist_flux_new.Fit('hybrid_flux_2d')

    E_el_vis = hybrid_flux_2d.GetParameter(0)
    CR_efficiency = hybrid_flux_2d.GetParameter(1)
    CR_efficiency_err = hybrid_flux_2d.GetParError(1)
    print ('2D fit E_el_vis = %0.3f e-16'%(E_el_vis*1e16))
    print ('2D fit CR_efficiency = %0.3f +/- %0.3f'%(CR_efficiency,CR_efficiency_err))

    D_ism = hybrid_flux_2d.GetParameter(2)
    print ('2D fit D_ism = %0.3f e28 cm2/s'%(D_ism/1e28))

    RA_PSR_initial = hybrid_flux_2d.GetParameter(3)
    Dec_PSR_initial = hybrid_flux_2d.GetParameter(4)
    print ('2D fit RA_PSR_initial = %0.2f'%(RA_PSR_initial))
    print ('2D fit Dec_PSR_initial = %0.2f'%(Dec_PSR_initial))


    return E_el_vis, D_ism, CR_efficiency, RA_PSR_initial, Dec_PSR_initial

def FitLeptonicModel2D(Hist_flux):

    MapEdge_left = Hist_flux.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_flux.GetXaxis().GetBinLowEdge(Hist_flux.GetNbinsX()+1)
    MapEdge_lower = Hist_flux.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_flux.GetYaxis().GetBinLowEdge(Hist_flux.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = Hist_flux.GetXaxis().GetBinLowEdge(2)-Hist_flux.GetXaxis().GetBinLowEdge(1)

    RA_PSR_initial, Dec_PSR_initial, RA_PSR, Dec_PSR, d_PSR, t_PSR, n_0 = GetPulsarParameters()

    diffusion_radius = 0.2
    deg_to_cm = 3.14/180.*d_PSR*1000.*3.086e18
    vol_diff = 4.*3.14/3.*pow(diffusion_radius*deg_to_cm,3)
    total_flux = 0.
    total_area = 0.
    for binx in range(0,Hist_flux.GetNbinsX()):
        for biny in range(0,Hist_flux.GetNbinsY()):
            cell_x = Hist_flux.GetXaxis().GetBinCenter(binx+1)
            cell_y = Hist_flux.GetYaxis().GetBinCenter(biny+1)
            distance = pow(pow(cell_x-RA_PSR,2)+pow(cell_y-Dec_PSR,2),0.5)
            if distance>diffusion_radius: continue
            total_flux += Hist_flux.GetBinContent(binx+1,biny+1)
            total_area += 1.
    E_el_vis = total_flux/100.

    npar = 4
    D_ism = 0.082*1e28
    D_ism_upper = (0.082+0.209)*1e28
    D_ism_lower = (0.082-0.059)*1e28
    print ('Initial E_el_vis = %0.3f e-16'%(E_el_vis*1e16))
    print ('Initial D_ism = %0.3f e28 cm2/s'%(D_ism/1e28))
    fix_r_diff = True
    leptonic_flux_2d = ROOT.TF2('leptonic_flux_2d',leptonic_model_2d,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    leptonic_flux_2d.SetParameter(0,E_el_vis)
    leptonic_flux_2d.SetParameter(1,D_ism)
    leptonic_flux_2d.SetParLimits(1,D_ism_lower,D_ism_upper)
    if fix_r_diff: leptonic_flux_2d.FixParameter(1,D_ism)
    leptonic_flux_2d.SetParameter(2,RA_PSR_initial)
    leptonic_flux_2d.SetParameter(3,Dec_PSR_initial)

    Hist_flux_new = Hist_flux.Clone()

    Hist_flux_new.Fit('leptonic_flux_2d')
    E_el_vis = leptonic_flux_2d.GetParameter(0)
    D_ism = leptonic_flux_2d.GetParameter(1)
    print ('2D fit E_el_vis = %0.3f e-16'%(E_el_vis*1e16))
    print ('2D fit D_ism = %0.3f e28 cm2/s'%(D_ism/1e28))

    return E_el_vis, D_ism

def FitHadronicModel2D(Hist_flux):

    MapEdge_left = Hist_flux.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_flux.GetXaxis().GetBinLowEdge(Hist_flux.GetNbinsX()+1)
    MapEdge_lower = Hist_flux.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_flux.GetYaxis().GetBinLowEdge(Hist_flux.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = Hist_flux.GetXaxis().GetBinLowEdge(2)-Hist_flux.GetXaxis().GetBinLowEdge(1)

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit
    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year

    #diffusion_radius = 0.5
    #diffusion_radius_up = DeriveSupernovaDiffusionLength(3.,t_SN,d_SN,n_0,current_gamma_energy)
    #diffusion_radius_dw = DeriveSupernovaDiffusionLength(1.,t_SN,d_SN,n_0,current_gamma_energy)
    #Syst_err = 10.
    #if CR_brightness>0.: Syst_err = CR_brightness_err/CR_brightness

    CR_efficiency = 3.14/5.*(1.-pow(t_Sedov/t_SN,0.4))

    diffusion_radius = 0.5
    esc_param = DeriveSupernovaEscapeParameter(diffusion_radius,t_SN,d_SN,n_0)
    #CR_brightness, CR_brightness_err = EstimateCosmicRayBrightness(Hist_flux,RA_SN,Dec_SN,0.5,diffusion_radius)
    #CR_efficiency = CR_brightness
    D_ism = 0.082*1e28
    D_ism_upper = (0.082+0.209)*1e28
    D_ism_lower = (0.082-0.059)*1e28
    #esc_param = 2.48
    #diffusion_radius = DeriveSupernovaDiffusionLength(esc_param,t_SN,d_SN,n_0,current_gamma_energy)
    #CR_efficiency = 3.14/5.*(1.-pow(t_Sedov/t_SN,0.4))
    print ('Initial CR_efficiency = %0.3f'%(CR_efficiency))
    print ('Initial esc_param = %0.3f'%(esc_param))
    print ('Initial diffusion_radius = %0.3f'%(diffusion_radius))
    npar = 4
    hadronic_flux_2d = ROOT.TF2('hadronic_flux_2d',hadronic_model_2d,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    hadronic_flux_2d.SetParameter(0,CR_efficiency)
    hadronic_flux_2d.SetParLimits(0,0.01*CR_efficiency,2.*CR_efficiency)
    hadronic_flux_2d.SetParameter(1,D_ism)
    hadronic_flux_2d.SetParLimits(1,D_ism_lower,D_ism_upper)
    hadronic_flux_2d.FixParameter(1,D_ism)
    hadronic_flux_2d.SetParameter(2,RA_SN)
    hadronic_flux_2d.FixParameter(2,RA_SN)
    hadronic_flux_2d.SetParameter(3,Dec_SN)
    hadronic_flux_2d.FixParameter(3,Dec_SN)

    Hist_flux_new = Hist_flux.Clone()

    Hist_flux_new.Fit('hadronic_flux_2d')

    CR_efficiency = hadronic_flux_2d.GetParameter(0)
    D_ism = hadronic_flux_2d.GetParameter(1)
    CR_efficiency_err = hadronic_flux_2d.GetParError(0)
    D_ism_err = hadronic_flux_2d.GetParError(1)
    print ('2D fit CR_efficiency = %0.3f +/- %0.3f'%(CR_efficiency,CR_efficiency_err))
    print ('2D fit D_ism = %0.3f'%(D_ism))

    return CR_efficiency, D_ism

def GetExtention2D(Hist_data_input, Hist_bkgd_input, Hist_syst_input, highlight_threshold, init_x, init_y, tag):

    excess_center_x_init = init_x
    excess_center_y_init = init_y
    excess_radius_init = 0.1

    excess_center_x_all = [0.,0.,0.]
    excess_center_y_all = [0.,0.,0.]
    excess_radius_all = [0.,0.,0.]
    excess_norm_all = [0.,0.,0.]
    excess_disk_all = [0.,0.,0.]

    print ("Hist_data_input.Integral() = %s"%(Hist_data_input.Integral()))
    print ("Hist_bkgd_input.Integral() = %s"%(Hist_bkgd_input.Integral()))
    print ("Hist_syst_input.Integral() = %s"%(Hist_syst_input.Integral()))
    if Hist_data_input.Integral()==0.:
        return excess_center_x_all, excess_center_y_all, excess_radius_all, excess_disk_all, 0.

    Hist_Excess = Hist_data_input.Clone()
    Hist_Excess_up = Hist_data_input.Clone()
    Hist_Excess_dw = Hist_data_input.Clone()

    Hist_Excess.Reset()
    Hist_Excess_up.Reset()
    Hist_Excess_dw.Reset()
    Hist_Excess.Add(Hist_data_input)
    Hist_Excess.Add(Hist_bkgd_input,-1.)
    Hist_Excess.Divide(Hist_bkgd_input)
    Hist_Excess_up.Add(Hist_data_input)
    Hist_Excess_up.Add(Hist_bkgd_input,-1.)
    Hist_Excess_up.Add(Hist_syst_input,-1.)
    Hist_Excess_up.Divide(Hist_bkgd_input)
    Hist_Excess_dw.Add(Hist_data_input)
    Hist_Excess_dw.Add(Hist_bkgd_input,-1.)
    Hist_Excess_dw.Add(Hist_syst_input,1.)
    Hist_Excess_dw.Divide(Hist_bkgd_input)

    Hist_Syst = Hist_syst_input.Clone()
    Hist_Syst.Divide(Hist_bkgd_input)

    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bkgd_content = Hist_bkgd_input.GetBinContent(bx+1,by+1)
            bkgd_error = Hist_bkgd_input.GetBinError(bx+1,by+1)
            if bkgd_error==0.:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
                Hist_Syst.SetBinContent(bx+1,by+1,0.)
                Hist_Excess.SetBinError(bx+1,by+1,0.)
                Hist_Excess_up.SetBinContent(bx+1,by+1,0.)
                Hist_Excess_dw.SetBinContent(bx+1,by+1,0.)
                continue
            if bkgd_content/bkgd_error<5.:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
                Hist_Syst.SetBinContent(bx+1,by+1,0.)
                Hist_Excess.SetBinError(bx+1,by+1,0.)
                Hist_Excess_up.SetBinContent(bx+1,by+1,0.)
                Hist_Excess_dw.SetBinContent(bx+1,by+1,0.)
                continue

    MapEdge_left = Hist_Excess.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Excess.GetXaxis().GetBinLowEdge(Hist_Excess.GetNbinsX()+1)
    MapEdge_lower = Hist_Excess.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Excess.GetYaxis().GetBinLowEdge(Hist_Excess.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = Hist_Excess.GetXaxis().GetBinLowEdge(2)-Hist_Excess.GetXaxis().GetBinLowEdge(1)

    map_x, map_y = np.linspace(MapEdge_left, MapEdge_right, Hist_Excess.GetNbinsX()), np.linspace(MapEdge_lower, MapEdge_upper, Hist_Excess.GetNbinsY())
    grid_x, grid_y = np.meshgrid(map_x, map_y)

    # The function to be fit is grid_z.
    grid_z = np.zeros(grid_x.shape)
    grid_z_up = np.zeros(grid_x.shape)
    grid_z_dw = np.zeros(grid_x.shape)
    grid_error = np.zeros(grid_x.shape)
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            grid_z[bx,by] = Hist_Excess.GetBinContent(bx+1,by+1)
            grid_z_up[bx,by] = Hist_Excess_up.GetBinContent(bx+1,by+1)
            grid_z_dw[bx,by] = Hist_Excess_dw.GetBinContent(bx+1,by+1)
            stat_err = Hist_Excess.GetBinError(bx+1,by+1)
            syst_err = Hist_Syst.GetBinContent(bx+1,by+1)
            grid_error[bx,by] = max(0.01,pow(stat_err*stat_err+syst_err*syst_err,0.5))

    # We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
    xydata = np.vstack((grid_x.ravel(), grid_y.ravel()))

    A_init = Hist_Excess.Integral()
    guess_prms = [(excess_center_x_init,excess_center_y_init,excess_radius_init,A_init,0.)]
    # Flatten the initial guess parameter list.
    start = [p for prms in guess_prms for p in prms]

    peak_location_range = 1.*excess_radius_init
    popt, pcov = curve_fit(_gaussian_2d, xydata, grid_z.ravel(), p0=start, sigma=grid_error.ravel(),bounds=((excess_center_x_init-peak_location_range, excess_center_y_init-peak_location_range, 0.1, 0, 0), (excess_center_x_init+peak_location_range, excess_center_y_init+peak_location_range, 2.0, A_init, A_init)))
    perr = np.sqrt(np.diag(pcov))
    model_fit = gaussian_2d(grid_x,grid_y, *popt)
    residual = grid_z - model_fit
    chisq = np.sum((residual/grid_error)**2)
    dof = Hist_Excess.GetNbinsX()*Hist_Excess.GetNbinsY()-3
    excess_center_x = popt[0]
    excess_center_y = popt[1]
    excess_radius = popt[2]
    excess_norm = popt[3]
    excess_disk = popt[4]
    excess_radius = pow(max(0.,excess_radius*excess_radius-smooth_size_spectroscopy*smooth_size_spectroscopy),0.5)
    excess_center_x_err = perr[0]
    excess_center_y_err = perr[1]
    excess_radius_err = perr[2]
    excess_norm_err = perr[3]
    excess_disk_err = perr[4]
    print ('excess_radius = %0.3f'%(excess_radius))
    print ('excess_norm = %0.3f'%(excess_norm))
    print ('excess_disk = %0.3f'%(excess_disk))

    ax.cla()
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.imshow(grid_z, origin='lower', cmap='plasma', extent=(map_x.min(), map_x.max(), map_y.min(), map_y.max()))
    ax.contour(grid_x, grid_y, model_fit, colors='w')
    plotname = 'Fit_source_extent_%s'%(tag)
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    ax.cla()
    #ax_3d = fig.gca(projection='3d')
    #ax_3d.plot_surface(grid_x, grid_y, grid_z, cmap='plasma')
    #cset = ax_3d.contour(grid_x, grid_y, model_fit, zdir='z', offset=-4, colors='b')
    #ax_3d.set_zlim(-4,np.max(grid_z))
    #plotname = 'Fit_source_extent_%s'%(tag)
    #fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')

    popt, pcov = curve_fit(_gaussian_2d, xydata, grid_z_up.ravel(), p0=start, sigma=grid_error.ravel(), bounds=((excess_center_x_init-peak_location_range, excess_center_y_init-peak_location_range, 0.1, 0, 0), (excess_center_x_init+peak_location_range, excess_center_y_init+peak_location_range, 2.0, A_init, A_init)))
    excess_center_x_up = popt[0]
    excess_center_y_up = popt[1]
    excess_radius_up = popt[2]
    excess_norm_up = popt[3]
    excess_disk_up = popt[4]
    excess_radius_up = pow(max(0.,excess_radius_up*excess_radius_up-smooth_size_spectroscopy*smooth_size_spectroscopy),0.5)
    print ('excess_radius_up = %0.3f'%(excess_radius_up))
    print ('excess_norm_up = %0.3f'%(excess_norm_up))

    popt, pcov = curve_fit(_gaussian_2d, xydata, grid_z_dw.ravel(), p0=start, sigma=grid_error.ravel(), bounds=((excess_center_x_init-peak_location_range, excess_center_y_init-peak_location_range, 0.1, 0, 0), (excess_center_x_init+peak_location_range, excess_center_y_init+peak_location_range, 2.0, A_init, A_init)))
    excess_center_x_dw = popt[0]
    excess_center_y_dw = popt[1]
    excess_radius_dw = popt[2]
    excess_norm_dw = popt[3]
    excess_disk_dw = popt[4]
    excess_radius_dw = pow(max(0.,excess_radius_dw*excess_radius_dw-smooth_size_spectroscopy*smooth_size_spectroscopy),0.5)
    print ('excess_radius_dw = %0.3f'%(excess_radius_dw))
    print ('excess_norm_dw = %0.3f'%(excess_norm_dw))

    excess_center_x_all = [excess_center_x,excess_center_x_up,excess_center_x_dw]
    excess_center_y_all = [excess_center_y,excess_center_y_up,excess_center_y_dw]
    excess_radius_all = [excess_radius,excess_radius_up,excess_radius_dw]
    excess_norm_all = [excess_norm,excess_norm_up,excess_norm_dw]
    excess_disk_all = [excess_disk,excess_disk_up,excess_disk_dw]

    return excess_center_x_all, excess_center_y_all, excess_radius_all, excess_disk_all, chisq/dof

def Make2DProjectionPlot(Hist_Data,xtitle,ytitle,name,doProj):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    #pad1.SetGrid()
    pad1.Draw()
    pad1.cd()
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.Draw("COL4Z")
    #Hist_Data.Draw("CONT3 same")
    bins = []
    for b in range(0,Hist_Data.GetNbinsX()+1):
        bins += [Hist_Data.GetXaxis().GetBinLowEdge(b+1)]
    Hist_1D = ROOT.TH1D("Hist_1D","",len(bins)-1,array('d',bins))
    for b in range(0,Hist_Data.GetNbinsX()):
        hist_temp = Hist_Data.ProjectionY("hist_temp",b+1,b+1)
        Hist_1D.SetBinContent(b+1,hist_temp.GetMean())
        Hist_1D.SetBinError(b+1,hist_temp.GetRMS())
    Hist_1D.SetLineColor(2)
    if doProj: Hist_1D.Draw("E same")

    canvas.SaveAs('output_plots/%s.png'%(name))

def MatrixDecompositionDemo(name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.90,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.90)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.05)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    line1 = ROOT.TLine(-MSCL_blind_cut,MSCW_blind_cut,MSCL_blind_cut,MSCW_blind_cut)
    line1.SetLineStyle(10)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line2 = ROOT.TLine(MSCL_blind_cut,-MSCW_blind_cut,MSCL_blind_cut,MSCW_blind_cut)
    line2.SetLineStyle(10)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line3 = ROOT.TLine(-MSCL_blind_cut,-MSCW_blind_cut,MSCL_blind_cut,-MSCW_blind_cut)
    line3.SetLineStyle(10)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line4 = ROOT.TLine(-MSCL_blind_cut,-MSCW_blind_cut,-MSCL_blind_cut,MSCW_blind_cut)
    line4.SetLineStyle(10)
    line4.SetLineColor(2)
    line4.SetLineWidth(2)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.25,'M^{Full} = #sum_{k=1}^{8} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Sum.GetYaxis().SetTitle('scaled width')
    Hist2D_OnData_Sum.GetXaxis().SetTitle('scaled length')
    Hist2D_OnData_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    #line3.Draw("same")
    #line4.Draw("same")
    canvas.SaveAs('output_plots/OnData_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.25,'M^{Full} = #sum_{k=1}^{8} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnDark_Sum.GetYaxis().SetTitle('scaled width')
    Hist2D_OnDark_Sum.GetXaxis().SetTitle('scaled length')
    Hist2D_OnDark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    #line3.Draw("same")
    #line4.Draw("same")
    canvas.SaveAs('output_plots/OnDark_%s_%s.png'%(name,selection_tag))

    entry = 1
    for entry in range(0,3):
        Hists = []
        legends = []
        colors = []
        Hists += [Hist2D_Coeff_Data_Sum.ProjectionX("hist_temp_x",Hist2D_Coeff_Data_Sum.GetNbinsY()-entry,Hist2D_Coeff_Data_Sum.GetNbinsY()-entry)]
        legends += ['C_{i,1}']
        colors += [1]
        Hists += [Hist2D_Coeff_Data_Sum.ProjectionY("hist_temp_y",Hist2D_Coeff_Data_Sum.GetNbinsX()-entry,Hist2D_Coeff_Data_Sum.GetNbinsX()-entry)]
        legends += ['C_{1,j}']
        colors += [2]
        MakeMultiplePlot(Hists,legends,colors,'entry','C_{ij}','Coeff_C_entry%s_%s'%(entry,name),0,0,False,False)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'t_{ij} coefficents' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.20,'(Zenith = %0.1f #pm %0.1f, NSB = %0.1f #pm %0.1f)'%(Zenith_mean_data,Zenith_RMS_data,NSB_mean_data,NSB_RMS_data) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.2)
    lumilab2.Draw()
    pad1.cd()
    #pad1.SetLogz()
    Hist2D_Coeff_Bkgd_Sum.GetYaxis().SetTitle('n')
    Hist2D_Coeff_Bkgd_Sum.GetXaxis().SetTitle('k')
    Hist2D_Coeff_Bkgd_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_C_bkgd_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'t_{ij} coefficents')
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.20,'(Zenith = %0.1f #pm %0.1f, NSB = %0.1f #pm %0.1f)'%(Zenith_mean_data,Zenith_RMS_data,NSB_mean_data,NSB_RMS_data) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.2)
    lumilab2.Draw()
    pad1.cd()
    #pad1.SetLogz()
    #Hist2D_Coeff_Data_Sum.SetMaximum(Hist2D_Coeff_Bkgd_Sum.GetMaximum());
    #Hist2D_Coeff_Data_Sum.SetMinimum(Hist2D_Coeff_Bkgd_Sum.GetMinimum());
    Hist2D_Coeff_Data_Sum.GetYaxis().SetTitle('n')
    Hist2D_Coeff_Data_Sum.GetXaxis().SetTitle('k')
    Hist2D_Coeff_Data_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_C_data_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.20,'#delta H_{ij} coefficents' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #pad1.SetLogz()
    ##Hist2D_U_Proj_Sum.SetMaximum(1e-1);
    ##Hist2D_U_Proj_Sum.SetMinimum(1e-4);
    #Hist2D_U_Proj_Sum.GetYaxis().SetTitle('y')
    #Hist2D_U_Proj_Sum.GetXaxis().SetTitle('x')
    #Hist2D_U_Proj_Sum.Draw("COL4Z")
    #canvas.SaveAs('output_plots/Coeff_U_proj_%s_%s.png'%(name,selection_tag))
    #pad1.SetLogz(0)

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.20,'#delta H_{ij} coefficents' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #pad1.SetLogz()
    ##Hist2D_V_Proj_Sum.SetMaximum(1e-1);
    ##Hist2D_V_Proj_Sum.SetMinimum(1e-4);
    #Hist2D_V_Proj_Sum.GetYaxis().SetTitle('y')
    #Hist2D_V_Proj_Sum.GetXaxis().SetTitle('x')
    #Hist2D_V_Proj_Sum.Draw("COL4Z")
    #canvas.SaveAs('output_plots/Coeff_V_proj_%s_%s.png'%(name,selection_tag))
    #pad1.SetLogz(0)

    #pad3.Clear()
    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sum_{k=1}^{1} #sigma_{k} u_{k} v_{k}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank0_Data_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank0_Data_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank0_Data_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank0_Data_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sigma_{1} u_{1} v_{1}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank0_Dark_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank0_Dark_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank0_Dark_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank0_Dark_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sum_{k=1}^{2} #sigma_{k} u_{k} v_{k}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank1_Data_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank1_Data_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank1_Data_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank1_Data_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sigma_{2} u_{2} v_{2}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank1_Dark_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank1_Dark_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank1_Dark_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank1_Dark_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sum_{k=1}^{3} #sigma_{k} u_{k} v_{k}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank2_Data_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank2_Data_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank2_Data_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank2_Data_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sigma_{3} u_{3} v_{3}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank2_Dark_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank2_Dark_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank2_Dark_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank2_Dark_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sum_{k=1}^{4} #sigma_{k} u_{k} v_{k}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank3_Data_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank3_Data_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank3_Data_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank3_Data_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sigma_{4} u_{4} v_{4}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank3_Dark_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank3_Dark_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank3_Dark_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank3_Dark_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'M^{ON} - #sum_{k=1}^{5} #sigma_{k} u_{k} v_{k}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank4_Data_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank4_Data_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank4_Data_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank4_Data_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'#sigma_{5} u_{5} v_{5}^{T}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_Rank4_Dark_Sum.GetYaxis().SetTitle('scaled width')
    #Hist2D_Rank4_Dark_Sum.GetXaxis().SetTitle('scaled length')
    #Hist2D_Rank4_Dark_Sum.Draw("COL4Z")
    #line1.Draw("same")
    #line2.Draw("same")
    #canvas.SaveAs('output_plots/Rank4_Dark_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrDark = ROOT.TH2D("Hist2D_ErrDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrDark.Add(Hist2D_OnData_Sum)
    Hist2D_ErrDark.Add(Hist2D_OnDark_Sum,-1)
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrDark.GetNbinsX()):
        for biny in range (0,Hist2D_ErrDark.GetNbinsY()):
            if Hist2D_ErrDark.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrDark.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnDark_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    error_avg = sum_error/sum_bins
    sum_var = 0.
    for binx in range (0,Hist2D_ErrDark.GetNbinsX()):
        for biny in range (0,Hist2D_ErrDark.GetNbinsY()):
            if Hist2D_ErrDark.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrDark.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnDark_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd-error_avg,2)
    error_var = sum_var/sum_bins
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.25,'')
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrDark.GetYaxis().SetTitle('scaled width')
    Hist2D_ErrDark.GetXaxis().SetTitle('scaled length')
    Hist2D_ErrDark.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")
    line4.Draw("same")
    canvas.SaveAs('output_plots/ErrorDark_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrBkgd = ROOT.TH2D("Hist2D_ErrBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrBkgd.Add(Hist2D_OnData_Sum)
    Hist2D_ErrBkgd.Add(Hist2D_OnBkgd_Sum,-1)
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    error_avg = sum_error/sum_bins
    sum_var = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd-error_avg,2)
    error_var = sum_var/sum_bins
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)>MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)>MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    sum_var = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)>MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)>MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd,2)
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.25,'')
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.25)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrBkgd.SetMaximum(Hist2D_ErrDark.GetMaximum());
    Hist2D_ErrBkgd.SetMinimum(Hist2D_ErrDark.GetMinimum());
    Hist2D_ErrBkgd.GetYaxis().SetTitle('scaled width')
    Hist2D_ErrBkgd.GetXaxis().SetTitle('scaled length')
    Hist2D_ErrBkgd.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")
    line4.Draw("same")
    canvas.SaveAs('output_plots/ErrorBkgd_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrRing = ROOT.TH2D("Hist2D_ErrRing","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrRing.Add(Hist2D_OnData_Point_Sum)
    Hist2D_ErrRing.Add(Hist2D_OnData_Ring_Sum,-1)
    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'M^{ON}-M^{Ring}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrRing.SetMaximum(Hist2D_ErrDark.GetMaximum());
    Hist2D_ErrRing.SetMinimum(Hist2D_ErrDark.GetMinimum());
    Hist2D_ErrRing.GetYaxis().SetTitle('scaled width')
    Hist2D_ErrRing.GetXaxis().SetTitle('scaled length')
    Hist2D_ErrRing.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")
    line4.Draw("same")
    canvas.SaveAs('output_plots/ErrorRing_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'M^{ON}-M^{Ring}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Point_Sum.GetYaxis().SetTitle('scaled width')
    Hist2D_OnData_Point_Sum.GetXaxis().SetTitle('scaled length')
    Hist2D_OnData_Point_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")
    line4.Draw("same")
    canvas.SaveAs('output_plots/MatrixPoint_%s_%s.png'%(name,selection_tag))

    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.25,'M^{ON}-M^{Ring}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.25)
    #lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Ring_Sum.GetYaxis().SetTitle('scaled width')
    Hist2D_OnData_Ring_Sum.GetXaxis().SetTitle('scaled length')
    Hist2D_OnData_Ring_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")
    line4.Draw("same")
    canvas.SaveAs('output_plots/MatrixRing_%s_%s.png'%(name,selection_tag))

    #norm_data_sr = Hist_OnData_SR_XYoff_Sum.Integral()
    #Hist_OnData_SR_XYoff_Sum.Scale(1./norm_data_sr)
    #norm_dark_sr = Hist_OnDark_SR_XYoff_Sum.Integral()
    #Hist_OnDark_SR_XYoff_Sum.Scale(1./norm_dark_sr)
    #norm_data_cr = Hist_OnData_CR_XYoff_Sum.Integral()
    #Hist_OnData_CR_XYoff_Sum.Scale(1./norm_data_cr)
    #Hist2D_ErrXYoff_DarkSR = ROOT.TH2D("Hist2D_ErrXYoff_DarkSR","",12,-3,3,12,-3,3)
    #Hist2D_ErrXYoff_DarkSR.Add(Hist_OnData_SR_XYoff_Sum)
    #Hist2D_ErrXYoff_DarkSR.Add(Hist_OnDark_SR_XYoff_Sum,-1)
    #Hist2D_ErrXYoff_DarkSR.Divide(Hist_OnData_SR_XYoff_Sum)
    #Hist2D_ErrXYoff_DataCR = ROOT.TH2D("Hist2D_ErrXYoff_DataCR","",12,-3,3,12,-3,3)
    #Hist2D_ErrXYoff_DataCR.Add(Hist_OnData_SR_XYoff_Sum)
    #Hist2D_ErrXYoff_DataCR.Add(Hist_OnData_CR_XYoff_Sum,-1)
    #Hist2D_ErrXYoff_DataCR.Divide(Hist_OnData_SR_XYoff_Sum)
    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.20,'P^{ON}-P^{OFF}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.3)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_ErrXYoff_DarkSR.SetMaximum(0.2);
    #Hist2D_ErrXYoff_DarkSR.SetMinimum(-0.2);
    #Hist2D_ErrXYoff_DarkSR.GetYaxis().SetTitle('Xoff')
    #Hist2D_ErrXYoff_DarkSR.GetXaxis().SetTitle('Yoff')
    #Hist2D_ErrXYoff_DarkSR.Draw("COL4Z")
    #canvas.SaveAs('output_plots/ErrXYoff_DarkSR_%s_%s.png'%(name,selection_tag))
    #pad3.cd()
    #lumilab1 = ROOT.TLatex(0.15,0.20,'P^{ON}-P^{CR}' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.3)
    #lumilab1.Draw()
    #pad1.cd()
    #Hist2D_ErrXYoff_DataCR.SetMaximum(0.2);
    #Hist2D_ErrXYoff_DataCR.SetMinimum(-0.2);
    #Hist2D_ErrXYoff_DataCR.GetYaxis().SetTitle('Xoff')
    #Hist2D_ErrXYoff_DataCR.GetXaxis().SetTitle('Yoff')
    #Hist2D_ErrXYoff_DataCR.Draw("COL4Z")
    #canvas.SaveAs('output_plots/ErrXYoff_DataCR_%s_%s.png'%(name,selection_tag))

def MakeOneHistPlot(Hist,title_x,title_y,name,logy):
    
    global max_chi2_diff2_position_this_energy

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.2)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.2)
    pad1.SetBorderMode(0)
    pad1.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

    Hist.GetXaxis().SetTitleOffset(1.0)
    Hist.GetXaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetTitleOffset(1.0)
    Hist.GetYaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetTitle(title_x)
    Hist.GetYaxis().SetTitle(title_y)
    Hist.Draw("E")

    if Hist.GetXaxis().GetLabelOffset()==999:
        x1 = Hist.GetXaxis().GetXmin()
        x2 = Hist.GetXaxis().GetXmax()
        y1 = Hist.GetYaxis().GetXmin()
        y2 = Hist.GetYaxis().GetXmax()
        IncValues = ROOT.TF1( "IncValues", "x", 0, 256 )
        raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
        raLowerAxis.SetLabelSize(Hist.GetXaxis().GetLabelSize())
        raLowerAxis.Draw()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeAccumulatedSignificancePlot(hist_data,hist_bkgd,hist_syst,direction,title_x,title_y,name):

    hist_sig = hist_data.Clone()
    hist_syst_content = hist_data.Clone()
    for binx in range(0,hist_data.GetNbinsX()):
        hist_syst_content.SetBinContent(binx+1,hist_syst.GetBinError(binx+1))
    for binx in range(0,hist_data.GetNbinsX()):
        data = 0.
        bkgd = 0.
        syst = 0.
        if direction==0:
            data = hist_data.Integral(1,binx+1)
            bkgd = hist_bkgd.Integral(1,binx+1)
            syst = hist_syst_content.Integral(1,binx+1)
        else:
            data = hist_data.Integral(binx+1,hist_data.GetNbinsX())
            bkgd = hist_bkgd.Integral(binx+1,hist_data.GetNbinsX())
            syst = hist_syst_content.Integral(binx+1,hist_data.GetNbinsX())
        sig = CalculateSignificance(data-bkgd,bkgd,syst)
        print ('data = %0.2f, bkgd = %0.2f, syst = %0.2f, sig = %0.2f'%(data,bkgd,syst,sig))
        hist_sig.SetBinContent(binx+1,sig)
        hist_sig.SetBinError(binx+1,0.)
    MakeOneHistPlot(hist_sig,title_x,title_y,'%s_%s'%(name,selection_tag),False)

def MakeRankResidualPlots(name):

    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    data_integral = Hist2D_OnData_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)
    data_rank_integral = []
    data_rank_integral += [Hist2D_Rank0_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank1_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank2_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank3_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank4_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]

    Hist_Residual_Gamma = ROOT.TH1D("Hist_Residual_Gamma","",5,0,5)
    for binx in range(0,len(data_rank_integral)):
        bkgd_integral = 0.
        #for rank in range (0,binx+1):
        #    bkgd_integral += data_rank_integral[rank]
        bkgd_integral = data_rank_integral[binx]
        #Hist_Residual_Gamma.SetBinContent(binx+1,abs(1.-bkgd_integral/data_integral))
        Hist_Residual_Gamma.SetBinContent(binx+1,abs(bkgd_integral/data_integral))
        ratio_err = bkgd_integral/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Gamma.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Gamma,'n ranks','Residual in gamma region','Residual_Gamma_%s_%s'%(name,selection_tag),True)

    data_integral = Hist2D_OnData_Sum.Integral()
    data_rank_integral = []
    data_rank_integral += [Hist2D_Rank0_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank1_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank2_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank3_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank4_Data_Sum.Integral()]

    Hist_Residual_Full = ROOT.TH1D("Hist_Residual_Full","",5,0,5)
    for binx in range(0,len(data_rank_integral)):
        bkgd_integral = 0.
        #for rank in range (0,binx+1):
        #    bkgd_integral += data_rank_integral[rank]
        bkgd_integral = data_rank_integral[binx]
        #Hist_Residual_Full.SetBinContent(binx+1,abs(1.-bkgd_integral/data_integral))
        Hist_Residual_Full.SetBinContent(binx+1,abs(bkgd_integral/data_integral))
        ratio_err = bkgd_integral/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Full.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Full,'n ranks','Residual in gamma region','Residual_Full_%s_%s'%(name,selection_tag),True)

def MakeSystematicPlots(hists, colors, fillstyles, legends, exposure, name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'Total exposure = %0.1f h'%(exposure) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(hists)):
        legend.AddEntry(hists[h],legends[h])
    legend.Draw("same")

    pad1.cd()

    hists[0].GetXaxis().SetTitleOffset(0.8)
    hists[0].GetXaxis().SetTitleSize(0.06)
    hists[0].GetXaxis().SetLabelSize(0.06)
    hists[0].GetYaxis().SetLabelSize(0.06)
    hists[0].GetYaxis().SetTitleOffset(1.2)
    hists[0].GetYaxis().SetTitleSize(0.06)
    hists[0].GetXaxis().SetTitle('relative error')
    hists[0].GetYaxis().SetTitle('number of measurements')
    for h in range(0,len(hists)):
        hists[h].SetFillStyle(fillstyles[h])
        hists[h].SetFillColor(colors[h])
        hists[h].SetLineColor(colors[h])
        hists[h].SetLineWidth(2)
        hists[h].Draw("hist same")
    pad3.cd()
    c_both.SaveAs('output_plots/%s.png'%(name))

def SystematicAnalysis():

    FilePath_List = []
    for elev in range(0,len(root_file_tags)):
        FilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_"%(folder_path)+"Syst_%s"%(root_file_tags[elev])+".root"
        FilePath_List += [FilePath]
        print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
        InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
        bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
        bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
        bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

        Hist2D_OnData_AllGroups = []
        Hist2D_OnDark_AllGroups = []
        Hist2D_OnBkgd_AllGroups = []
        for nth_group in range(0,InfoTree.number_groups):
            Hist2D_OnData_AllGroups += [ROOT.TH2D("Hist2D_OnData_%s"%(nth_group),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnDark_AllGroups += [ROOT.TH2D("Hist2D_OnDark_%s"%(nth_group),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnBkgd_AllIterations = []
            for nth_iteration in range(0,20):
                Hist2D_OnBkgd_AllIterations += [ROOT.TH2D("Hist2D_OnBkgd_%s_%s"%(nth_group,nth_iteration),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnBkgd_AllGroups += [Hist2D_OnBkgd_AllIterations]

        for e in range(0,len(energy_bin)-1):
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            if ErecS_upper_cut<=energy_bin[energy_bin_cut_low]: continue
            if ErecS_lower_cut>=energy_bin[energy_bin_cut_up]: continue
            ErecS_lower_cut_int = int(ErecS_lower_cut)
            ErecS_upper_cut_int = int(ErecS_upper_cut)

            for nth_group in range(0,InfoTree.number_groups):
                HistName = "Hist_NthGroup_Data_MSCLW_V%s_ErecS%sto%s"%(nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                Hist2D_OnData_AllGroups[nth_group].Add(InputFile.Get(HistName))
                for nth_iteration in range(0,20):
                    HistName = "Hist_NthGroup_Bkgd_MSCLW_I%s_V%s_ErecS%sto%s"%(nth_iteration,nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                    Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration].Add(InputFile.Get(HistName))
                HistName = "Hist_NthGroup_Dark_MSCLW_V%s_ErecS%sto%s"%(nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                Hist2D_OnDark_AllGroups[nth_group].Add(InputFile.Get(HistName))

        Hist2D_Converge = ROOT.TH2D("Hist_Converge","",20,0,20,20,0.,0.01)
        for nth_iteration in range(0,20):
            bkgd_chi2 = 0.
            bkgd_syst_err = 0.
            dark_syst_err = 0.
            Hist_Syst_Bkgd = ROOT.TH1D("Hist_Syst_Bkgd","",19,-0.2,0.2)
            Hist_Syst_Dark = ROOT.TH1D("Hist_Syst_Dark","",19,-0.2,0.2)
            for nth_group in range(0,InfoTree.number_groups):
                Total_Data = Hist2D_OnData_AllGroups[nth_group].Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)
                if Total_Data == 0.: continue

                Hist2D_Diff = Hist2D_OnData_AllGroups[nth_group].Clone()
                Hist2D_Diff.Add(Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration],-1.)
                Hist2D_Chi2Diff = Hist2D_OnBkgd_AllGroups[nth_group][19].Clone()
                Hist2D_Chi2Diff.Add(Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration],-1.)
                syst_this = Hist2D_Diff.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)/Total_Data
                Hist_Syst_Bkgd.Fill(syst_this)
                bkgd_syst_err += pow(syst_this,2)/InfoTree.number_groups
                chi2_this = 0.
                for binx in range(0,Hist2D_Chi2Diff.GetNbinsX()):
                    for biny in range(0,Hist2D_Chi2Diff.GetNbinsY()):
                        if binx>bin_upper_x : continue
                        if biny>bin_upper_y : continue
                        chi2_this += pow(Hist2D_Chi2Diff.GetBinContent(binx+1,biny+1),2)
                chi2_this = pow(chi2_this,0.5)/Total_Data
                bkgd_chi2 += pow(chi2_this,2)/InfoTree.number_groups
                Hist2D_Converge.Fill(nth_iteration,chi2_this)

                Hist2D_Diff = Hist2D_OnData_AllGroups[nth_group].Clone()
                Hist2D_Diff.Add(Hist2D_OnDark_AllGroups[nth_group],-1.)
                syst_this = Hist2D_Diff.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)/Total_Data
                Hist_Syst_Dark.Fill(syst_this)
                dark_syst_err += pow(syst_this,2)/InfoTree.number_groups

            bkgd_chi2 = pow(bkgd_chi2,0.5)
            bkgd_syst_err = pow(bkgd_syst_err,0.5)
            dark_syst_err = pow(dark_syst_err,0.5)

            hists = []
            colors = []
            fillstyles = []
            legends = []
            hists += [Hist_Syst_Bkgd]
            colors += [38]
            fillstyles += [3325]
            legends += ['MDM syst. = %0.3f'%(bkgd_syst_err)]
            hists += [Hist_Syst_Dark]
            colors += [46]
            fillstyles += [3352]
            legends += ['Initial syst. = %0.3f'%(dark_syst_err)]
            name = "SystAnalysis_I%s"%(nth_iteration)
            MakeSystematicPlots(hists, colors, fillstyles, legends, InfoTree.exposure_hours_usable, name)

        canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
        pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
        pad3.SetBottomMargin(0.0)
        pad3.SetTopMargin(0.03)
        pad3.SetBorderMode(1)
        pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
        pad1.SetBottomMargin(0.15)
        pad1.SetRightMargin(0.15)
        pad1.SetLeftMargin(0.15)
        pad1.SetTopMargin(0.0)
        pad1.SetBorderMode(0)
        pad1.Draw()
        pad3.Draw()
        pad3.cd()
        lumilab1 = ROOT.TLatex(0.15,0.20,'#chi^{2} = #sum_{ij} (N^{nth}_{ij}-N^{20th}_{ij})^{2} / (#sum_{ij} N^{data}_{ij})^{2}' )
        lumilab1.SetNDC()
        lumilab1.SetTextSize(0.3)
        lumilab1.Draw()
        pad1.cd()
        Hist2D_Converge.GetYaxis().SetTitle('#sqrt{#chi^{2}}')
        Hist2D_Converge.GetXaxis().SetTitle('iterations')
        Hist2D_Converge.Draw("COL4Z")
        canvas.SaveAs('output_plots/Converge_%s.png'%(selection_tag))

        InputFile.Close()

def GetPulsarParameters():

    RA_PSR_initial = source_ra
    Dec_PSR_initial = source_dec
    RA_PSR_final = source_ra
    Dec_PSR_final = source_dec
    d_PSR = 3.2 #kpc
    t_PSR = 19.5*1000. #yr
    n_0 = 3.0 # cm^{-3} ambient density

    if sys.argv[1]=='MGRO_J1908_ON':

        if J1908_model=='SNR_7p9kpc_PSR_stationary' or J1908_model=='SNR_3p4kpc_PSR_stationary':
            #PSR J1907+0602
            RA_PSR_initial = 286.975
            Dec_PSR_initial = 6.03777777778
            RA_PSR_final = 286.975
            Dec_PSR_final = 6.03777777778
            d_PSR = 3.2 #kpc
            t_PSR = 19.5*1000. #yr
            n_0 = 3.0 # cm^{-3} ambient density

        else:
            #PSR J1907+0602
            RA_PSR_initial = 287.2
            Dec_PSR_initial = 6.6
            RA_PSR_final = 286.975
            Dec_PSR_final = 6.03777777778
            d_PSR = 1.884 #kpc
            #d_PSR = 3.2 #kpc
            t_PSR = 19.5*1000. #yr
            n_0 = 3.0 # cm^{-3} ambient density

    elif sys.argv[1]=='Geminga_ON':

        RA_PSR_initial = 98.476
        Dec_PSR_initial = 17.770
        RA_PSR_final = 98.476
        Dec_PSR_final = 17.770
        d_PSR = 0.25 #kpc
        t_PSR = 342.*1000. #yr
        n_0 = 3.0 # cm^{-3} ambient density


    return RA_PSR_initial, Dec_PSR_initial, RA_PSR_final, Dec_PSR_final, d_PSR, t_PSR, n_0

def GetSupernovaParameters():

    RA_SN = 286.786
    Dec_SN= 6.498
    d_SN = 7.9 #kpc
    #t_SN = 20.*1000. # lower limit by Downes et al. (1980)
    t_SN = 11.*1000. # age of PSR J1907+0631 in year
    r_SN = 0.43/2. # degree, Yang et al. 2006, ChJAA, 6, 210.
    n_0 = 3.0 # cm^{-3} ambient density

    if sys.argv[1]=='MGRO_J1908_ON':

        if J1908_model=='SNR_7p9kpc_PSR_moving' or J1908_model=='SNR_7p9kpc_PSR_stationary':
            #RA_SN = 286.786
            #Dec_SN= 6.498
            #d_SN = 8.7 #kpc
            #t_SN = 11.*1000. # age of PSR J1907+0631 in year
            #r_SN = 0.5*0.43 # degree, Yang et al. 2006, ChJAA, 6, 210.
            #n_0 = 3.0 # cm^{-3} ambient density
            RA_SN = 286.9 # G41.1-0.3
            Dec_SN= 7.1
            d_SN = 8.5 #kpc
            t_SN = 5.*1000. 
            r_SN = 0.5*0.43 
            n_0 = 3.0 # cm^{-3} ambient density

        if J1908_model=='SNR_3p4kpc_PSR_moving' or J1908_model=='SNR_3p4kpc_PSR_stationary':
            RA_SN = 286.786
            Dec_SN= 6.498
            d_SN = 3.4 #kpc, Yang et al. 2006
            t_SN = 11.*1000.
            r_SN = 0.5*0.43 # degree, Yang et al. 2006, ChJAA, 6, 210.
            n_0 = 3.0 # cm^{-3} ambient density

        if J1908_model=='plerion':
            RA_SN = 287.2
            Dec_SN= 6.6
            d_SN = 1.884 #kpc
            #d_SN = 3.2 #kpc
            t_SN = 19.5*1000.
            r_SN = 0.5*1.0 
            n_0 = 3.0 # cm^{-3} ambient density

    elif sys.argv[1]=='GammaCygni_ON':

        RA_SN = 305.208
        Dec_SN= 40.43
        d_SN = 1.7 #kpc
        t_SN = 7.*1000. # Mavromatakis (2003)
        r_SN = 0.5*1.0 # Yang et al. 2006, ChJAA, 6, 210.
        n_0 = 3.0 # cm^{-3} ambient density

    elif sys.argv[1]=='IC443HotSpot_ON':

        RA_SN = 94.213
        Dec_SN= 22.503
        d_SN = 1.9 #kpc
        t_SN = 30.*1000. # year
        #t_SN = 3.*1000. # year
        r_SN = 0.5*45./60.
        n_0 = 3. # cm^{-3} ambient density

    return RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0

def DeriveSupernovaDiffusionLength(esc_param,t_SN,d_SN,n_0,gamma_energy):

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit

    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
    E_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit

    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    E_vis = (1./0.17)*2.*gamma_energy  # visible CR energy threshold
    t_esc = t_Sedov*pow(E_vis/E_max,-1./esc_param)

    diffusion_ism = 1e28*pow(E_vis/10.,0.5)*pow(B_ISM/3.,-0.5) # cm2/s
    t_diff = max(0.,t_SN-t_esc) # year
    diffusion_radius = pow(4.*diffusion_ism*t_diff*365.*24.*60.*60.,0.5)/deg_to_cm
    diffusion_radius = max(smooth_size_spectroscopy,diffusion_radius)

    E_min = pow(t_SN/t_Sedov,-esc_param)*E_max
    total_energy_in_CR = 1./(Gamma_CR-2.)*(pow(1.,-Gamma_CR+2.)-pow(E_max,-Gamma_CR+2.))
    visible_energy_in_CR = 1./(Gamma_CR-2.)*(pow(E_vis,-Gamma_CR+2.)-pow(E_max,-Gamma_CR+2.))
    visible_energy_frac = visible_energy_in_CR/total_energy_in_CR

    return diffusion_radius

def DeriveSupernovaEscapeParameter(diffusion_radius,t_SN,d_SN,n_0):

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit

    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
    E_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit

    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    E_vis = (1./0.17)*2.*current_gamma_energy  # visible CR energy threshold

    diffusion_ism = 1e28*pow(E_vis/10.,0.5)*pow(B_ISM/3.,-0.5) # cm2/s

    t_diff = pow(diffusion_radius*deg_to_cm,2)/(4.*diffusion_ism*365.*24.*60.*60.)
    t_esc = t_SN-t_diff
    if t_esc>0.:
        esc_param = -1.*math.log(E_vis/E_max)/math.log(t_esc/t_Sedov)
        R_esc = 2.5*vel_init*t_Sedov*365.*24.*60.*60.*(pow(t_esc/t_Sedov,0.4)-0.6)/deg_to_cm
        R_esc = max(smooth_size_spectroscopy,R_esc)
    else:
        print ("!!! CR diffusion time %0.3f kyr is longer than SNR age %0.3f kyr !!!"%(t_diff/1e3,t_SN/1e3))
        esc_param = 0.
        R_esc = 0.

    return esc_param

def DeriveSupernovaParameters(CR_efficiency, esc_param):

    RA_SN, Dec_SN, r_SN, d_SN, t_SN, n_0 = GetSupernovaParameters()

    Gamma_CR = 2.1
    f_Gamma = 0.9
    #Gamma_CR = 2.2
    #f_Gamma = 0.43
    #Gamma_CR = 2.3
    #f_Gamma = 0.19

    E_SN = 1. # 10^{51} erg
    B_SNR = 100. # micro G
    B_ISM = 6. # micro G
    M_ej = 1. # ejecta mass in solar mass unit

    t_Sedov = 423.*pow(E_SN,-0.5)*pow(M_ej,5./6.)*pow(n_0,-1./3.) # year
    vel_init = 7090.*1e5*pow(E_SN,-0.5)*pow(M_ej,-0.5) # cm/s 
    print ('t_Sedov = %0.2f year'%(t_Sedov))
    print ('vel_init = %0.2f km/s'%(vel_init/(100.*1000.)))
    E_max = 1e6*vel_init*vel_init*t_Sedov*365.*24.*60.*60./(3.4*1e28)*B_SNR # GeV, Bohm limit
    print ('E_max = %0.2f TeV'%(E_max/1000.))

    deg_to_cm = 3.14/180.*d_SN*1000.*3.086e18
    # E_pion = k*E_proton, k = 0.17 from https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA9.pdf
    E_vis = (1./0.17)*2.*current_gamma_energy  # visible CR energy threshold
    print ('E_gamma = %0.2f TeV'%(current_gamma_energy/1000.))
    print ('E_CR = %0.2f TeV'%(E_vis/1000.))
    diffusion_ism = 1e28*pow(E_vis/10.,0.5)*pow(B_ISM/3.,-0.5) # cm2/s

    if esc_param>0.:
        t_esc = t_Sedov*pow(E_vis/E_max,-1./esc_param)
        print ('t_esc = %0.2f kyr'%(t_esc/1000.))
        t_diff = t_SN-t_esc
        print ('t_diff = %0.2f kyr'%(t_diff/1000.))
        diffusion_radius = pow(t_diff*(4.*diffusion_ism*365.*24.*60.*60.),0.5)/deg_to_cm
    else:
        print ("!!! CR diffusion time is longer than SNR age !!!")
        diffusion_radius = 0.1

    E_min = pow(t_SN/t_Sedov,-esc_param)*E_max
    print ('E_min = %0.2f TeV'%(E_min/1000.))
    total_energy_in_CR = 1./(Gamma_CR-2.)*(pow(1.,-Gamma_CR+2.)-pow(E_max,-Gamma_CR+2.))
    visible_energy_in_CR = 1./(Gamma_CR-2.)*(pow(E_vis,-Gamma_CR+2.)-pow(E_max,-Gamma_CR+2.))
    visible_energy_frac = visible_energy_in_CR/total_energy_in_CR
    print ('visible_energy_frac = %0.5f'%(visible_energy_frac))

    R_SN = 2.5*vel_init*t_Sedov*365.*24.*60.*60.*(pow(t_SN/t_Sedov,0.4)-0.6)/deg_to_cm
    print ('expected R_SN = %0.2f deg'%(R_SN))

def tick_function_deg_to_pc(theta_deg,distance_kpc):
    radius_pc = theta_deg*3.14/180.*distance_kpc*1000.
    return ["%.2f" % z for z in radius_pc]

def FindProjection(Hist_Data_input,Hist_Syst_input,roi_x1,roi_y1,roi_x2,roi_y2,roi_d,roi_width,isTransverse,fraction):

    return CommonPlotFunctions.FindProjection(Hist_Data_input,Hist_Syst_input,roi_x1,roi_y1,roi_x2,roi_y2,roi_d,roi_width,isTransverse,fraction)

def MakeIntegratedSignificance(Hist_Data_input,Hist_Bkgd_input,Hist_Syst_input,roi_x,roi_y):

    Hist_Zscore_Theta2 = ROOT.TH1D("Hist_Zscore_Theta2","",20,0,2.)
    for br in range(0,Hist_Zscore_Theta2.GetNbinsX()):
        range_limit = Hist_Zscore_Theta2.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Zscore_Theta2.GetBinLowEdge(br+1)
        total_data = 0.
        total_bkgd = 0.
        total_syst = 0.
        total_stat = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                bkgd_content = Hist_Bkgd_input.GetBinContent(bx+1,by+1)
                bkgd_error = Hist_Bkgd_input.GetBinError(bx+1,by+1)
                if distance_sq>pow(range_limit,2): continue
                xy_bin_up = 1
                bin_size_up = integration_radii[1]
                bin_size_low = integration_radii[0]
                for xy_bin in range(1,len(integration_radii)):
                    current_bin_size = integration_radii[xy_bin]
                    if current_bin_size>range_limit:
                        xy_bin_up = xy_bin
                        bin_size_up = integration_radii[xy_bin]
                        bin_size_low = integration_radii[xy_bin-1]
                        break
                syst_up = Hist_Syst_input[xy_bin_up].GetBinContent(bx+1,by+1)
                syst_low = Hist_Syst_input[xy_bin_up-1].GetBinContent(bx+1,by+1)
                syst_content = syst_low + (syst_up-syst_low)/(bin_size_up-bin_size_low)*(range_limit-bin_size_low)
                total_data += data_content
                total_bkgd += bkgd_content
                total_syst += syst_content
                total_stat += bkgd_error*bkgd_error
        total_stat = pow(total_stat,0.5)
        zscore = CalculateSignificance(total_data-total_bkgd,total_bkgd,pow(total_stat*total_stat+total_syst*total_syst,0.5))
        Hist_Zscore_Theta2.SetBinContent(br+1,zscore)

    zscore = []
    theta2 = []
    for binx in range(0,Hist_Zscore_Theta2.GetNbinsX()):
        center = Hist_Zscore_Theta2.GetBinCenter(binx+1)
        content = Hist_Zscore_Theta2.GetBinContent(binx+1)
        theta2 += [center]
        zscore += [content]

    return zscore, theta2

def MakeSignificanceDistribution(Hist_Data_input,Hist_Bkgd_input,Hist_Syst_input):

    Hist_Zscore = GetSignificanceMap(Hist_Data_input,Hist_Bkgd_input,Hist_Syst_input,False)
    Hist_Zscore_Dist = ROOT.TH1D("Hist_Zscore_Dist","",65,-5,8)
    for bx in range(0,Hist_Zscore.GetNbinsX()):
        for by in range(0,Hist_Zscore.GetNbinsY()):
            if Hist_Data_input.GetBinContent(bx+1,by+1)==0: continue
            content = Hist_Zscore.GetBinContent(bx+1,by+1)
            Hist_Zscore_Dist.Fill(content)

    zscore = []
    count = []
    count_err = []
    for binx in range(0,Hist_Zscore_Dist.GetNbinsX()):
        center = Hist_Zscore_Dist.GetBinCenter(binx+1)
        content = Hist_Zscore_Dist.GetBinContent(binx+1)
        zscore += [center]
        count += [content]
        count_err += [max(0.9,pow(content,0.5))]

    return count, count_err, zscore

def GetCRcounts(name):

    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
    for bx in range(1,Hist2D_OnData_Sum.GetNbinsX()+1):
        for by in range(1,Hist2D_OnData_Sum.GetNbinsY()+1):
            if bx>bin_upper_x and by<=bin_upper_y:
                counts = Hist2D_OnData_Sum.GetBinContent(bx,by)
                counts = max(counts,1.)
                Hist_CR_Counts_MSCL.Fill(min(5,math.log10(counts)))
            if bx<=bin_upper_x and by>bin_upper_y:
                counts = Hist2D_OnData_Sum.GetBinContent(bx,by)
                counts = max(counts,1.)
                Hist_CR_Counts_MSCW.Fill(min(5,math.log10(counts)))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_CR_Counts_MSCW]
    legends += ['CR MSCW']
    colors += [1]
    Hists += [Hist_CR_Counts_MSCL]
    legends += ['CR MSCL']
    colors += [2]
    MakeMultiplePlot(Hists,legends,colors,'log10 counts per bin','N bins','CR_Counts_%s'%(name),0,0,False,False)

def SingleSourceAnalysis(source_list,e_low,e_up):

    global ErecS_lower_cut
    global ErecS_upper_cut
    global Syst_MDM
    global energy_bin_cut_low
    global energy_bin_cut_up
    global selection_tag
    global max_chi2_diff2_position_this_energy

    energy_bin_cut_low = e_low
    energy_bin_cut_up = e_up
    selection_tag = root_file_tags[0]
    selection_tag += '_%s'%(folder_path)
    selection_tag += '_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)

    FilePath_List = []
    ResetStackedShowerHistograms()
    for source in range(0,len(source_list)):
        source_name = source_list[source]
        for elev in range(0,len(root_file_tags)):
            for x_idx in range(0,Skymap_nzones_x):
                for y_idx in range(0,Skymap_nzones_y):
                    n_groups = 0
                    file_exists = True
                    while file_exists:
                        SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(n_groups,x_idx,y_idx)+".root"
                        print ('Read file: %s'%(SourceFilePath))
                        if os.path.exists(SourceFilePath):
                            n_groups += 1
                            print ('file exists.')
                        else:
                            file_exists = False
                            print ('file does not exist.')
                    for g_idx in range(0,n_groups):
                        FilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_"%(folder_path)+source_list[source]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(g_idx,x_idx,y_idx)+".root"
                        GetBrightStarInfo(FilePath)
                        print ('len(bright_star_ra) = %s'%(len(bright_star_ra)))
                        FilePath_List += [FilePath]
                        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
                        print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
                        InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
                        InfoTree = InputFile.Get("InfoTree")
                        InfoTree.GetEntry(0)

                        MSCW_blind_cut = InfoTree.MSCW_cut_blind
                        MSCL_blind_cut = InfoTree.MSCL_cut_blind
                        bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
                        bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
                        bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
                        bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
                        for e in range(0,len(energy_bin)-1):
                            print ('energy_bin = %s'%(energy_bin))
                            max_chi2_diff2_position_this_energy = max_chi2_diff2_position[e]
                            ErecS_lower_cut = energy_bin[e]
                            ErecS_upper_cut = energy_bin[e+1]
                            print ('max_chi2_diff2_position_this_energy = %0.3f'%(max_chi2_diff2_position_this_energy))
                            if ErecS_upper_cut<=energy_bin[energy_bin_cut_low]: continue
                            if ErecS_lower_cut>=energy_bin[energy_bin_cut_up]: continue
                            Syst_MDM = energy_syst[e]
                            GetShowerHistogramsFromFile(FilePath_List[len(FilePath_List)-1])
                            StackShowerHistograms()
                            NormalizeEnergyHistograms(FilePath_List[len(FilePath_List)-1])
                            StackEnergyHistograms()
                            NormalizeTheta2Histograms(FilePath_List[len(FilePath_List)-1],e)
                            StackTheta2Histograms()
                            NormalizeSkyMapHistograms(FilePath_List[len(FilePath_List)-1],e)
                            StackSkymapHistograms(e)

                        InputFile.Close()

    for ebin in range(0,len(energy_bin)-1):
        for binx in range(0,Hist_SumSyst_Energy_Skymap[ebin].GetNbinsX()):
            for biny in range(0,Hist_SumSyst_Energy_Skymap[ebin].GetNbinsY()):
                norm_syst = Hist_NormSyst_Energy_Skymap[ebin].GetBinContent(binx+1,biny+1)
                xy_bin_up = 1
                bin_size_up = integration_radii[1]
                bin_size_low = integration_radii[0]
                for xy_bin in range(1,len(integration_radii)):
                    current_bin_size = integration_radii[xy_bin]
                    if current_bin_size>smooth_size_spectroscopy:
                        xy_bin_up = xy_bin
                        bin_size_up = integration_radii[xy_bin]
                        bin_size_low = integration_radii[xy_bin-1]
                        break
                shape_syst_up = Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin_up].GetBinContent(binx+1,biny+1)
                shape_syst_low = Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin_up-1].GetBinContent(binx+1,biny+1)
                shape_syst = shape_syst_low + (shape_syst_up-shape_syst_low)/(bin_size_up-bin_size_low)*(smooth_size_spectroscopy-bin_size_low)
                Hist_SumSyst_Energy_Skymap[ebin].SetBinContent(binx+1,biny+1,pow(norm_syst*norm_syst+shape_syst*shape_syst,0.5))


    slice_center_x = roi_ra[0]
    slice_center_y = roi_dec[0]
    slice_radius = roi_radius_outer[0]
    first_bin_x = 0 
    last_bin_x = -1
    first_bin_y = 0 
    last_bin_y = -1
    first_bin_x = Hist_OnData_Skymap_Sum.GetXaxis().FindBin(slice_center_x-slice_radius)
    last_bin_x = Hist_OnData_Skymap_Sum.GetXaxis().FindBin(slice_center_x+slice_radius)
    first_bin_y = Hist_OnData_Skymap_Sum.GetYaxis().FindBin(slice_center_y-slice_radius)
    last_bin_y = Hist_OnData_Skymap_Sum.GetYaxis().FindBin(slice_center_y+slice_radius)
    if 'OFF' in ONOFF_tag:
        first_bin_x = 0 
        last_bin_x = -1
        first_bin_y = 0 
        last_bin_y = -1
    Hist_OnData_Skymap_ProjX_Sum.Add(Hist_OnData_Skymap_Sum.ProjectionX("",first_bin_x,last_bin_x,"eo"))
    Hist_OnDark_Skymap_ProjX_Sum.Add(Hist_OnDark_Skymap_Sum.ProjectionX("",first_bin_x,last_bin_x,"eo"))
    Hist_OnBkgd_Skymap_ProjX_Sum.Add(Hist_OnBkgd_Skymap_Sum.ProjectionX("",first_bin_x,last_bin_x,"eo"))
    Hist_OnData_Skymap_ProjY_Sum.Add(Hist_OnData_Skymap_Sum.ProjectionY("",first_bin_y,last_bin_y,"eo"))
    Hist_OnDark_Skymap_ProjY_Sum.Add(Hist_OnDark_Skymap_Sum.ProjectionY("",first_bin_y,last_bin_y,"eo"))
    Hist_OnBkgd_Skymap_ProjY_Sum.Add(Hist_OnBkgd_Skymap_Sum.ProjectionY("",first_bin_y,last_bin_y,"eo"))
    Hist_OnBkgd_Skymap_Syst_ProjX_Sum.Add(Hist_OnBkgd_Skymap_Syst_Sum.ProjectionX("",first_bin_x,last_bin_x,"eo"))
    Hist_OnBkgd_Skymap_Syst_ProjY_Sum.Add(Hist_OnBkgd_Skymap_Syst_Sum.ProjectionY("",first_bin_x,last_bin_x,"eo"))
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_ProjX_Sum.GetNbinsX()):
        bin_content = Hist_OnBkgd_Skymap_Syst_ProjX_Sum.GetBinContent(binx+1)
        Hist_OnBkgd_Skymap_Syst_ProjX_Sum.SetBinContent(binx+1,0.)
        Hist_OnBkgd_Skymap_Syst_ProjX_Sum.SetBinError(binx+1,bin_content)
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_ProjY_Sum.GetNbinsX()):
        bin_content = Hist_OnBkgd_Skymap_Syst_ProjY_Sum.GetBinContent(binx+1)
        Hist_OnBkgd_Skymap_Syst_ProjY_Sum.SetBinContent(binx+1,0.)
        Hist_OnBkgd_Skymap_Syst_ProjY_Sum.SetBinError(binx+1,bin_content)

    ErecS_lower_cut = energy_bin[energy_bin_cut_low]
    ErecS_upper_cut = energy_bin[energy_bin_cut_up]

    source_name = sys.argv[1]
    MatrixDecompositionDemo(source_name)
    GetCRcounts(source_name)

    GetSourceInfo(FilePath_List)

    Syst_MDM = energy_syst[energy_bin_cut_low]

    #############################################################
    #PlotsStackedHistograms('%s%s'%(source_list[0],selection_tag))
    print ('finish stacked plots.')
    print ('selection_tag = %s'%(selection_tag))
    #############################################################

    #MakeRankResidualPlots('%s%s'%(source_list[0],PercentCrab))

    if not drawMap: return

    event_rate = Hist_OnBkgd_Energy_CamCenter_Sum.Integral()/exposure_hours*(smooth_size_spectroscopy*smooth_size_spectroscopy)/(3.14*1.0*1.0)

    Hist_IndexMap = MakeSpectrumIndexSkymap(exposure_hours,Hist_Data_Energy_Skymap,Hist_Bkgd_Energy_Skymap,Hist_Rfov_Energy_Skymap,Hist_NormSyst_Energy_Skymap,Hist_SumSyst_Energy_Skymap,Hist_Expo_Energy_Skymap,Hist_EffArea_Energy_Skymap,'RA','Dec','%s%s_RaDec'%(source_name,PercentCrab),skymap_zoomin_scale)

    fig.clf()
    axbig = fig.add_subplot()
    cycol = cycle('krgbcmy')
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        if Hist_Data_Energy_Skymap[ebin].Integral()<100.: continue
        next_color = next(cycol)
        list_maxsig, list_binsize = VariableSkymapBins(Hist_Data_Energy_Skymap[ebin],Hist_Bkgd_Energy_Skymap[ebin],Hist_SumSyst_Energy_Skymap[ebin])
        axbig.plot(list_binsize, list_maxsig,color=next_color,label='%s-%s GeV'%(energy_bin[ebin],energy_bin[ebin+1]))
    axbig.set_xlabel('bin size [degree]')
    axbig.set_ylabel('max z score')
    axbig.legend(loc='best')
    plotname = 'MaxZscore_wSyst_%s'%(source_list[0])
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    axbig.remove()

    Hist_Syst_Energy_Skymap_Sum = []
    for xy_bin in range(0,len(integration_radii)):
        Hist_Syst_Energy_Skymap_Sum += [ROOT.TH2D("Hist_Syst_Energy_Skymap_Sum_%s"%(xy_bin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_Data_Energy_Skymap_Sum = ROOT.TH2D("Hist_Data_Energy_Skymap_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
    Hist_Bkgd_Energy_Skymap_Sum = ROOT.TH2D("Hist_Bkgd_Energy_Skymap_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
    for ebin in range(0,len(energy_bin)-1):
        Hist_Data_Energy_Skymap_Sum.Add(Hist_Data_Energy_Skymap[ebin])
        Hist_Bkgd_Energy_Skymap_Sum.Add(Hist_Bkgd_Energy_Skymap[ebin])
        for xy_bin in range(0,len(integration_radii)):
            for bx in range(0,Hist_Syst_Energy_Skymap_Sum[xy_bin].GetNbinsX()):
                for by in range(0,Hist_Syst_Energy_Skymap_Sum[xy_bin].GetNbinsY()):
                    norm_syst = Hist_NormSyst_Energy_Skymap[ebin].GetBinContent(bx+1,by+1)
                    shape_syst = Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin].GetBinContent(bx+1,by+1)
                    syst_new_content = pow(norm_syst*norm_syst+shape_syst*shape_syst,0.5)
                    syst_old_content = Hist_Syst_Energy_Skymap_Sum[xy_bin].GetBinContent(bx+1,by+1)
                    Hist_Syst_Energy_Skymap_Sum[xy_bin].SetBinContent(bx+1,by+1,pow(syst_new_content*syst_new_content+syst_old_content*syst_old_content,0.5))

    fig.clf()
    axbig = fig.add_subplot()
    cycol = cycle('cmy')
    #for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    #    if Hist_Data_Energy_Skymap[ebin].Integral()<100.: continue
    #    next_color = next(cycol)
    #    count, count_err, zscore = MakeSignificanceDistribution(Hist_Data_Energy_Skymap[ebin],Hist_Bkgd_Energy_Skymap[ebin],Hist_SumSyst_Energy_Skymap[ebin])
    #    axbig.errorbar(zscore, count, count_err,color=next_color,marker='_',ls='none',label='%s-%s GeV'%(energy_bin[ebin],energy_bin[ebin+1]))
    count, count_err, zscore = MakeSignificanceDistribution(Hist_Data_Energy_Skymap_Sum,Hist_Bkgd_Energy_Skymap_Sum,Hist_Syst_Energy_Skymap_Sum[0])
    axbig.errorbar(zscore, count, count_err,color='k',marker='_',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up-1]))
    start = (100., 0., 1.)
    popt, pcov = curve_fit(gaussian_func, np.array(zscore), np.array(count), p0=start, sigma=np.array(count_err))
    axbig.plot(np.array(zscore), gaussian_func(np.array(zscore), *popt),color='r',label='$\mu = %0.2f$, $\sigma = %0.2f$'%(popt[1],popt[2]))
    axbig.set_xlabel('significance z score')
    axbig.set_ylabel('bin count')
    axbig.legend(loc='best')
    plt.ylim(0.9, 200)
    axbig.set_yscale('log')
    plotname = 'ZscoreDist_wSyst'
    fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
    axbig.remove()


    fig.clf()
    axbig = fig.add_subplot()
    cycol = cycle('krgbcmy')
    next_color = next(cycol)
    zscore, theta2 = MakeIntegratedSignificance(Hist_Data_Energy_Skymap_Sum,Hist_Bkgd_Energy_Skymap_Sum,Hist_Syst_Energy_Skymap_Sum,roi_ra[0],roi_dec[0])
    axbig.plot(theta2,zscore,color=next_color,label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for xy_bin in range(0,len(integration_radii)):
            Hist_Syst_Energy_Skymap_Sum[xy_bin].Reset()
        for xy_bin in range(0,len(integration_radii)):
            for bx in range(0,Hist_Syst_Energy_Skymap_Sum[xy_bin].GetNbinsX()):
                for by in range(0,Hist_Syst_Energy_Skymap_Sum[xy_bin].GetNbinsY()):
                    norm_syst = Hist_NormSyst_Energy_Skymap[ebin].GetBinContent(bx+1,by+1)
                    shape_syst = Hist_ShapeSyst_Energy_Skymap[ebin][xy_bin].GetBinContent(bx+1,by+1)
                    syst_new_content = pow(norm_syst*norm_syst+shape_syst*shape_syst,0.5)
                    syst_old_content = Hist_Syst_Energy_Skymap_Sum[xy_bin].GetBinContent(bx+1,by+1)
                    Hist_Syst_Energy_Skymap_Sum[xy_bin].SetBinContent(bx+1,by+1,pow(syst_new_content*syst_new_content+syst_old_content*syst_old_content,0.5))
        next_color = next(cycol)
        zscore, theta2 = MakeIntegratedSignificance(Hist_Data_Energy_Skymap[ebin],Hist_Bkgd_Energy_Skymap[ebin],Hist_Syst_Energy_Skymap_Sum,roi_ra[0],roi_dec[0])
        axbig.plot(theta2,zscore,color=next_color,label='%s-%s GeV'%(energy_bin[ebin],energy_bin[ebin+1]))
    axbig.set_ylabel('significance z score')
    axbig.set_xlabel('radius of integration region [degree]')
    plt.ylim(-5.0, 12.0)
    axbig.legend(loc='best')
    plotname = 'ZscoreVsTheta2_%s'%(source_list[0])
    fig.savefig("output_plots/%s_%s.png"%(plotname,selection_tag),bbox_inches='tight')
    axbig.remove()


def FindSourceIndex(source_name):
    for source in range(0,len(sample_list)):
        if source_name==sample_list[source]:
            return source
    return 0

Hist_NSB = ROOT.TH1D("Hist_NSB","",20,0,8)
Hist_L3Rate = ROOT.TH1D("Hist_L3Rate","",30,0,600)
Hist_EffArea = ROOT.TH1D("Hist_EffArea","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_EffArea_Sum = ROOT.TH1D("Hist_EffArea_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
print ('len(sample_list) = %s'%(len(sample_list)))
print ('len(root_file_tags) = %s'%(len(root_file_tags)))
print ('Skymap_nzones_x = %s'%(Skymap_nzones_x))
print ('Skymap_nzones_y = %s'%(Skymap_nzones_y))
for source in range(0,len(sample_list)):
    source_idx = FindSourceIndex(sample_list[source])
    FilePath_Folder = []
    for elev in range(0,len(root_file_tags)):
        for x_idx in range(0,Skymap_nzones_x):
            for y_idx in range(0,Skymap_nzones_y):
                n_groups = 0
                file_exists = True
                while file_exists:
                    SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(n_groups,x_idx,y_idx)+".root"
                    print ('Read file: %s'%(SourceFilePath))
                    if os.path.exists(SourceFilePath):
                        n_groups += 1
                        print ('file exists.')
                    else:
                        file_exists = False
                        print ('file does not exist.')
                for g_idx in range(0,n_groups):
                    SourceFilePath = "/gamma_raid/userspace/rshang/SMI_output/%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(g_idx,x_idx,y_idx)+".root"
                    FilePath_Folder += [SourceFilePath]
                    print ('Get %s...'%(FilePath_Folder[len(FilePath_Folder)-1]))
                    if not os.path.isfile(FilePath_Folder[len(FilePath_Folder)-1]): 
                        print ('Found no file!!')
                        continue
                    else:
                        print ('Found a file.')
    GetSourceInfo(FilePath_Folder)
    MakeOneHistPlot(Hist_NSB,'NSB','number of runs','NSB_%s'%(sample_list[source]),False)
    MakeOneHistPlot(Hist_L3Rate,'L3 rate','number of runs','L3Rate_%s'%(sample_list[source]),False)

print ('analysis cut: MSCL = %s, MSCW = %s'%(MSCL_blind_cut,MSCW_blind_cut))
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))+MSCL_blind_cut
MSCW_plot_lower = -gamma_hadron_low_end*gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))-MSCW_blind_cut
MSCL_plot_lower = -gamma_hadron_low_end*gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))-MSCL_blind_cut
#MSCW_plot_lower = -0.5
#MSCL_plot_lower = -0.5
print ('plot range: MSCL = %s, MSCW = %s'%(MSCL_plot_upper,MSCW_plot_upper))
print ('N_bins_for_deconv = %d'%(N_bins_for_deconv))

Hist2D_OnData_Point_Sum = ROOT.TH2D("Hist2D_OnData_Point_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData_Ring_Sum = ROOT.TH2D("Hist2D_OnData_Ring_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData_Sum = ROOT.TH2D("Hist2D_OnData_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd_Sum = ROOT.TH2D("Hist2D_OnBkgd_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnGamma_Sum = ROOT.TH2D("Hist2D_OnGamma_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark_Sum = ROOT.TH2D("Hist2D_OnDark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData_Point = ROOT.TH2D("Hist2D_OnData_Point","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData_Ring = ROOT.TH2D("Hist2D_OnData_Ring","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData = ROOT.TH2D("Hist2D_OnData","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd = ROOT.TH2D("Hist2D_OnBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnGamma = ROOT.TH2D("Hist2D_OnGamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark = ROOT.TH2D("Hist2D_OnDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCL_Sum = ROOT.TH1D("Hist_OnData_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL_Sum = ROOT.TH1D("Hist_OnDark_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCL = ROOT.TH1D("Hist_OnData_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL = ROOT.TH1D("Hist_OnBkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL = ROOT.TH1D("Hist_OnDark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCW_Sum = ROOT.TH1D("Hist_OnData_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnGamma_MSCW_Sum = ROOT.TH1D("Hist_OnGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnGamma_MSCL_Sum = ROOT.TH1D("Hist_OnGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCW_Sum = ROOT.TH1D("Hist_OnDark_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCW = ROOT.TH1D("Hist_OnData_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW = ROOT.TH1D("Hist_OnBkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCW = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCW = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnGamma_MSCW = ROOT.TH1D("Hist_OnGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCL = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCL = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnGamma_MSCL = ROOT.TH1D("Hist_OnGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCW = ROOT.TH1D("Hist_OnDark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist2D_U_Proj = ROOT.TH2D("Hist2D_U_Proj","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_V_Proj = ROOT.TH2D("Hist2D_V_Proj","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_U_Proj_Sum = ROOT.TH2D("Hist2D_U_Proj_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_V_Proj_Sum = ROOT.TH2D("Hist2D_V_Proj_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Data = ROOT.TH2D("Hist2D_Coeff_Data","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Bkgd = ROOT.TH2D("Hist2D_Coeff_Bkgd","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Data_Sum = ROOT.TH2D("Hist2D_Coeff_Data_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Bkgd_Sum = ROOT.TH2D("Hist2D_Coeff_Bkgd_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)

Hist2D_Rank0_Data = ROOT.TH2D("Hist2D_Rank0_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Data = ROOT.TH2D("Hist2D_Rank1_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Data = ROOT.TH2D("Hist2D_Rank2_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Data = ROOT.TH2D("Hist2D_Rank3_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Data = ROOT.TH2D("Hist2D_Rank4_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Data_Sum = ROOT.TH2D("Hist2D_Rank0_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Data_Sum = ROOT.TH2D("Hist2D_Rank1_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Data_Sum = ROOT.TH2D("Hist2D_Rank2_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Data_Sum = ROOT.TH2D("Hist2D_Rank3_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Data_Sum = ROOT.TH2D("Hist2D_Rank4_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Dark = ROOT.TH2D("Hist2D_Rank0_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Dark = ROOT.TH2D("Hist2D_Rank1_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Dark = ROOT.TH2D("Hist2D_Rank2_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Dark = ROOT.TH2D("Hist2D_Rank3_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Dark = ROOT.TH2D("Hist2D_Rank4_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Dark_Sum = ROOT.TH2D("Hist2D_Rank0_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Dark_Sum = ROOT.TH2D("Hist2D_Rank1_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Dark_Sum = ROOT.TH2D("Hist2D_Rank2_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Dark_Sum = ROOT.TH2D("Hist2D_Rank3_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Dark_Sum = ROOT.TH2D("Hist2D_Rank4_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

n_iterations = 100
optimiz_lower = -6.
optimiz_upper = -1.
Hist_Bkgd_Chi2 = ROOT.TH1D("Hist_Bkgd_Chi2","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Chi2_Diff = ROOT.TH1D("Hist_Bkgd_Chi2_Diff","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Chi2_Diff2 = ROOT.TH1D("Hist_Bkgd_Chi2_Diff2","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Converge_Blind = ROOT.TH1D("Hist_Bkgd_Converge_Blind","",n_iterations,0,n_iterations)
Hist_Bkgd_Converge_Unblind = ROOT.TH1D("Hist_Bkgd_Converge_Unblind","",n_iterations,0,n_iterations)
Hist_VVV_Eigenvalues = ROOT.TH1D("Hist_VVV_Eigenvalues","",N_bins_for_deconv*N_bins_for_deconv,0,N_bins_for_deconv*N_bins_for_deconv)

Hist1D_Data_Rank0_LeftVector = ROOT.TH1D("Hist1D_Data_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank1_LeftVector = ROOT.TH1D("Hist1D_Data_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank2_LeftVector = ROOT.TH1D("Hist1D_Data_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank3_LeftVector = ROOT.TH1D("Hist1D_Data_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank0_RightVector = ROOT.TH1D("Hist1D_Data_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank1_RightVector = ROOT.TH1D("Hist1D_Data_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank2_RightVector = ROOT.TH1D("Hist1D_Data_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank3_RightVector = ROOT.TH1D("Hist1D_Data_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank0_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank1_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank2_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank3_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank0_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank1_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank2_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank3_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank0_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank1_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank2_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank3_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank0_RightVector = ROOT.TH1D("Hist1D_Dark_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank1_RightVector = ROOT.TH1D("Hist1D_Dark_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank2_RightVector = ROOT.TH1D("Hist1D_Dark_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank3_RightVector = ROOT.TH1D("Hist1D_Dark_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)

Hist_OnData_SR_XYoff_Sum = ROOT.TH2D("Hist_OnData_SR_XYoff_Sum","",12,-3,3,12,-3,3)
Hist_OnData_SR_XYoff = ROOT.TH2D("Hist_OnData_SR_XYoff","",12,-3,3,12,-3,3)
Hist_OnDark_SR_XYoff_Sum = ROOT.TH2D("Hist_OnDark_SR_XYoff_Sum","",30,-3,3,30,-3,3)
Hist_OnDark_SR_XYoff = ROOT.TH2D("Hist_OnDark_SR_XYoff","",30,-3,3,30,-3,3)
Hist_OnData_CR_XYoff_Sum = ROOT.TH2D("Hist_OnData_CR_XYoff_Sum","",30,-3,3,30,-3,3)
Hist_OnData_CR_XYoff = ROOT.TH2D("Hist_OnData_CR_XYoff","",30,-3,3,30,-3,3)
Hist_OnData_R2off_Sum = ROOT.TH1D("Hist_OnData_R2off_Sum","",36,0,9)
Hist_OnData_ISR_R2off_Sum = ROOT.TH1D("Hist_OnData_ISR_R2off_Sum","",36,0,9)
Hist_OnData_Yoff_Sum = ROOT.TH1D("Hist_OnData_Yoff_Sum","",30,-3,3)
Hist_OnBkgd_R2off_Sum = ROOT.TH1D("Hist_OnBkgd_R2off_Sum","",36,0,9)
Hist_OnBkgd_Yoff_Sum = ROOT.TH1D("Hist_OnBkgd_Yoff_Sum","",30,-3,3)
Hist_OnBkgd_Yoff_Raw_Sum = ROOT.TH1D("Hist_OnBkgd_Yoff_Raw_Sum","",30,-3,3)
Hist_OnData_R2off = ROOT.TH1D("Hist_OnData_R2off","",36,0,9)
Hist_OnData_ISR_R2off = ROOT.TH1D("Hist_OnData_ISR_R2off","",36,0,9)
Hist_OnData_Yoff = ROOT.TH1D("Hist_OnData_Yoff","",30,-3,3)
Hist_OnBkgd_R2off = ROOT.TH1D("Hist_OnBkgd_R2off","",36,0,9)
Hist_OnBkgd_Yoff = ROOT.TH1D("Hist_OnBkgd_Yoff","",30,-3,3)
Hist_OnBkgd_Yoff_Raw = ROOT.TH1D("Hist_OnBkgd_Yoff_Raw","",30,-3,3)
Hist_OnData_Theta2_Sum = ROOT.TH1D("Hist_OnData_Theta2_Sum","",50,0,10)
Hist_OnBkgd_Theta2_Sum = ROOT.TH1D("Hist_OnBkgd_Theta2_Sum","",50,0,10)
Hist_OnRFoV_Theta2_Sum = ROOT.TH1D("Hist_OnRFoV_Theta2_Sum","",50,0,10)
Hist_OnDark_Theta2_Sum = ROOT.TH1D("Hist_OnDark_Theta2_Sum","",50,0,10)
Hist_OnData_Theta2 = ROOT.TH1D("Hist_OnData_Theta2","",50,0,10)
Hist_OnBkgd_Theta2 = ROOT.TH1D("Hist_OnBkgd_Theta2","",50,0,10)
Hist_OnRFoV_Theta2 = ROOT.TH1D("Hist_OnRFoV_Theta2","",50,0,10)
Hist_OnDark_Theta2 = ROOT.TH1D("Hist_OnDark_Theta2","",50,0,10)
Hist_OnData_Energy_Sum = ROOT.TH1D("Hist_OnData_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnRFoV_Energy_Sum = ROOT.TH1D("Hist_OnRFoV_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy_CamCenter_Sum = ROOT.TH1D("Hist_OnData_Energy_CamCenter_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_CamCenter_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_CamCenter_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnDark_Energy_Sum = ROOT.TH1D("Hist_OnDark_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy = ROOT.TH1D("Hist_OnData_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy = ROOT.TH1D("Hist_OnBkgd_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnRFoV_Energy = ROOT.TH1D("Hist_OnRFoV_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy_CamCenter = ROOT.TH1D("Hist_OnData_Energy_CamCenter","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_CamCenter = ROOT.TH1D("Hist_OnBkgd_Energy_CamCenter","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnDark_Energy = ROOT.TH1D("Hist_OnDark_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Zenith_Sum = ROOT.TH1D("Hist_OnData_Zenith_Sum","",45,0,90)
Hist_OnBkgd_Zenith_Sum = ROOT.TH1D("Hist_OnBkgd_Zenith_Sum","",45,0,90)
Hist_OnData_Zenith = ROOT.TH1D("Hist_OnData_Zenith","",45,0,90)
Hist_OnBkgd_Zenith = ROOT.TH1D("Hist_OnBkgd_Zenith","",45,0,90)
Hist_OnData_Height_Sum = ROOT.TH1D("Hist_OnData_Height_Sum","",25,0,25.)
Hist_OnBkgd_Height_Sum = ROOT.TH1D("Hist_OnBkgd_Height_Sum","",25,0,25.)
Hist_OnData_Height = ROOT.TH1D("Hist_OnData_Height","",25,0,25.)
Hist_OnBkgd_Height = ROOT.TH1D("Hist_OnBkgd_Height","",25,0,25.)
Hist_OnData_Depth_Sum = ROOT.TH1D("Hist_OnData_Depth_Sum","",25,0,25.)
Hist_OnBkgd_Depth_Sum = ROOT.TH1D("Hist_OnBkgd_Depth_Sum","",25,0,25.)
Hist_OnData_Depth = ROOT.TH1D("Hist_OnData_Depth","",25,0,25.)
Hist_OnBkgd_Depth = ROOT.TH1D("Hist_OnBkgd_Depth","",25,0,25.)
Hist_OnData_Rcore_Sum = ROOT.TH1D("Hist_OnData_Rcore_Sum","",25,0,500.)
Hist_OnBkgd_Rcore_Sum = ROOT.TH1D("Hist_OnBkgd_Rcore_Sum","",25,0,500.)
Hist_OnData_Rcore = ROOT.TH1D("Hist_OnData_Rcore","",25,0,500.)
Hist_OnBkgd_Rcore = ROOT.TH1D("Hist_OnBkgd_Rcore","",25,0,500.)

Hist_Data_Skymap = ROOT.TH2D("Hist_Data_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_Data_Elev_Skymap = ROOT.TH2D("Hist_Data_Elev_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_Data_Azim_Skymap = ROOT.TH2D("Hist_Data_Azim_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_Data_NSB_Skymap = ROOT.TH2D("Hist_Data_NSB_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_Data_MJD_Skymap = ROOT.TH2D("Hist_Data_MJD_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_Data_ElevNSB = ROOT.TH2D("Hist_Data_ElevNSB","",20,0,10,18,0,90)
Hist_Data_ElevAzim = ROOT.TH2D("Hist_Data_ElevAzim","",18,0,360,18,0,90)
Hist_Dark_ElevNSB = ROOT.TH2D("Hist_Dark_ElevNSB","",20,0,10,18,0,90)
Hist_Dark_ElevAzim = ROOT.TH2D("Hist_Dark_ElevAzim","",18,0,360,18,0,90)

Hist_SumSyst_Energy_Skymap = []
Hist_NormSyst_Energy_Skymap = []
Hist_ShapeSyst_Energy_Skymap = []
Hist_EffArea_Energy_Skymap = []
Hist_Expo_Energy_Skymap = []
Hist_Data_Energy_Skymap = []
Hist_Bkgd_Energy_Skymap = []
Hist_Rfov_Energy_Skymap = []
Hist_Data_CR_XYoff = []
for ebin in range(0,len(energy_bin)-1):
    Hist_SumSyst_Energy_Skymap += [ROOT.TH2D("Hist_SumSyst_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_NormSyst_Energy_Skymap += [ROOT.TH2D("Hist_NormSyst_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_EffArea_Energy_Skymap += [ROOT.TH2D("Hist_EffArea_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_Expo_Energy_Skymap += [ROOT.TH2D("Hist_Expo_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_Data_Energy_Skymap += [ROOT.TH2D("Hist_Data_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_Bkgd_Energy_Skymap += [ROOT.TH2D("Hist_Bkgd_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_Data_CR_XYoff += [ROOT.TH2D("Hist_OnData_CR_XYoff_%s"%(ebin),"",30,-3,3,30,-3,3)]
    Hist_Rfov_Energy_Skymap += [ROOT.TH2D("Hist_Rfov_Energy_Skymap_%s"%(ebin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_ShapeSyst_Energy_Skymap_ThisBin = []
    for xy_bin in range(0,len(integration_radii)):
        Hist_ShapeSyst_Energy_Skymap_ThisBin += [ROOT.TH2D("Hist_ShapeSyst_Energy_Skymap_%s_%s"%(ebin,xy_bin),"",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)]
    Hist_ShapeSyst_Energy_Skymap += [Hist_ShapeSyst_Energy_Skymap_ThisBin]

Hist_bkgd_global = ROOT.TH2D("Hist_bkgd_global","",int(Skymap_nbins_x/skymap_zoomin_scale),source_ra-Skymap_size_x/skymap_zoomin_scale,source_ra+Skymap_size_x/skymap_zoomin_scale,int(Skymap_nbins_y/skymap_zoomin_scale),source_dec-Skymap_size_y/skymap_zoomin_scale,source_dec+Skymap_size_y/skymap_zoomin_scale)
Hist_MWL_global = ROOT.TH2D("Hist_MWL_global","",int(Skymap_nbins_x/skymap_zoomin_scale),source_ra-Skymap_size_x/skymap_zoomin_scale,source_ra+Skymap_size_x/skymap_zoomin_scale,int(Skymap_nbins_y/skymap_zoomin_scale),source_dec-Skymap_size_y/skymap_zoomin_scale,source_dec+Skymap_size_y/skymap_zoomin_scale)

Hist_OnData_Skymap = ROOT.TH2D("Hist_OnData_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnDark_Skymap = ROOT.TH2D("Hist_OnDark_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap = ROOT.TH2D("Hist_OnBkgd_Skymap","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnData_Skymap_Sum = ROOT.TH2D("Hist_OnData_Skymap_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnDark_Skymap_Sum = ROOT.TH2D("Hist_OnDark_Skymap_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_Norm = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_Norm","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_Shape = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_Shape","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_RBM = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_RBM","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x,Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnData_Skymap_Galactic = ROOT.TH2D("Hist_OnData_Skymap_Galactic","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)
Hist_OnBkgd_Skymap_Galactic = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)
Hist_OnData_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnData_Skymap_Galactic_Sum","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)
Hist_OnBkgd_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Sum","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)
Hist_OnBkgd_Skymap_Galactic_Syst_MDM = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst_MDM","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)
Hist_OnBkgd_Skymap_Galactic_Syst_RBM = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst_RBM","",Skymap_nbins_x,source_l-Skymap_size_x,source_l+Skymap_size_x,Skymap_nbins_y,source_b-Skymap_size_y,source_b+Skymap_size_y)

Hist_OnData_Skymap_ProjX_Sum = ROOT.TH1D("Hist_OnData_Skymap_ProjX_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x)
Hist_OnDark_Skymap_ProjX_Sum = ROOT.TH1D("Hist_OnDark_Skymap_ProjX_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x)
Hist_OnBkgd_Skymap_ProjX_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_ProjX_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x)
Hist_OnBkgd_Skymap_Syst_ProjX = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjX","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x)
Hist_OnBkgd_Skymap_Syst_ProjX_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjX_Sum","",Skymap_nbins_x,source_ra-Skymap_size_x,source_ra+Skymap_size_x)
Hist_OnData_Skymap_ProjY_Sum = ROOT.TH1D("Hist_OnData_Skymap_ProjY_Sum","",Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnDark_Skymap_ProjY_Sum = ROOT.TH1D("Hist_OnDark_Skymap_ProjY_Sum","",Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_ProjY_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_ProjY_Sum","",Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_ProjY = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjY","",Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)
Hist_OnBkgd_Skymap_Syst_ProjY_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjY_Sum","",Skymap_nbins_y,source_dec-Skymap_size_y,source_dec+Skymap_size_y)

Hist_SystErr_MSCL = ROOT.TH1D("Hist_SystErr_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_SystErr_MSCW = ROOT.TH1D("Hist_SystErr_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_SystErr_Energy = ROOT.TH1D("Hist_SystErr_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_SystErr_Theta2_Sum = ROOT.TH1D("Hist_SystErr_Theta2_Sum","",50,0,10)
Hist_NormSyst_Theta2 = ROOT.TH1D("Hist_NormSyst_Theta2","",50,0,10)
Hist_ShapeSyst_Theta2 = ROOT.TH1D("Hist_ShapeSyst_Theta2","",50,0,10)

Hist_CR_Counts_MSCW = ROOT.TH1D("Hist_CR_Counts_MSCW","",6,0,6)
Hist_CR_Counts_MSCL = ROOT.TH1D("Hist_CR_Counts_MSCL","",6,0,6)

print ('MJD_Start = %s'%(MJD_Start))
print ('MJD_End = %s'%(MJD_End))
#time = Time(MJD_Start, format='mjd')
#time.format = 'decimalyear'
#year_start = time.value
#print ('MJD_End = %s'%(MJD_End))
#time = Time(MJD_End, format='mjd')
#time.format = 'decimalyear'
#year_end = time.value
#print ('year_start = %s'%(year_start))
#print ('year_end = %s'%(year_end))
print ('roi_ra = %s'%(roi_ra[0]))
print ('roi_dec = %s'%(roi_dec[0]))
Hist_OnData_RoI_Energy_Sum = [] 
Hist_OnBkgd_RoI_Energy_Sum = [] 
Hist_OnData_RoI_Energy = [] 
Hist_OnBkgd_RoI_Energy = [] 
Hist_SumSyst_RoI_Energy = [] 
Hist_NormSyst_RoI_Energy = [] 
Hist_ShapeSyst_RoI_Energy = [] 
Hist_OnData_RoI_Theta2_Sum = []
Hist_OnBkgd_RoI_Theta2_Sum = []
Hist_OnData_RoI_Theta2 = []
Hist_OnBkgd_RoI_Theta2 = []
Hist_SumSyst_RoI_Theta2 = [] 
Hist_NormSyst_RoI_Theta2 = [] 
Hist_ShapeSyst_RoI_Theta2 = [] 
Hist_OnData_RoI_X_Sum = []
Hist_OnBkgd_RoI_X_Sum = []
Hist_OnData_RoI_X = []
Hist_OnBkgd_RoI_X = []
Hist_OnData_RoI_Y_Sum = []
Hist_OnBkgd_RoI_Y_Sum = []
Hist_OnData_RoI_Y = []
Hist_OnBkgd_RoI_Y = []
Hist_OnData_RoI_MJD = [] 
Hist_OnBkgd_RoI_MJD = [] 
Hist_OnData_RoI_MJD_Sum = [] 
Hist_OnBkgd_RoI_MJD_Sum = [] 
Hist_RoI_Flux = []
for nth_roi in range(0,len(roi_ra)):
    roi_range = roi_radius_outer[nth_roi]
    Hist_OnData_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnData_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Energy += [ROOT.TH1D("Hist_OnData_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_SumSyst_RoI_Energy += [ROOT.TH1D("Hist_SumSyst_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_NormSyst_RoI_Energy += [ROOT.TH1D("Hist_NormSyst_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_ShapeSyst_RoI_Energy += [ROOT.TH1D("Hist_ShapeSyst_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnData_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_Theta2 += [ROOT.TH1D("Hist_OnData_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2 += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_SumSyst_RoI_Theta2 += [ROOT.TH1D("Hist_SumSyst_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_NormSyst_RoI_Theta2 += [ROOT.TH1D("Hist_NormSyst_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_ShapeSyst_RoI_Theta2 += [ROOT.TH1D("Hist_ShapeSyst_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_X_Sum += [ROOT.TH1D("Hist_OnData_RoI_X_Sum_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnBkgd_RoI_X_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_X_Sum_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnData_RoI_X += [ROOT.TH1D("Hist_OnData_RoI_X_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnBkgd_RoI_X += [ROOT.TH1D("Hist_OnBkgd_RoI_X_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnData_RoI_Y_Sum += [ROOT.TH1D("Hist_OnData_RoI_Y_Sum_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnBkgd_RoI_Y_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Y_Sum_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnData_RoI_Y += [ROOT.TH1D("Hist_OnData_RoI_Y_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnBkgd_RoI_Y += [ROOT.TH1D("Hist_OnBkgd_RoI_Y_%s"%(nth_roi),"",60,-roi_range,roi_range)]
    Hist_OnData_RoI_MJD += [ROOT.TH1D("Hist_OnData_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnData_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnData_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_RoI_Flux += [ROOT.TH1D("Hist_RoI_Flux_%s"%(nth_roi),"",len(energy_bin)-1,array('d',energy_bin))]


Hist2D_OffData = []
Hist2D_OffData_Sum = []
Hist_OffData_MSCL = []
Hist_OffData_MSCL_Sum = []
Hist_OffBkgd_MSCL = []
Hist_OffBkgd_MSCL_Sum = []
Hist_OffData_MSCW = []
Hist_OffData_MSCW_Sum = []
Hist_OffBkgd_MSCW = []
Hist_OffBkgd_MSCW_Sum = []
Hist_OffData_Energy = []
Hist_OffData_Energy_Sum = []
Hist_OffBkgd_Energy = []
Hist_OffBkgd_Energy_Sum = []
Hist_OffData_CameraFoV_Theta2 = []
Hist_OffData_CameraFoV_Theta2_Sum = []
Hist_OffBkgd_CameraFoV_Theta2 = []
Hist_OffBkgd_CameraFoV_Theta2_Sum = []
for nth_sample in range(0,n_control_samples):
    Hist2D_OffData += [ROOT.TH2D("Hist2D_OffData_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffData_Sum += [ROOT.TH2D("Hist2D_OffData_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCL += [ROOT.TH1D("Hist_OnData_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCL_Sum += [ROOT.TH1D("Hist_OnData_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCW += [ROOT.TH1D("Hist_OnData_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCW_Sum += [ROOT.TH1D("Hist_OnData_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCL += [ROOT.TH1D("Hist_OnBkgd_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCL_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCW += [ROOT.TH1D("Hist_OnBkgd_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCW_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_Energy += [ROOT.TH1D("Hist_OffData_Energy_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffData_Energy_Sum += [ROOT.TH1D("Hist_OffData_Energy_Sum_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffBkgd_Energy += [ROOT.TH1D("Hist_OffBkgd_Energy_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffBkgd_Energy_Sum += [ROOT.TH1D("Hist_OffBkgd_Energy_Sum_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffData_CameraFoV_Theta2 += [ROOT.TH1D("Hist_OffData_CameraFoV_Theta2_%s"%(nth_sample),"",50,0,10)]
    Hist_OffData_CameraFoV_Theta2_Sum += [ROOT.TH1D("Hist_OffData_CameraFoV_Theta2_Sum_%s"%(nth_sample),"",50,0,10)]
    Hist_OffBkgd_CameraFoV_Theta2 += [ROOT.TH1D("Hist_OffBkgd_CameraFoV_Theta2_%s"%(nth_sample),"",50,0,10)]
    Hist_OffBkgd_CameraFoV_Theta2_Sum += [ROOT.TH1D("Hist_OffBkgd_CameraFoV_Theta2_Sum_%s"%(nth_sample),"",50,0,10)]

GetGammaSourceInfo()

#SystematicAnalysis()

drawMap = True
doMWLMap = False
doUpperLimit = False
doReferenceFlux = True
doExtentFit = False

if sys.argv[1]=='MGRO_J1908_ON': doMWLMap = False
if sys.argv[1]=='GammaCygni_ON': doMWLMap = False
if sys.argv[1]=='IC443HotSpot_ON': doMWLMap = False

#J1908_model = 'SNR_7p9kpc_PSR_stationary'
#J1908_model = 'SNR_3p4kpc_PSR_stationary'
#J1908_model = 'SNR_7p9kpc_PSR_moving'
#J1908_model = 'SNR_3p4kpc_PSR_moving'
J1908_model = 'plerion'

#set_palette('default')
#set_palette('gray')

exclude_roi = []
#exclude_roi += ['HAWC region (VERITAS)']
exclude_roi += ['VHE region']
exclude_roi += ['G40.5-0.5']
exclude_roi += ['PSR tail 1']
exclude_roi += ['PSR tail 2']
exclude_roi += ['PSR region']
exclude_roi += ['Hot spot north']
exclude_roi += ['Hot spot west']
exclude_roi += ['Hot spot east']
exclude_roi += ['Ring 1']
exclude_roi += ['Ring 2']
exclude_roi += ['Ring 3']

#exclude_roi += ['Crab int. r=0.3 deg']
#exclude_roi += ['Crab int. r=0.5 deg']
#exclude_roi += ['Crab int. r=1.0 deg']
exclude_roi += ['Crab int. r=1.5 deg']

exclude_roi += ['Control region r=0.5 deg']
exclude_roi += ['Control region r=1.0 deg']
exclude_roi += ['Control region r=1.5 deg']
exclude_roi += ['Control region r=2.0 deg']
exclude_roi += ['Validation r=0.5 deg']
exclude_roi += ['Validation r=1.0 deg']
exclude_roi += ['Validation r=1.5 deg']
#exclude_roi += ['Validation r=2.0 deg']

exclude_roi += ['Unknown']

current_gamma_energy = 0.
next_gamma_energy = 0.
calibration = [1., 1., 1., 1., 1., 1.]
ref_rate = [1., 1., 1., 1., 1., 1.]
ref_rate_err = [1., 1., 1., 1., 1., 1.]

ref_rate = [45.22654248763829, 817.5534855109867, 807.4094245295635, 258.22163163239907, 54.18710411731964, 8.846775595651803]
ref_rate_err = [2.4996653879391553, 9.86773110448764, 9.161501194667895, 4.790992130002874, 2.1073161634064816, 0.8422544884101056]
calibration = [1.7171573302771546e-08, 3.815540141778941e-09, 5.205686613199803e-10, 6.111128565697704e-11, 6.175722680747704e-12, 5.443641190302156e-13]

SingleSourceAnalysis(sample_list,int(sys.argv[2]),int(sys.argv[3]))

Hist_EffArea_Sum.Print("All")

print ('exposure_hours = %0.1f hours'%(exposure_hours))
