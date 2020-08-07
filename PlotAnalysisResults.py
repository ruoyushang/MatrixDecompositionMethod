
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
#from scipy import special

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

#method_tag = '8bins_constrained'
#method_tag = '8bins_unconstrained'
method_tag = '16bins_constrained'
#method_tag = '16bins_unconstrained'

energy_bin_cut_low = 0
energy_bin_cut_up = 5

elev_bins = [45,85]
#theta2_bins = [0,1]
#theta2_bins = [0,1,4]
theta2_bins = [0,1,4,9]

ONOFF_tag = 'ON'
sample_list = []
sky_coord = []

##ONOFF_tag = 'OFF'
#sample_list += ['SgrAV6']
#sky_coord += ['17 45 39.6 -29 00 22']

#ONOFF_tag = 'OFF'
#sample_list += ['OJ287V6']
#sky_coord += ['08 54 49.1 +20 05 58.89']

#ONOFF_tag = 'OFF'
#sample_list += ['1ES0229V6']
#sky_coord += ['02 32 53.2 +20 16 21']
#sample_list += ['1ES0229V5']
#sky_coord += ['02 32 53.2 +20 16 21']

#ONOFF_tag = 'OFF'
#sample_list += ['H1426V6']
#sky_coord += ['14 28 32.609 +42 40 21.05']

#ONOFF_tag = 'OFF'
#sample_list += ['PKS1424V6']
#sky_coord += ['14 27 00 +23 47 00']

#ONOFF_tag = 'OFF'
#sample_list += ['3C264V6']
#sky_coord += ['11 45 5.009 +19 36 22.74']
# The VERITAS observations of 3C 264 were taken from February through May 2017, from February through April 2018, and from January through May 2019.

#ONOFF_tag = 'ON'
#sample_list += ['PG1553V5']
#sky_coord += ['15 55 44.7 +11 11 41']

#ONOFF_tag = 'OFF'
#sample_list += ['1ES1011V6']
#sky_coord += ['10 15 4.139 +49 26 0.71']

#ONOFF_tag = 'OFF'
#sample_list += ['RBS0413V6']
#sky_coord += ['03 19 47 +18 45 42']
#sample_list += ['RBS0413V5']
#sky_coord += ['03 19 47 +18 45 42']

#ONOFF_tag = 'OFF'
#sample_list += ['1ES0647V6']
#sky_coord += ['06 50 46.490 +25 02 59.62']

#ONOFF_tag = 'OFF'
#sample_list += ['RGBJ0710V5']
#sky_coord += ['07 10 26.4 +59 09 00']

#ONOFF_tag = 'OFF'
#sample_list += ['Segue1V6']
#sky_coord += ['10 07 04 +16 04 55']
#sample_list += ['Segue1V5']
#sky_coord += ['10 07 04 +16 04 55']
# only V5 data published

#ONOFF_tag = 'ON'
#sample_list += ['NGC1275V6']
#sky_coord += ['03 19 48.1 +41 30 42']

ONOFF_tag = 'OFF'
sample_list += ['M82V4']
sky_coord += ['09 55 52.7 +69 40 46']
sample_list += ['M82V5']
sky_coord += ['09 55 52.7 +69 40 46']

#ONOFF_tag = 'ON'
#sample_list += ['CrabV5']
#sky_coord += ['05 34 31.97 +22 00 52.1']
#sample_list += ['CrabV4']
#sky_coord += ['05 34 31.97 +22 00 52.1']

#ONOFF_tag = 'OFF'
#sample_list += ['Mrk421V5']
#sky_coord += ['11 04 19 +38 11 41']

#ONOFF_tag = 'ON'
#sample_list += ['2HWC_J1930V6']
#sky_coord += ['19 30 32 +18 52 12']

#ONOFF_tag = 'ON'
#sample_list += ['2HWC_J1953V6']
#sky_coord += ['19 53 02.4 +29 28 48']

#ONOFF_tag = 'ON'
#sample_list += ['WComaeV6']
#sky_coord += ['12 21 31.7 +28 13 59']
#sample_list += ['WComaeV5']
#sky_coord += ['12 21 31.7 +28 13 59']
#sample_list += ['WComaeV4']
#sky_coord += ['12 21 31.7 +28 13 59']
# https://arxiv.org/pdf/2002.04119.pdf VERITAS observations of 1ES 1215+303 from 2008 December to 2017 May.

#ONOFF_tag = 'ON'
#sample_list += ['IC443HotSpotV6']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#sample_list += ['IC443HotSpotV5']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#sample_list += ['IC443HotSpotV4']
#sky_coord += ['06 18 2.700 +22 39 36.00']

#ONOFF_tag = 'ON'
#sample_list += ['BoomerangV6']
#sky_coord += ['22 28 44 +61 10 00']
#sample_list += ['BoomerangV5']
#sky_coord += ['22 28 44 +61 10 00']
#sample_list += ['BoomerangV4']
#sky_coord += ['22 28 44 +61 10 00']

#ONOFF_tag = 'ON'
#sample_list += ['MGRO_J1908_V6']
#sky_coord += ['19 07 54 +06 16 07']
#sample_list += ['MGRO_J1908_V5']
#sky_coord += ['19 07 54 +06 16 07']
#sample_list += ['MGRO_J1908_V4']
#sky_coord += ['19 07 54 +06 16 07']

#ONOFF_tag = 'ON'
#sample_list += ['MAGIC_J1857_V6']
#sky_coord += ['18 57 35.6 +02 58 02']
#sample_list += ['MAGIC_J1857_V5']
#sky_coord += ['18 57 35.6 +02 58 02']
#sample_list += ['MAGIC_J1857_V4']
#sky_coord += ['18 57 35.6 +02 58 02']

#ONOFF_tag = 'ON'
#sample_list += ['MGRO_J2031_V6']
#sky_coord += ['20 28 43.2 +41 18 36']
#sample_list += ['MGRO_J2031_V5']
#sky_coord += ['20 28 43.2 +41 18 36']
#sample_list += ['MGRO_J2031_V4']
#sky_coord += ['20 28 43.2 +41 18 36']

#ONOFF_tag = 'ON'
#sample_list += ['CygnusV6']
#sky_coord += ['20 18 35.03 +36 50 00.0']
#sample_list += ['CygnusV5']
#sky_coord += ['20 18 35.03 +36 50 00.0']

#ONOFF_tag = 'ON'
#sample_list += ['GemingaV6']
#sky_coord += ['06 32 28 +17 22 00']
#sample_list += ['GemingaV5']
#sky_coord += ['06 32 28 +17 22 00']

#ONOFF_tag = 'ON'
#sample_list += ['ComaV4']
#sky_coord += ['12 59 48.7 +27 58 50']
#sample_list += ['ComaV6']
#sky_coord += ['12 59 48.7 +27 58 50']

root_file_tags = []
# all time
mjd_tag = []
mjd_tag += ['']

for elev in range(0,len(elev_bins)-1):
    elev_tag = '_TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

print 'Get %s'%(root_file_tags[0])

#for elev in range(0,len(elev_bins)-1):
#    # all time
#    root_file_tags += [method_tag+'_TelElev%sto%s_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    # 1ES 1215 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56372to56464_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    # WComae and 1ES 1218 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD54400to54700_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    # Geminga
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56000to56500_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58000to58700_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58700to59000_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    # MGRO J2031
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to57997_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57997to59000_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    # 3 slices of 1ES 0229
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to56500_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56501to57500_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57501to59500_%s'%(elev_bins[elev],elev_bins[elev+1],ONOFF_tag)]

selection_tag = root_file_tags[0]
selection_tag += '_E%s'%(energy_bin_cut_low)

folder_path = 'output_test'
#folder_path = 'output_root'
PercentCrab = ''

N_bins_for_deconv = 8
gamma_hadron_dim_ratio = 1.
MSCW_blind_cut = 0.5
MSCL_blind_cut = 0.5
MSCW_lower_cut = -1.
MSCL_lower_cut = -1.
MSCW_plot_lower = -1.
MSCL_plot_lower = -1.
MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
ErecS_lower_cut = 0
ErecS_upper_cut = 0

n_good_matches = 0
exposure_hours = 0.
exposure_hours_ref = 0.
Skymap_size = 3.
source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
n_control_samples = 3
MJD_Start = 2147483647
MJD_End = 0
roi_name = ROOT.std.vector("string")(10)
roi_ra = ROOT.std.vector("double")(10)
roi_dec = ROOT.std.vector("double")(10)
roi_radius = ROOT.std.vector("double")(10)

Syst_MDM = 0.02
Syst_Init = 0.02
Syst_Redu = 0.02
Syst_Corr = 0.02
Syst_Clos = 0.02

energy_bin = []
energy_bin += [int(pow(10,2.0))]
energy_bin += [int(pow(10,2.3))]
energy_bin += [int(pow(10,2.6))]
energy_bin += [int(pow(10,3.0))]
energy_bin += [int(pow(10,3.4))]
energy_bin += [int(pow(10,4.0))]

energy_fine_bin = []
energy_fine_bin += [pow(10,2.0)]
energy_fine_bin += [pow(10,2.1)]
energy_fine_bin += [pow(10,2.2)]
energy_fine_bin += [pow(10,2.3)]
energy_fine_bin += [pow(10,2.4)]
energy_fine_bin += [pow(10,2.5)]
energy_fine_bin += [pow(10,2.6)]
energy_fine_bin += [pow(10,2.7)]
energy_fine_bin += [pow(10,2.8)]
energy_fine_bin += [pow(10,2.9)]
energy_fine_bin += [pow(10,3.0)]
energy_fine_bin += [pow(10,3.1)]
energy_fine_bin += [pow(10,3.2)]
energy_fine_bin += [pow(10,3.3)]
energy_fine_bin += [pow(10,3.4)]
energy_fine_bin += [pow(10,3.5)]
energy_fine_bin += [pow(10,3.6)]
energy_fine_bin += [pow(10,3.7)]
energy_fine_bin += [pow(10,3.8)]
energy_fine_bin += [pow(10,3.9)]
energy_fine_bin += [pow(10,4.0)]

#sample_list += ['Crab']
#sky_coord += ['05 34 31.97 +22 00 52.1']
#sample_list += ['S3_1227_V6']
#sky_coord += ['12 30 14.1 +25 18 07']
#sample_list += ['MS1221V6']
#sky_coord += ['12 24 24.2 +24 36 24']
#sample_list += ['PKS1441V6']
#sky_coord += ['14 43 56.9 +25 01 44']
#sample_list += ['ComaV6']
#sky_coord += ['12 59 48.7 +27 58 50']
#sample_list += ['1ES1440V6']
#sky_coord += ['14 42 48.277 +12 00 40.37']
#sample_list += ['1ES1741V6']
#sky_coord += ['17 44 01.2 +19 32 47']
#sample_list += ['CasA']
#sky_coord += ['23 23 13.8 +58 48 26']
#sample_list += ['G079']
#sky_coord += ['20 32 28.56 +40 19 41.52']
#sample_list += ['1ES1218V6']
#sky_coord += ['12 21 26.3 +30 11 29']
#sample_list += ['TychoV6']
#sky_coord += ['00 25 21.6 +64 07 48']
#sample_list += ['CTA1V5']
#sky_coord += ['00 06 26 +72 59 01.0']

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

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

def GetBrightStarInfo(file_list):

    global bright_star_ra
    global bright_star_dec
    global bright_star_brightness
    global faint_star_ra
    global faint_star_dec
    global faint_star_brightness

    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        StarTree = InputFile.Get("StarTree")
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

    inputFile = open('TeVCat_RaDec_w_Names.txt')
    for line in inputFile:
        gamma_source_name = line.split(',')[0]
        gamma_source_ra = float(line.split(',')[1])
        gamma_source_dec = float(line.split(',')[2])
        other_stars += [gamma_source_name]
        other_star_coord += [[gamma_source_ra,gamma_source_dec]]

def ResetStackedShowerHistograms():

    Hist_EffArea_Sum.Reset()

    Hist2D_OnData_Sum.Reset()
    Hist2D_OnBkgd_Sum.Reset()
    Hist2D_OnDark_Sum.Reset()
    Hist_OnData_MSCL_Sum.Reset()
    Hist_OnBkgd_MSCL_Sum.Reset()
    Hist_OnDark_MSCL_Sum.Reset()
    Hist_OnData_MSCW_Sum.Reset()
    Hist_OnBkgd_MSCW_Sum.Reset()
    Hist_OnDark_MSCW_Sum.Reset()

    Hist2D_AllRanks_Sum.Reset()
    Hist2D_Rank0_Sum.Reset()
    Hist2D_Rank1_Sum.Reset()
    Hist2D_Rank2_Sum.Reset()
    Hist2D_Rank3_Sum.Reset()

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Reset()
        Hist2D_OffBkgd_Sum[nth_sample].Reset()
        Hist2D_OnSyst_Sum[nth_sample].Reset()
        Hist_OffData_MSCL_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCL_Sum[nth_sample].Reset()
        Hist_OffData_MSCW_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCW_Sum[nth_sample].Reset()

def GetSourceInfo(file_list):

    global N_bins_for_deconv
    global MSCW_blind_cut
    global MSCL_blind_cut
    global n_good_matches
    global exposure_hours
    global exposure_hours_ref
    global Skymap_size
    global source_ra
    global source_dec
    global source_l
    global source_b
    global n_control_samples
    global MJD_Start
    global MJD_End
    global roi_name
    global roi_ra
    global roi_dec
    global roi_radius

    n_good_matches = 0
    exposure_hours = 0.
    exposure_hours_ref = 0.
    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.SetBranchAddress('roi_name',ROOT.AddressOf(roi_name))
        InfoTree.SetBranchAddress('roi_ra',ROOT.AddressOf(roi_ra))
        InfoTree.SetBranchAddress('roi_dec',ROOT.AddressOf(roi_dec))
        InfoTree.SetBranchAddress('roi_radius',ROOT.AddressOf(roi_radius))
        InfoTree.GetEntry(0)
        N_bins_for_deconv = InfoTree.N_bins_for_deconv
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        Skymap_size = InfoTree.Skymap_size
        source_ra = InfoTree.mean_tele_point_ra
        source_dec = InfoTree.mean_tele_point_dec
        source_l = InfoTree.mean_tele_point_l
        source_b = InfoTree.mean_tele_point_b

        NewInfoTree = InputFile.Get("NewInfoTree")
        NewInfoTree.GetEntry(0)
        MJD_Start = min(NewInfoTree.MJD_Start,MJD_Start)
        MJD_End = max(NewInfoTree.MJD_End,MJD_End)
        if 'Theta20' in file_list[path]:
            exposure_hours += NewInfoTree.exposure_hours

        HistName = "Hist_EffArea"
        Hist_EffArea.Reset()
        Hist_EffArea.Add(InputFile.Get(HistName))
        Hist_EffArea.Scale(NewInfoTree.exposure_hours*3600.)
        Hist_EffArea_Sum.Add(Hist_EffArea)

        InputFile.Close()

def GetShowerHistogramsFromFile(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut

    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnData.Reset()
    Hist2D_OnData.Add(InputFile.Get(HistName))

    HistName = "Hist_OnDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnDark.Reset()
    Hist2D_OnDark.Add(InputFile.Get(HistName))

    HistName = "Hist_OnBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnBkgd.Reset()
    Hist2D_OnBkgd.Add(InputFile.Get(HistName))

    Hist2D_Rank0.Reset()
    Hist2D_Rank1.Reset()
    Hist2D_Rank2.Reset()
    Hist2D_Rank3.Reset()

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

    HistName = "Hist_Rank0_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank0.Add(InputFile.Get(HistName))
    #Hist2D_Rank1.Add(InputFile.Get(HistName))
    #Hist2D_Rank2.Add(InputFile.Get(HistName))
    #Hist2D_Rank3.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank1_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank1.Add(InputFile.Get(HistName))
    #Hist2D_Rank2.Add(InputFile.Get(HistName))
    #Hist2D_Rank3.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank2_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank2.Add(InputFile.Get(HistName))
    #Hist2D_Rank3.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank3_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank3.Add(InputFile.Get(HistName))

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

    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank0_LeftVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank0_LeftVector]
    legends += ['inital vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank0_LeftVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{1}}','Rank0_LeftVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank1_LeftVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank1_LeftVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank1_LeftVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{2}}','Rank1_LeftVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank2_LeftVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank2_LeftVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank2_LeftVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{3}}','Rank2_LeftVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank3_LeftVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank3_LeftVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank3_LeftVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{3}}','Rank3_LeftVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank0_RightVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank0_RightVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank0_RightVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{1}}','Rank0_RightVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank1_RightVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank1_RightVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank1_RightVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{2}}','Rank1_RightVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank2_RightVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank2_RightVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank2_RightVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{3}}','Rank2_RightVector_%s'%(ErecS_lower_cut),0,0,False,False)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist1D_Data_Rank3_RightVector]
    legends += ['data vector']
    colors += [1]
    Hists += [Hist1D_Dark_Rank3_RightVector]
    legends += ['initial vector']
    colors += [2]
    Hists += [Hist1D_Bkgd_Rank3_RightVector]
    legends += ['predict. vector']
    colors += [3]
    MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{3}}','Rank3_RightVector_%s'%(ErecS_lower_cut),0,0,False,False)

    if Hist2D_OnData.Integral()<1600. or Hist2D_OnDark.Integral()<1600.:
        Hist2D_OnData.Reset()
        Hist2D_OnDark.Reset()
        Hist2D_OnBkgd.Reset()
        Hist2D_Rank0.Reset()
        Hist2D_Rank1.Reset()
        Hist2D_Rank2.Reset()
        Hist2D_Rank3.Reset()

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

    for nth_sample in range(0,n_control_samples):

        HistName = "Hist_OffData2_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffData[nth_sample].Reset()
        Hist2D_OffData[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OffBkgd_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffBkgd[nth_sample].Reset()
        Hist2D_OffBkgd[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OnSyst_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OnSyst[nth_sample].Reset()
        Hist2D_OnSyst[nth_sample].Add(InputFile.Get(HistName))

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist2D_OffData[nth_sample].Reset()
            Hist2D_OffBkgd[nth_sample].Reset()
            Hist2D_OnSyst[nth_sample].Reset()
        if math.isnan(Hist2D_OffData[nth_sample].Integral()):
            Hist2D_OffData[nth_sample].Reset()
            Hist2D_OffBkgd[nth_sample].Reset()
            Hist2D_OnSyst[nth_sample].Reset()

        Hist_OffData_MSCL[nth_sample].Reset()
        Hist_OffData_MSCL[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionX("Hist1D_OffData_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffBkgd_MSCL[nth_sample].Reset()
        Hist_OffBkgd_MSCL[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionX("Hist1D_OffBkgd_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffData_MSCW[nth_sample].Reset()
        Hist_OffData_MSCW[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionY("Hist1D_OffData_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

        Hist_OffBkgd_MSCW[nth_sample].Reset()
        Hist_OffBkgd_MSCW[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionY("Hist1D_OffBkgd_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

def StackShowerHistograms():

    Hist2D_OnData_Sum.Add(Hist2D_OnData)
    Hist2D_OnDark_Sum.Add(Hist2D_OnDark)
    Hist2D_OnBkgd_Sum.Add(Hist2D_OnBkgd)

    Hist2D_Rank0_Sum.Add(Hist2D_Rank0)
    Hist2D_Rank1_Sum.Add(Hist2D_Rank1)
    Hist2D_Rank2_Sum.Add(Hist2D_Rank2)
    Hist2D_Rank3_Sum.Add(Hist2D_Rank3)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank0)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank1)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank2)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank3)

    Hist_OnData_MSCL_Sum.Add(Hist_OnData_MSCL)
    Hist_OnData_MSCW_Sum.Add(Hist_OnData_MSCW)
    Hist_OnDark_MSCL_Sum.Add(Hist_OnDark_MSCL)
    Hist_OnDark_MSCW_Sum.Add(Hist_OnDark_MSCW)
    Hist_OnBkgd_MSCL_Sum.Add(Hist_OnBkgd_MSCL)
    Hist_OnBkgd_MSCW_Sum.Add(Hist_OnBkgd_MSCW)

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Add(Hist2D_OffData[nth_sample])
        Hist2D_OffBkgd_Sum[nth_sample].Add(Hist2D_OffBkgd[nth_sample])
        Hist2D_OnSyst_Sum[nth_sample].Add(Hist2D_OnSyst[nth_sample])

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

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
    
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
        Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        if Hists[h]!=0:
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

    pad2.cd()
    Hist_Band = Hists[1].Clone()
    Hist_Band.Add(Hists[0],-1)
    Hist_Band.GetXaxis().SetTitle(title_x)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("error")
    Hist_Band.GetYaxis().SetTitleSize(0.10)
    Hist_Band.GetYaxis().SetNdivisions(505)
    plot_max = max(Hist_Band.GetMaximum(),-1.*Hist_Band.GetMinimum())
    Hist_Band.SetMaximum(plot_max*1.1)
    Hist_Band.SetMinimum(-1.*plot_max*1.1)
    Hist_Band.Draw("hist")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    for h in range(0,len(Hists)):
        Hist_Err = Hists[h].Clone()
        Hist_Err.Add(Hists[0],-1.)
        Hist_Err.SetLineColor(colors[h])
        Hist_Err.SetLineWidth(2)
        Hist_Err.Draw("hist same")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

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

def MakeSpectrumInNonCrabUnit(Hists,title,name,syst,nth_roi):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.1)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    Hist_EffArea_tmp = Hist_EffArea_Sum.Clone()

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 200, 4000)
    func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_MDM = Hists[0].Clone()
    Hist_MDM.Add(Hists[1],-1.)

    for binx in range(0,Hist_MDM.GetNbinsX()):
        if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
        Hist_MDM.SetBinContent(binx+1,Hist_MDM.GetBinContent(binx+1)*scale)
        Hist_MDM.SetBinError(binx+1,Hist_MDM.GetBinError(binx+1)*scale)

    Hist_MDM.GetYaxis().SetTitle("flux [counts/GeV/s/cm2]")
    Hist_MDM.GetXaxis().SetTitle(title)
    Hist_MDM.GetXaxis().SetTitleOffset(1.3)
    Hist_MDM.GetXaxis().SetTitleSize(0.05)
    Hist_MDM.Draw("E")

    if 'Crab' in name:
        func_crab.SetLineColor(4)
        func_crab.Draw("same")
    # MGRO J1908
    if 'MGRO_J1908' in name:
        func_1908 = ROOT.TF1("func_1908","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_1908.SetParameters(4.23,-2.2)
        func_1908.SetLineColor(4)
        func_1908.Draw("same")
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    if 'IC443' in name:
        print 'Compare to official IC 443 flux...'
        func_ic443 = ROOT.TF1("func_ic443","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_ic443.SetParameters(0.838,-2.99)
        func_ic443.SetLineColor(4)
        func_ic443.Draw("same")
    # 1ES 1218 https://arxiv.org/pdf/0810.0301.pdf
    if '1ES1218' in name:
        func_1218 = ROOT.TF1("func_1218","[0]*pow(10,-12)*pow(x/500.,[1])", 200, 4000)
        func_1218.SetParameters(7.5,-3.08)
        func_1218.SetLineColor(4)
        func_1218.Draw("same")


    pad3.cd()

    legend = ROOT.TLegend(0.15,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_MDM,"MDM flux","pl")
    if len(Hists)==3:
        legend.AddEntry(Hist_Ring,"Ring flux","pl")
    if 'IC443' in name:
        legend.AddEntry(func_ic443,"flux from arXiv:0905.3291","pl")
    if 'MGRO_J1908' in name:
        legend.AddEntry(func_1908,"flux from arXiv:1404.7185","pl")
    legend.Draw("SAME")

    #lumilab1 = ROOT.TLatex(0.15,0.80,'Exposure %.1f hrs'%(exposure_hours) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    pad1.SetLogy()
    pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))


def MakeLightCurvePlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    n_roi = len(Hist_data)
    pads = []

    #step_size = 1./float(n_roi)
    #for roi in range(0,n_roi):
    #    pads += [ROOT.TPad("pad%s"%(roi),"",0,1-(roi+1)*step_size,1,1-(roi)*step_size)]
    #    pads[roi].SetBorderMode(1)
    #    pads[roi].SetTopMargin(0.0)
    #    pads[roi].SetBottomMargin(0.0)
    #    if roi==0:
    #        pads[roi].SetTopMargin(0.1/step_size)
    #    if roi==n_roi-1:
    #        pads[roi].SetBottomMargin(0.1/step_size)
    #    pads[roi].Draw("same")

    step_size = 0.5
    for nth_pad in range(0,2):
        pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,1-(nth_pad+1)*step_size,1,1-(nth_pad)*step_size)]
        pads[nth_pad].SetBorderMode(1)
        pads[nth_pad].SetTopMargin(0.0)
        pads[nth_pad].SetBottomMargin(0.0)
        if nth_pad==0:
            pads[nth_pad].SetTopMargin(0.1/step_size)
        if nth_pad==1:
            pads[nth_pad].SetBottomMargin(0.1/step_size)
        pads[nth_pad].Draw("same")

    time = Time(MJD_Start, format='mjd')
    time.format = 'decimalyear'
    year_start = time.value
    time = Time(MJD_End, format='mjd')
    time.format = 'decimalyear'
    year_end = time.value

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    Func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    Func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    Hist_ratio_mjd = []
    IncValues = ROOT.TF1( "IncValues", "x", year_start , year_end )
    raLowerAxis = []
    days_per_bin = 90.
    for roi in range(0,n_roi):

        Hist_data_mjd += [ROOT.TH1D("Hist_data_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_data_mjd[roi].Fill(mjd,Hist_data[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
            Hist_data_mjd[roi].SetBinError(binx+1,pow(Hist_data_mjd[roi].GetBinContent(binx+1),0.5))

        Hist_bkgd_mjd += [ROOT.TH1D("Hist_bkgd_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_bkgd_mjd[roi].Fill(mjd,Hist_bkgd[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_mjd[roi].GetNbinsX()):
            Hist_bkgd_mjd[roi].SetBinError(binx+1,Hist_bkgd_mjd[roi].GetBinContent(binx+1)*Syst_MDM)


        Hist_ratio_mjd += [ROOT.TH1D("Hist_ratio_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        Hist_ratio_mjd[roi].Add(Hist_data_mjd[roi])
        Hist_ratio_mjd[roi].Add(Hist_bkgd_mjd[roi],-1.)
        Hist_ratio_mjd[roi].Divide(Hist_bkgd_mjd[roi])

        ratio_threshold = 1.0
        for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
            if Hist_ratio_mjd[roi].GetBinError(binx+1)==0.: continue
            significance = Hist_ratio_mjd[roi].GetBinContent(binx+1)/Hist_ratio_mjd[roi].GetBinError(binx+1)
            if significance<4.0: continue
            if Hist_ratio_mjd[roi].GetBinContent(binx+1)>ratio_threshold:
                print '>%s Crab, bin edge: %s, %s'%(ratio_threshold,Hist_ratio_mjd[roi].GetBinLowEdge(binx+1),Hist_ratio_mjd[roi].GetBinLowEdge(binx+2))

        pads[0].cd()
        Hist_ratio_mjd[roi].SetLineColor(1)
        Hist_ratio_mjd[roi].SetLineWidth(3)
        Hist_ratio_mjd[roi].Draw("E")
        c_both.Update()
        x1 = ROOT.gPad.GetUxmin()
        x2 = ROOT.gPad.GetUxmax()
        y1 = ROOT.gPad.GetUymin()
        y2 = ROOT.gPad.GetUymax()
        raLowerAxis += [ROOT.TGaxis( x1, y2, x2, y2,"IncValues", 510, "-")]
        raLowerAxis[roi].SetLabelSize(Hist_ratio_mjd[roi].GetXaxis().GetLabelSize())
        raLowerAxis[roi].SetTitle("Year")
        raLowerAxis[roi].Draw()

        pads[1].cd()
        Hist_data_mjd[roi].GetXaxis().SetTitle("MJD")
        Hist_data_mjd[roi].SetLineColor(1)
        Hist_data_mjd[roi].SetLineWidth(3)
        Hist_data_mjd[roi].Draw("E")
        fill_color = 38
        stack = ROOT.THStack("stack", "")
        set_histStyle( Hist_bkgd_mjd[roi] , fill_color)
        stack.Add( Hist_bkgd_mjd[roi] )
        stack.Draw("hist same")
        Hist_data_mjd[roi].Draw("E same")

        c_both.SaveAs('output_plots/%s_RoI%s_%s.png'%(name,roi,selection_tag))

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
            found_syst_hist = True
        if colors[h]!=0 and colors[h]!=1:
            Hist_Bkgd.Add(Hists[h].Clone())
    if not found_syst_hist:
        for binx in range(0,Hist_Syst.GetNbinsX()):
            syst_err = Syst_MDM*Hist_Bkgd.GetBinContent(binx+1)
            stat_err = Hist_Bkgd.GetBinError(binx+1)
            Hist_Syst.SetBinError(binx+1,pow(syst_err*syst_err+stat_err*stat_err,0.5))

    fill_color = [0,0,46,0,38,30]
    if doSum:
        stack = ROOT.THStack("stack", "")
        for h in range(1,len(Hists)):
            if stack_it[h]:
                set_histStyle( Hists[h] , fill_color[colors[h]])
                stack.Add( Hists[h] )
        stack.Draw("hist same")
        Hist_Syst.SetFillColor(1)
        Hist_Syst.SetFillStyle(3004)
        Hist_Syst.SetMarkerSize(0)
        Hist_Syst.Draw("e2 same")

    for h in range(0,len(Hists)):
        if colors[h]==0: continue
        Hists[h].SetLineWidth(3)
        Hists[h].Draw("E same")
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
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
        if doSum:
            legend.AddEntry(Hists[h],legends[h],"f")
        else:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()

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
    lumilab2 = ROOT.TLatex(0.15,0.60,'Excess = %0.1f#pm%0.1f (%0.1f#sigma)'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5),Sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.40,'Bkg = %0.1f#pm%0.1f'%(predict_bkg,err_bkg) )
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
    lumilab4 = ROOT.TLatex(0.15,0.20,'S/B = %0.3f#pm%0.3f'%(sbratio,sbratio_err) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()

    pad2.cd()
    Hist_Band = Hist_Syst.Clone()
    Hist_Band.Divide(Hist_Bkgd)
    Hist_Band.SetFillColor(1)
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
    Hist_Ratio = []
    for h in range(1,len(Hists)):
        Hist_Bkg = Hists[h].Clone()
        Hist_Ratio += [Hists[0].Clone()]
        Hist_Ratio[h-1].Divide(Hist_Bkg)
        Hist_Ratio[h-1].SetLineWidth(2)
        Hist_Ratio[h-1].SetLineColor(colors[h])
        Hist_Ratio[h-1].Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def PlotsStackedHistograms(tag):

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
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_lower_cut,MSCW_blind_cut,-1)

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
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_lower_cut,MSCW_blind_cut,-1)

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
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_lower_cut,MSCL_blind_cut,-1)

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
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_lower_cut,MSCL_blind_cut,-1)

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
    plotname = 'Stack_Theta2_MDM_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    for nth_roi in range(0,len(roi_ra)):
        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Theta2_Sum[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
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

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        legends += ['obs. data (%s)'%(roi_name[nth_roi])]
        colors += [1]
        Hists += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        legends += ['MDM bkg.']
        colors += [4]
        title = 'energy [GeV]'
        plotname = 'CrabUnit_RoI%s_Energy_MDM_%s'%(nth_roi,tag)
        MakeSpectrumInNonCrabUnit(Hists,title,plotname,-1,nth_roi)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    legends = []
    colors = []
    for nth_roi in range(0,len(roi_ra)):
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [nth_roi+1]
        Hist_data_mjd += [Hist_OnData_RoI_MJD_Sum[nth_roi]]
        Hist_bkgd_mjd += [Hist_OnBkgd_RoI_MJD_Sum[nth_roi]]
    plotname = 'LightCurve_RoI_Year_MDM_%s'%(tag)
    title = 'Year'
    MakeLightCurvePlot(Hist_data_mjd,Hist_bkgd_mjd,legends,colors,title,plotname)

    Hist_data_energy = []
    Hist_bkgd_energy = []
    legends = []
    colors = []
    for nth_roi in range(0,len(roi_ra)):
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [nth_roi+1]
        Hist_data_energy += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        Hist_bkgd_energy += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
    plotname = 'Spectrum_RoI%s_CrabUnit_MDM_%s'%(nth_roi,tag)
    title = 'GeV'
    MakeCrabUnitSpectrumPlot(Hist_data_energy,Hist_bkgd_energy,legends,colors,title,plotname)

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
    plotname = 'Stack_Energy_MDM_%s'%(tag)
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


    for nth_sample in range(0,n_control_samples):

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OffData_MSCW_Sum[nth_sample]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OffBkgd_MSCW_Sum[nth_sample]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [False]
        plotname = 'Stack_MSCW_MDM_%s_Control%s'%(tag,nth_sample)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCW_lower_cut,MSCW_blind_cut,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OffData_MSCL_Sum[nth_sample]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OffBkgd_MSCL_Sum[nth_sample]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [False]
        plotname = 'Stack_MSCL_MDM_%s_Control%s'%(tag,nth_sample)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCL_lower_cut,MSCL_blind_cut,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OffData_Energy_Sum[nth_sample]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OffBkgd_Energy_Sum[nth_sample]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [False]
        plotname = 'Stack_Energy_MDM_%s_Control%s'%(tag,nth_sample)
        title = 'energy [GeV]'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,pow(10,4.0),-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        #Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
        #Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
        Hists += [Hist_OffData_CameraFoV_Theta2_Sum[nth_sample]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [False]
        plotname = 'Stack_CameraFoV_Theta2_MDM_%s_Control%s'%(tag,nth_sample)
        title = 'squared angle from camera center #theta^{2}'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,1.,-1)

def CalculateSystError():

    global Syst_MDM
    global Syst_Init
    global Syst_Redu
    global Syst_Corr
    global Syst_Clos
    Syst_MDM = 0.
    Syst_Init = 0.
    Syst_Redu = 0.
    Syst_Corr = 0.
    Syst_Clos = 0.
    Hist_SystErr_MSCL.Reset()
    Hist_SystErr_MSCW.Reset()

    binx_low_target = Hist_OnData_MSCL_Sum.FindBin(MSCL_lower_cut)
    binx_up_target = Hist_OnData_MSCL_Sum.FindBin(MSCL_plot_upper)-1
    binx_blind_target = Hist_OnData_MSCL_Sum.FindBin(MSCL_blind_cut)-1
    biny_low_target = Hist_OnData_MSCW_Sum.FindBin(MSCW_lower_cut)
    biny_up_target = Hist_OnData_MSCW_Sum.FindBin(MSCW_plot_upper)-1
    biny_blind_target = Hist_OnData_MSCW_Sum.FindBin(MSCW_blind_cut)-1

    Hist2D_Diff = Hist2D_OnBkgd_Sum.Clone()
    Hist2D_Diff.Add(Hist2D_OnData_Sum,-1.)
    for binx in range(0,Hist2D_Diff.GetNbinsX()):
        for biny in range(0,Hist2D_Diff.GetNbinsY()):
            if binx<=binx_blind_target and biny<=biny_blind_target: continue
            Syst_Clos += pow(Hist2D_Diff.GetBinContent(binx+1,biny+1),2)
    Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
    Syst_Clos = pow(Syst_Clos/3.,0.5)/Total_OnBkgd

    n_samples_used = 0.
    for nth_sample in range(0,n_control_samples):

        Total_OffData = Hist2D_OffBkgd_Sum[nth_sample].Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
        if Total_OffData == 0.: continue
        Hist2D_Diff = Hist2D_OffData_Sum[nth_sample].Clone()
        Hist2D_Diff.Add(Hist2D_OffBkgd_Sum[nth_sample],-1.)
        syst_this = Hist2D_Diff.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)/Total_OffData
        Syst_Corr += syst_this
        Syst_Redu += pow(syst_this,2)

        Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
        if Total_OffData == 0.: continue
        Hist2D_Diff = Hist2D_OnBkgd_Sum.Clone()
        Hist2D_Diff.Add(Hist2D_OnSyst_Sum[nth_sample],-1.)
        syst_this = Hist2D_Diff.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)/Total_OnBkgd
        Syst_Init += pow(syst_this,2)

        print '%s th sample syst = %s'%(nth_sample,syst_this)
        n_samples_used += 1.

        Hist_Diff_MSCL = Hist_OffData_MSCL_Sum[nth_sample].Clone()
        Hist_Diff_MSCL.Add(Hist_OffBkgd_MSCL_Sum[nth_sample],-1.)
        for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
            syst_err = max(0.,pow(Hist_Diff_MSCL.GetBinContent(binx+1),2)-pow(Hist_Diff_MSCL.GetBinError(binx+1),2))
            old_syst = Hist_SystErr_MSCL.GetBinError(binx+1)
            new_syst = old_syst+syst_err
            Hist_SystErr_MSCL.SetBinError(binx+1,new_syst)

        Hist_Diff_MSCW = Hist_OffData_MSCW_Sum[nth_sample].Clone()
        Hist_Diff_MSCW.Add(Hist_OffBkgd_MSCW_Sum[nth_sample],-1.)
        for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
            syst_err = max(0.,pow(Hist_Diff_MSCW.GetBinContent(binx+1),2)-pow(Hist_Diff_MSCW.GetBinError(binx+1),2))
            old_syst = Hist_SystErr_MSCW.GetBinError(binx+1)
            new_syst = old_syst+syst_err
            Hist_SystErr_MSCW.SetBinError(binx+1,new_syst)

    Syst_MDM = Syst_Init + Syst_Redu
    Syst_MDM = pow(Syst_MDM/n_samples_used,0.5)
    Syst_MDM = pow(Syst_MDM*Syst_MDM + Syst_Clos*Syst_Clos,0.5)
    Syst_Init = pow(Syst_Init/n_samples_used,0.5)
    Syst_Redu = pow(Syst_Redu/n_samples_used,0.5)
    Syst_Corr = Syst_Corr/n_samples_used

    for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
            old_syst = Hist_SystErr_MSCL.GetBinError(binx+1)
            new_syst = pow(old_syst/(n_control_samples),0.5)
            Hist_SystErr_MSCL.SetBinError(binx+1,new_syst)

    for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
            old_syst = Hist_SystErr_MSCW.GetBinError(binx+1)
            new_syst = pow(old_syst/(n_control_samples),0.5)
            Hist_SystErr_MSCW.SetBinError(binx+1,new_syst)

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
    bin_lower = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Energy.Reset()
    Hist_OnData_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Energy.Reset()
    Hist_OnBkgd_Energy.Add(InputFile.Get(HistName))

    HistName = "Hist_OnData_SR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Zenith.Reset()
    Hist_OnData_Zenith.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Zenith.Reset()
    Hist_OnBkgd_Zenith.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Energy[nth_roi].Reset()
        Hist_OnData_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
        Hist_OnBkgd_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_SR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_MJD[nth_roi].Reset()
        Hist_OnData_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_MJD[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Energy.Reset()
        Hist_OnBkgd_Energy.Reset()
        Hist_OnData_Zenith.Reset()
        Hist_OnBkgd_Zenith.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Energy[nth_roi].Reset()
            Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
            Hist_OnData_RoI_MJD[nth_roi].Reset()
            Hist_OnBkgd_RoI_MJD[nth_roi].Reset()

def NormalizeTheta2Histograms(FilePath):

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

    HistName = "Hist_OnData_SR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Theta2.Reset()
    Hist_OnData_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Theta2.Reset()
    Hist_OnBkgd_Theta2.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Theta2[nth_roi].Reset()
        Hist_OnData_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()
        Hist_OnBkgd_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))


    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Theta2.Reset()
        Hist_OnBkgd_Theta2.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Theta2[nth_roi].Reset()
            Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()

def CorrectEnergyHistograms():

    Hist_OnBkgd_Energy_Sum.Scale(1.+Syst_Corr)
    Hist_OnBkgd_Zenith_Sum.Scale(1.+Syst_Corr)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Scale(1.+Syst_Corr)
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Scale(1.+Syst_Corr)

def StackEnergyHistograms():

    Hist_OnData_Energy_Sum.Add(Hist_OnData_Energy)
    Hist_OnBkgd_Energy_Sum.Add(Hist_OnBkgd_Energy)
    Hist_OnData_Zenith_Sum.Add(Hist_OnData_Zenith)
    Hist_OnBkgd_Zenith_Sum.Add(Hist_OnBkgd_Zenith)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Add(Hist_OnData_RoI_Energy[nth_roi])
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Energy[nth_roi])
        Hist_OnData_RoI_MJD_Sum[nth_roi].Add(Hist_OnData_RoI_MJD[nth_roi])
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Add(Hist_OnBkgd_RoI_MJD[nth_roi])

def CorrectTheta2Histograms():

    Hist_OnBkgd_Theta2_Sum.Scale(1.+Syst_Corr)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Scale(1.+Syst_Corr)
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Scale(1.+Syst_Corr)

def StackTheta2Histograms():

    Hist_OnData_Theta2_Sum.Add(Hist_OnData_Theta2)
    Hist_OnBkgd_Theta2_Sum.Add(Hist_OnBkgd_Theta2)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Add(Hist_OnData_RoI_Theta2[nth_roi])
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Theta2[nth_roi])
        Hist_OnData_RoI_MJD_Sum[nth_roi].Add(Hist_OnData_RoI_MJD[nth_roi])
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Add(Hist_OnBkgd_RoI_MJD[nth_roi])

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

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()

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

def NormalizeSkyMapHistograms(FilePath):

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

    HistName = "Hist_OnData_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap.Reset()
    Hist_OnData_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap_Galactic.Reset()
    Hist_OnData_Skymap_Galactic.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_VR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_VR_Skymap.Reset()
    Hist_OnData_VR_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_VR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_VR_Skymap_Galactic.Reset()
    Hist_OnData_VR_Skymap_Galactic.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap.Reset()
    Hist_OnBkgd_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Galactic.Reset()
    Hist_OnBkgd_Skymap_Galactic.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_Skymap_Syst_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Syst.Reset()
    Hist_OnBkgd_Skymap_Syst.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_Skymap_Galactic_Syst_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Galactic_Syst.Reset()
    Hist_OnBkgd_Skymap_Galactic_Syst.Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Skymap.Reset()
        Hist_OnData_Skymap_Galactic.Reset()
        Hist_OnData_VR_Skymap.Reset()
        Hist_OnData_VR_Skymap_Galactic.Reset()
        Hist_OnBkgd_Skymap.Reset()
        Hist_OnBkgd_Skymap_Galactic.Reset()
        Hist_OnBkgd_Skymap_Syst.Reset()
        Hist_OnBkgd_Skymap_Galactic_Syst.Reset()

def CorrectSkymapHistograms():

    Hist_OnBkgd_Skymap_Sum.Scale(1.+Syst_Corr)
    Hist_OnBkgd_Skymap_Galactic_Sum.Scale(1.+Syst_Corr)

def StackSkymapHistograms():

    Hist_OnData_Skymap_Sum.Add(Hist_OnData_Skymap)
    Hist_OnData_Skymap_Galactic_Sum.Add(Hist_OnData_Skymap_Galactic)
    Hist_OnData_VR_Skymap_Sum.Add(Hist_OnData_VR_Skymap)
    Hist_OnData_VR_Skymap_Galactic_Sum.Add(Hist_OnData_VR_Skymap_Galactic)
    Hist_OnBkgd_Skymap_Sum.Add(Hist_OnBkgd_Skymap)
    Hist_OnBkgd_Skymap_Galactic_Sum.Add(Hist_OnBkgd_Skymap_Galactic)
    Hist_OnBkgd_Skymap_Syst_Sum.Add(Hist_OnBkgd_Skymap_Syst)
    Hist_OnBkgd_Skymap_Galactic_Syst_Sum.Add(Hist_OnBkgd_Skymap_Galactic_Syst)

def Smooth2DMap(Hist_Old,smooth_size,addLinearly):

    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(4*smooth_size/bin_size) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            bin_content = 0
            bin_error = 0
            bin_norm = 0
            locationx1 = Hist_Old.GetXaxis().GetBinCenter(bx1)
            locationy1 = Hist_Old.GetYaxis().GetBinCenter(by1)
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2>=1 and bx2<=Hist_Old.GetNbinsX():
                        if by2>=1 and by2<=Hist_Old.GetNbinsY():
                            locationx2 = Hist_Old.GetXaxis().GetBinCenter(bx2)
                            locationy2 = Hist_Old.GetYaxis().GetBinCenter(by2)
                            distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5)
                            bin_content += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinContent(bx2,by2)
                            bin_norm += ROOT.TMath.Gaus(distance,0,smooth_size)
                            if not addLinearly:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*pow(Hist_Old.GetBinError(bx2,by2),2)
                            else:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    return Hist_Smooth

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_VR,Hist_Syst,syst_method):

    Hist_Skymap_Syst = Hist_VR.Clone()
    Hist_Skymap_Syst.Add(Hist_Syst,-1.)
    Hist_Skymap_Syst.Divide(Hist_Syst)

    Hist_Skymap = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            sky_x = Hist_SR.GetXaxis().GetBinCenter(bx+1)
            sky_y = Hist_SR.GetYaxis().GetBinCenter(by+1)
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            bx2 = Hist_Skymap_Syst.GetXaxis().FindBin(sky_x)
            by2 = Hist_Skymap_Syst.GetYaxis().FindBin(sky_y)
            Shape_Err = Hist_Skymap_Syst.GetBinContent(bx2+1,by2+1)*NBkg
            #Shape_Err = 0.
            Norm_Err = syst_method*NBkg
            Stat_Err = Hist_Bkg.GetBinError(bx+1,by+1)
            NBkg_Err = pow(pow(Stat_Err,2)+pow(Norm_Err,2)+pow(Shape_Err,2),0.5)
            if syst_method==0.: NBkg_Err = Stat_Err
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
    return Hist_Skymap

def reflectXaxis(hist):

    # taken from VPlotAnasumHistograms.cpp
	
    # temporary histogram
    hT = ROOT.TH2D( "%s_REFLECTED"%(hist.GetName()), "", hist.GetNbinsX(), -1.*hist.GetXaxis().GetXmax(), -1.*hist.GetXaxis().GetXmin(), hist.GetNbinsY(), hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax() )
    hT.SetStats( 0 )
    hT.SetXTitle( hist.GetXaxis().GetTitle() )
    hT.SetYTitle( hist.GetYaxis().GetTitle() )
	
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsX()+1):
            hT.SetBinContent( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinContent( binx, biny ) )
    return hT

def Make2DSignificancePlot(syst_method,Hist_SR,Hist_Bkg,Hist_VR,Hist_Syst,xtitle,ytitle,name):

    Hist_Skymap = GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_VR,Hist_Syst,syst_method)
    Hist_Skymap = reflectXaxis(Hist_Skymap)
    max_sig = Hist_Skymap.GetMaximum()

    Hist_Skymap_Syst = Hist_VR.Clone()
    Hist_Skymap_Syst.Add(Hist_Syst,-1.)
    Hist_Skymap_Syst.Divide(Hist_Syst)
    Hist_Skymap_Syst = reflectXaxis(Hist_Skymap_Syst)

    other_star_labels = []
    other_star_markers = []
    other_star_names = []
    other_star_significance = []
    bright_star_labels = []
    bright_star_markers = []
    faint_star_labels = []
    faint_star_markers = []
    star_range = 3.0
    if xtitle=="RA":
        for star in range(0,len(other_stars)):
            if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
            other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
            other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.02)
            other_star_names += [other_stars[star]]
            #binx = Hist_Skymap.GetXaxis().FindBin(-other_star_coord[star][0])
            #biny = Hist_Skymap.GetYaxis().FindBin(other_star_coord[star][1])
            #other_star_significance += [Hist_Skymap.GetBinContent(binx,biny)]
            other_star_significance += [FindLocalMaximum(Hist_Skymap, -other_star_coord[star][0], other_star_coord[star][1])]
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
    else:
        for star in range(0,len(other_stars)):
            gal_l, gal_b = ConvertRaDecToGalactic(other_star_coord[star][0],other_star_coord[star][1])
            if pow(source_l-gal_l,2)+pow(source_b-gal_b,2)>star_range*star_range: continue
            other_star_markers += [ROOT.TMarker(-gal_l,gal_b,2)]
            other_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.02)
            other_star_names += [other_stars[star]]
            #binx = Hist_Skymap.GetXaxis().FindBin(-gal_l)
            #biny = Hist_Skymap.GetYaxis().FindBin(gal_b)
            #other_star_significance += [Hist_Skymap.GetBinContent(binx,biny)]
            other_star_significance += [FindLocalMaximum(Hist_Skymap, -gal_l, gal_b)]
        for star in range(0,len(bright_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(bright_star_ra[star],bright_star_dec[star])
            bright_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
        for star in range(0,len(faint_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(faint_star_ra[star],faint_star_dec[star])
            faint_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
            faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
            faint_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(faint_star_brightness[star]))]
            faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
            faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)

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

    pad1.cd()

    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for star in range(0,len(other_star_names)):
        legend.AddEntry(other_stars[star],'%s (%0.1f#sigma)'%(other_star_names[star],other_star_significance[star]))

    Hist_Contour = Hist_Skymap.Clone()
    Hist_Contour.Reset()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            Hist_Contour.SetBinContent(bx+1,by+1,(Hist_Skymap.GetBinContent(bx+1,by+1)))
    Hist_Contour.SetContour(3)
    Hist_Contour.SetContourLevel(0,3)
    Hist_Contour.SetContourLevel(1,4)
    Hist_Contour.SetContourLevel(2,5)

    Hist_Skymap.GetYaxis().SetTitle(ytitle)
    Hist_Skymap.GetXaxis().SetTitle(xtitle)
    Hist_Skymap.SetMaximum(5)
    Hist_Skymap.SetMinimum(-5)

    pad1.cd()
    Hist_Skymap.Draw("COL4Z")
    Hist_Contour.Draw("CONT3 same")
    Hist_Skymap.GetXaxis().SetLabelOffset(999)
    Hist_Skymap.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap.GetXaxis().GetXmin()
    x2 = Hist_Skymap.GetXaxis().GetXmax()
    y1 = Hist_Skymap.GetYaxis().GetXmin()
    y2 = Hist_Skymap.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'max. %0.1f#sigma (syst = %0.1f%%)'%(max_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.50,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.30,'MJD %s-%s'%(MJD_Start,MJD_End) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/SkymapSig_%s_%s.png'%(name,selection_tag))

    pad1.cd()
    #Hist_Skymap_Syst.SetMaximum(0.1)
    #Hist_Skymap_Syst.SetMinimum(-0.1)
    pad1.SetLogz(1)
    Hist_Skymap_Syst.Draw("COL4Z")
    Hist_Contour.Draw("CONT3 same")
    Hist_Skymap_Syst.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Syst.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_Syst.GetXaxis().GetXmin()
    x2 = Hist_Skymap_Syst.GetXaxis().GetXmax()
    y1 = Hist_Skymap_Syst.GetYaxis().GetXmin()
    y2 = Hist_Skymap_Syst.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_Syst.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'max. %0.1f#sigma (syst = %0.1f%%)'%(max_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.50,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.30,'MJD %s-%s'%(MJD_Start,MJD_End) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/SkymapSyst_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)


    Hist_Skymap_Excess = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            Hist_Skymap_Excess.SetBinContent(bx+1,by+1,NSR-NBkg)
    Hist_Skymap_Excess = reflectXaxis(Hist_Skymap_Excess)

    pad1.cd()
    Hist_Skymap_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_Excess.Draw("COL4Z")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    Hist_Skymap_Excess.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Excess.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_Excess.GetXaxis().GetXmin()
    x2 = Hist_Skymap_Excess.GetXaxis().GetXmax()
    y1 = Hist_Skymap_Excess.GetYaxis().GetXmin()
    y2 = Hist_Skymap_Excess.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_Excess.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapExcess_%s_%s.png'%(name,selection_tag))

    Hist_Skymap_Ratio = Hist_SR.Clone()
    Hist_Skymap_Ratio.Add(Hist_Bkg,-1.)
    Hist_Skymap_Ratio.Divide(Hist_Bkg)
    Hist_Skymap_Ratio = reflectXaxis(Hist_Skymap_Ratio)
    for bx in range(0,Hist_Skymap_Ratio.GetNbinsX()):
        for by in range(0,Hist_Skymap_Ratio.GetNbinsY()):
            if Hist_SR.GetBinContent(bx+1,by+1)<100.:
                Hist_Skymap_Ratio.SetBinContent(bx+1,by+1,0.)

    pad1.cd()
    Hist_Skymap_Ratio.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_Ratio.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_Ratio.Draw("COL4Z")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    Hist_Skymap_Ratio.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Ratio.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_Ratio.GetXaxis().GetXmin()
    x2 = Hist_Skymap_Ratio.GetXaxis().GetXmax()
    y1 = Hist_Skymap_Ratio.GetYaxis().GetXmin()
    y2 = Hist_Skymap_Ratio.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_Ratio.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapRatio_%s_%s.png'%(name,selection_tag))

    MapEdge_left = Hist_Skymap_Excess.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Skymap_Excess.GetXaxis().GetBinLowEdge(Hist_Skymap_Excess.GetNbinsX()+1)
    MapEdge_lower = Hist_Skymap_Excess.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Skymap_Excess.GetYaxis().GetBinLowEdge(Hist_Skymap_Excess.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.
    Hist_Skymap_zoomin = ROOT.TH2D("Hist_Skymap_zoomin","",25,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,25,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_Contour_zoomin = ROOT.TH2D("Hist_Contour_zoomin","",25,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,25,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_Skymap_zoomin.Rebin2D(n_rebin,n_rebin)
    Hist_Contour_zoomin.Rebin2D(n_rebin,n_rebin)
    for bx in range(0,Hist_Skymap_Excess.GetNbinsX()):
        for by in range(0,Hist_Skymap_Excess.GetNbinsY()):
            bx_center = Hist_Skymap_Excess.GetXaxis().GetBinCenter(bx+1)
            by_center = Hist_Skymap_Excess.GetYaxis().GetBinCenter(by+1)
            bx2 = Hist_Skymap_zoomin.GetXaxis().FindBin(bx_center)
            by2 = Hist_Skymap_zoomin.GetYaxis().FindBin(by_center)
            #Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Skymap_Excess.GetBinContent(bx+1,by+1))
            Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Skymap_Ratio.GetBinContent(bx+1,by+1))
            Hist_Contour_zoomin.SetBinContent(bx2,by2,Hist_Contour.GetBinContent(bx+1,by+1))
    Hist_Contour_zoomin.SetContour(3)
    Hist_Contour_zoomin.SetContourLevel(0,3)
    Hist_Contour_zoomin.SetContourLevel(1,4)
    Hist_Contour_zoomin.SetContourLevel(2,5)

    pad1.cd()
    Hist_Skymap_zoomin.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_zoomin.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_zoomin.Draw("COL4Z")
    Hist_Contour_zoomin.Draw("CONT3 same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    Hist_Skymap_zoomin.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_zoomin.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_zoomin.GetXaxis().GetXmin()
    x2 = Hist_Skymap_zoomin.GetXaxis().GetXmax()
    y1 = Hist_Skymap_zoomin.GetYaxis().GetXmin()
    y2 = Hist_Skymap_zoomin.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_zoomin.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapZoomin_%s_%s.png'%(name,selection_tag))

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

def GetExtention(Hist_data, Hist_bkgd, Hist_sig, highlight_threshold, init_x, init_y):

    Hist_Excess = Hist_data.Clone()
    Hist_Excess.Add(Hist_bkgd,-1.)

    xx, yy, zz = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
    maxbin = Hist_sig.GetMaximumBin()
    Hist_sig.GetBinXYZ(maxbin, xx, yy, zz)
    print 'max excess at %s, %s'%(Hist_sig.GetXaxis().GetBinCenter(xx),Hist_sig.GetYaxis().GetBinCenter(yy))
    minbin = Hist_sig.GetMinimumBin()
    Hist_sig.GetBinXYZ(minbin, xx, yy, zz)
    print 'min excess at %s, %s'%(Hist_sig.GetXaxis().GetBinCenter(xx),Hist_sig.GetYaxis().GetBinCenter(yy))

    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_x = Hist_Excess.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_Excess.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_x-init_x,2)+pow(bin_y-init_y,2),0.5)
            if distance > 0.5:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
            if not Hist_sig.GetBinContent(bx+1,by+1)>=highlight_threshold: 
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
    excess_center_x_init = Hist_Excess.GetMean(1)
    excess_center_y_init = Hist_Excess.GetMean(2)
    excess_radius_init = pow(pow(Hist_Excess.GetRMS(1),2)+pow(Hist_Excess.GetRMS(2),2),0.5)

    Hist_Excess.Reset()
    Hist_Excess.Add(Hist_data)
    Hist_Excess.Add(Hist_bkgd,-1.)
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_delta_ra = Hist_Excess.GetXaxis().GetBinCenter(bx+1)-excess_center_x_init
            bin_delta_dec = Hist_Excess.GetYaxis().GetBinCenter(by+1)-excess_center_y_init
            if pow(bin_delta_ra*bin_delta_ra+bin_delta_dec*bin_delta_dec,0.5)>(2.*excess_radius_init):
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
    excess_center_x = Hist_Excess.GetMean(1)
    excess_center_y = Hist_Excess.GetMean(2)
    excess_radius = pow(pow(Hist_Excess.GetRMS(1),2)+pow(Hist_Excess.GetRMS(2),2),0.5)

    #return excess_center_x_init, excess_center_y_init, excess_radius_init
    return excess_center_x, excess_center_y, excess_radius

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

    line1 = ROOT.TLine(MSCL_plot_lower,MSCW_blind_cut,MSCL_blind_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line2 = ROOT.TLine(MSCL_blind_cut,MSCW_plot_lower,MSCL_blind_cut,MSCW_blind_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M_{Data}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_OnData_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_OnData_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/OnData_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M_{Initial}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnDark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_OnDark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_OnDark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/OnDark_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#lambda_{1} q_{1} p_{1}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank0_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank0_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank0_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank0_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#lambda_{2} q_{2} p_{2}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank1_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank1_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank1_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank1_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#lambda_{3} q_{3} p_{3}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank2_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank2_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank2_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank2_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#lambda_{4} q_{4} p_{4}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank3_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank3_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank3_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank3_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrDark = ROOT.TH2D("Hist2D_ErrDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrDark.Add(Hist2D_OnData_Sum)
    Hist2D_ErrDark.Add(Hist2D_OnDark_Sum,-1)
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'Error map of initial matrix' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrDark.GetYaxis().SetTitle('MSCW')
    Hist2D_ErrDark.GetXaxis().SetTitle('MSCL')
    Hist2D_ErrDark.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/ErrorDark_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrBkgd = ROOT.TH2D("Hist2D_ErrBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrBkgd.Add(Hist2D_OnData_Sum)
    Hist2D_ErrBkgd.Add(Hist2D_OnBkgd_Sum,-1)
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'Error map of MDM prediction' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrBkgd.SetMaximum(Hist2D_ErrDark.GetMaximum());
    Hist2D_ErrBkgd.SetMinimum(Hist2D_ErrDark.GetMinimum());
    Hist2D_ErrBkgd.GetYaxis().SetTitle('MSCW')
    Hist2D_ErrBkgd.GetXaxis().SetTitle('MSCL')
    Hist2D_ErrBkgd.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/ErrorBkgd_%s_%s.png'%(name,selection_tag))

def MakeOneHistPlot(Hist,title_x,title_y,name,logy):
    
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
    #if logy: pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad3.cd()
    lumilab2 = ROOT.TLatex(0.15,0.50,'#epsilon = (#sum_{#gamma region} M_{data} - #sum_{i=1}^{i=n} #lambda_{i} q_{i} p_{i}^{T})  / #sum_{#gamma region} M_{data}' )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.2)
    lumilab2.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

    Hist.GetXaxis().SetTitleOffset(0.8)
    Hist.GetXaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetTitleOffset(1.2)
    Hist.GetYaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetTitle(title_x)
    Hist.GetYaxis().SetTitle(title_y)
    Hist.Draw("E")

    pad3.cd()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeRankResidualPlots(name):

    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    data_integral = Hist2D_OnData_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)
    bkgd_rank_integral = []
    bkgd_rank_integral += [Hist2D_Rank0_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    bkgd_rank_integral += [Hist2D_Rank1_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    bkgd_rank_integral += [Hist2D_Rank2_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    bkgd_rank_integral += [Hist2D_Rank3_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]

    Hist_Residual_Rank = ROOT.TH1D("Hist_Residual_Rank","",4,0,4)
    for binx in range(0,len(bkgd_rank_integral)):
        bkgd_integral = 0.
        for rank in range (0,binx+1):
            bkgd_integral += bkgd_rank_integral[rank]
        Hist_Residual_Rank.SetBinContent(binx+1,abs(1.-bkgd_integral/data_integral))
        ratio_err = bkgd_integral/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Rank.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Rank,'n ranks','Residual in gamma region','Residual_Rank_%s'%(name),True)

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
        FilePath = "%s/Netflix_"%(folder_path)+"Syst_%s"%(root_file_tags[elev])+".root"
        FilePath_List += [FilePath]
        print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
        InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
        bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
        bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
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
            MakeSystematicPlots(hists, colors, fillstyles, legends, InfoTree.exposure_hours, name)

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
        lumilab1 = ROOT.TLatex(0.15,0.50,'#chi^{2} = #sum_{ij} (N^{nth}_{ij}-N^{20th}_{ij})^{2} / (#sum_{ij} N^{data}_{ij})^{2}' )
        lumilab1.SetNDC()
        lumilab1.SetTextSize(0.3)
        lumilab1.Draw()
        pad1.cd()
        Hist2D_Converge.GetYaxis().SetTitle('#sqrt{#chi^{2}}')
        Hist2D_Converge.GetXaxis().SetTitle('iterations')
        Hist2D_Converge.Draw("COL4Z")
        canvas.SaveAs('output_plots/Converge_%s.png'%(selection_tag))

def MakeGaussComparisonPlot(Hists,legends,colors,title,name):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.15)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

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
            Hists[h].GetXaxis().SetTitle(title)
            #Hists[h].GetXaxis().SetRangeUser(0,8)
            Hists[h].GetXaxis().SetRangeUser(-5,8)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #Hists[max_hist].SetMinimum(0)
    #Hists[max_hist].SetMaximum(1)
    Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            if h==0:
                Hists[h].Draw("hist same")
            else:
                Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s, mean = %.2f, RMS = %.2f'%(legends[h],Hists[h].GetMean(),Hists[h].GetRMS()),"pl")
            #legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    pad1.SetLogy()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeSignificanceDistribution(Hist2D_Sig,Hist2D_SR,Hist2D_Bkg,name):

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    for bx in range(0,Hist2D_Sig.GetNbinsX()):
        for by in range(0,Hist2D_Sig.GetNbinsY()):
            if Hist2D_SR.GetBinContent(bx+1,by+1)==0: continue
            content = Hist2D_Sig.GetBinContent(bx+1,by+1)
            Hist_Sig.Fill(content)
    Hist_Model = ROOT.TH1D("Hist_Model","",65,-5,8)
    Hist_Model.FillRandom("func",10000*int(Hist_Sig.GetEntries()))
    Hist_Model.Scale(1./10000.)
    Hist_Model.SetMinimum(0.5)
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Ref. Gaussian']
    color_list += [2]
    Hist_list += [Hist_Sig]
    legend_list += ['Data']
    color_list += [4]
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))

def SingleSourceAnalysis(source_list,doMap):

    global ErecS_lower_cut
    global ErecS_upper_cut

    FilePath_List = []
    ResetStackedShowerHistograms()
    for source in range(0,len(source_list)):
        source_name = source_list[source]
        for elev in range(0,len(root_file_tags)):
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
            InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
            InfoTree = InputFile.Get("InfoTree")
            InfoTree.GetEntry(0)

            MSCW_blind_cut = InfoTree.MSCW_cut_blind
            MSCL_blind_cut = InfoTree.MSCL_cut_blind
            bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
            bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
            bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
            bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
            for e in range(0,len(energy_bin)-1):
                ErecS_lower_cut = energy_bin[e]
                ErecS_upper_cut = energy_bin[e+1]
                if ErecS_upper_cut<=energy_bin[energy_bin_cut_low]: continue
                if ErecS_lower_cut>=energy_bin[energy_bin_cut_up]: continue
                GetShowerHistogramsFromFile(FilePath_List[len(FilePath_List)-1])
                StackShowerHistograms()
                NormalizeEnergyHistograms(FilePath_List[len(FilePath_List)-1])
                StackEnergyHistograms()
                NormalizeTheta2Histograms(FilePath_List[len(FilePath_List)-1])
                StackTheta2Histograms()
                NormalizeSkyMapHistograms(FilePath_List[len(FilePath_List)-1])
                StackSkymapHistograms()

    MatrixDecompositionDemo(source_name)

    GetSourceInfo(FilePath_List)
    GetBrightStarInfo(FilePath_List)

    CalculateSystError()
    #CorrectEnergyHistograms()
    #CorrectTheta2Histograms()
    #CorrectSkymapHistograms()

    PlotsStackedHistograms('%s%s'%(source_list[0],PercentCrab))

    MakeRankResidualPlots('%s%s'%(source_list[0],PercentCrab))

    if not doMap: return

    Hist_OnData_Skymap_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnData_Skymap_Galactic_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Galactic_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnData_VR_Skymap_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnData_VR_Skymap_Galactic_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Syst_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Galactic_Syst_Sum.Rebin2D(n_rebin,n_rebin)

    Hist_OnData_Skymap_smooth = Smooth2DMap(Hist_OnData_Skymap_Sum,smooth_size,False)
    Hist_OnData_VR_Skymap_smooth = Smooth2DMap(Hist_OnData_VR_Skymap_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_Syst_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Syst_Sum,smooth_size,False)
    Hist_OnData_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnData_Skymap_Galactic_Sum,smooth_size,False)
    Hist_OnData_VR_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnData_VR_Skymap_Galactic_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_Galactic_Syst_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Syst_Sum,smooth_size,False)

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnData_VR_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_smooth,'RA','Dec','Skymap_RaDec_MDM_%s%s'%(source_name,PercentCrab))

    Hist_Significance_Skymap = GetSignificanceMap(Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnData_VR_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_smooth,Syst_MDM)
    MakeSignificanceDistribution(Hist_Significance_Skymap,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,'SigDist_MDM_%s%s'%(source_name,PercentCrab))
    Hist_Significance_Skymap_NoSyst = GetSignificanceMap(Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnData_VR_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_smooth,0.)
    MakeSignificanceDistribution(Hist_Significance_Skymap_NoSyst,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,'SigDist_NoSyst_MDM_%s%s'%(source_name,PercentCrab))

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_smooth,Hist_OnBkgd_Skymap_smooth,Hist_OnData_VR_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_smooth,'RA','Dec','Skymap_Smooth_RaDec_MDM_%s%s'%(source_name,PercentCrab))
    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Galactic_smooth,Hist_OnBkgd_Skymap_Galactic_smooth,Hist_OnData_VR_Skymap_Galactic_smooth,Hist_OnBkgd_Skymap_Galactic_Syst_smooth,'gal. l.','gal. b.','Skymap_Smooth_Galactic_MDM_%s%s'%(source_name,PercentCrab))

    Hist_Significance_Skymap_smooth = GetSignificanceMap(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth,Hist_OnData_VR_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_smooth,Syst_MDM)


    init_x = source_ra
    init_y = source_dec
    init_x = 284
    init_y = 4.2
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_Sum, Hist_OnBkgd_Skymap_Sum, Hist_Significance_Skymap_smooth,5,init_x,init_y)
    print 'Excess (5 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (5 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (5 sigma) radius = %0.3f'%(excess_radius)
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_Sum, Hist_OnBkgd_Skymap_Sum, Hist_Significance_Skymap_smooth,4,init_x,init_y)
    print 'Excess (4 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (4 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (4 sigma) radius = %0.3f'%(excess_radius)
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_Sum, Hist_OnBkgd_Skymap_Sum, Hist_Significance_Skymap_smooth,3,init_x,init_y)
    print 'Excess (3 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (3 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (3 sigma) radius = %0.3f'%(excess_radius)

def FindSourceIndex(source_name):
    for source in range(0,len(sample_list)):
        if source_name==sample_list[source]:
            return source
    return 0

Hist_EffArea = ROOT.TH1D("Hist_EffArea","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_EffArea_Sum = ROOT.TH1D("Hist_EffArea_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
for source in range(0,len(sample_list)):
    source_idx = FindSourceIndex(sample_list[source])
    FilePath_Folder = []
    for elev in range(0,len(root_file_tags)):
        SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+".root"
        FilePath_Folder += [SourceFilePath]
        if not os.path.isfile(FilePath_Folder[elev]): 
            continue
        else:
            GetSourceInfo(FilePath_Folder)
print 'analysis cut: MSCL = %s, MSCW = %s'%(MSCL_blind_cut,MSCW_blind_cut)
MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
print 'plot range: MSCL = %s, MSCW = %s'%(MSCL_plot_upper,MSCW_plot_upper)

Hist2D_OnData_Sum = ROOT.TH2D("Hist2D_OnData_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd_Sum = ROOT.TH2D("Hist2D_OnBkgd_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark_Sum = ROOT.TH2D("Hist2D_OnDark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData = ROOT.TH2D("Hist2D_OnData","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd = ROOT.TH2D("Hist2D_OnBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark = ROOT.TH2D("Hist2D_OnDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCL_Sum = ROOT.TH1D("Hist_OnData_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL_Sum = ROOT.TH1D("Hist_OnDark_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCL = ROOT.TH1D("Hist_OnData_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL = ROOT.TH1D("Hist_OnBkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL = ROOT.TH1D("Hist_OnDark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCW_Sum = ROOT.TH1D("Hist_OnData_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnDark_MSCW_Sum = ROOT.TH1D("Hist_OnDark_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCW = ROOT.TH1D("Hist_OnData_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW = ROOT.TH1D("Hist_OnBkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnDark_MSCW = ROOT.TH1D("Hist_OnDark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist2D_Rank0 = ROOT.TH2D("Hist2D_Rank0","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1 = ROOT.TH2D("Hist2D_Rank1","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2 = ROOT.TH2D("Hist2D_Rank2","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3 = ROOT.TH2D("Hist2D_Rank3","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_AllRanks_Sum = ROOT.TH2D("Hist2D_AllRanks_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Sum = ROOT.TH2D("Hist2D_Rank0_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Sum = ROOT.TH2D("Hist2D_Rank1_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Sum = ROOT.TH2D("Hist2D_Rank2_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Sum = ROOT.TH2D("Hist2D_Rank3_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

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

Hist_OnData_Theta2_Sum = ROOT.TH1D("Hist_OnData_Theta2_Sum","",50,0,10)
Hist_OnBkgd_Theta2_Sum = ROOT.TH1D("Hist_OnBkgd_Theta2_Sum","",50,0,10)
Hist_OnData_Theta2 = ROOT.TH1D("Hist_OnData_Theta2","",50,0,10)
Hist_OnBkgd_Theta2 = ROOT.TH1D("Hist_OnBkgd_Theta2","",50,0,10)
Hist_OnData_Energy_Sum = ROOT.TH1D("Hist_OnData_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy = ROOT.TH1D("Hist_OnData_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy = ROOT.TH1D("Hist_OnBkgd_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Zenith_Sum = ROOT.TH1D("Hist_OnData_Zenith_Sum","",45,0,90)
Hist_OnBkgd_Zenith_Sum = ROOT.TH1D("Hist_OnBkgd_Zenith_Sum","",45,0,90)
Hist_OnData_Zenith = ROOT.TH1D("Hist_OnData_Zenith","",45,0,90)
Hist_OnBkgd_Zenith = ROOT.TH1D("Hist_OnBkgd_Zenith","",45,0,90)

Hist_OnData_Skymap = ROOT.TH2D("Hist_OnData_Skymap","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_VR_Skymap = ROOT.TH2D("Hist_OnData_VR_Skymap","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap = ROOT.TH2D("Hist_OnBkgd_Skymap","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_Skymap_Sum = ROOT.TH2D("Hist_OnData_Skymap_Sum","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_VR_Skymap_Sum = ROOT.TH2D("Hist_OnData_VR_Skymap_Sum","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Sum","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_Sum","",75,source_ra-Skymap_size,source_ra+Skymap_size,75,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_Skymap_Galactic = ROOT.TH2D("Hist_OnData_Skymap_Galactic","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnData_VR_Skymap_Galactic = ROOT.TH2D("Hist_OnData_VR_Skymap_Galactic","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Syst = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnData_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnData_Skymap_Galactic_Sum","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnData_VR_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnData_VR_Skymap_Galactic_Sum","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Sum","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Syst_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst_Sum","",75,source_l-Skymap_size,source_l+Skymap_size,75,source_b-Skymap_size,source_b+Skymap_size)

Hist_SystErr_MSCL = ROOT.TH1D("Hist_SystErr_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_SystErr_MSCW = ROOT.TH1D("Hist_SystErr_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

print 'MJD_Start = %s'%(MJD_Start)
print 'MJD_End = %s'%(MJD_End)
time = Time(MJD_Start, format='mjd')
time.format = 'decimalyear'
year_start = time.value
time = Time(MJD_End, format='mjd')
time.format = 'decimalyear'
year_end = time.value
print 'year_start = %s'%(year_start)
print 'year_end = %s'%(year_end)
print 'roi_ra = %s'%(roi_ra[0])
print 'roi_dec = %s'%(roi_dec[0])
Hist_OnData_RoI_Energy_Sum = [] 
Hist_OnBkgd_RoI_Energy_Sum = [] 
Hist_OnData_RoI_Energy = [] 
Hist_OnBkgd_RoI_Energy = [] 
Hist_OnData_RoI_Theta2_Sum = []
Hist_OnBkgd_RoI_Theta2_Sum = []
Hist_OnData_RoI_Theta2 = []
Hist_OnBkgd_RoI_Theta2 = []
Hist_OnData_RoI_MJD = [] 
Hist_OnBkgd_RoI_MJD = [] 
Hist_OnData_RoI_MJD_Sum = [] 
Hist_OnBkgd_RoI_MJD_Sum = [] 
for nth_roi in range(0,len(roi_ra)):
    Hist_OnData_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnData_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Energy += [ROOT.TH1D("Hist_OnData_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnData_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_Theta2 += [ROOT.TH1D("Hist_OnData_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2 += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_MJD += [ROOT.TH1D("Hist_OnData_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnData_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnData_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]


Hist2D_OffData = []
Hist2D_OffData_Sum = []
Hist2D_OffBkgd = []
Hist2D_OffBkgd_Sum = []
Hist2D_OnSyst = []
Hist2D_OnSyst_Sum = []
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
    Hist2D_OffBkgd += [ROOT.TH2D("Hist2D_OffBkgd_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffBkgd_Sum += [ROOT.TH2D("Hist2D_OffBkgd_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OnSyst += [ROOT.TH2D("Hist2D_OnSyst_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OnSyst_Sum += [ROOT.TH2D("Hist2D_OnSyst_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
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

#n_rebin = 1
#smooth_size = 0.05
n_rebin = 1
smooth_size = 0.1
#n_rebin = 1
#smooth_size = 0.15
if energy_bin[energy_bin_cut_low]>500.:
    #n_rebin = 2
    smooth_size = 0.15
if energy_bin[energy_bin_cut_low]>1000.:
    #n_rebin = 4
    smooth_size = 0.2

GetGammaSourceInfo()

#SystematicAnalysis()

SingleSourceAnalysis(sample_list,True)
print 'n_good_matches = %s'%(n_good_matches)

print "Syst_MDM = %s"%(Syst_MDM) 
print "Syst_Init = %s"%(Syst_Init) 
print "Syst_Redu = %s"%(Syst_Redu) 
print "Syst_Clos = %s"%(Syst_Clos) 
