
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

root_file_tags = []
root_file_tags += ['16bins_constrained']
#root_file_tags += ['16bins_unconstrained']

selection_tag = root_file_tags[0]

folder_path = 'output_root'
PercentCrab = ''

energy_fine_bin_cut_low = 0
energy_fine_bin_cut_up = 20
selection_tag += '_E%s'%(energy_fine_bin_cut_low)

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
source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
n_control_samples = 6
MJD_Start = 2147483647
MJD_End = 0
roi_ra = ROOT.std.vector("double")(10)
roi_dec = ROOT.std.vector("double")(10)
roi_radius = ROOT.std.vector("double")(10)

Syst_MDM = 0.02

energy_list = []
energy_list += [int(pow(10,2.0))]
energy_list += [int(pow(10,4.0))]

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

sample_list = []
sky_coord = []

#sample_list += ['OJ287V6']
#sky_coord += ['08 54 49.1 +20 05 58.89']

#sample_list += ['WComaeV6']
#sky_coord += ['12 21 31.7 +28 13 59']
#sample_list += ['WComaeV5']
#sky_coord += ['12 21 31.7 +28 13 59']
#sample_list += ['WComaeV4']
#sky_coord += ['12 21 31.7 +28 13 59']

#sample_list += ['IC443HotSpotV6']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#sample_list += ['IC443HotSpotV5']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#sample_list += ['IC443HotSpotV4']
#sky_coord += ['06 18 2.700 +22 39 36.00']

sample_list += ['MGRO_J1908_V6']
sky_coord += ['19 07 54 +06 16 07']
sample_list += ['MGRO_J1908_V5']
sky_coord += ['19 07 54 +06 16 07']
sample_list += ['MGRO_J1908_V4']
sky_coord += ['19 07 54 +06 16 07']

#sample_list += ['Segue1V6']
#sky_coord += ['10 07 04 +16 04 55']
#sample_list += ['Segue1V5']
#sky_coord += ['10 07 04 +16 04 55']

#sample_list += ['Crab']
#sky_coord += ['05 34 31.97 +22 00 52.1']
#sample_list += ['CrabV5']
#sky_coord += ['05 34 31.97 +22 00 52.1']
#sample_list += ['Mrk421']
#sky_coord += ['11 04 19 +38 11 41']
#sample_list += ['H1426']
#sky_coord += ['14 28 32.609 +42 40 21.05']
#sample_list += ['PKS1424']
#sky_coord += ['14 27 00 +23 47 00']
#sample_list += ['3C264']
#sky_coord += ['11 45 5.009 +19 36 22.74']
#sample_list += ['1ES0229']
#sky_coord += ['02 32 53.2 +20 16 21']
#sample_list += ['S3_1227_V6']
#sky_coord += ['12 30 14.1 +25 18 07']
#sample_list += ['MS1221V6']
#sky_coord += ['12 24 24.2 +24 36 24']
#sample_list += ['PKS1441V6']
#sky_coord += ['14 43 56.9 +25 01 44']
#sample_list += ['RBS0413V6']
#sky_coord += ['03 19 47 +18 45 42']
#sample_list += ['PG1553V6']
#sky_coord += ['15 55 44.7 +11 11 41']
#sample_list += ['ComaV6']
#sky_coord += ['12 59 48.7 +27 58 50']
#sample_list += ['1ES1011V6']
#sky_coord += ['10 15 4.139 +49 26 0.71']
#sample_list += ['NGC1275V6']
#sky_coord += ['03 19 48.1 +41 30 42']
#sample_list += ['1ES0647V6']
#sky_coord += ['06 50 46.490 +25 02 59.62']
#sample_list += ['1ES1440V6']
#sky_coord += ['14 42 48.277 +12 00 40.37']
#sample_list += ['1ES1741V6']
#sky_coord += ['17 44 01.2 +19 32 47']
#sample_list += ['RGBJ0710']
#sky_coord += ['07 10 26.4 +59 09 00']
#sample_list += ['CasA']
#sky_coord += ['23 23 13.8 +58 48 26']
#sample_list += ['M82']
#sky_coord += ['09 55 52.7 +69 40 46']
#sample_list += ['G079']
#sky_coord += ['20 32 28.56 +40 19 41.52']
#sample_list += ['1ES1218V6']
#sky_coord += ['12 21 26.3 +30 11 29']
#sample_list += ['CygnusV6']
#sky_coord += ['20 18 35.03 +36 50 00.0']
#sample_list += ['GemingaV6']
#sky_coord += ['06 32 28 +17 22 00']
#sample_list += ['GemingaV5']
#sky_coord += ['06 32 28 +17 22 00']
#sample_list += ['SgrAV6']
#sky_coord += ['17 45 39.6 -29 00 22']
#sample_list += ['TychoV6']
#sky_coord += ['00 25 21.6 +64 07 48']
#sample_list += ['CTA1V5']
#sky_coord += ['00 06 26 +72 59 01.0']
#sample_list += ['2HWC_J1953V6_new']
#sky_coord += ['19 53 02.4 +29 28 48']
#sample_list += ['2HWC_J1930V6']
#sky_coord += ['19 30 32 +18 52 12']

other_stars = []
other_star_coord = []
#other_stars += ['PSR J1907']
#other_star_coord += [[286.978,6.038]]
#other_stars += ['SNR G40.5']
#other_star_coord += [[286.786,6.498]]
#other_stars += ['ARGO J1910']
#other_star_coord += [[287.650,7.350]]
#other_stars += ['PSR J1954']
#other_star_coord += [[298.830,28.590]]
#other_stars += ["4FGL J2017.9"]
#other_star_coord += [[304.490833333,36.4277777778]]
#other_stars += ["4FGL J2021.1"]
#other_star_coord += [[305.279583333,36.8561111111]]
#other_stars += ["4FGL J2015.5"]
#other_star_coord += [[303.8925,37.1761111111]]
#other_stars += ["4FGL J2021.9"]
#other_star_coord += [[305.48625,36.1577777778]]
#other_stars += ["4FGL J2013.5"]
#other_star_coord += [[303.384583333,36.2252777778]]
#other_stars += ["4FGL J2017.3"]
#other_star_coord += [[304.349166667,35.4225]]
#other_stars += ["4FGL J2022.3"]
#other_star_coord += [[305.584583333,38.6711111111]]

bright_star_ra = []
bright_star_dec = []
bright_star_brightness = []

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b;
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

def GetBrightStarInfo(file_list):

    global bright_star_ra
    global bright_star_dec
    global bright_star_brightness

    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        StarTree = InputFile.Get("StarTree")
        for entry in range(0,StarTree.GetEntries()):
            StarTree.GetEntry(entry)
            if StarTree.star_brightness>5.0: continue
            distance_to_center = pow(pow(source_ra-StarTree.star_ra,2)+pow(source_dec-StarTree.star_dec,2),0.5)
            if distance_to_center>2.: continue
            bright_star_ra += [StarTree.star_ra]
            bright_star_dec += [StarTree.star_dec]
            bright_star_brightness += [StarTree.star_brightness]
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

    for nth_sample in range(0,n_control_samples-1):

        Hist2D_OffData_Sum[nth_sample].Reset()
        Hist2D_OffBkgd_Sum[nth_sample].Reset()
        Hist2D_OffDark_Sum[nth_sample].Reset()
        Hist_OffData_MSCL_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCL_Sum[nth_sample].Reset()
        Hist_OffDark_MSCL_Sum[nth_sample].Reset()
        Hist_OffData_MSCW_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCW_Sum[nth_sample].Reset()
        Hist_OffDark_MSCW_Sum[nth_sample].Reset()

def GetSourceInfo(file_list):

    global N_bins_for_deconv
    global MSCW_blind_cut
    global MSCL_blind_cut
    global n_good_matches
    global exposure_hours
    global source_ra
    global source_dec
    global source_l
    global source_b
    global n_control_samples
    global MJD_Start
    global MJD_End
    global roi_ra
    global roi_dec
    global roi_radius

    n_good_matches = 0
    exposure_hours = 0.
    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.SetBranchAddress('roi_ra',ROOT.AddressOf(roi_ra))
        InfoTree.SetBranchAddress('roi_dec',ROOT.AddressOf(roi_dec))
        InfoTree.SetBranchAddress('roi_radius',ROOT.AddressOf(roi_radius))
        InfoTree.GetEntry(0)
        n_good_matches += InfoTree.n_good_matches
        exposure_hours += InfoTree.exposure_hours
        N_bins_for_deconv = InfoTree.N_bins_for_deconv
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        source_ra = InfoTree.mean_tele_point_ra
        source_dec = InfoTree.mean_tele_point_dec
        source_l = InfoTree.mean_tele_point_l
        source_b = InfoTree.mean_tele_point_b
        n_control_samples = InfoTree.n_control_samples
        MJD_Start = min(InfoTree.MJD_Start,MJD_Start)
        MJD_End = max(InfoTree.MJD_End,MJD_End)
        HistName = "Hist_EffArea"
        Hist_EffArea.Reset()
        Hist_EffArea.Add(InputFile.Get(HistName))
        Hist_EffArea.Scale(exposure_hours*3600.)
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

    HistName = "Hist_Rank0_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_Rank0.Add(InputFile.Get(HistName))
    Hist2D_Rank1.Add(InputFile.Get(HistName))
    Hist2D_Rank2.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank1_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_Rank1.Add(InputFile.Get(HistName))
    Hist2D_Rank2.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank2_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_Rank2.Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist2D_OnData.Reset()
        Hist2D_OnDark.Reset()
        Hist2D_OnBkgd.Reset()
        Hist2D_Rank0.Reset()
        Hist2D_Rank1.Reset()
        Hist2D_Rank2.Reset()

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

    for nth_sample in range(0,n_control_samples-1):

        HistName = "Hist_OffData_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffData[nth_sample].Reset()
        Hist2D_OffData[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OffDark_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffDark[nth_sample].Reset()
        Hist2D_OffDark[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OffBkgd_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffBkgd[nth_sample].Reset()
        Hist2D_OffBkgd[nth_sample].Add(InputFile.Get(HistName))

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist2D_OffData[nth_sample].Reset()
            Hist2D_OffDark[nth_sample].Reset()
            Hist2D_OffBkgd[nth_sample].Reset()

        Hist_OffData_MSCL[nth_sample].Reset()
        Hist_OffData_MSCL[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionX("Hist1D_OffData_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffDark_MSCL[nth_sample].Reset()
        Hist_OffDark_MSCL[nth_sample].Add(Hist2D_OffDark[nth_sample].ProjectionX("Hist1D_OffDark_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffBkgd_MSCL[nth_sample].Reset()
        Hist_OffBkgd_MSCL[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionX("Hist1D_OffBkgd_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffData_MSCW[nth_sample].Reset()
        Hist_OffData_MSCW[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionY("Hist1D_OffData_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

        Hist_OffDark_MSCW[nth_sample].Reset()
        Hist_OffDark_MSCW[nth_sample].Add(Hist2D_OffDark[nth_sample].ProjectionY("Hist1D_OffDark_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

        Hist_OffBkgd_MSCW[nth_sample].Reset()
        Hist_OffBkgd_MSCW[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionY("Hist1D_OffBkgd_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

def StackShowerHistograms():

    Hist2D_OnData_Sum.Add(Hist2D_OnData)
    Hist2D_OnDark_Sum.Add(Hist2D_OnDark)
    Hist2D_OnBkgd_Sum.Add(Hist2D_OnBkgd)

    Hist2D_Rank0_Sum.Add(Hist2D_Rank0)
    Hist2D_Rank1_Sum.Add(Hist2D_Rank1)
    Hist2D_Rank2_Sum.Add(Hist2D_Rank2)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank0)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank1)
    Hist2D_AllRanks_Sum.Add(Hist2D_Rank2)

    Hist_OnData_MSCL_Sum.Add(Hist_OnData_MSCL)
    Hist_OnData_MSCW_Sum.Add(Hist_OnData_MSCW)
    Hist_OnDark_MSCL_Sum.Add(Hist_OnDark_MSCL)
    Hist_OnDark_MSCW_Sum.Add(Hist_OnDark_MSCW)
    Hist_OnBkgd_MSCL_Sum.Add(Hist_OnBkgd_MSCL)
    Hist_OnBkgd_MSCW_Sum.Add(Hist_OnBkgd_MSCW)

    for nth_sample in range(0,n_control_samples-1):

        Hist2D_OffData_Sum[nth_sample].Add(Hist2D_OffData[nth_sample])
        Hist2D_OffDark_Sum[nth_sample].Add(Hist2D_OffDark[nth_sample])
        Hist2D_OffBkgd_Sum[nth_sample].Add(Hist2D_OffBkgd[nth_sample])

        Hist_OffData_MSCL_Sum[nth_sample].Add(Hist_OffData_MSCL[nth_sample])
        Hist_OffData_MSCW_Sum[nth_sample].Add(Hist_OffData_MSCW[nth_sample])
        Hist_OffDark_MSCL_Sum[nth_sample].Add(Hist_OffDark_MSCL[nth_sample])
        Hist_OffDark_MSCW_Sum[nth_sample].Add(Hist_OffDark_MSCW[nth_sample])
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
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetBorderMode(0)
    if logy: pad1.SetGrid()
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
    legend = ROOT.TLegend(0.1,0.1,0.9,0.9)
    legend.SetNColumns(2)
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

def MakeLightCurvePlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.0,1,0.8)
    pad1.SetBottomMargin(0.39)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    time = Time(MJD_Start, format='mjd')
    time.format = 'decimalyear'
    year_start = int(time.value)
    time = Time(MJD_End, format='mjd')
    time.format = 'decimalyear'
    year_end = int(time.value)

    Hist_data_year = []
    for hist in range(0,len(Hist_data)):
        Hist_data_year += [ROOT.TH1D("Hist_data_year_%s"%(hist),"",year_end-year_start+1,year_start,year_end+1)]
        for binx in range(0,Hist_data[hist].GetNbinsX()):
            mjd = Hist_data[hist].GetBinCenter(binx+1)
            time = Time(mjd, format='mjd')
            time.format = 'decimalyear'
            year = int(time.value)
            Hist_data_year[hist].Fill(year,Hist_data[hist].GetBinContent(binx+1))
        for binx in range(0,Hist_data_year[hist].GetNbinsX()):
            Hist_data_year[hist].SetBinError(binx+1,pow(Hist_data_year[hist].GetBinContent(binx+1),0.5))
    Hist_bkgd_year = []
    for hist in range(0,len(Hist_data)):
        Hist_bkgd_year += [ROOT.TH1D("Hist_bkgd_year_%s"%(hist),"",year_end-year_start+1,year_start,year_end+1)]
        for binx in range(0,Hist_data[hist].GetNbinsX()):
            mjd = Hist_data[hist].GetBinCenter(binx+1)
            time = Time(mjd, format='mjd')
            time.format = 'decimalyear'
            year = int(time.value)
            Hist_bkgd_year[hist].Fill(year,Hist_bkgd[hist].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_year[hist].GetNbinsX()):
            Hist_bkgd_year[hist].SetBinError(binx+1,Hist_bkgd_year[hist].GetBinContent(binx+1)*Syst_MDM)
    Hist_ratio_year = []
    for hist in range(0,len(Hist_data_year)):
        Hist_ratio_year += [ROOT.TH1D("Hist_ratio_year_%s"%(hist),"",year_end-year_start+1,year_start,year_end+1)]
        Hist_ratio_year[hist].Add(Hist_data_year[hist])
        Hist_ratio_year[hist].Add(Hist_bkgd_year[hist],-1.)
        #Hist_ratio_year[hist].Divide(Hist_bkgd_year[hist])

    max_heigh = 0
    max_hist = 0
    for h in range(0,len(Hist_ratio_year)):
        if Hist_ratio_year[h]!=0:
            Hist_ratio_year[h].SetLineColor(colors[h])
            Hist_ratio_year[h].GetXaxis().SetTitle(title)
            if max_heigh < Hist_ratio_year[h].GetMaximum(): 
                max_heigh = Hist_ratio_year[h].GetMaximum()
                max_hist = h

    Hist_ratio_year[max_hist].Draw("E")
    for h in range(0,len(Hist_ratio_year)):
        Hist_ratio_year[h].SetLineWidth(3)
        Hist_ratio_year[h].Draw("E same")

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
    for h in range(0,len(Hist_ratio_year)):
        legend.AddEntry(Hist_ratio_year[h],legends[h],"pl")
    legend.Draw("SAME")

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
    Hists[max_hist].Draw("E")

    Hist_Bkgd = Hists[1].Clone()
    Hist_Bkgd.Reset()
    Hist_Syst = Hists[1].Clone()
    Hist_Syst.Reset()
    found_syst_hist = False
    for h in range(0,len(Hists)):
        if stack_it[h] or 'syst' in legends[h]:
            Hist_Syst.Add(Hists[h])
            found_syst_hist = True
        if colors[h]==4:
            Hist_Bkgd.Add(Hists[h].Clone())
    if not found_syst_hist:
        for binx in range(0,Hist_Syst.GetNbinsX()):
            syst_err = Syst_MDM*Hist_Bkgd.GetBinContent(binx+1)
            stat_err = Hist_Bkgd.GetBinError(binx+1)
            Hist_Syst.SetBinError(binx+1,pow(syst_err*syst_err+stat_err*stat_err,0.5))

    fill_color = [0,0,0,30,38,46]
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
    lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(energy_fine_bin[energy_fine_bin_cut_low],exposure_hours) )
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
    Hists += [Hist_OnDark_MSCW_Sum]
    legends += ['init. bkg.']
    colors += [2]
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
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_MSCL_Sum]
    legends += ['init. bkg.']
    colors += [2]
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
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI_Theta2_MDM_V%s_%s'%(nth_roi,tag)
        title = 'squared angle from source location #theta^{2}'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI_Energy_MDM_V%s_%s'%(nth_roi,tag)
        title = 'energy [GeV]'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    legends = []
    colors = []
    for nth_roi in range(0,len(roi_ra)):
        legends += ['RoI %s'%(nth_roi)]
        colors += [nth_roi+1]
        Hist_data_mjd += [Hist_OnData_RoI_MJD_Sum[nth_roi]]
        Hist_bkgd_mjd += [Hist_OnBkgd_RoI_MJD_Sum[nth_roi]]
    plotname = 'LightCurve_RoI_Year_MDM_V%s_%s'%(nth_roi,tag)
    title = 'Year'
    MakeLightCurvePlot(Hist_data_mjd,Hist_bkgd_mjd,legends,colors,title,plotname)

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


    for nth_sample in range(0,n_control_samples-1):

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OffData_MSCW_Sum[nth_sample]]
        legends += ['obs. data']
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OffDark_MSCW_Sum[nth_sample]]
        legends += ['init. bkg.']
        colors += [2]
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
        Hists += [Hist_OffDark_MSCL_Sum[nth_sample]]
        legends += ['init. bkg.']
        colors += [2]
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
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,1.,-1)

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
    Syst_MDM = 0.
    Hist_SystErr_MSCL.Reset()
    Hist_SystErr_MSCW.Reset()

    n_samples_used = 0.
    for nth_sample in range(0,n_control_samples-1):

        norm_bin_low_target = Hist_OffData_MSCL_Sum[nth_sample].FindBin(MSCL_lower_cut)
        norm_bin_up_target = Hist_OffData_MSCL_Sum[nth_sample].FindBin(MSCL_blind_cut)-1
        Total_Bkgd = Hist_OffData_MSCL_Sum[nth_sample].Integral(norm_bin_low_target,norm_bin_up_target)
        Hist_Diff_MSCL = Hist_OffData_MSCL_Sum[nth_sample].Clone()
        Hist_Diff_MSCL.Add(Hist_OffBkgd_MSCL_Sum[nth_sample],-1.)
        if Total_Bkgd == 0.: continue
        syst_this = Hist_Diff_MSCL.Integral(norm_bin_low_target,norm_bin_up_target)/Total_Bkgd
        print '%s th sample syst = %s'%(nth_sample,syst_this)
        n_samples_used += 1.
        Syst_MDM += pow(syst_this,2)
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

    Syst_MDM = pow(Syst_MDM/n_samples_used,0.5)
    print "Syst_MDM = %s"%(Syst_MDM) 

    for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
            old_syst = Hist_SystErr_MSCL.GetBinError(binx+1)
            new_syst = pow(old_syst/(n_control_samples-1.),0.5)
            Hist_SystErr_MSCL.SetBinError(binx+1,new_syst)

    for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
            old_syst = Hist_SystErr_MSCW.GetBinError(binx+1)
            new_syst = pow(old_syst/(n_control_samples-1.),0.5)
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

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Energy[nth_roi].Reset()
        Hist_OnData_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
        Hist_OnBkgd_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Energy.Reset()
        Hist_OnBkgd_Energy.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Energy[nth_roi].Reset()
            Hist_OnBkgd_RoI_Energy[nth_roi].Reset()

    bkg_total, bkg_err = IntegralAndError(Hist_OnBkgd_MSCW,bin_lower,bin_upper)
    old_integral = Hist_OnBkgd_Energy.Integral()
    bkgd_scale = 0
    bkgd_scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        bkgd_scale = bkg_total/old_integral
        bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
    else:
        bkgd_scale = 0
        bkgd_scale_err = 0
    if not bkg_total==0:
        Theta2HistScale(Hist_OnBkgd_Energy,bkgd_scale,bkgd_scale_err)
        for nth_roi in range(0,len(roi_ra)):
            Theta2HistScale(Hist_OnBkgd_RoI_Energy[nth_roi],bkgd_scale,bkgd_scale_err)
    else:
        Hist_OnBkgd_Energy.Scale(0)
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnBkgd_RoI_Energy[nth_roi].Scale(0)

    for nth_sample in range(0,n_control_samples-1):

        HistName = "Hist_OffData_SR_Energy_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffData_Energy[nth_sample].Reset()
        Hist_OffData_Energy[nth_sample].Add(InputFile.Get(HistName))
        HistName = "Hist_OffData_CR_Energy_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffBkgd_Energy[nth_sample].Reset()
        Hist_OffBkgd_Energy[nth_sample].Add(InputFile.Get(HistName))

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist_OffData_Energy[nth_sample].Reset()
            Hist_OffBkgd_Energy[nth_sample].Reset()

        bkg_total, bkg_err = IntegralAndError(Hist_OffBkgd_MSCW[nth_sample],bin_lower,bin_upper)
        old_integral = Hist_OffBkgd_Energy[nth_sample].Integral()
        bkgd_scale = 0
        bkgd_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            bkgd_scale = bkg_total/old_integral
            bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        else:
            bkgd_scale = 0
            bkgd_scale_err = 0
        if not bkg_total==0:
            Theta2HistScale(Hist_OffBkgd_Energy[nth_sample],bkgd_scale,bkgd_scale_err)
        else:
            Hist_OffBkgd_Energy[nth_sample].Scale(0)

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
    HistName = "Hist_OnData_CR_CameraFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_R2off.Reset()
    Hist_OnBkgd_R2off.Add(InputFile.Get(HistName))

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

    bkg_total, bkg_err = IntegralAndError(Hist_OnBkgd_Energy,bin_lower,bin_upper)
    if bkg_total==0.:
        Hist_OnData_Theta2.Reset()
        Hist_OnBkgd_Theta2.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Theta2[nth_roi].Reset()
            Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()
        return

    old_integral = Hist_OnBkgd_R2off.Integral()
    bkgd_scale = 0
    bkgd_scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        bkgd_scale = bkg_total/old_integral
        bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
    else:
        bkgd_scale = 0
        bkgd_scale_err = 0
    if not bkg_total==0:
        Theta2HistScale(Hist_OnBkgd_Theta2,bkgd_scale,bkgd_scale_err)
        for nth_roi in range(0,len(roi_ra)):
            Theta2HistScale(Hist_OnBkgd_RoI_Theta2[nth_roi],bkgd_scale,bkgd_scale_err)
    else:
        Hist_OnBkgd_Theta2.Scale(0)
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnBkgd_RoI_Theta2[nth_roi].Scale(0)

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_MJD[nth_roi].Reset()
        Hist_OnData_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_MJD[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_MJD[nth_roi].Reset()
            Hist_OnBkgd_RoI_MJD[nth_roi].Reset()

    for nth_roi in range(0,len(roi_ra)):
        bkg_total, bkg_err = IntegralAndError(Hist_OnBkgd_RoI_Energy[nth_roi],bin_lower,bin_upper)
        if bkg_total==0.:
            Hist_OnData_RoI_MJD[nth_roi].Reset()
            Hist_OnBkgd_RoI_MJD[nth_roi].Reset()
            return

        old_integral = Hist_OnBkgd_RoI_MJD[nth_roi].Integral()
        bkgd_scale = 0
        bkgd_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            bkgd_scale = bkg_total/old_integral
            bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        else:
            bkgd_scale = 0
            bkgd_scale_err = 0
        if not bkg_total==0:
            Hist_OnBkgd_RoI_MJD[nth_roi].Scale(bkgd_scale)
        else:
            Hist_OnBkgd_RoI_MJD[nth_roi].Scale(0)


def StackEnergyHistograms():

    Hist_OnData_Energy_Sum.Add(Hist_OnData_Energy)
    Hist_OnBkgd_Energy_Sum.Add(Hist_OnBkgd_Energy)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Add(Hist_OnData_RoI_Energy[nth_roi])
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Energy[nth_roi])

    for nth_sample in range(0,n_control_samples-1):

        Hist_OffData_Energy_Sum[nth_sample].Add(Hist_OffData_Energy[nth_sample])
        Hist_OffBkgd_Energy_Sum[nth_sample].Add(Hist_OffBkgd_Energy[nth_sample])

def StackTheta2Histograms():

    Hist_OnData_Theta2_Sum.Add(Hist_OnData_Theta2)
    Hist_OnBkgd_Theta2_Sum.Add(Hist_OnBkgd_Theta2)
    Hist_OnBkgd_R2off_Sum.Add(Hist_OnBkgd_R2off)

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

    for nth_sample in range(0,n_control_samples-1):

        HistName = "Hist_OffData_SR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffData_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))
        HistName = "Hist_OffData_CR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()

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

    for nth_sample in range(0,n_control_samples-1):
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
    HistName = "Hist_OnData_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap.Reset()
    Hist_OnBkgd_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Galactic.Reset()
    Hist_OnBkgd_Skymap_Galactic.Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Skymap.Reset()
        Hist_OnData_Skymap_Galactic.Reset()
        Hist_OnBkgd_Skymap.Reset()
        Hist_OnBkgd_Skymap_Galactic.Reset()

    bkg_total, bkg_err = IntegralAndError(Hist_OnBkgd_Energy,bin_lower,bin_upper)
    if bkg_total==0.:
        Hist_OnData_Skymap.Reset()
        Hist_OnData_Skymap_Galactic.Reset()
        Hist_OnBkgd_Skymap.Reset()
        Hist_OnBkgd_Skymap_Galactic.Reset()
        return

    old_integral = Hist_OnBkgd_R2off.Integral()
    scale = 0
    scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        scale = bkg_total/old_integral
        scale_err = scale*(bkg_err/bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not bkg_total==0:
        RaDecHistScale(Hist_OnBkgd_Skymap,scale,scale_err)
        RaDecHistScale(Hist_OnBkgd_Skymap_Galactic,scale,scale_err)
    else:
        Hist_OnBkgd_Skymap.Scale(0)
        Hist_OnBkgd_Skymap_Galactic.Scale(0)

def StackSkymapHistograms():

    Hist_OnData_Skymap_Sum.Add(Hist_OnData_Skymap)
    Hist_OnData_Skymap_Galactic_Sum.Add(Hist_OnData_Skymap_Galactic)
    Hist_OnBkgd_Skymap_Sum.Add(Hist_OnBkgd_Skymap)
    Hist_OnBkgd_Skymap_Galactic_Sum.Add(Hist_OnBkgd_Skymap_Galactic)

def Smooth2DMap(Hist_Old,smooth_size,addLinearly):

    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(4*smooth_size/bin_size) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            bin_content = 0
            bin_error = 0
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
                            if not addLinearly:
                                bin_error += pow(ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2),2)
                            else:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    return Hist_Smooth

def GetSignificanceMap(Hist_SR,Hist_Bkg,syst_method):

    Hist_Skymap = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+pow(syst_method*Hist_Bkg.GetBinContent(bx+1,by+1),2),0.5)
            if syst_method==0.: NBkg_Err = pow(Hist_Bkg.GetBinContent(bx+1,by+1),0.5)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
    return Hist_Skymap

def reflectXaxis(hist):

    # taken from VPlotAnasumHistograms.cpp
	
    # temporary histogram
    hT = ROOT.TH2D( "%s_REFLECTED"%(hist.GetName()), "", hist.GetNbinsX(), -1.*hist.GetXaxis().GetXmax(), -1.*hist.GetXaxis().GetXmin(), hist.GetNbinsY(), hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax() )
    hT.SetStats( 0 )
    hT.SetXTitle( hist.GetXaxis().GetTitle() );
    hT.SetYTitle( hist.GetYaxis().GetTitle() );
	
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsX()+1):
            hT.SetBinContent( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinContent( binx, biny ) )
    return hT

def Make2DSignificancePlot(syst_method,Hist_SR,Hist_Bkg,xtitle,ytitle,name):

    max_sig = 0.
    Hist_Skymap = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+pow(syst_method*Hist_Bkg.GetBinContent(bx+1,by+1),2),0.5)
            if syst_method==0.: NBkg_Err = pow(Hist_Bkg.GetBinContent(bx+1,by+1),0.5)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            if Sig>max_sig: max_sig = Sig
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
    Hist_Skymap = reflectXaxis(Hist_Skymap)

    other_star_labels = []
    other_star_markers = []
    other_star_names = []
    other_star_significance = []
    bright_star_labels = []
    bright_star_markers = []
    if xtitle=="RA":
        for star in range(0,len(other_stars)):
            if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>3.0*3.0: continue
            other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
            other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.02)
            other_star_names += [other_stars[star]]
            binx = Hist_Skymap.GetXaxis().FindBin(-other_star_coord[star][0])
            biny = Hist_Skymap.GetYaxis().FindBin(other_star_coord[star][1])
            other_star_significance += [Hist_Skymap.GetBinContent(binx,biny)]
        for star in range(0,len(bright_star_ra)):
            bright_star_markers += [ROOT.TMarker(-bright_star_ra[star],bright_star_dec[star],30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-bright_star_ra[star]-0.15,bright_star_dec[star]+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
    else:
        for star in range(0,len(other_stars)):
            gal_l, gal_b = ConvertRaDecToGalactic(other_star_coord[star][0],other_star_coord[star][1])
            if pow(source_l-gal_l,2)+pow(source_b-gal_b,2)>3.0*3.0: continue
            other_star_markers += [ROOT.TMarker(-gal_l,gal_b,2)]
            other_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.02)
            other_star_names += [other_stars[star]]
            binx = Hist_Skymap.GetXaxis().FindBin(-gal_l)
            biny = Hist_Skymap.GetYaxis().FindBin(gal_b)
            other_star_significance += [Hist_Skymap.GetBinContent(binx,biny)]
        for star in range(0,len(bright_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(bright_star_ra[star],bright_star_dec[star])
            bright_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)

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
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'max. %0.1f#sigma (syst = %0.1f%%)'%(max_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.50,'exposure %0.1f hrs'%(exposure_hours) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.30,'E >%0.1f GeV (%.1f hrs)'%(energy_fine_bin[energy_fine_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/SkymapSig_%s_%s.png'%(name,selection_tag))

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
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)<30: 
                Hist_Skymap_Ratio.SetBinContent(bx+1,by+1,0.)
    Hist_Skymap_Ratio = reflectXaxis(Hist_Skymap_Ratio)

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
    Hist_Skymap_zoomin = ROOT.TH2D("Hist_Skymap_zoomin","",50,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,50,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_Contour_zoomin = ROOT.TH2D("Hist_Contour_zoomin","",50,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,50,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
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

def GetExtention(Hist_data, Hist_bkgd, Hist_sig, highlight_threshold):

    Hist_Excess = Hist_data.Clone()
    Hist_Excess.Add(Hist_bkgd,-1.)

    maxbin = Hist_Excess.GetMaximumBin()
    xx, yy, zz = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
    Hist_Excess.GetBinXYZ(maxbin, xx, yy, zz)
    print 'max excess at %s, %s'%(Hist_Excess.GetXaxis().GetBinCenter(xx),Hist_Excess.GetYaxis().GetBinCenter(yy))

    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
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

    return excess_center_x_init, excess_center_y_init, excess_radius_init
    #return excess_center_x, excess_center_y, excess_radius

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

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M_{Data}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_OnData_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_OnData_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/OnData_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M = #sum_{i=1}^{i=1} #lambda_{i} q_{i} p_{i}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank0_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank0_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank0_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Rank0_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M = #sum_{i=1}^{i=2} #lambda_{i} q_{i} p_{i}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank1_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank1_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank1_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Rank1_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M = #sum_{i=1}^{i=3} #lambda_{i} q_{i} p_{i}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank2_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank2_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank2_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Rank2_%s_%s.png'%(name,selection_tag))

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
    lumilab2 = ROOT.TLatex(0.15,0.50,'#epsilon = 1/#sum_{#gamma region} M_{data} (#sum_{#gamma region} M_{data} - #sum_{i=1}^{i=n} #lambda_{i} q_{i} p_{i}^{T})' )
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

    Hist_Residual_Rank = ROOT.TH1D("Hist_Residual_Rank","",3,0,3)
    for binx in range(0,len(bkgd_rank_integral)):
        Hist_Residual_Rank.SetBinContent(binx+1,abs(1.-bkgd_rank_integral[binx]/data_integral))
        ratio_err = bkgd_rank_integral[binx]/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Rank.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Rank,'n ranks','Residual in gamma region','Residual_Rank_%s'%(name),True)

def SingleSourceAnalysis(source_list,doMap):

    global ErecS_lower_cut
    global ErecS_upper_cut

    FilePath_List = []
    ResetStackedShowerHistograms()
    Hist_Data_ShowerElevNSB_Sum.Reset()
    Hist_Dark_ShowerElevNSB_Sum.Reset()
    Hist_Data_Unbiased_Energy_Sum.Reset()
    for source in range(0,len(source_list)):
        source_name = source_list[source]
        for elev in range(0,len(root_file_tags)):
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+"_%s"%(root_file_tags[elev])+".root";
            FilePath_List += [FilePath]
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
            InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
            InfoTree = InputFile.Get("InfoTree")
            InfoTree.GetEntry(0)

            HistName = "Hist_Data_ElevNSB"
            Hist_Data_ShowerElevNSB = InputFile.Get(HistName)
            Hist_Data_ShowerElevNSB_Sum.Add(Hist_Data_ShowerElevNSB)
            HistName = "Hist_Dark_ElevNSB"
            Hist_Dark_ShowerElevNSB = InputFile.Get(HistName)
            Hist_Dark_ShowerElevNSB_Sum.Add(Hist_Dark_ShowerElevNSB)
            HistName = "Hist_Data_Unbiased_Energy"
            Hist_Data_Unbiased_Energy = InputFile.Get(HistName)
            Hist_Data_Unbiased_Energy_Sum.Add(Hist_Data_Unbiased_Energy)

            MSCW_blind_cut = InfoTree.MSCW_cut_blind
            MSCL_blind_cut = InfoTree.MSCL_cut_blind
            bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
            bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
            bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
            bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                if ErecS_upper_cut<=energy_fine_bin[energy_fine_bin_cut_low]: continue
                if ErecS_lower_cut>=energy_fine_bin[energy_fine_bin_cut_up]: continue
                GetShowerHistogramsFromFile(FilePath_List[len(FilePath_List)-1])
                StackShowerHistograms()
                NormalizeEnergyHistograms(FilePath_List[len(FilePath_List)-1])
                StackEnergyHistograms()
                for e2 in range(energy_fine_bin_cut_low,energy_fine_bin_cut_up):
                    ErecS_lower_cut = energy_fine_bin[e2]
                    ErecS_upper_cut = energy_fine_bin[e2+1]
                    NormalizeTheta2Histograms(FilePath_List[len(FilePath_List)-1])
                    StackTheta2Histograms()
                    NormalizeCameraFoVHistograms(FilePath_List[len(FilePath_List)-1])
                    StackCameraFoVHistograms()
                    NormalizeSkyMapHistograms(FilePath_List[len(FilePath_List)-1])
                    StackSkymapHistograms()

    Make2DProjectionPlot(Hist_Data_ShowerElevNSB_Sum,'Ped. Var.','Zenith','Data_ShowerElevNSB_%s%s'%(source_name,PercentCrab),False)
    Make2DProjectionPlot(Hist_Dark_ShowerElevNSB_Sum,'Ped. Var.','Zenith','Dark_ShowerElevNSB_%s%s'%(source_name,PercentCrab),False)

    hists = [Hist_Data_Unbiased_Energy_Sum]
    legends = ['unfiltered cosmic-ray events']
    colors = [1]
    MakeComparisonPlot(hists,legends,colors,'E [GeV]','counts','Data_Unbiased_Energy_%s%s'%(source_name,PercentCrab),0,0,True,True)

    MatrixDecompositionDemo(source_name)

    GetSourceInfo(FilePath_List)
    GetBrightStarInfo(FilePath_List)

    CalculateSystError()
    PlotsStackedHistograms('%s%s'%(source_list[0],PercentCrab))

    MakeRankResidualPlots('%s%s'%(source_list[0],PercentCrab))

    if not doMap: return

    Hist_OnData_Skymap_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnData_Skymap_Galactic_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Sum.Rebin2D(n_rebin,n_rebin)
    Hist_OnBkgd_Skymap_Galactic_Sum.Rebin2D(n_rebin,n_rebin)

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,'RA','Dec','Skymap_RaDec_MDM_%s%s'%(source_name,PercentCrab))

    Hist_OnData_Skymap_smooth = Smooth2DMap(Hist_OnData_Skymap_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Sum,smooth_size,False)
    Hist_OnData_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnData_Skymap_Galactic_Sum,smooth_size,False)
    Hist_OnBkgd_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Sum,smooth_size,False)

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_smooth,Hist_OnBkgd_Skymap_smooth,'RA','Dec','Skymap_Smooth_RaDec_MDM_%s%s'%(source_name,PercentCrab))
    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Galactic_smooth,Hist_OnBkgd_Skymap_Galactic_smooth,'gal. l.','gal. b.','Skymap_Smooth_Galactic_MDM_%s%s'%(source_name,PercentCrab))

    Hist_Significance_Skymap_smooth = GetSignificanceMap(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth,Syst_MDM)

    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth, Hist_Significance_Skymap_smooth,5)
    print 'Excess (5 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (5 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (5 sigma) radius = %0.3f'%(excess_radius)
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth, Hist_Significance_Skymap_smooth,4)
    print 'Excess (4 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (4 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (4 sigma) radius = %0.3f'%(excess_radius)
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth, Hist_Significance_Skymap_smooth,3)
    print 'Excess (3 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (3 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (3 sigma) radius = %0.3f'%(excess_radius)

def FindSourceIndex(source_name):
    for source in range(0,len(sample_list)):
        if source_name==sample_list[source]:
            return source
    return 0

Hist_Data_ShowerElevNSB_Sum = ROOT.TH2D("Hist_Data_ShowerElevNSB_Sum","",20,0,10,18,0,90)
Hist_Dark_ShowerElevNSB_Sum = ROOT.TH2D("Hist_Dark_ShowerElevNSB_Sum","",20,0,10,18,0,90)
Hist_EffArea = ROOT.TH1D("Hist_EffArea","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_EffArea_Sum = ROOT.TH1D("Hist_EffArea_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_Data_Unbiased_Energy_Sum = ROOT.TH1D("Hist_Data_Unbiased_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
source_idx = FindSourceIndex(sample_list[0])
FilePath_Folder = []
for elev in range(0,len(root_file_tags)):
    SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+".root";
    FilePath_Folder += [SourceFilePath]
    if not os.path.isfile(FilePath_Folder[0]): 
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
Hist2D_AllRanks_Sum = ROOT.TH2D("Hist2D_AllRanks_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Sum = ROOT.TH2D("Hist2D_Rank0_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Sum = ROOT.TH2D("Hist2D_Rank1_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Sum = ROOT.TH2D("Hist2D_Rank2_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist_OnBkgd_R2off_Sum = ROOT.TH1D("Hist_OnBkgd_R2off_Sum","",50,0,10)
Hist_OnBkgd_R2off = ROOT.TH1D("Hist_OnBkgd_R2off","",50,0,10)
Hist_OnData_Theta2_Sum = ROOT.TH1D("Hist_OnData_Theta2_Sum","",50,0,10)
Hist_OnBkgd_Theta2_Sum = ROOT.TH1D("Hist_OnBkgd_Theta2_Sum","",50,0,10)
Hist_OnData_Theta2 = ROOT.TH1D("Hist_OnData_Theta2","",50,0,10)
Hist_OnBkgd_Theta2 = ROOT.TH1D("Hist_OnBkgd_Theta2","",50,0,10)
Hist_OnData_Energy_Sum = ROOT.TH1D("Hist_OnData_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy = ROOT.TH1D("Hist_OnData_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy = ROOT.TH1D("Hist_OnBkgd_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

Hist_OnData_Skymap = ROOT.TH2D("Hist_OnData_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_OnBkgd_Skymap = ROOT.TH2D("Hist_OnBkgd_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_OnData_Skymap_Sum = ROOT.TH2D("Hist_OnData_Skymap_Sum","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_OnBkgd_Skymap_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Sum","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_OnData_Skymap_Galactic = ROOT.TH2D("Hist_OnData_Skymap_Galactic","",150,source_l-3,source_l+3,150,source_b-3,source_b+3)
Hist_OnBkgd_Skymap_Galactic = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic","",150,source_l-3,source_l+3,150,source_b-3,source_b+3)
Hist_OnData_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnData_Skymap_Galactic_Sum","",150,source_l-3,source_l+3,150,source_b-3,source_b+3)
Hist_OnBkgd_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Sum","",150,source_l-3,source_l+3,150,source_b-3,source_b+3)

Hist_SystErr_MSCL = ROOT.TH1D("Hist_SystErr_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_SystErr_MSCW = ROOT.TH1D("Hist_SystErr_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

print 'MJD_Start = %s'%(MJD_Start)
print 'MJD_End = %s'%(MJD_End)
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
Hist2D_OffDark = []
Hist2D_OffDark_Sum = []
Hist2D_OffBkgd = []
Hist2D_OffBkgd_Sum = []
Hist_OffData_MSCL = []
Hist_OffData_MSCL_Sum = []
Hist_OffDark_MSCL = []
Hist_OffDark_MSCL_Sum = []
Hist_OffBkgd_MSCL = []
Hist_OffBkgd_MSCL_Sum = []
Hist_OffData_MSCW = []
Hist_OffData_MSCW_Sum = []
Hist_OffDark_MSCW = []
Hist_OffDark_MSCW_Sum = []
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
for nth_sample in range(0,n_control_samples-1):
    Hist2D_OffData += [ROOT.TH2D("Hist2D_OffData_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffData_Sum += [ROOT.TH2D("Hist2D_OffData_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffDark += [ROOT.TH2D("Hist2D_OffDark_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffDark_Sum += [ROOT.TH2D("Hist2D_OffDark_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffBkgd += [ROOT.TH2D("Hist2D_OffBkgd_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffBkgd_Sum += [ROOT.TH2D("Hist2D_OffBkgd_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCL += [ROOT.TH1D("Hist_OnData_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCL_Sum += [ROOT.TH1D("Hist_OnData_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCW += [ROOT.TH1D("Hist_OnData_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCW_Sum += [ROOT.TH1D("Hist_OnData_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffDark_MSCL += [ROOT.TH1D("Hist_OnDark_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffDark_MSCL_Sum += [ROOT.TH1D("Hist_OnDark_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffDark_MSCW += [ROOT.TH1D("Hist_OnDark_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffDark_MSCW_Sum += [ROOT.TH1D("Hist_OnDark_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
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
n_rebin = 2
smooth_size = 0.1
#n_rebin = 2
#smooth_size = 0.15
if energy_fine_bin[energy_fine_bin_cut_low]>500.:
    n_rebin = 2
    smooth_size = 0.15
if energy_fine_bin[energy_fine_bin_cut_low]>1000.:
    n_rebin = 4
    smooth_size = 0.2

GetGammaSourceInfo()

SingleSourceAnalysis(sample_list,True)
#SingleSourceAnalysis(sample_list,False)
print 'n_good_matches = %s'%(n_good_matches)
