
import os
import sys,ROOT
import array
import math
import csv
from array import *
from ROOT import *
from astropy import units as my_unit
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
from scipy import special
from scipy import interpolate
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import cycle
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import NullFormatter
from operator import itemgetter, attrgetter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Great examples of matplotlib plots: https://atmamani.github.io/cheatsheets/matplotlib/matplotlib_2/


#folder_path = 'output_test'
folder_path = 'output_default'

#folder_path = 'output_loose'
#folder_path = 'output_medium'
#folder_path = 'output_tight'

#folder_path = 'output_1hrs'
#folder_path = 'output_2hrs'
#folder_path = 'output_5hrs'
#folder_path = 'output_10hrs'
#folder_path = 'output_20hrs'

#folder_path = 'output_FreeElev'
#folder_path = 'output_FreeAzim'
#folder_path = 'output_FreeNSB'
#folder_path = 'output_FreeMJD'

#folder_path = 'output_elbow1p2'
#folder_path = 'output_elbow1p5'
#folder_path = 'output_elbow2p0'
#folder_path = 'output_elbow2p5'
#folder_path = 'output_elbow3p0'

#folder_path = 'output_6x6'
#folder_path = 'output_8x8'
#folder_path = 'output_12x12'
#folder_path = 'output_16x16'
#folder_path = 'output_20x20'

#folder_path = 'output_dAz_05'
#folder_path = 'output_dAz_10'
#folder_path = 'output_dAz_20'
#folder_path = 'output_dAz_40'
#folder_path = 'output_dAz_80'
#folder_path = 'output_dAz_180'

#folder_path = 'output_dNSB_0p5'
#folder_path = 'output_dNSB_1p0'
#folder_path = 'output_dNSB_1p5'
#folder_path = 'output_dNSB_2p0'
#folder_path = 'output_dNSB_2p5'

#folder_path = 'output_LogAlpha_m1p5'
#folder_path = 'output_LogAlpha_m1p0'
#folder_path = 'output_LogAlpha_m0p5'
#folder_path = 'output_LogAlpha_p0p0'
#folder_path = 'output_LogAlpha_p0p5'
#folder_path = 'output_LogAlpha_p1p0'
#folder_path = 'output_LogAlpha_p1p5'

#N_bins_for_deconv = 6
N_bins_for_deconv = 12
if '4x4' in folder_path:
    N_bins_for_deconv = 4
if '6x6' in folder_path:
    N_bins_for_deconv = 6
if '8x8' in folder_path:
    N_bins_for_deconv = 8
if '12x12' in folder_path:
    N_bins_for_deconv = 12
if '16x16' in folder_path:
    N_bins_for_deconv = 16
if '20x20' in folder_path:
    N_bins_for_deconv = 20
gamma_hadron_dim_ratio_w = 1.
gamma_hadron_dim_ratio_l = 1.
gamma_hadron_low_end = 0.

MSCW_blind_cut = 0.5
MSCL_blind_cut = 0.7
if 'loose' in folder_path:
    MSCW_blind_cut = 0.6
if 'tight' in folder_path:
    MSCW_blind_cut = 0.4

skymap_zoomin_scale = 1
#skymap_zoomin_scale = 1.5
#skymap_zoomin_scale = 2
smooth_size_spectroscopy = 0.1
#smooth_size_spectroscopy = 0.14
#smooth_size_spectroscopy = 0.2
#smooth_size_spectroscopy = 0.3

Smoothing = True
#Smoothing = False

additional_tag = ''

#UseEffectiveArea = True
#additional_tag = ''
UseEffectiveArea = False
additional_tag = '_DD'

#calibration_radius = 0.2 # need to be larger than the PSF and smaller than the integration radius
calibration_radius = 0.3 # need to be larger than the PSF and smaller than the integration radius
#calibration_radius = 1.0

#energy_index_scale = 0
energy_index_scale = 2

#doGalacticCoord = True
doGalacticCoord = False

Skymap_nzones_x = 1
Skymap_nzones_y = 1
Skymap_size_x = 2.0
Skymap_size_y = 2.0
Skymap_nbins_x = 50
Skymap_nbins_y = 50
if doGalacticCoord:
    Skymap_nzones_x = 3
    Skymap_nzones_y = 3
    Skymap_size_x = 6.
    Skymap_nbins_x = 60
    Skymap_size_y = 6.
    Skymap_nbins_y = 60
#Skymap_nbins = 15 
#Skymap_nbins = 9 # Crab calibration 
#Skymap_nbins = 5
#Skymap_nbins = 3

#Skymap_normalization_nbins = 1

#target_max_dist_cut = 3.
target_max_dist_cut = 100.

elev_range = [40,90]
#elev_range = [30,90]
#elev_range = [35,45]

#energy_bin = [100.,316.,1000.,3162.,10000.]
energy_bin = [200.,398.,794.,1585.,3162.,6310.,12589.]
energy_fine_bin = energy_bin

def Hist2DIntegralAndError(Hist):

    integral = 0.
    error = 0.
    for bx in range(0,Hist.GetNbinsX()):
        for by in range(0,Hist.GetNbinsY()):
            integral += Hist.GetBinContent(bx+1,by+1)
            error += Hist.GetBinError(bx+1,by+1)
    return integral, error

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
            hT.SetBinError( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinError( binx, biny ) )
    return hT

def Smooth2DMap(Hist_Old,smooth_size,addLinearly,normalized):

    nbins = Hist_Old.GetNbinsX()
    MapEdge_left = Hist_Old.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Old.GetXaxis().GetBinLowEdge(Hist_Old.GetNbinsX()+1)
    map_size = (MapEdge_right-MapEdge_left)/2.
    Hist_Kernel = ROOT.TH2D("Hist_Kernel","",nbins,-map_size,map_size,nbins,-map_size,map_size)
    Hist_Kernel.Reset()
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            cell_x = Hist_Kernel.GetXaxis().GetBinCenter(bx1)
            cell_y = Hist_Kernel.GetYaxis().GetBinCenter(by1)
            distance = pow(cell_x*cell_x+cell_y*cell_y,0.5)
            bin_content = ROOT.TMath.Gaus(distance,0,smooth_size)
            Hist_Kernel.SetBinContent(bx1,by1,bin_content)
    #print ('Hist_Kernel.Integral() = %s'%(Hist_Kernel.Integral()))

    Hist_Smooth = Hist_Old.Clone()

    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(2*smooth_size/bin_size) + 1
    central_bin = int(nbins/2) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            old_content = Hist_Old.GetBinContent(bx1,by1)
            old_error = Hist_Old.GetBinError(bx1,by1)
            bin_content = 0
            bin_error = 0
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2<1: 
                        continue
                    if by2<1: 
                        continue
                    if bx2>Hist_Old.GetNbinsX(): 
                        continue
                    if by2>Hist_Old.GetNbinsY(): 
                        continue
                    bin_content += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*Hist_Old.GetBinContent(bx2,by2)
                    bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*pow(Hist_Old.GetBinError(bx2,by2),2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))

    if normalized:
        Hist_Smooth.Scale(1./Hist_Kernel.Integral())
        #for bx1 in range(1,Hist_Smooth.GetNbinsX()+1):
        #    for by1 in range(1,Hist_Smooth.GetNbinsY()+1):
        #        new_error = pow(Hist_Smooth.GetBinContent(bx1,by1),0.5)
        #        Hist_Smooth.SetBinError(bx1,by1,max(new_error,1.))

    return Hist_Smooth

def FindProjection(Hist_Data_input,Hist_Syst_input,roi_x1,roi_y1,roi_x2,roi_y2,roi_d,roi_width,isTransverse,fraction):

    integration_range = pow(pow(roi_x1-roi_x2,2)+pow(roi_y1-roi_y2,2),0.5)
    roi_x0 = roi_x1
    roi_y0 = roi_y1
    extension_low = 1.5
    extension_up = 1.5
    if isTransverse:
        integration_range = 0.
        extension_low = 2.
        extension_up = 2.
        roi_x0 = roi_x1+(roi_x2-roi_x1)*fraction
        roi_y0 = roi_y1+(roi_y2-roi_y1)*fraction
    n_bins = int((integration_range+extension_low+extension_up)/0.1)
    Hist_Profile_ProjectedX = ROOT.TH1D("Hist_Profile_ProjectedX","",10,-1.5,1.5)

    # ax + by + c = 0
    b = 1.
    a = -b*(roi_y1-roi_y2)/(roi_x1-roi_x2)
    c = -(b*roi_y0 + a*roi_x0)
    if isTransverse:
        a = -1./a
        c = -(b*roi_y0 + a*roi_x0)

    for br in range(0,Hist_Profile_ProjectedX.GetNbinsX()):
        range_limit = Hist_Profile_ProjectedX.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_ProjectedX.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_syst_err = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_cell_to_line = abs(a*cell_x+b*cell_y+c)/(a*a+b*b)
                if distance_cell_to_line>roi_width: continue
                closest_x_on_line = (b*(b*cell_x-a*cell_y)-a*c)/(a*a+b*b)
                closest_y_on_line = (a*(-b*cell_x+a*cell_y)-b*c)/(a*a+b*b)
                projected_x = pow(pow(closest_x_on_line-roi_x0,2)+pow(closest_y_on_line-roi_y0,2),0.5)
                if not isTransverse:
                    if (closest_y_on_line-roi_y0)<0.:
                        projected_x = -projected_x
                else:
                    if (closest_x_on_line-roi_x0)>0.:
                        projected_x = -projected_x
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                syst_error = 0.
                if Hist_Syst_input!=None:
                    syst_error = Hist_Syst_input.GetBinContent(bx+1,by+1)
                if projected_x>=range_limit_previous and projected_x<range_limit:
                    slice_data += data_content
                    slice_data_err += data_error*data_error
                    slice_syst_err += syst_error
        slice_data_err = pow(slice_data_err,0.5)
        #if slice_data==0.: slice_data_err = 1.
        Hist_Profile_ProjectedX.SetBinContent(br+1,slice_data)
        #Hist_Profile_ProjectedX.SetBinError(br+1,slice_data_err)
        Hist_Profile_ProjectedX.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))

    profile = []
    profile_err = []
    theta2 = []
    for binx in range(0,Hist_Profile_ProjectedX.GetNbinsX()):
        center = Hist_Profile_ProjectedX.GetBinCenter(binx+1)
        range_limit = Hist_Profile_ProjectedX.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_ProjectedX.GetBinLowEdge(binx+1)
        solid_angle = 2.*roi_width*(range_limit-range_limit_previous)
        profile_content = Hist_Profile_ProjectedX.GetBinContent(binx+1)/solid_angle
        profile_error = Hist_Profile_ProjectedX.GetBinError(binx+1)/solid_angle
        theta2 += [center]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2

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

def FindCountProjection(Hist_Data_input,proj_type="Y"):

    n_bins_y = Hist_Data_input.GetNbinsY()
    n_bins_x = Hist_Data_input.GetNbinsX()
    #n_bins_y = 15
    #n_bins_x = 15
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)

    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    if proj_type=="X":
        MapCenter_x = (MapEdge_upper+MapEdge_lower)/2.

    hist_nbins = n_bins_y
    hist_edge_lower = MapEdge_lower
    hist_edge_upper = MapEdge_upper
    if proj_type=="X":
        hist_nbins = n_bins_x
        hist_edge_lower = MapEdge_left
        hist_edge_upper = MapEdge_right
    Hist_Profile_Y = ROOT.TH1D("Hist_Profile_Y","",hist_nbins,hist_edge_lower,hist_edge_upper)
    for br in range(0,Hist_Profile_Y.GetNbinsX()):
        range_limit = Hist_Profile_Y.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                if proj_type=="X":
                    cell_y = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                    cell_x = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                if abs(MapCenter_x-cell_x)>2.0: continue
                if cell_y>=range_limit_previous and cell_y<range_limit:
                    if not data_error==0.:
                        total_cell_size += 1.
                        slice_data += data_content
                        slice_data_err += data_error*data_error
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data
            slice_data_err = pow(slice_data_err,0.5)
        Hist_Profile_Y.SetBinContent(br+1,slice_data)
        Hist_Profile_Y.SetBinError(br+1,slice_data_err)

    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Profile_Y.GetNbinsX()):
        center = Hist_Profile_Y.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Y.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Y.GetBinContent(binx+1)
        profile_error = Hist_Profile_Y.GetBinError(binx+1)
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2, theta2_err

def FindGalacticProjection_v2(Hist_Data_input,proj_type="Y"):

    #n_bins_y = Hist_Data_input.GetNbinsY()
    #n_bins_x = Hist_Data_input.GetNbinsX()
    n_bins_y = 15
    n_bins_x = 15
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)

    cell_size = (MapEdge_right-MapEdge_left)/float(n_bins_x)*(MapEdge_upper-MapEdge_lower)/float(n_bins_y)
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    if proj_type=="X":
        MapCenter_x = (MapEdge_upper+MapEdge_lower)/2.

    hist_nbins = n_bins_y
    hist_edge_lower = MapEdge_lower
    hist_edge_upper = MapEdge_upper
    if proj_type=="X":
        hist_nbins = n_bins_x
        hist_edge_lower = MapEdge_left
        hist_edge_upper = MapEdge_right
    Hist_Profile_Y = ROOT.TH1D("Hist_Profile_Y","",hist_nbins,hist_edge_lower,hist_edge_upper)
    for br in range(0,Hist_Profile_Y.GetNbinsX()):
        range_limit = Hist_Profile_Y.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                if proj_type=="X":
                    cell_y = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                    cell_x = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                if abs(MapCenter_x-cell_x)>0.7: continue
                if cell_y>=range_limit_previous and cell_y<range_limit:
                    if not data_error==0.:
                        total_cell_size += cell_size
                        slice_data += data_content
                        slice_data_err += data_error*data_error
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data/total_cell_size
            slice_data_err = pow(slice_data_err,0.5)/total_cell_size
        Hist_Profile_Y.SetBinContent(br+1,slice_data)
        Hist_Profile_Y.SetBinError(br+1,slice_data_err)

    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Profile_Y.GetNbinsX()):
        center = Hist_Profile_Y.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Y.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(binx+1)
        sr_to_deg2 = 3282.80635
        #sr_to_deg2 = 1.
        profile_content = Hist_Profile_Y.GetBinContent(binx+1)*sr_to_deg2
        profile_error = Hist_Profile_Y.GetBinError(binx+1)*sr_to_deg2
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2, theta2_err

def FindExtension_v1(Hist_Data_input,Hist_Bkgd_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins_2d = Hist_Data_input.GetNbinsX()
    integration_range = min(integration_range,1.8)
    n_bins_1d = int(10.*integration_range)
    #n_bins_1d = int(5.*integration_range)

    n_bins_y = Hist_Data_input.GetNbinsY()
    n_bins_x = Hist_Data_input.GetNbinsX()
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)
    cell_size = (MapEdge_right-MapEdge_left)/float(n_bins_x)*(MapEdge_upper-MapEdge_lower)/float(n_bins_y)

    Hist_Bkgd_Theta2 = ROOT.TH1D("Hist_Bkgd_Theta2","",n_bins_1d,0,integration_range)
    Hist_Data_Theta2 = ROOT.TH1D("Hist_Data_Theta2","",n_bins_1d,0,integration_range)
    for br in range(0,Hist_Data_Theta2.GetNbinsX()):
        range_limit = Hist_Data_Theta2.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Data_Theta2.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_bkgd = 0.
        slice_bkgd_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                #data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                data_error = pow(Hist_Data_input.GetBinContent(bx+1,by+1),0.5)
                bkgd_content = Hist_Bkgd_input.GetBinContent(bx+1,by+1)
                #bkgd_error = Hist_Bkgd_input.GetBinError(bx+1,by+1)
                bkgd_error = pow(Hist_Bkgd_input.GetBinContent(bx+1,by+1),0.5)
                if distance_sq>=pow(range_limit_previous,2) and distance_sq<pow(range_limit,2):
                    if not data_error==0.:
                        slice_data += data_content
                        slice_data_err += data_error*data_error
                        slice_bkgd += bkgd_content
                        slice_bkgd_err += bkgd_error*bkgd_error
                        total_cell_size += cell_size
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
            slice_bkgd = 0.
            slice_bkgd_err = 0.
        else:
            slice_data = slice_data/total_cell_size
            slice_data_err = pow(slice_data_err,0.5)/total_cell_size
            slice_bkgd = slice_bkgd/total_cell_size
            slice_bkgd_err = pow(slice_bkgd_err,0.5)/total_cell_size
        Hist_Data_Theta2.SetBinContent(br+1,slice_data)
        Hist_Data_Theta2.SetBinError(br+1,slice_data_err)
        Hist_Bkgd_Theta2.SetBinContent(br+1,slice_bkgd)
        Hist_Bkgd_Theta2.SetBinError(br+1,slice_bkgd_err)

    data_cnt = []
    data_cnt_err = []
    bkgd_cnt = []
    bkgd_cnt_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Data_Theta2.GetNbinsX()):
        center = Hist_Data_Theta2.GetBinCenter(binx+1)
        range_limit = Hist_Data_Theta2.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Data_Theta2.GetBinLowEdge(binx+1)
        data_content = Hist_Data_Theta2.GetBinContent(binx+1)
        data_error = Hist_Data_Theta2.GetBinError(binx+1)
        bkgd_content = Hist_Bkgd_Theta2.GetBinContent(binx+1)
        bkgd_error = Hist_Bkgd_Theta2.GetBinError(binx+1)
        if center>integration_range: continue
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        data_cnt += [data_content]
        data_cnt_err += [data_error]
        bkgd_cnt += [bkgd_content]
        bkgd_cnt_err += [bkgd_error]

    return data_cnt, data_cnt_err, bkgd_cnt, bkgd_cnt_err, theta2, theta2_err


def FindExtension_v2(Hist_Data_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins_2d = Hist_Data_input.GetNbinsX()
    integration_range = min(integration_range,2.0)
    #n_bins_1d = 5
    n_bins_1d = 10

    n_bins_y = Hist_Data_input.GetNbinsY()
    n_bins_x = Hist_Data_input.GetNbinsX()
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)
    cell_size = (MapEdge_right-MapEdge_left)/float(n_bins_x)*(MapEdge_upper-MapEdge_lower)/float(n_bins_y)

    Hist_Profile_Theta2 = ROOT.TH1D("Hist_Profile_Theta2","",n_bins_1d,0,integration_range)
    for br in range(0,Hist_Profile_Theta2.GetNbinsX()):
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        total_error_weight = 0.
        total_cell_size = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2)
                #data_content = Hist_Data_input.GetBinContent(bx+1,by+1)/cell_size
                #data_error = Hist_Data_input.GetBinError(bx+1,by+1)/cell_size
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                if distance_sq>=pow(range_limit_previous,2) and distance_sq<pow(range_limit,2):
                    if not data_error==0.:
                        #error_weight = 1./(data_error*data_error)
                        #slice_data += data_content*error_weight
                        #total_error_weight += error_weight
                        #slice_data_err += 1.
                        slice_data += data_content
                        slice_data_err += data_error*data_error
                        total_cell_size += cell_size
        if total_cell_size==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data/total_cell_size
            slice_data_err = pow(slice_data_err,0.5)/total_cell_size
        Hist_Profile_Theta2.SetBinContent(br+1,slice_data)
        Hist_Profile_Theta2.SetBinError(br+1,slice_data_err)

    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
    for binx in range(0,Hist_Profile_Theta2.GetNbinsX()):
        center = Hist_Profile_Theta2.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Theta2.GetBinContent(binx+1)
        profile_error = Hist_Profile_Theta2.GetBinError(binx+1)
        if center>integration_range: continue
        theta2 += [center]
        theta2_err += [range_limit-range_limit_previous]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2, theta2_err


def SaveAsFITS(hist_map,fits_name):

    map_nbins = hist_map.GetNbinsX()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)
    grid_z = np.zeros((map_nbins, map_nbins))
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins)
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins)
    max_z = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            grid_z[ybin,xbin] = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if max_z<hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                max_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)

    print ("Constructing WCS from a dictionary header:")
    bin_ref_x = 2
    bin_ref_y = 2
    bin_ref_center_x = -1.*hist_map.GetXaxis().GetBinCenter(bin_ref_x)
    bin_ref_center_y = hist_map.GetYaxis().GetBinCenter(bin_ref_y)
    bin_delta_x = -1.*(hist_map.GetXaxis().GetBinCenter(2)-hist_map.GetXaxis().GetBinCenter(1))
    bin_delta_y = hist_map.GetYaxis().GetBinCenter(2)-hist_map.GetYaxis().GetBinCenter(1)
    header = {'NAXIS': 2,
              'NAXIS1': map_nbins,
              'CTYPE1': 'RA---TAN',
              'CRVAL1': bin_ref_center_x,
              'CRPIX1': bin_ref_x,
              'CDELT1': bin_delta_x,
              'NAXIS2': map_nbins,
              'CTYPE2': 'DEC--TAN',
              'CRVAL2': bin_ref_center_y,
              'CRPIX2': bin_ref_y,
              'CDELT2': bin_delta_y}
    wcs = WCS(header)

    outfile = 'output_plots/%s.fits'%(fits_name)
    hdu = fits.PrimaryHDU(grid_z)
    hdu.header.update(wcs.to_header())
    hdu.writeto(outfile, overwrite=True)

def FindExtension(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins = 8
    integration_range = 1.6
    if sys.argv[1]=='Geminga_ON':
        n_bins = 5
        integration_range = 1.6
    if sys.argv[1]=='LHAASO_J2108_ON':
        n_bins = 5
        integration_range = 1.6

    Hist_Profile_Theta2 = ROOT.TH1D("Hist_Profile_Theta2","",n_bins,0,integration_range)
    for br in range(0,Hist_Profile_Theta2.GetNbinsX()):
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_syst_err = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                distance_sq = pow(cell_x-roi_x,2)+pow(cell_y-roi_y,2)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                syst_error = 0.
                if Hist_Syst_input!=None:
                    syst_error = Hist_Syst_input.GetBinContent(bx+1,by+1)
                if distance_sq>=pow(range_limit_previous,2) and distance_sq<pow(range_limit,2):
                    slice_data += data_content
                    slice_data_err += data_error*data_error
                    slice_syst_err += syst_error
        slice_data_err = pow(slice_data_err,0.5)
        if slice_data==0.: slice_data_err = 1.
        Hist_Profile_Theta2.SetBinContent(br+1,slice_data)
        Hist_Profile_Theta2.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))


    profile = []
    profile_err = []
    theta2 = []
    for binx in range(0,Hist_Profile_Theta2.GetNbinsX()):
        center = Hist_Profile_Theta2.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Theta2.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Theta2.GetBinLowEdge(binx+1)
        #solid_angle = 3.14*(range_limit*range_limit-range_limit_previous*range_limit_previous)
        solid_angle = 2.*3.14*center*(range_limit-range_limit_previous)
        profile_content = Hist_Profile_Theta2.GetBinContent(binx+1)/solid_angle
        #profile_content = max(0.,profile_content)
        profile_error = Hist_Profile_Theta2.GetBinError(binx+1)/solid_angle
        if center>integration_range: continue
        theta2 += [center]
        profile += [profile_content]
        profile_err += [profile_error]

    return profile, profile_err, theta2

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def FillSkymapHoles(hist_map, map_resolution):

    hist_map_new = hist_map.Clone()
    bin_size = hist_map.GetXaxis().GetBinCenter(2)-hist_map.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(map_resolution/bin_size) + 1
    for bx1 in range(1,hist_map.GetNbinsX()+1):
        for by1 in range(1,hist_map.GetNbinsY()+1):
            bin_content_1 = hist_map.GetBinContent(bx1,by1)
            if bin_content_1!=0.: continue
            min_distance = 1e10
            min_distance_content = 0.
            locationx1 = hist_map.GetXaxis().GetBinCenter(bx1)
            locationy1 = hist_map.GetYaxis().GetBinCenter(by1)
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2>=1 and bx2<=hist_map.GetNbinsX():
                        if by2>=1 and by2<=hist_map.GetNbinsY():
                            bin_content_2 = hist_map.GetBinContent(bx2,by2)
                            if bin_content_2==0.: continue
                            locationx2 = hist_map.GetXaxis().GetBinCenter(bx2)
                            locationy2 = hist_map.GetYaxis().GetBinCenter(by2)
                            distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5)
                            if min_distance>distance:
                                min_distance = distance
                                min_distance_content = bin_content_2
            hist_map_new.SetBinContent(bx1,by1,min_distance_content)
    return hist_map_new

def GetFITSMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    nbins_x = hist_map.GetNbinsX()
    nbins_y = hist_map.GetNbinsY()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.

    hdu = fits.open(map_file)[0]
    wcs = WCS(hdu.header)
    image_data = hdu.data

    print ("wcs")
    print (wcs)

    image_shape = np.shape(image_data)
    for binx in range(0,nbins_x):
        for biny in range(0,nbins_y):
            ra = hist_map.GetXaxis().GetBinCenter(binx+1)
            dec = hist_map.GetYaxis().GetBinCenter(biny+1)
            pixs = wcs.all_world2pix(ra,dec,1)
            if int(pixs[0])<0: continue
            if int(pixs[1])<0: continue
            if int(pixs[0])>=image_shape[0]: continue
            if int(pixs[1])>=image_shape[1]: continue
            fits_data = image_data[int(pixs[1]),int(pixs[0])]
            hist_map.SetBinContent(binx+1,biny+1,fits_data)

    return hist_map

def GetGalacticCoordMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open(map_file)
    for line in inputFile:
        sig = float(line.split(' ')[2])
        l = float(line.split(' ')[0])
        b = float(line.split(' ')[1])
        if isRaDec: 
            l, b = ConvertGalacticToRaDec(l,b)
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,old_sig))

    map_resolution = 0.1
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split(':')]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split(':')]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)           
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def ReadBrightStarListFromFile():
    star_name = []
    star_ra = []
    star_dec = []
    inputFile = open('Hipparcos_MAG8_1997.dat')
    for line in inputFile:
        if line[0]=="#": continue
        if line[0]=="*": continue
        if len(line.split())<5: continue
        target_ra = float(line.split()[0])
        target_dec = float(line.split()[1])
        target_brightness = float(line.split()[3])+float(line.split()[4])
        #if target_brightness>7.: continue
        if target_brightness>6.: continue
        #if target_brightness>5.: continue
        star_ra += [target_ra]
        star_dec += [target_dec]
        star_name += ['bmag = %0.1f'%(target_brightness)]
    #print ('Found %s bright stars.'%(len(star_name)))
    return star_name, star_ra, star_dec

def ReadFermiCatelog():
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open('gll_psc_v26.xml')
    target_name = ''
    target_type = ''
    target_info = ''
    target_flux = ''
    target_ra = ''
    target_dec = ''
    for line in inputFile:
        if line.split(' ')[0]=='<source':
            for block in range(0,len(line.split(' '))):
                if 'Unc_' in line.split(' ')[block]: continue
                if 'name=' in line.split(' ')[block]:
                    target_name = line.split('name="')[1].split('"')[0]
                if 'type=' in line.split(' ')[block]:
                    target_type = line.split('type="')[1].split('"')[0]
                if 'Flux1000=' in line.split(' ')[block]:
                    target_flux = line.split('Flux1000="')[1].split('"')[0]
                if 'Energy_Flux100=' in line.split(' ')[block]:
                    target_info = line.split(' ')[block]
                    target_info = target_info.strip('>\n')
                    target_info = target_info.strip('"')
                    target_info = target_info.lstrip('Energy_Flux100="')
        if '<parameter' in line and 'name="RA"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_ra = line.split(' ')[block].split('"')[1]
        if '<parameter' in line and 'name="DEC"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_dec = line.split(' ')[block].split('"')[1]
        if 'source>' in line:
            if target_ra=='': 
                target_name = ''
                target_type = ''
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #if target_type=='PointSource': 
            #if float(target_info)<1e-5: 
            if float(target_flux)<1e-9: 
                target_name = ''
                target_type = ''
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #print ('target_name = %s'%(target_name))
            #print ('target_type = %s'%(target_type))
            #print ('target_ra = %s'%(target_ra))
            #print ('target_dec = %s'%(target_dec))
            source_name += [target_name]
            #source_name += ['%0.2e'%(float(target_info))]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_type = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def ReadHAWCTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        if '- name:' in line:
            target_name = line.lstrip('   - name: ')
            target_name = target_name.strip('\n')
        if 'RA:' in line:
            target_ra = line.lstrip('     RA: ')
        if 'Dec:' in line:
            target_dec = line.lstrip('     Dec: ')
        if 'flux measurements:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def ReadATNFTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        if target_name=="\n": continue
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        target_dist = line.split(',')[3].strip(" ")
        target_age = line.split(',')[4].strip(" ")
        target_edot = line.split(',')[5].strip(" ")
        if target_dist=='*': continue
        if target_age=='*': continue
        if target_edot=='*': continue
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        #if float(target_age)/1000.>pow(10,4.5): continue
        if float(target_edot)<1e35: continue
        if float(target_dist)>6.: continue
        #if target_dist>target_max_dist_cut: continue
        #if target_dist>2.0: continue

        #ra_deg = float(HMS2deg(target_ra,target_dec)[0])
        #dec_deg = float(HMS2deg(target_ra,target_dec)[1])
        #gal_l, gal_b = ConvertRaDecToGalactic(ra_deg,dec_deg)
        #if abs(gal_b)<5.: continue

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [float(target_dist)]
        source_age += [float(target_age)]
        source_edot += [float(target_edot)]
    return source_name, source_ra, source_dec, source_dist, source_age

def ReadSNRTargetListFromCSVFile():
    source_name = []
    source_ra = []
    source_dec = []
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '1000'
            if float(target_min_dist)>6.: continue
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            #print(row)
    return source_name, source_ra, source_dec

def GetGammaSourceInfo(hist_contour,prime_psr_name=None,prime_psr_ra=None,prime_psr_dec=None):

    other_stars = []
    other_stars_type = []
    other_star_coord = []

    near_source_cut = 0.1

    drawBrightStar = True
    drawPulsar = True
    drawSNR = True
    drawFermi = False
    drawHAWC = False
    drawTeV = False

    #if not prime_psr_name==None:
    #    for src in range(0,len(prime_psr_name)):
    #        src_ra = prime_psr_ra[src]
    #        src_dec = prime_psr_dec[src]
    #        if doGalacticCoord:
    #            src_ra, src_dec = ConvertRaDecToGalactic(src_ra,src_dec)
    #        near_a_source = False
    #        for entry in range(0,len(other_stars)):
    #            distance = pow(src_ra-other_star_coord[entry][0],2)+pow(src_dec-other_star_coord[entry][1],2)
    #            if distance<near_source_cut*near_source_cut:
    #                near_a_source = True
    #        if not near_a_source:
    #            other_stars += [prime_psr_name[src]]
    #            other_star_coord += [[src_ra,src_dec]]

    if drawBrightStar:
        star_name, star_ra, star_dec = ReadBrightStarListFromFile()
        for src in range(0,len(star_name)):
            src_ra = star_ra[src]
            src_dec = star_dec[src]
            if doGalacticCoord:
                src_ra, src_dec = ConvertRaDecToGalactic(src_ra,src_dec)
            other_stars += [star_name[src]]
            other_stars_type += ['Star']
            other_star_coord += [[src_ra,src_dec]]

    if drawSNR:
        target_snr_name, target_snr_ra, target_snr_dec = ReadSNRTargetListFromCSVFile()
        for src in range(0,len(target_snr_name)):
            gamma_source_name = target_snr_name[src]
            gamma_source_ra = target_snr_ra[src]
            gamma_source_dec = target_snr_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['SNR']
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    if drawPulsar:
        target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
        for src in range(0,len(target_psr_name)):
            gamma_source_name = target_psr_name[src]
            gamma_source_ra = target_psr_ra[src]
            gamma_source_dec = target_psr_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                if target_psr_age[src]/1000.>1e6:
                        other_stars_type += ['MSP']
                else:
                        other_stars_type += ['PSR']
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    if drawFermi:
        fermi_name, fermi_ra, fermi_dec = ReadFermiCatelog()
        for src in range(0,len(fermi_name)):
            gamma_source_name = fermi_name[src]
            gamma_source_ra = fermi_ra[src]
            gamma_source_dec = fermi_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['Fermi']
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    if drawHAWC:
        target_hwc_name, target_hwc_ra, target_hwc_dec = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
        for src in range(0,len(target_hwc_name)):
            gamma_source_name = target_hwc_name[src]
            gamma_source_ra = target_hwc_ra[src]
            gamma_source_dec = target_hwc_dec[src]
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source:
                other_stars += [gamma_source_name]
                other_stars_type += ['HAWC']
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    if drawTeV:
        inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            if doGalacticCoord:
                gamma_source_ra, gamma_source_dec = ConvertRaDecToGalactic(gamma_source_ra,gamma_source_dec)
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<near_source_cut*near_source_cut:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_stars_type += ['TeV']
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    return other_stars, other_stars_type, other_star_coord

def GetRegionIntegral(hist_data_skymap,hist_syst_skymap,hist_mask_skymap,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r):

    flux_sum = 0.
    flux_stat_err = 0.
    flux_syst_err = 0.
    for bx in range(0,hist_data_skymap.GetNbinsX()):
        for by in range(0,hist_data_skymap.GetNbinsY()):
            bin_ra = hist_data_skymap.GetXaxis().GetBinCenter(bx+1)
            bin_dec = hist_data_skymap.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_ra-roi_x,2) + pow(bin_dec-roi_y,2),0.5)
            excl_distance = pow(pow(bin_ra-excl_roi_x,2) + pow(bin_dec-excl_roi_y,2),0.5)
            if distance>roi_r: 
                continue
            if excl_distance<excl_roi_r: 
                continue
            if not hist_mask_skymap==None:
                mask = hist_mask_skymap.GetBinContent(bx+1,by+1)
                if mask!=0.:
                    continue
            flux_sum += hist_data_skymap.GetBinContent(bx+1,by+1)
            if hist_syst_skymap!=None:
                flux_syst_err += hist_syst_skymap.GetBinContent(bx+1,by+1)
            flux_stat_err += pow(hist_data_skymap.GetBinError(bx+1,by+1),2)
    flux_stat_err = pow(flux_stat_err,0.5)
    return flux_sum, flux_stat_err, flux_syst_err

def GetRegionSpectrum_v2(hist_data_skymap,hist_syst_skymap,hist_mask_skymap,ebin_low,ebin_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r):

    x_axis = []
    x_error = []
    y_axis = []
    y_error = []
    for ebin in range(ebin_low,ebin_up):
        flux_sum = 0.
        flux_stat_err = 0.
        flux_syst_err = 0.
        flux_sum, flux_stat_err, flux_syst_err = GetRegionIntegral(hist_data_skymap[ebin],hist_syst_skymap[ebin],hist_mask_skymap,roi_x,roi_y,roi_r[ebin],excl_roi_x,excl_roi_y,excl_roi_r[ebin])
        x_axis += [0.5*(energy_bin[ebin]+energy_bin[ebin+1])]
        x_error += [0.5*(energy_bin[ebin+1]-energy_bin[ebin])]
        y_axis += [flux_sum]
        #y_error += [pow(flux_stat_err*flux_stat_err+flux_syst_err*flux_syst_err,0.5)]
        y_error += [flux_stat_err]

    return x_axis, x_error, y_axis, y_error

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn):

    Hist_Skymap = Hist_SR.Clone()
    MapEdge_left = Hist_Skymap.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Skymap.GetXaxis().GetBinLowEdge(Hist_Skymap.GetNbinsX()+1)
    MapEdge_lower = Hist_Skymap.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Skymap.GetYaxis().GetBinLowEdge(Hist_Skymap.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.
    Hist_Skymap_zoomin = ROOT.TH2D("Hist_Skymap_zoomin","",int(Skymap_nbins_x/3),MapCenter_x-MapSize_x/2,MapCenter_x+MapSize_x/2,int(Skymap_nbins_y/3),MapCenter_y-MapSize_y/2,MapCenter_y+MapSize_y/2)
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            Shape_Err = Hist_Syst.GetBinContent(bx+1,by+1)
            Data_Stat_Err = Hist_SR.GetBinError(bx+1,by+1)
            Bkgd_Stat_Err = Hist_Bkg.GetBinError(bx+1,by+1)
            NBkg_Err = pow(pow(Bkgd_Stat_Err,2)+pow(Shape_Err,2),0.5)
            #Sig = CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Sig = 0.
            if pow(pow(NBkg_Err,2)+pow(Data_Stat_Err,2),0.5)>0.:
                Sig = (NSR-NBkg)/pow(pow(NBkg_Err,2)+pow(Data_Stat_Err,2),0.5)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
            bx_center = Hist_Skymap.GetXaxis().GetBinCenter(bx+1)
            by_center = Hist_Skymap.GetYaxis().GetBinCenter(by+1)
            bx2 = Hist_Skymap_zoomin.GetXaxis().FindBin(bx_center)
            by2 = Hist_Skymap_zoomin.GetYaxis().FindBin(by_center)
            Hist_Skymap_zoomin.SetBinContent(bx2,by2,Sig)
    if isZoomIn:
        return Hist_Skymap_zoomin
    else:
        return Hist_Skymap

def GetHawcSkymap(hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open('MWL_maps/hawc_map.txt')
    for line in inputFile:
        sig = float(line.split(' ')[0])
        l = float(line.split(' ')[1])
        b = float(line.split(' ')[2])
        if not isRaDec: 
            l, b = ConvertRaDecToGalactic(l,b)
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,0.))

    map_resolution = 0.1
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def MatplotlibHist2D(hist_map,fig,label_x,label_y,label_z,plotname,zmax=0,zmin=0):

    top = cm.get_cmap('Blues_r', 128) # r means reversed version
    bottom = cm.get_cmap('Oranges', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')
    colormap = 'coolwarm'

    map_nbins = hist_map.GetNbinsX()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins
    nbins_per_deg = map_nbins/(MapEdge_right-MapEdge_left)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins)
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins)
    grid_z = np.zeros((map_nbins, map_nbins))
    max_z = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            grid_z[ybin,xbin] = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if max_z<hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                max_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)

    fig.clf()
    figsize_x = 6.
    figsize_y = 6.
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    if zmax==0 and zmin==0:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0,vmin=zmin,vmax=zmax)

    if 'Matrix' in plotname:
        x_cut = 0.7
        y_cut = 0.5
        x1, y1 = [-x_cut, -x_cut], [-y_cut, y_cut]
        x2, y2 = [-x_cut, x_cut],  [y_cut, y_cut]
        x3, y3 = [x_cut, x_cut],  [y_cut, -y_cut]
        x4, y4 = [x_cut, -x_cut],  [-y_cut, -y_cut]
        plt.plot(x1, y1, x2, y2, x3, y3, x4, y4, color='k')

    #divider = make_axes_locatable(axbig)
    #cax = divider.append_axes("bottom", size="5%", pad=0.7)
    #cbar = fig.colorbar(im,orientation="horizontal",cax=cax)
    #cbar.set_label(label_z)
    
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

def BackgroundSubtractMap(fig,hist_data,hist_bkgd,label_x,label_y,label_z,plotname):

    colormap = 'coolwarm'

    Old_map_nbins_x = hist_data.GetNbinsX()
    Old_map_nbins_y = hist_data.GetNbinsY()
    Old_MapEdge_left = hist_data.GetXaxis().GetBinLowEdge(1)
    Old_MapEdge_right = hist_data.GetXaxis().GetBinLowEdge(hist_data.GetNbinsX()+1)
    Old_MapEdge_lower = hist_data.GetYaxis().GetBinLowEdge(1)
    Old_MapEdge_upper = hist_data.GetYaxis().GetBinLowEdge(hist_data.GetNbinsY()+1)
    Old_MapEdge_center_x = 0.5*(Old_MapEdge_left+Old_MapEdge_right)
    Old_MapEdge_center_y = 0.5*(Old_MapEdge_lower+Old_MapEdge_upper)
    Old_MapEdge_size_x = Old_MapEdge_right-Old_MapEdge_center_x
    Old_MapEdge_size_y = Old_MapEdge_upper-Old_MapEdge_center_y

    map_bin_size = 2.*Old_MapEdge_size_x/float(Old_map_nbins_x)
    map_nbins_x = 50
    map_nbins_y = 50
    MapEdge_left = Old_MapEdge_center_x-int(map_nbins_x/2)*map_bin_size
    MapEdge_right = Old_MapEdge_center_x+int(map_nbins_x/2)*map_bin_size
    MapEdge_lower = Old_MapEdge_center_y-int(map_nbins_y/2)*map_bin_size
    MapEdge_upper = Old_MapEdge_center_y+int(map_nbins_y/2)*map_bin_size

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins_x
    nbins_per_deg = map_nbins_x/(MapEdge_right-MapEdge_left)
    old_x_axis = np.linspace(Old_MapEdge_left,Old_MapEdge_right,Old_map_nbins_x)
    old_y_axis = np.linspace(Old_MapEdge_lower,Old_MapEdge_upper,Old_map_nbins_y)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins_x)
    x_axis_sparse = np.linspace(MapEdge_left,MapEdge_right,5)
    x_axis_reflect = ["{:6.1f}".format(-1.*i) for i in x_axis_sparse]
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins_y)

    grid_data = np.zeros((map_nbins_y, map_nbins_x))
    grid_bkgd = np.zeros((map_nbins_y, map_nbins_x))
    grid_excs = np.zeros((map_nbins_y, map_nbins_x))
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            #hist_bin_x = hist_data.GetXaxis().FindBin(x_axis[xbin])
            #hist_bin_y = hist_data.GetYaxis().FindBin(y_axis[ybin])
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_data.GetNbinsX(): continue
            if hist_bin_y>hist_data.GetNbinsY(): continue
            grid_data[ybin,xbin] = hist_data.GetBinContent(hist_bin_x,hist_bin_y)
            grid_bkgd[ybin,xbin] = hist_bkgd.GetBinContent(hist_bin_x,hist_bin_y)
            grid_excs[ybin,xbin] = hist_data.GetBinContent(hist_bin_x,hist_bin_y)-hist_bkgd.GetBinContent(hist_bin_x,hist_bin_y)

    # Define the locations for the axes
    left, width = 0.12, 0.6
    bottom, height = 0.12, 0.6
    bottom_h = bottom+height+0.03
    left_h = left+width+0.03
     
    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.20] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.20, height] # dimensions of y-histogram

    # Set up the size of the figure
    #fig = plt.figure(1, figsize=(9.5,9))
    fig = plt.figure(1)
    fig.set_figheight(8)
    fig.set_figwidth(8)

    fig.clf()
    # Make the three plots
    axTemperature = plt.axes(rect_temperature) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Plot the temperature data
    cax = (axTemperature.imshow(grid_data, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()), origin='lower', cmap=colormap))

    #Plot the axes labels
    axTemperature.set_xlabel(label_x)
    axTemperature.set_ylabel(label_y)
    axTemperature.set_xticks(x_axis_sparse)
    axTemperature.set_xticklabels(x_axis_reflect)

    x_profile, x_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_data,proj_type='X')
    y_profile, y_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_data,proj_type='Y')
    x_bkgd_profile, x_bkgd_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_bkgd,proj_type='X')
    y_bkgd_profile, y_bkgd_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_bkgd,proj_type='Y')
    #Plot the histograms
    axHistx.errorbar(x_ax,x_profile,yerr=x_profile_stat_err,color='k',marker='.',ls='none')
    axHisty.errorbar(y_profile,y_ax,xerr=y_profile_stat_err,color='k',marker='.',ls='none')
    axHistx.plot(x_ax,x_bkgd_profile,color='r',ls='solid')
    axHisty.plot(y_bkgd_profile,y_ax,color='r',ls='solid')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    axHistx.yaxis.set_major_formatter(formatter)
    axHisty.xaxis.set_major_formatter(formatter)

    fig.savefig("output_plots/%s_Before.png"%(plotname),bbox_inches='tight')

    fig.clf()
    # Make the three plots
    axTemperature = plt.axes(rect_temperature) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Plot the temperature data
    cax = (axTemperature.imshow(grid_excs, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()), origin='lower', cmap=colormap))

    #Plot the axes labels
    axTemperature.set_xlabel(label_x)
    axTemperature.set_ylabel(label_y)
    axTemperature.set_xticks(x_axis_sparse)
    axTemperature.set_xticklabels(x_axis_reflect)

    nbins = hist_data.GetNbinsX()
    MapEdge_left = hist_data.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_data.GetXaxis().GetBinLowEdge(hist_data.GetNbinsX()+1)
    MapEdge_bottom = hist_data.GetYaxis().GetBinLowEdge(1)
    MapEdge_top = hist_data.GetYaxis().GetBinLowEdge(hist_data.GetNbinsY()+1)
    hist_excs = ROOT.TH2D("hist_excs","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_bottom,MapEdge_top)
    hist_excs.Reset()
    hist_excs.Add(hist_data)
    hist_excs.Add(hist_bkgd,-1.)

    x_profile, x_profile_stat_err, x_ax, x_ax_err = FindCountProjection(hist_excs,proj_type='X')
    y_profile, y_profile_stat_err, y_ax, y_ax_err = FindCountProjection(hist_excs,proj_type='Y')
    #Plot the histograms
    axHistx.errorbar(x_ax,x_profile,yerr=x_profile_stat_err,color='k',marker='.',ls='none')
    axHisty.errorbar(y_profile,y_ax,xerr=y_profile_stat_err,color='k',marker='.',ls='none')

    fig.savefig("output_plots/%s_After.png"%(plotname),bbox_inches='tight')


def MatplotlibMap2D(hist_map,hist_contour,fig,label_x,label_y,label_z,plotname,roi_x=0.,roi_y=0.,roi_r=0.,rotation_angle=0.,fill_gaps=False,prime_psr_name=None,prime_psr_ra=None,prime_psr_dec=None):

    print ('Making plot %s...'%(plotname))
    isGalactic = False
    if label_x=='gal. l':
        isGalactic = True

    top = cm.get_cmap('Blues_r', 128) # r means reversed version
    bottom = cm.get_cmap('Oranges', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')
    colormap = 'coolwarm'

    map_nbins_x = hist_map.GetNbinsX()
    map_nbins_y = hist_map.GetNbinsY()
    MapEdge_left = hist_map.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_map.GetXaxis().GetBinLowEdge(hist_map.GetNbinsX()+1)
    MapEdge_lower = hist_map.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_map.GetYaxis().GetBinLowEdge(hist_map.GetNbinsY()+1)

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins_x
    nbins_per_deg = map_nbins_x/(MapEdge_right-MapEdge_left)
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins_x)
    x_axis_sparse = np.linspace(MapEdge_left,MapEdge_right,5)
    x_axis_reflect = ["{:6.2f}".format(-1.*i) for i in x_axis_sparse]
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins_y)

    MapWidth = (MapEdge_right-MapEdge_left)/2.
    MapHeight = (MapEdge_upper-MapEdge_lower)/2.
    New_MapEdge_left = MapEdge_left+0.5*MapWidth
    New_MapEdge_right = MapEdge_right-0.5*MapWidth
    New_MapEdge_lower = MapEdge_lower+0.5*MapHeight
    New_MapEdge_upper = MapEdge_upper-0.5*MapHeight
    new_x_axis = np.linspace(New_MapEdge_left,New_MapEdge_right,int(0.5*map_nbins_x))
    new_y_axis = np.linspace(New_MapEdge_lower,New_MapEdge_upper,int(0.5*map_nbins_y))

    mean_z = 0.
    rms_z = 0.
    n_bins = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            #hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
            #hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            if hist_map.GetBinContent(hist_bin_x,hist_bin_y)==0.: continue
            mean_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            n_bins += 1.
    if n_bins>0.:
        mean_z = mean_z/n_bins
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            #hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
            #hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            if hist_map.GetBinContent(hist_bin_x,hist_bin_y)==0.: continue
            rms_z += pow(hist_map.GetBinContent(hist_bin_x,hist_bin_y)-mean_z,2)
            n_bins += 1.
    if n_bins>0.:
        rms_z = rms_z/n_bins
        rms_z = pow(rms_z,0.5)

    grid_z = np.zeros((map_nbins_y, map_nbins_x))
    max_z = 0.
    min_z = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            #hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
            #hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            grid_z[ybin,xbin] = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if max_z<hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                max_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            if min_z>hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                min_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)

    if not hist_contour==None and 'SkymapFlux' in plotname:
        max_zscore = hist_contour.GetMaximum()
        mid_class_z = 0.
        mid_class_z_bins = 0.
        low_class_z = 0.
        low_class_z_bins = 0.
        high_class_z = 0.
        high_class_z_bins = 0.
        for ybin in range(0,len(y_axis)):
            for xbin in range(0,len(x_axis)):
                #hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
                #hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
                hist_bin_x = xbin+1
                hist_bin_y = ybin+1
                if hist_bin_x<1: continue
                if hist_bin_y<1: continue
                if hist_bin_x>hist_map.GetNbinsX(): continue
                if hist_bin_y>hist_map.GetNbinsY(): continue
                zscore = hist_contour.GetBinContent(hist_bin_x,hist_bin_y)
                if zscore>1.-0.5 and zscore<1.+0.5:
                    mid_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    mid_class_z_bins += 1.
                if zscore>0.-0.5 and zscore<0.+0.5:
                    low_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    low_class_z_bins += 1.
                if zscore>max_zscore-0.5:
                    high_class_z += hist_map.GetBinContent(hist_bin_x,hist_bin_y)
                    high_class_z_bins += 1.
        if mid_class_z_bins>0.:
            mid_class_z = mid_class_z/mid_class_z_bins
        if low_class_z_bins>0.:
            low_class_z = low_class_z/low_class_z_bins
        if high_class_z_bins>0.:
            high_class_z = high_class_z/high_class_z_bins
        min_z = low_class_z-2.*(mid_class_z-low_class_z)
        max_z = max(high_class_z*1.1,mid_class_z+2.*(mid_class_z-low_class_z))

    #levels = np.arange(3.0, 6.0, 1.0)
    levels = np.arange(2.0, 6.0, 1.0)

    grid_contour = np.zeros((map_nbins_y, map_nbins_x))
    if not hist_contour==None:
        max_z_contour = 5.0
        if label_z!='Z score' and max_z>0.:
            levels = np.arange(2.0*max_z/max_z_contour, 6.0*max_z/max_z_contour, 1.0*max_z/max_z_contour)
        if 'SkymapFlux' in plotname and max_z>0.:
            levels = np.arange(2.0*max_z/max_z_contour, 6.0*max_z/max_z_contour, 1.0*max_z/max_z_contour)
        for ybin in range(0,len(y_axis)):
            for xbin in range(0,len(x_axis)):
                #hist_bin_x = hist_contour.GetXaxis().FindBin(x_axis[xbin])
                #hist_bin_y = hist_contour.GetYaxis().FindBin(y_axis[ybin])
                hist_bin_x = xbin+1
                hist_bin_y = ybin+1
                if hist_bin_x<1: continue
                if hist_bin_y<1: continue
                if hist_bin_x>hist_contour.GetNbinsX(): continue
                if hist_bin_y>hist_contour.GetNbinsY(): continue
                grid_contour[ybin,xbin] = hist_contour.GetBinContent(hist_bin_x,hist_bin_y)
                if label_z!='Z score' and max_z_contour>0.:
                    grid_contour[ybin,xbin] = hist_contour.GetBinContent(hist_bin_x,hist_bin_y)*max_z/max_z_contour

    other_stars, other_star_type, other_star_coord = GetGammaSourceInfo(hist_contour,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec) 

    other_star_labels = []
    other_star_types = []
    other_star_markers = []
    star_range = 0.9*(MapEdge_right-MapEdge_left)/2.
    source_ra = (MapEdge_left+MapEdge_right)/2.
    source_dec = (MapEdge_lower+MapEdge_upper)/2.
    n_stars = 0
    for star in range(0,len(other_stars)):
        if doGalacticCoord:
            if abs(source_ra+other_star_coord[star][0])>star_range: continue
            if abs(source_dec-other_star_coord[star][1])>0.5*star_range: continue
        else:
            #if pow(source_ra+other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
            if abs(source_ra+other_star_coord[star][0])>star_range: continue
            if abs(source_dec-other_star_coord[star][1])>star_range: continue
        #if '#' in other_stars[star]: continue
        other_star_markers += [[-other_star_coord[star][0],other_star_coord[star][1]]]
        other_star_labels += ['(%s) %s'%(n_stars,other_stars[star])]
        other_star_types += [other_star_type[star]]
        n_stars += 1
    #print ('len(other_stars) = %s'%len(other_stars))
    #print ('len(other_star_markers) = %s'%len(other_star_markers))

    fig.clf()
    figsize_x = 14
    figsize_y = 7
    if doGalacticCoord:
        figsize_x = 10
        figsize_y = 5
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    if doGalacticCoord:
        label_x = 'gal. l'
        label_y = 'gal. b'
    else:
        label_x = 'RA'
        label_y = 'Dec'
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    if label_z=='Z score':
        max_z = 5.
        #max_z = 8.
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=-max_z,vmax=max_z,zorder=0)
    elif fill_gaps:
        max_z = mean_z+3.*rms_z
        min_z = mean_z-3.*rms_z
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=min_z,vmax=max_z,zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=min_z,vmax=max_z,zorder=0)
    if not doGalacticCoord and not 'SkymapZscore_Sum' in plotname:
        axbig.contour(grid_contour, levels, colors='w', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=1)
    #cbar = fig.colorbar(im)
    #im_ratio = grid_z.shape[1]/grid_z.shape[0]
    #cbar = fig.colorbar(im,orientation="horizontal",fraction=0.047*im_ratio)
    for star in range(0,len(other_star_markers)):
        marker_size = 60
        if other_star_types[star]=='PSR':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='lime', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='SNR':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='yellow', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='HAWC':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='violet', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='Fermi':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='cyan', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='MSP':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='tomato', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='TeV':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='k', marker='^', label=other_star_labels[star])
        if other_star_types[star]=='Star':
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], s=marker_size, c='k', marker='o', label=other_star_labels[star])
        text_offset_x = 0.
        text_offset_y = 0.
        text_offset_x = 0.05
        #plt.annotate(other_star_labels[star], (other_star_markers[star][0]+text_offset_x, other_star_markers[star][1]+text_offset_y), fontsize=10, color='k', rotation = rotation_angle)
        plt.annotate('%s'%(star), (other_star_markers[star][0]+text_offset_x, other_star_markers[star][1]+text_offset_y), fontsize=12, color='k')
    divider = make_axes_locatable(axbig)
    cax = divider.append_axes("bottom", size="5%", pad=0.7)
    cbar = fig.colorbar(im,orientation="horizontal",cax=cax)
    cbar.set_label(label_z)
    axbig.set_xticks(x_axis_sparse)
    axbig.set_xticklabels(x_axis_reflect)
    if 'SNR_G189_p03' in plotname:
        mycircle = plt.Circle((-94.25, 22.57), 0.35, color='b', fill=False)
        axbig.add_patch(mycircle)
        mycircle = plt.Circle((-94.62, 22.17), 1.0, color='r', fill=False)
        axbig.add_patch(mycircle)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    if 'SkymapZscore_Sum' in plotname:
        if not doGalacticCoord:
            if roi_r>0.:
                mycircle = plt.Circle((-roi_x, roi_y), roi_r, color='b', fill=False)
                axbig.add_patch(mycircle)
        axbig.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0, fontsize=7)
        fig.savefig("output_plots/%s_legend.png"%(plotname),bbox_inches='tight')
    axbig.remove()

    if 'SkymapZscore_Sum' in plotname:
        for star in range(0,len(other_star_markers)):
            print ('%s, (%0.2f,%0.2f)'%(other_star_labels[star],-other_star_markers[star][0],other_star_markers[star][1]))

    figsize_x = 6.4
    figsize_y = 4.8
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)

def GetSkyViewMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open(map_file)
    for line in inputFile:
        sig = float(line.split(' ')[2])
        l = float(line.split(' ')[0])
        b = float(line.split(' ')[1])
        if not isRaDec: 
            l, b = ConvertRaDecToGalactic(l,b)
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        if binx<1: continue
        if binx>=hist_map.GetNbinsX()-1: continue
        if biny<1: continue
        if biny>=hist_map.GetNbinsY()-1: continue
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,old_sig))
        #hist_map.SetBinContent(binx+1,biny+1,sig+old_sig)

    map_resolution = 0.1
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)
    return hist_map_new

