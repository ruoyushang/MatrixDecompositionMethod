
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

folder_path = 'output_test'
#folder_path = 'output_2hrs'
#folder_path = 'output_5hrs'
#folder_path = 'output_10hrs'
#folder_path = 'output_20hrs'
#folder_path = 'output_4x4'
#folder_path = 'output_8x8'
#folder_path = 'output_12x12'
#folder_path = 'output_FreeNSB'
#folder_path = 'output_FreeAzim'
#folder_path = 'output_FreeElev'
#folder_path = 'output_radialcorrect'

#N_bins_for_deconv = 4
#N_bins_for_deconv = 8
N_bins_for_deconv = 12
gamma_hadron_dim_ratio_w = 1.
gamma_hadron_dim_ratio_l = 1.
gamma_hadron_low_end = 0.

MSCW_blind_cut = 0.6
MSCL_blind_cut = 0.6

skymap_zoomin_scale = 1
#skymap_zoomin_scale = 1.5
#skymap_zoomin_scale = 2
smooth_size_spectroscopy = 0.1
#smooth_size_spectroscopy = 0.15
#smooth_size_spectroscopy = 0.2
#smooth_size_spectroscopy = 0.3

#Smoothing = True
Smoothing = False


calibration_radius = 0.5 # need to be larger than the PSF and smaller than the integration radius

energy_index_scale = 2

#doGalacticCoord = True
#Skymap_size_x = 4.
#Skymap_nbins_x = 90
#Skymap_size_y = 2.
#Skymap_nbins_y = 45
doGalacticCoord = False
Skymap_size_x = 2.
Skymap_nbins_x = 45
Skymap_size_y = 2.
Skymap_nbins_y = 45
#Skymap_nbins = 15 
#Skymap_nbins = 9 # Crab calibration 
#Skymap_nbins = 5
#Skymap_nbins = 3

Skymap_normalization_nbins = 1

elev_range = [45,90]
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

def Smooth2DMap_v2(Hist_Old,Hist_Smooth,smooth_size,addLinearly,normalized):

    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    if bin_size>smooth_size: 
        Hist_Smooth.Reset()
        Hist_Smooth.Add(Hist_Old)
        return

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
                    if not addLinearly:
                        bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*pow(Hist_Old.GetBinError(bx2,by2),2)
                    else:
                        bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    if normalized:
        Hist_Smooth.Scale(1./Hist_Kernel.Integral())
        #for bx1 in range(1,Hist_Smooth.GetNbinsX()+1):
        #    for by1 in range(1,Hist_Smooth.GetNbinsY()+1):
        #        old_error = Hist_Old.GetBinError(bx1,by1)
        #        Hist_Smooth.SetBinError(bx1,by1,old_error)

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
                    if not addLinearly:
                        bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*pow(Hist_Old.GetBinError(bx2,by2),2)
                    else:
                        bin_error += Hist_Kernel.GetBinContent(bx2-bx1+central_bin,by2-by1+central_bin)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    if normalized:
        Hist_Smooth.Scale(1./Hist_Kernel.Integral())
        #for bx1 in range(1,Hist_Smooth.GetNbinsX()+1):
        #    for by1 in range(1,Hist_Smooth.GetNbinsY()+1):
        #        old_error = Hist_Old.GetBinError(bx1,by1)
        #        Hist_Smooth.SetBinError(bx1,by1,old_error)

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

def FindGalacticProjection(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range):

    global calibration_radius

    roi_l, roi_b = ConvertRaDecToGalactic(roi_x,roi_y)

    n_bins = 10
    Hist_Profile_Gal_l = ROOT.TH1D("Hist_Profile_Gal_l","",n_bins,roi_l-2.,roi_l+2.)
    n_bins = 10
    Hist_Profile_Gal_b = ROOT.TH1D("Hist_Profile_Gal_b","",n_bins,roi_b-2.,roi_b+2.)

    for br in range(0,Hist_Profile_Gal_l.GetNbinsX()):
        range_limit = Hist_Profile_Gal_l.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Gal_l.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_syst_err = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                cell_l, cell_b = ConvertRaDecToGalactic(cell_x,cell_y)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                syst_error = 0.
                if Hist_Syst_input!=None:
                    syst_error = Hist_Syst_input.GetBinContent(bx+1,by+1)
                cell_err = pow(syst_error*syst_error+data_error*data_error,0.5)
                if cell_err>0.9: continue
                delta_b = abs(cell_b-roi_b)
                if delta_b>integration_range: continue
                if cell_l>=range_limit_previous and cell_l<range_limit:
                    slice_data += data_content
                    slice_data_err += data_error*data_error
                    slice_syst_err += syst_error
        slice_data_err = pow(slice_data_err,0.5)
        if slice_data==0.: slice_data_err = 1.
        Hist_Profile_Gal_l.SetBinContent(br+1,slice_data)
        Hist_Profile_Gal_l.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))

    for br in range(0,Hist_Profile_Gal_b.GetNbinsX()):
        range_limit = Hist_Profile_Gal_b.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Gal_b.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_syst_err = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                cell_l, cell_b = ConvertRaDecToGalactic(cell_x,cell_y)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                syst_error = 0.
                if Hist_Syst_input!=None:
                    syst_error = Hist_Syst_input.GetBinContent(bx+1,by+1)
                cell_err = pow(syst_error*syst_error+data_error*data_error,0.5)
                if cell_err>0.9: continue
                delta_l = abs(cell_l-roi_l)
                if delta_l>integration_range: continue
                if cell_b>=range_limit_previous and cell_b<range_limit:
                    slice_data += data_content
                    slice_data_err += data_error*data_error
                    slice_syst_err += syst_error
        slice_data_err = pow(slice_data_err,0.5)
        if slice_data==0.: slice_data_err = 1.
        Hist_Profile_Gal_b.SetBinContent(br+1,slice_data)
        Hist_Profile_Gal_b.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))


    profile_l = []
    profile_err_l = []
    theta2_l = []
    for binx in range(0,Hist_Profile_Gal_l.GetNbinsX()):
        center = Hist_Profile_Gal_l.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Gal_l.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Gal_l.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Gal_l.GetBinContent(binx+1)
        profile_error = Hist_Profile_Gal_l.GetBinError(binx+1)
        if profile_error>0.9: continue
        theta2_l += [center]
        profile_l += [profile_content]
        profile_err_l += [profile_error]
    profile_b = []
    profile_err_b = []
    theta2_b = []
    for binx in range(0,Hist_Profile_Gal_b.GetNbinsX()):
        center = Hist_Profile_Gal_b.GetBinCenter(binx+1)
        range_limit = Hist_Profile_Gal_b.GetBinLowEdge(binx+2)
        range_limit_previous = Hist_Profile_Gal_b.GetBinLowEdge(binx+1)
        profile_content = Hist_Profile_Gal_b.GetBinContent(binx+1)
        profile_error = Hist_Profile_Gal_b.GetBinError(binx+1)
        if profile_error>0.9: continue
        theta2_b += [center]
        profile_b += [profile_content]
        profile_err_b += [profile_error]

    profile = []
    profile_err = []
    theta2 = []
    profile += [profile_l]
    profile_err += [profile_err_l]
    theta2 += [theta2_l]
    profile += [profile_b]
    profile_err += [profile_err_b]
    theta2 += [theta2_b]

    return profile, profile_err, theta2

def FindGalacticProjection_v2(Hist_Data_input,Hist_Syst_input):

    n_bins_y = Hist_Data_input.GetNbinsY()
    n_bins_x = Hist_Data_input.GetNbinsX()
    MapEdge_left = Hist_Data_input.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Data_input.GetXaxis().GetBinLowEdge(Hist_Data_input.GetNbinsX()+1)
    MapEdge_lower = Hist_Data_input.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Data_input.GetYaxis().GetBinLowEdge(Hist_Data_input.GetNbinsY()+1)
    cell_size = (MapEdge_right-MapEdge_left)/float(n_bins_x)*(MapEdge_upper-MapEdge_lower)/float(n_bins_y)

    Hist_Profile_Y = ROOT.TH1D("Hist_Profile_Y","",n_bins_y,MapEdge_lower,MapEdge_upper)
    for br in range(0,Hist_Profile_Y.GetNbinsX()):
        range_limit = Hist_Profile_Y.GetBinLowEdge(br+2)
        range_limit_previous = Hist_Profile_Y.GetBinLowEdge(br+1)
        slice_data = 0.
        slice_data_err = 0.
        slice_syst_err = 0.
        solid_angle = 0.
        for bx in range(0,Hist_Data_input.GetNbinsX()):
            for by in range(0,Hist_Data_input.GetNbinsY()):
                cell_x = Hist_Data_input.GetXaxis().GetBinCenter(bx+1)
                cell_y = Hist_Data_input.GetYaxis().GetBinCenter(by+1)
                data_content = Hist_Data_input.GetBinContent(bx+1,by+1)
                data_error = Hist_Data_input.GetBinError(bx+1,by+1)
                syst_error = 0.
                if Hist_Syst_input!=None:
                    syst_error = Hist_Syst_input.GetBinContent(bx+1,by+1)
                if cell_y>=range_limit_previous and cell_y<range_limit:
                    if not data_content==0.:
                        slice_data += data_content
                        slice_data_err += data_error*data_error
                        slice_syst_err += syst_error
                        solid_angle += cell_size
        if solid_angle==0.: 
            slice_data = 0.
            slice_data_err = 0.
        else:
            slice_data = slice_data/solid_angle
            slice_data_err = pow(slice_data_err,0.5)/solid_angle
        Hist_Profile_Y.SetBinContent(br+1,slice_data)
        Hist_Profile_Y.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))

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


def FindExtension_v2(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins_2d = Hist_Data_input.GetNbinsX()
    #n_bins_1d = 4
    #n_bins_1d = min(8,int(float(n_bins_2d)/2.))
    n_bins_1d = min(16,int(float(n_bins_2d)/2.))
    #integration_range = 0.8
    integration_range = 1.6

    Hist_Profile_Theta2 = ROOT.TH1D("Hist_Profile_Theta2","",n_bins_1d,0,integration_range)
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
        if slice_data==0.: 
            slice_data = 0.
            slice_data_err = 0.
        Hist_Profile_Theta2.SetBinContent(br+1,slice_data)
        Hist_Profile_Theta2.SetBinError(br+1,pow(slice_data_err*slice_data_err+slice_syst_err*slice_syst_err,0.5))


    profile = []
    profile_err = []
    theta2 = []
    theta2_err = []
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

def GetGammaSourceInfo():

    other_stars = []
    other_star_coord = []
    if 'MGRO_J1908' in sys.argv[1]:
        inputFile = open('J1908_sources_RaDec_w_Names.txt')
        #inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<0.*0.:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]
    else:

        inputFile = open('PSR_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<0.3*0.3:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]
        inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
                if distance<0.3*0.3:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    return other_stars, other_star_coord

def GetRegionIntegral(hist_data_skymap,hist_syst_skymap,hist_mask_skymap,roi_x,roi_y,roi_r):

    flux_sum = 0.
    flux_stat_err = 0.
    flux_syst_err = 0.
    for bx in range(0,hist_data_skymap.GetNbinsX()):
        for by in range(0,hist_data_skymap.GetNbinsY()):
            bin_ra = hist_data_skymap.GetXaxis().GetBinCenter(bx+1)
            bin_dec = hist_data_skymap.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_ra-roi_x,2) + pow(bin_dec-roi_y,2),0.5)
            if distance>roi_r: 
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

def GetRegionSpectrum_v2(hist_data_skymap,hist_syst_skymap,hist_mask_skymap,ebin_low,ebin_up,roi_x,roi_y,roi_r):

    x_axis = []
    x_error = []
    y_axis = []
    y_error = []
    for ebin in range(ebin_low,ebin_up):
        flux_sum = 0.
        flux_stat_err = 0.
        flux_syst_err = 0.
        flux_sum, flux_stat_err, flux_syst_err = GetRegionIntegral(hist_data_skymap[ebin],hist_syst_skymap[ebin],hist_mask_skymap,roi_x,roi_y,roi_r)
        x_axis += [0.5*(energy_bin[ebin]+energy_bin[ebin+1])]
        x_error += [0.5*(energy_bin[ebin+1]-energy_bin[ebin])]
        y_axis += [flux_sum]
        #y_error += [pow(flux_stat_err*flux_stat_err+flux_syst_err*flux_syst_err,0.5)]
        y_error += [flux_stat_err]

    return x_axis, x_error, y_axis, y_error

def GetRegionSpectrum(hist_data_skymap,hist_syst_skymap,hist_mask_skymap,ebin_low,ebin_up,roi_x,roi_y,roi_r):

    x_axis = []
    y_axis = []
    y_error = []
    for ebin in range(ebin_low,ebin_up):
        flux_sum = 0.
        flux_stat_err = 0.
        flux_syst_err = 0.
        flux_sum, flux_stat_err, flux_syst_err = GetRegionIntegral(hist_data_skymap[ebin],hist_syst_skymap[ebin],hist_mask_skymap,roi_x,roi_y,roi_r)
        x_axis += [energy_bin[ebin]]
        y_axis += [flux_sum]
        #y_error += [pow(flux_stat_err*flux_stat_err+flux_syst_err*flux_syst_err,0.5)]
        y_error += [flux_stat_err]

    return x_axis, y_axis, y_error

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

def MatplotlibMultiMaps(hist_maps,fig,label_x,label_y,label_z,plotname):

    top = cm.get_cmap('Blues_r', 5) # r means reversed version
    bottom = cm.get_cmap('Reds', 5)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 5)),bottom(np.linspace(0, 1, 5))))# create a new colormaps with a name of OrangeBlue
    redblue = ListedColormap(newcolors, name='RedBlue')

    N = 256
    red = np.ones((N, 4))
    red[:, 0] = np.linspace(255/256, 1, N) # R = 255
    red[:, 1] = np.linspace(0/256, 1, N)   # G = 0
    red[:, 2] = np.linspace(0/256, 1, N)   # B = 0
    red_cmp = ListedColormap(red[::-1])
    blue = np.ones((N, 4))
    blue[:, 0] = np.linspace(0/256, 1, N)   # R = 0
    blue[:, 1] = np.linspace(0/256, 1, N)   # G = 0
    blue[:, 2] = np.linspace(255/256, 1, N) # B = 255
    blue_cmp = ListedColormap(blue[::-1])
    map_colors = [blue_cmp,red_cmp]

    map_nbins = hist_maps[0].GetNbinsX()
    MapEdge_left = hist_maps[0].GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_maps[0].GetXaxis().GetBinLowEdge(hist_maps[0].GetNbinsX()+1)
    MapEdge_lower = hist_maps[0].GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_maps[0].GetYaxis().GetBinLowEdge(hist_maps[0].GetNbinsY()+1)

    deg_per_bin = (MapEdge_right-MapEdge_left)/map_nbins
    nbins_per_deg = map_nbins/(MapEdge_right-MapEdge_left)
    nbins_extend = int(0.*nbins_per_deg)
    deg_extend = nbins_extend/nbins_per_deg
    x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins)
    x_axis_extend = np.linspace(MapEdge_left-deg_extend,MapEdge_right,map_nbins+nbins_extend)
    x_axis_sparse = np.linspace(MapEdge_left-deg_extend,MapEdge_right,5)
    x_axis_reflect = ["{:6.2f}".format(-1.*i) for i in x_axis_sparse]
    y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins)

    grid_z_list = []
    zmax_list = []
    for nmap in range(0,len(hist_maps)):
        grid_z = np.zeros((map_nbins, map_nbins+nbins_extend))
        zmax = 0.
        for ybin in range(0,len(y_axis)):
            for xbin in range(0,len(x_axis_extend)):
                hist_bin_x = xbin+1
                hist_bin_y = ybin+1
                if hist_bin_x<1: continue
                if hist_bin_y<1: continue
                if hist_bin_x>hist_maps[nmap].GetNbinsX(): continue
                if hist_bin_y>hist_maps[nmap].GetNbinsY(): continue
                grid_z[ybin,xbin] = hist_maps[nmap].GetBinContent(hist_bin_x,hist_bin_y)
                if zmax<grid_z[ybin,xbin]:
                    zmax = grid_z[ybin,xbin]
        grid_z_list += [grid_z]
        zmax_list += [zmax]

    other_stars, other_star_coord = GetGammaSourceInfo() 
    other_star_labels = []
    other_star_markers = []
    star_range = 2.0
    source_ra = (MapEdge_left+MapEdge_right)/2.
    source_dec = (MapEdge_lower+MapEdge_upper)/2.
    for star in range(0,len(other_stars)):
        if pow(source_ra+other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
        #if '#' in other_stars[star]: continue
        other_star_markers += [[-other_star_coord[star][0],other_star_coord[star][1]]]
        other_star_labels += ['%s'%(other_stars[star])]

    fig.clf()
    axbig = fig.add_subplot()
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    #level_low = 2.0
    #level_high = 7.0
    for nmap in range(0,len(hist_maps)):
        mycolor = cm.get_cmap(map_colors[nmap], 3) # r means reversed version
        im = axbig.imshow(grid_z_list[nmap], origin='lower', cmap=mycolor, alpha=0.5, extent=(x_axis_extend.min(),x_axis_extend.max(),y_axis.min(),y_axis.max()), vmin=0.5*zmax_list[nmap])
        cbar = fig.colorbar(im)
        cbar.set_label(label_z)
        #axbig.contour(grid_z_list[nmap], levels, cmap=mycolor, extent=(x_axis_extend.min(),x_axis_extend.max(),y_axis.min(),y_axis.max()), vmin=0.)
    for star in range(0,len(other_star_markers)):
        axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], c='k', marker='+', label=other_star_labels[star])
        text_offset_x = 0.
        text_offset_y = 0.
        text_length = len(other_star_labels[star])/10.
        if source_ra+other_star_coord[star][0]>0.:
            text_offset_x = -text_length
        else: 
            text_offset_x = 0.1
        plt.annotate(other_star_labels[star], (other_star_markers[star][0]+text_offset_x, other_star_markers[star][1]+text_offset_y), fontsize=10, color='k')
    axbig.set_xticks(x_axis_sparse)
    axbig.set_xticklabels(x_axis_reflect)
    #axbig.legend(fontsize=7)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

def MatplotlibHist2D(hist_map,fig,label_x,label_y,label_z,plotname):

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
    axbig = fig.add_subplot()
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0)
    cbar = fig.colorbar(im)
    cbar.set_label(label_z)
    
    if 'Matrix' in plotname:
        x1, y1 = [-0.6, -0.6], [-0.6, 0.6]
        x2, y2 = [-0.6, 0.6],  [0.6, 0.6]
        x3, y3 = [0.6, 0.6],  [0.6, -0.6]
        x4, y4 = [0.6, -0.6],  [-0.6, -0.6]
        plt.plot(x1, y1, x2, y2, x3, y3, x4, y4, color='k')

    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

def MatplotlibMap2D(hist_map,hist_contour,fig,label_x,label_y,label_z,plotname,roi_x=0.,roi_y=0.,roi_r=0.,rotation_angle=0.,fill_gaps=False):

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

    mean_z = 0.
    rms_z = 0.
    n_bins = 0.
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
            hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
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
            hist_bin_x = hist_map.GetXaxis().FindBin(x_axis[xbin])
            hist_bin_y = hist_map.GetYaxis().FindBin(y_axis[ybin])
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
    for ybin in range(0,len(y_axis)):
        for xbin in range(0,len(x_axis)):
            hist_bin_x = xbin+1
            hist_bin_y = ybin+1
            if hist_bin_x<1: continue
            if hist_bin_y<1: continue
            if hist_bin_x>hist_map.GetNbinsX(): continue
            if hist_bin_y>hist_map.GetNbinsY(): continue
            grid_z[ybin,xbin] = hist_map.GetBinContent(hist_bin_x,hist_bin_y)
            #if fill_gaps and grid_z[ybin,xbin]==0.:
            #    grid_z[ybin,xbin] = mean_z
            if max_z<hist_map.GetBinContent(hist_bin_x,hist_bin_y):
                max_z = hist_map.GetBinContent(hist_bin_x,hist_bin_y)

    levels = np.arange(3.0, 6.0, 1.0)

    grid_contour = np.zeros((map_nbins_y, map_nbins_x))
    if not hist_contour==None:
        max_z_contour = 5.0
        if label_z!='Z score' and max_z>0.:
            levels = np.arange(3.0*max_z/max_z_contour, 6.0*max_z/max_z_contour, 1.0*max_z/max_z_contour)
        for ybin in range(0,len(y_axis)):
            for xbin in range(0,len(x_axis)):
                hist_bin_x = hist_contour.GetXaxis().FindBin(x_axis[xbin])
                hist_bin_y = hist_contour.GetYaxis().FindBin(y_axis[ybin])
                if hist_bin_x<1: continue
                if hist_bin_y<1: continue
                if hist_bin_x>hist_contour.GetNbinsX(): continue
                if hist_bin_y>hist_contour.GetNbinsY(): continue
                grid_contour[ybin,xbin] = hist_contour.GetBinContent(hist_bin_x,hist_bin_y)
                if label_z!='Z score' and max_z_contour>0.:
                    grid_contour[ybin,xbin] = hist_contour.GetBinContent(hist_bin_x,hist_bin_y)*max_z/max_z_contour
    other_stars, other_star_coord = GetGammaSourceInfo() 
    if isGalactic:
        for star in range(0,len(other_stars)):
            other_star_coord[star][0], other_star_coord[star][1] = ConvertRaDecToGalactic(other_star_coord[star][0], other_star_coord[star][1]) 

    other_star_labels = []
    other_star_markers = []
    star_range = 2.0
    source_ra = (MapEdge_left+MapEdge_right)/2.
    source_dec = (MapEdge_lower+MapEdge_upper)/2.
    for star in range(0,len(other_stars)):
        if pow(source_ra+other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
        #if '#' in other_stars[star]: continue
        other_star_markers += [[-other_star_coord[star][0],other_star_coord[star][1]]]
        other_star_labels += ['%s'%(other_stars[star])]
    #print ('len(other_stars) = %s'%len(other_stars))
    #print ('len(other_star_markers) = %s'%len(other_star_markers))

    fig.clf()
    axbig = fig.add_subplot()
    axbig.set_xlabel(label_x)
    axbig.set_ylabel(label_y)
    if label_z=='Z score':
        max_z = 5.
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=-max_z,vmax=max_z,zorder=0)
    elif fill_gaps:
        max_z = mean_z+3.*rms_z
        min_z = mean_z-3.*rms_z
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),vmin=min_z,vmax=max_z,zorder=0)
    else:
        im = axbig.imshow(grid_z, origin='lower', cmap=colormap, extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0)
    axbig.contour(grid_contour, levels, colors='w', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=1)
    #cbar = fig.colorbar(im,orientation="horizontal")
    cbar = fig.colorbar(im)
    cbar.set_label(label_z)
    if not doGalacticCoord:
        for star in range(0,len(other_star_markers)):
            axbig.scatter(other_star_markers[star][0], other_star_markers[star][1], c='k', marker='+', label=other_star_labels[star])
            text_offset_x = 0.
            text_offset_y = 0.
            text_length = len(other_star_labels[star])/10.
            text_offset_x = 0.1
            #if source_ra+other_star_coord[star][0]>0.:
            #    text_offset_x = -text_length
            #else: 
            #    text_offset_x = 0.1
            plt.annotate(other_star_labels[star], (other_star_markers[star][0]+text_offset_x, other_star_markers[star][1]+text_offset_y), fontsize=10, color='k', rotation = rotation_angle)
        if roi_r>0.:
            mycircle = plt.Circle((-roi_x, roi_y), roi_r, color='b', fill=False)
            axbig.add_patch(mycircle)
        if 'J1908' in plotname:
            mycircle = plt.Circle((-286.975, 6.038), 0.43, color='g', fill=False) # PSR J1907+0602
            axbig.add_patch(mycircle)
            mycircle = plt.Circle((-286.786, 6.498), 0.39, color='r', fill=False) # SNR GG40.5-0.5
            #mycircle = plt.Circle((-286.786, 6.498), 0.12, color='r', fill=False) # SNR GG40.5-0.5
            axbig.add_patch(mycircle)
    axbig.set_xticks(x_axis_sparse)
    axbig.set_xticklabels(x_axis_reflect)
    #axbig.legend(fontsize=7)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

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

