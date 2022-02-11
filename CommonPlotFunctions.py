
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


skymap_zoomin_scale = 1
#skymap_zoomin_scale = 1.5
#skymap_zoomin_scale = 2
#n_rebins = 1
n_rebins = 2
#n_rebins = 4
#n_rebins = 8
#smooth_size_spectroscopy = 0.1
smooth_size_spectroscopy = 0.15
#smooth_size_spectroscopy = 0.2
#smooth_size_spectroscopy = 0.3

calibration_radius = 0.2

Skymap_size = 3.
Skymap_nbins = 120

energy_bin = []
energy_bin += [100]
energy_bin += [251]
energy_bin += [631]
energy_bin += [1585]
energy_bin += [3981]
energy_bin += [10000]
energy_bin += [25118]

energy_bin_big = []
energy_bin_big += [100]
energy_bin_big += [631]
energy_bin_big += [3981]
energy_bin_big += [25118]

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

    Hist_Kernel = ROOT.TH2D("Hist_Kernel","",Skymap_nbins,-Skymap_size,Skymap_size,Skymap_nbins,-Skymap_size,Skymap_size)
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
    central_bin = int(Skymap_nbins/2)
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
    n_bins = int((integration_range+extension_low+extension_up)/(0.1*n_rebins))
    #Hist_Profile_ProjectedX = ROOT.TH1D("Hist_Profile_ProjectedX","",n_bins,-extension_low,extension_up+integration_range)
    #Hist_Profile_ProjectedX = ROOT.TH1D("Hist_Profile_ProjectedX","",int(20/n_rebins),-1.5,1.5)
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


def FindExtension(Hist_Data_input,Hist_Syst_input,roi_x,roi_y,integration_range):

    global calibration_radius

    n_bins = 10
    #integration_range = 1.8
    if sys.argv[1]=='Geminga_ON':
        n_bins = 5
        integration_range = 1.8
    if sys.argv[1]=='LHAASO_J2108_ON':
        n_bins = 5
        integration_range = 1.8

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
    if sys.argv[1]=='MGRO_J1908_ON':
        inputFile = open('J1908_sources_RaDec_w_Names.txt')
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

