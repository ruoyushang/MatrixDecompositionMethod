
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

def simple_diffusion_model_PSR_flux(x,par):

    travel_interval = 0.15
    travel_angle = par[0]
    head_width = par[1]
    head_norm = par[2]
    core_width = par[3]
    core_norm = par[4]
    tail_width = par[5]
    tail_norm = par[6]

    RA_PSR = 286.975
    Dec_PSR = 6.03777777778

    PSR_head_x = 1.*travel_interval*np.cos(travel_angle*3.14/180.)+RA_PSR
    PSR_head_y = 1.*travel_interval*np.sin(travel_angle*3.14/180.)+Dec_PSR
    PSR_core_x = 2.*travel_interval*np.cos(travel_angle*3.14/180.)+PSR_head_x
    PSR_core_y = 2.*travel_interval*np.sin(travel_angle*3.14/180.)+PSR_head_y
    PSR_tail_x = 3.*travel_interval*np.cos(travel_angle*3.14/180.)+PSR_core_x
    PSR_tail_y = 3.*travel_interval*np.sin(travel_angle*3.14/180.)+PSR_core_y

    ## ax + by + c = 0
    #b = 1.
    #a = -b*(PSR_tail_y-PSR_head_y)/(PSR_tail_x-PSR_head_x)
    #c = -(b*PSR_head_y + a*PSR_head_x)

    #projected_x = abs(a*x[0]+b*x[1]+c)/pow(a*a+b*b,0.5)  # transverse distance to line
    #closest_x_on_line = (b*(b*x[0]-a*x[1])-a*c)/(a*a+b*b)
    #closest_y_on_line = (a*(-b*x[0]+a*x[1])-b*c)/(a*a+b*b)
    #projected_y_core = pow(pow(closest_x_on_line-PSR_core_x,2)+pow(closest_y_on_line-PSR_core_y,2),0.5)
    #if closest_y_on_line-PSR_core_y<0.:
    #    projected_y_core = -projected_y_core
    #projected_y_head = pow(pow(closest_x_on_line-PSR_head_x,2)+pow(closest_y_on_line-PSR_head_y,2),0.5)
    #if closest_y_on_line-PSR_head_y<0.:
    #    projected_y_head = -projected_y_head

    distance_to_head = pow(pow(x[0]-PSR_head_x,2)+pow(x[1]-PSR_head_y,2),0.5)
    distance_to_core = pow(pow(x[0]-PSR_core_x,2)+pow(x[1]-PSR_core_y,2),0.5)
    distance_to_tail = pow(pow(x[0]-PSR_tail_x,2)+pow(x[1]-PSR_tail_y,2),0.5)

    func_head = head_norm*np.exp(-distance_to_head*distance_to_head/(2.*head_width*head_width))
    func_core = core_norm*np.exp(-distance_to_core*distance_to_core/(2.*core_width*core_width))
    func_tail = tail_norm*np.exp(-distance_to_tail*distance_to_tail/(2.*tail_width*tail_width))

    return (func_head+func_core+func_tail)

def simple_diffusion_model_PSR(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    mask = Hist_bkgd_global.GetBinContent(binx,biny)
    bkgd_scale = 1.
    if mask!=0.:
        bkgd_scale = 0.

    travel_angle = par[0]
    travel_interval = par[1]
    head_width = par[2]
    head_norm = par[3]
    tail_width = par[4]
    tail_norm = par[5]

    RA_PSR = 286.975
    Dec_PSR = 6.03777777778

    PSR_head_x = 0.*travel_interval*np.cos(travel_angle*3.14/180.)+RA_PSR
    PSR_head_y = 0.*travel_interval*np.sin(travel_angle*3.14/180.)+Dec_PSR
    PSR_tail_x = 1.*travel_interval*np.cos(travel_angle*3.14/180.)+PSR_head_x
    PSR_tail_y = 1.*travel_interval*np.sin(travel_angle*3.14/180.)+PSR_head_y

    ## ax + by + c = 0
    #b = 1.
    #a = -b*(PSR_tail_y-PSR_head_y)/(PSR_tail_x-PSR_head_x)
    #c = -(b*PSR_head_y + a*PSR_head_x)

    #projected_x = abs(a*x[0]+b*x[1]+c)/pow(a*a+b*b,0.5)  # transverse distance to line
    #closest_x_on_line = (b*(b*x[0]-a*x[1])-a*c)/(a*a+b*b)
    #closest_y_on_line = (a*(-b*x[0]+a*x[1])-b*c)/(a*a+b*b)
    #projected_y_core = pow(pow(closest_x_on_line-PSR_core_x,2)+pow(closest_y_on_line-PSR_core_y,2),0.5)
    #if closest_y_on_line-PSR_core_y<0.:
    #    projected_y_core = -projected_y_core
    #projected_y_head = pow(pow(closest_x_on_line-PSR_head_x,2)+pow(closest_y_on_line-PSR_head_y,2),0.5)
    #if closest_y_on_line-PSR_head_y<0.:
    #    projected_y_head = -projected_y_head

    distance_to_head = pow(pow(x[0]-PSR_head_x,2)+pow(x[1]-PSR_head_y,2),0.5)
    distance_to_tail = pow(pow(x[0]-PSR_tail_x,2)+pow(x[1]-PSR_tail_y,2),0.5)

    func_head = head_norm*np.exp(-distance_to_head*distance_to_head/(2.*head_width*head_width))
    func_tail = tail_norm*np.exp(-distance_to_tail*distance_to_tail/(2.*tail_width*tail_width))

    return bkgd_scale*(func_head+func_tail)

def simple_diffusion_model_J1908(x,par):

    par_PWN = [par[0],par[1],par[2],par[3]]
    par_PSR_diffusion = [par[4],par[5],par[6],par[7],par[8],par[9]]
    par_CO_north = [par[10],par[11],par[12],par[13]]
    total = 0.
    total += symmetric_gauss_model(x,par_PWN)
    #total += disk_model(x,par_PWN)
    total += simple_diffusion_model_PSR(x,par_PSR_diffusion)
    total += symmetric_gauss_model(x,par_CO_north)
    return total

def simple_geometry_model_J1908(x,par):

    par_PWN = [par[0],par[1],par[2],par[3]]
    par_PSR_ellipse = [par[4],par[5],par[6],par[7],par[8],par[9]]
    total = 0.
    total += symmetric_gauss_model(x,par_PWN)
    #total += disk_model(x,par_PWN)
    total += asymmetric_gauss_model(x,par_PSR_ellipse)
    return total

def symmetric_gauss_model(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    mask = Hist_bkgd_global.GetBinContent(binx,biny)
    bkgd_scale = 1.
    if mask!=0.:
        bkgd_scale = 0.

    disk_center_x = par[0]
    disk_center_y = par[1]
    disk_radius = par[2]
    disk_norm = par[3]
    distance_sq = pow(x[0]-disk_center_x,2) + pow(x[1]-disk_center_y,2)
    return bkgd_scale*disk_norm*exp(-distance_sq/(2.*disk_radius*disk_radius))

def disk_model(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    mask = Hist_bkgd_global.GetBinContent(binx,biny)
    bkgd_scale = 1.
    if mask!=0.:
        bkgd_scale = 0.

    disk_center_x = par[0]
    disk_center_y = par[1]
    disk_radius = par[2]
    disk_norm = par[3]
    distance_sq = pow(x[0]-disk_center_x,2) + pow(x[1]-disk_center_y,2)
    if distance_sq>disk_radius*disk_radius:
        disk_norm = 0.
    return bkgd_scale*disk_norm

def asymmetric_gauss_model(x,par):

    binx = Hist_bkgd_global.GetXaxis().FindBin(x[0])
    biny = Hist_bkgd_global.GetYaxis().FindBin(x[1])
    mask = Hist_bkgd_global.GetBinContent(binx,biny)
    bkgd_scale = 1.
    if mask!=0.:
        bkgd_scale = 0.

    travel_distance = par[0]
    travel_angle = par[1]
    ellipse_center_x2 = par[2]
    ellipse_center_y2 = par[3]
    ellipse_center_x1 = travel_distance*np.cos(travel_angle*3.14/180.)+ellipse_center_x2
    ellipse_center_y1 = travel_distance*np.sin(travel_angle*3.14/180.)+ellipse_center_y2
    ellipse_semi_major_ax = par[4]
    ellipse_norm = par[5]
    gauss_center_x = 0.5*(ellipse_center_x1+ellipse_center_x2)
    gauss_center_y = 0.5*(ellipse_center_y1+ellipse_center_y2)
    eccentricity = 0.5*pow(pow(ellipse_center_x2-ellipse_center_x1,2)+pow(ellipse_center_y2-ellipse_center_y1,2),0.5)
    if ellipse_semi_major_ax*ellipse_semi_major_ax-eccentricity*eccentricity<0.:
        return 0.
    ellipse_semi_minor_ax = pow(ellipse_semi_major_ax*ellipse_semi_major_ax-eccentricity*eccentricity,0.5)
    theta_ellipse = ROOT.TMath.ATan2(ellipse_center_y1-ellipse_center_y2,ellipse_center_x1-ellipse_center_x2)
    x_rotated = (x[0]-gauss_center_x)*ROOT.TMath.Cos(theta_ellipse) - (x[1]-gauss_center_y)*ROOT.TMath.Sin(theta_ellipse)
    y_rotated = (x[0]-gauss_center_x)*ROOT.TMath.Sin(theta_ellipse) + (x[1]-gauss_center_y)*ROOT.TMath.Cos(theta_ellipse)
    return bkgd_scale*ellipse_norm*exp(-x_rotated*x_rotated/(2.*ellipse_semi_major_ax*ellipse_semi_major_ax))*exp(-y_rotated*y_rotated/(2.*ellipse_semi_minor_ax*ellipse_semi_minor_ax))

def FitSimpleGeometryModel1D_J1908(hist_mask_skymap,hist_data_skymap,hist_syst_skymap):

    global Hist_bkgd_global
    Hist_bkgd_global.Reset()

    global PWN_radius
    global PSR_travel_angle
    global PSR_travel_distance
    global CO_north_center_x
    global CO_north_center_y

    print ('====================================================================')
    print ('FitSimpleGeometryModel1D_J1908')

    nbins = hist_data_skymap.GetNbinsX()
    MapEdge_left = hist_data_skymap.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_data_skymap.GetXaxis().GetBinLowEdge(hist_data_skymap.GetNbinsX()+1)
    MapEdge_lower = hist_data_skymap.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_data_skymap.GetYaxis().GetBinLowEdge(hist_data_skymap.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = hist_data_skymap.GetXaxis().GetBinLowEdge(2)-hist_data_skymap.GetXaxis().GetBinLowEdge(1)

    hist_excess_skymap = hist_data_skymap.Clone()
    hist_excess_skymap.Reset()

    hist_syst_err_skymap = hist_syst_skymap.Clone() # convert syst. content map to syst. error map
    for bx in range(0,hist_syst_skymap.GetNbinsX()):
        for by in range(0,hist_syst_skymap.GetNbinsY()):
            old_content = hist_syst_skymap.GetBinError(bx+1,by+1)
            hist_syst_err_skymap.SetBinContent(bx+1,by+1,0.)
            hist_syst_err_skymap.SetBinError(bx+1,by+1,old_content)

    RA_PSR = 286.975
    Dec_PSR = 6.03777777778

    Hist_bkgd_global.Reset()
    Hist_bkgd_global.Add(hist_mask_skymap)
    hist_excess_skymap.Reset()
    hist_excess_skymap.Add(hist_data_skymap)
    hist_excess_skymap.Add(hist_syst_err_skymap)

    avg_count = hist_excess_skymap.Integral()
    ra_bin = hist_excess_skymap.GetXaxis().FindBin(RA_PSR)
    dec_bin = hist_excess_skymap.GetYaxis().FindBin(Dec_PSR)
    Source_RA = 286.975;
    Source_Dec = 6.269;
    model_PWN = [[Source_RA,0.],[Source_Dec,0.],[PWN_radius,0.],[avg_count,0.]]
    PSR_initial_RA = 287.09
    PSR_initial_Dec = 6.86
    model_PSR_diffusion = [[PSR_travel_angle,0.],[PSR_travel_distance,0.],[0.1,0.],[avg_count,0.], [0.1,0.], [avg_count,0.]]
    model_CO_north = [[CO_north_center_x,0.],[CO_north_center_y,0.],[0.2,0.],[avg_count,0.]]

    npar = 4
    simple_model_2d = ROOT.TF2('simple_model_2d',symmetric_gauss_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    #simple_model_2d = ROOT.TF2('simple_model_2d',disk_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    simple_model_2d.SetParameter(0,model_PWN[0][0])
    simple_model_2d.SetParameter(1,model_PWN[1][0])
    simple_model_2d.SetParameter(2,model_PWN[2][0])
    simple_model_2d.SetParLimits(2,0.1,2.0)
    simple_model_2d.SetParameter(3,model_PWN[3][0])
    simple_model_2d.SetParLimits(3,0.,10.0*(model_PWN[3][0]))
    simple_model_2d.FixParameter(0,model_PWN[0][0])
    simple_model_2d.FixParameter(1,model_PWN[1][0])
    print ('Initial PWN RA = %0.2f'%(model_PWN[0][0]))
    print ('Initial PWN Dec = %0.2f'%(model_PWN[1][0]))
    print ('Initial PWN radius = %0.2f deg'%(model_PWN[2][0]))
    print ('Initial PWN norm = %0.3e'%(model_PWN[3][0]))

    hist_excess_skymap.Fit('simple_model_2d')

    model_PWN[0][0] = simple_model_2d.GetParameter(0)
    model_PWN[1][0] = simple_model_2d.GetParameter(1)
    model_PWN[2][0] = simple_model_2d.GetParameter(2)
    model_PWN[3][0] = simple_model_2d.GetParameter(3)
    print ('2D fit PWN RA = %0.2f'%(model_PWN[0][0]))
    print ('2D fit PWN Dec = %0.2f'%(model_PWN[1][0]))
    print ('2D fit PWN radius = %0.2f deg'%(model_PWN[2][0]))
    print ('2D fit PWN norm = %0.3e'%(model_PWN[3][0]))
    PWN_radius = model_PWN[2][0]

    Hist_bkgd_global.Reset()

    npar = 14
    diffusion_model_2d = ROOT.TF2('diffusion_model_2d',simple_diffusion_model_J1908,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    diffusion_model_2d.SetParameter(0,model_PWN[0][0])
    diffusion_model_2d.SetParameter(1,model_PWN[1][0])
    diffusion_model_2d.SetParameter(2,model_PWN[2][0])
    diffusion_model_2d.SetParameter(3,model_PWN[3][0])
    diffusion_model_2d.FixParameter(0,model_PWN[0][0])
    diffusion_model_2d.FixParameter(1,model_PWN[1][0])
    diffusion_model_2d.FixParameter(2,model_PWN[2][0])
    diffusion_model_2d.FixParameter(3,model_PWN[3][0])
    print ('Initial PWN RA = %0.2f'%(model_PWN[0][0]))
    print ('Initial PWN Dec = %0.2f'%(model_PWN[1][0]))
    print ('Initial PWN radius = %0.2f deg'%(model_PWN[2][0]))
    print ('Initial PWN norm = %0.3e'%(model_PWN[3][0]))

    diffusion_model_2d.SetParameter(4,model_PSR_diffusion[0][0])
    diffusion_model_2d.SetParLimits(4,46.,89.)
    diffusion_model_2d.SetParameter(5,model_PSR_diffusion[1][0])
    diffusion_model_2d.SetParLimits(5,0.1,0.5)
    diffusion_model_2d.SetParameter(6,model_PSR_diffusion[2][0])
    diffusion_model_2d.SetParLimits(6,0.1,1.0)
    diffusion_model_2d.SetParameter(7,model_PSR_diffusion[3][0])
    diffusion_model_2d.SetParLimits(7,0.,10.*model_PSR_diffusion[3][0])
    diffusion_model_2d.SetParameter(8,model_PSR_diffusion[4][0])
    diffusion_model_2d.SetParLimits(8,0.1,1.0)
    diffusion_model_2d.SetParameter(9,model_PSR_diffusion[5][0])
    diffusion_model_2d.SetParLimits(9,0.,10.*model_PSR_diffusion[5][0])
    print ('Initial model_PSR_diffusion travel angle = %0.3e'%(model_PSR_diffusion[0][0]))
    print ('Initial model_PSR_diffusion travel distance = %0.3e'%(model_PSR_diffusion[1][0]))
    print ('Initial model_PSR_diffusion head width = %0.3e'%(model_PSR_diffusion[2][0]))
    print ('Initial model_PSR_diffusion head norm = %0.3e'%(model_PSR_diffusion[3][0]))
    print ('Initial model_PSR_diffusion tail width = %0.3e'%(model_PSR_diffusion[4][0]))
    print ('Initial model_PSR_diffusion tail norm = %0.3e'%(model_PSR_diffusion[5][0]))

    diffusion_model_2d.SetParameter(10,model_CO_north[0][0])
    diffusion_model_2d.SetParLimits(10,model_CO_north[0][0]-0.2,model_CO_north[0][0]+0.2)
    diffusion_model_2d.SetParameter(11,model_CO_north[1][0])
    diffusion_model_2d.SetParLimits(11,model_CO_north[1][0]-0.2,model_CO_north[1][0]+0.2)
    diffusion_model_2d.SetParameter(12,model_CO_north[2][0])
    diffusion_model_2d.SetParLimits(12,0.1,0.5)
    diffusion_model_2d.SetParameter(13,model_CO_north[3][0])
    diffusion_model_2d.SetParLimits(13,0.,10.0*(model_CO_north[3][0]))
    print ('Initial CO north RA = %0.2f'%(model_CO_north[0][0]))
    print ('Initial CO north Dec = %0.2f'%(model_CO_north[1][0]))
    print ('Initial CO north radius = %0.2f deg'%(model_CO_north[2][0]))
    print ('Initial CO north norm = %0.3e'%(model_CO_north[3][0]))


    hist_excess_skymap.Fit('diffusion_model_2d')

    model_PWN[0][0] = diffusion_model_2d.GetParameter(0)
    model_PWN[1][0] = diffusion_model_2d.GetParameter(1)
    model_PWN[2][0] = diffusion_model_2d.GetParameter(2)
    model_PWN[3][0] = diffusion_model_2d.GetParameter(3)
    print ('2D fit PWN RA = %0.2f'%(model_PWN[0][0]))
    print ('2D fit PWN Dec = %0.2f'%(model_PWN[1][0]))
    print ('2D fit PWN radius = %0.2f deg'%(model_PWN[2][0]))
    print ('2D fit PWN norm = %0.3e'%(model_PWN[3][0]))

    model_PSR_diffusion[0][0] = diffusion_model_2d.GetParameter(4)
    model_PSR_diffusion[1][0] = diffusion_model_2d.GetParameter(5)
    model_PSR_diffusion[2][0] = diffusion_model_2d.GetParameter(6)
    model_PSR_diffusion[3][0] = diffusion_model_2d.GetParameter(7)
    model_PSR_diffusion[4][0] = diffusion_model_2d.GetParameter(8)
    model_PSR_diffusion[5][0] = diffusion_model_2d.GetParameter(9)
    print ('2D fit model_PSR_diffusion travel angle = %0.3e'%(model_PSR_diffusion[0][0]))
    print ('2D fit model_PSR_diffusion travel distance = %0.3e'%(model_PSR_diffusion[1][0]))
    print ('2D fit model_PSR_diffusion head width = %0.3e'%(model_PSR_diffusion[2][0]))
    print ('2D fit model_PSR_diffusion head norm = %0.3e'%(model_PSR_diffusion[3][0]))
    print ('2D fit model_PSR_diffusion tail width = %0.3e'%(model_PSR_diffusion[4][0]))
    print ('2D fit model_PSR_diffusion tail norm = %0.3e'%(model_PSR_diffusion[5][0]))
    PSR_travel_angle = model_PSR_diffusion[0][0]
    PSR_travel_distance = model_PSR_diffusion[1][0]

    model_CO_north[0][0] = diffusion_model_2d.GetParameter(10)
    model_CO_north[1][0] = diffusion_model_2d.GetParameter(11)
    model_CO_north[2][0] = diffusion_model_2d.GetParameter(12)
    model_CO_north[3][0] = diffusion_model_2d.GetParameter(13)
    print ('2D fit CO north RA = %0.2f'%(model_CO_north[0][0]))
    print ('2D fit CO north Dec = %0.2f'%(model_CO_north[1][0]))
    print ('2D fit CO north radius = %0.2f deg'%(model_CO_north[2][0]))
    print ('2D fit CO north norm = %0.3e'%(model_CO_north[3][0]))
    CO_north_center_x = model_CO_north[0][0]
    CO_north_center_y = model_CO_north[1][0]


    npar = 4
    model_PWN_func = ROOT.TF2('model_PWN_func',symmetric_gauss_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    #model_PWN_func = ROOT.TF2('model_PWN_func',disk_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    model_PWN_func.SetParameter(0,model_PWN[0][0])
    model_PWN_func.SetParameter(1,model_PWN[1][0])
    model_PWN_func.SetParameter(2,model_PWN[2][0])
    model_PWN_func.SetParameter(3,model_PWN[3][0])

    npar = 6
    model_PSR_diffusion_func = ROOT.TF2('model_PSR_diffusion_func',simple_diffusion_model_PSR,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    model_PSR_diffusion_func.SetParameter(0,model_PSR_diffusion[0][0])
    model_PSR_diffusion_func.SetParameter(1,model_PSR_diffusion[1][0])
    model_PSR_diffusion_func.SetParameter(2,model_PSR_diffusion[2][0])
    model_PSR_diffusion_func.SetParameter(3,model_PSR_diffusion[3][0])
    model_PSR_diffusion_func.SetParameter(4,model_PSR_diffusion[4][0])
    model_PSR_diffusion_func.SetParameter(5,model_PSR_diffusion[5][0])

    npar = 4
    model_CO_north_func = ROOT.TF2('model_CO_north_func',symmetric_gauss_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    model_CO_north_func.SetParameter(0,model_CO_north[0][0])
    model_CO_north_func.SetParameter(1,model_CO_north[1][0])
    model_CO_north_func.SetParameter(2,model_CO_north[2][0])
    model_CO_north_func.SetParameter(3,model_CO_north[3][0])


    Hist_bkgd_global.Reset()
    hist_excess_skymap.Reset()
    hist_excess_skymap.Add(hist_data_skymap)
    hist_excess_skymap.Add(hist_syst_err_skymap)

    hist_PWN_skymap = hist_data_skymap.Clone()
    hist_PWN_skymap.Reset()
    hist_PSR_ellipse_skymap = hist_data_skymap.Clone()
    hist_PSR_ellipse_skymap.Reset()
    hist_CO_north_skymap = hist_data_skymap.Clone()
    hist_CO_north_skymap.Reset()
    skymap_bin_size_x = hist_data_skymap.GetXaxis().GetBinCenter(2)-hist_data_skymap.GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = hist_data_skymap.GetYaxis().GetBinCenter(2)-hist_data_skymap.GetYaxis().GetBinCenter(1)
    for bx in range(0,hist_data_skymap.GetNbinsX()):
        for by in range(0,hist_data_skymap.GetNbinsY()):
            map_x = hist_data_skymap.GetXaxis().GetBinCenter(bx+1)
            map_y = hist_data_skymap.GetYaxis().GetBinCenter(by+1)
            PWN_flux = model_PWN_func.Eval(map_x,map_y)
            hist_PWN_skymap.SetBinContent(bx+1,by+1,PWN_flux)
            PSR_ellipse_flux = model_PSR_diffusion_func.Eval(map_x,map_y)
            hist_PSR_ellipse_skymap.SetBinContent(bx+1,by+1,PSR_ellipse_flux)
            CO_north_flux = model_CO_north_func.Eval(map_x,map_y)
            hist_CO_north_skymap.SetBinContent(bx+1,by+1,CO_north_flux)
    print ('hist_PWN_skymap.Integral() = %0.2e'%(hist_PWN_skymap.Integral()))

    return model_PWN, hist_PWN_skymap, model_PSR_diffusion, hist_PSR_ellipse_skymap, model_CO_north, hist_CO_north_skymap

fig, ax = plt.subplots()
energy_bin = CommonPlotFunctions.energy_bin

energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])
print ('energy_bin_cut_low = %s'%(energy_bin_cut_low))
print ('energy_bin_cut_up = %s'%(energy_bin_cut_up))

PWN_radius = 1.5
PSR_travel_angle = 6.247e+01
PSR_travel_distance = 0.2
CO_north_center_x = 286.8
CO_north_center_y = 7.1

InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_data_skymap_sum"
nbins = InputFile.Get(HistName).GetNbinsX()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()

Hist_bkgd_global = ROOT.TH2D("Hist_bkgd_global","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_syst_skymap_sum = ROOT.TH2D("hist_flux_syst_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_data_skymap_sum = ROOT.TH2D("hist_data_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_bkgd_skymap_sum = ROOT.TH2D("hist_bkgd_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_expo_skymap_sum = ROOT.TH2D("hist_expo_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_mask_skymap_sum = ROOT.TH2D("hist_mask_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_syst_skymap_sum = ROOT.TH2D("hist_syst_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap = []
hist_flux_syst_skymap = []
hist_cali_skymap = []
hist_data_skymap = []
hist_bkgd_skymap = []
hist_expo_skymap = []
hist_syst_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_flux_skymap += [ROOT.TH2D("hist_flux_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_flux_syst_skymap += [ROOT.TH2D("hist_flux_syst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_cali_skymap += [ROOT.TH2D("hist_cali_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_data_skymap += [ROOT.TH2D("hist_data_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_bkgd_skymap += [ROOT.TH2D("hist_bkgd_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_expo_skymap += [ROOT.TH2D("hist_expo_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_syst_skymap += [ROOT.TH2D("hist_syst_skymap_v2_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_mask_skymap_sum"
hist_mask_skymap_sum.Add(InputFile.Get(HistName))
for ebin in range(0,len(energy_bin)-1):
    HistName = "hist_energy_flux_skymap_%s"%(ebin)
    hist_flux_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_energy_flux_syst_skymap_%s"%(ebin)
    hist_flux_syst_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_cali_skymap_%s"%(ebin)
    hist_cali_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_data_skymap_%s"%(ebin)
    hist_data_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    hist_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_expo_skymap_%s"%(ebin)
    hist_expo_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_syst_skymap_%s"%(ebin)
    hist_syst_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_flux_skymap_sum.Add(hist_flux_skymap[ebin])
    hist_flux_syst_skymap_sum.Add(hist_flux_syst_skymap[ebin])
    hist_data_skymap_sum.Add(hist_data_skymap[ebin])
    hist_bkgd_skymap_sum.Add(hist_bkgd_skymap[ebin])
    hist_syst_skymap_sum.Add(hist_syst_skymap[ebin])
    hist_expo_skymap_sum.Add(hist_expo_skymap[ebin])

hist_fit_PWN_skymap_sum = ROOT.TH2D("hist_fit_PWN_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_CO_north_skymap_sum = ROOT.TH2D("hist_fit_CO_north_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_PSR_ellipse_skymap_sum = ROOT.TH2D("hist_fit_PSR_ellipse_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_all_models_skymap_sum = ROOT.TH2D("hist_fit_all_models_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

model_PWN, hist_fit_PWN_skymap, model_PSR_ellipse, hist_fit_PSR_ellipse_skymap, model_CO_north, hist_fit_CO_north_skymap = FitSimpleGeometryModel1D_J1908(hist_mask_skymap_sum,hist_flux_skymap_sum,hist_flux_syst_skymap_sum)
hist_fit_PWN_skymap_sum.Add(hist_fit_PWN_skymap)
hist_fit_PSR_ellipse_skymap_sum.Add(hist_fit_PSR_ellipse_skymap)
hist_fit_CO_north_skymap_sum.Add(hist_fit_CO_north_skymap)
hist_fit_all_models_skymap_sum.Add(hist_fit_PWN_skymap_sum)
hist_fit_all_models_skymap_sum.Add(hist_fit_PSR_ellipse_skymap_sum)
hist_fit_all_models_skymap_sum.Add(hist_fit_CO_north_skymap_sum)
mycircle_PWN = ROOT.TEllipse(-1.*model_PWN[0][0],model_PWN[1][0],model_PWN[2][0])
mycircle_PWN.SetFillStyle(0)
mycircle_PWN.SetLineColor(2)
travel_angle = model_PSR_ellipse[0][0]
travel_interval = model_PSR_ellipse[1][0]
RA_PSR = 286.975
Dec_PSR = 6.03777777778
PSR_head_x = 0.*travel_interval*np.cos(travel_angle*3.14/180.)+RA_PSR
PSR_head_y = 0.*travel_interval*np.sin(travel_angle*3.14/180.)+Dec_PSR
PSR_tail_x = 1.*travel_interval*np.cos(travel_angle*3.14/180.)+PSR_head_x
PSR_tail_y = 1.*travel_interval*np.sin(travel_angle*3.14/180.)+PSR_head_y
print ('travel_angle = %0.2f deg'%(travel_angle))
print ('PSR_tail_x = %0.2f deg'%(PSR_tail_x))
print ('PSR_tail_y = %0.2f deg'%(PSR_tail_y))

mycircle_PSR_head = ROOT.TEllipse(-1.*PSR_head_x,PSR_head_y,model_PSR_ellipse[2][0])
mycircle_PSR_head.SetFillStyle(0)
mycircle_PSR_head.SetLineColor(2)
mycircle_PSR_tail = ROOT.TEllipse(-1.*PSR_tail_x,PSR_tail_y,model_PSR_ellipse[4][0])
mycircle_PSR_tail.SetFillStyle(0)
mycircle_PSR_tail.SetLineColor(2)

print ('hist_fit_all_models_skymap_sum.Integral() = %0.2e'%(hist_fit_all_models_skymap_sum.Integral()))
print ('hist_fit_PWN_skymap_sum.Integral() = %0.2e'%(hist_fit_PWN_skymap_sum.Integral()))
print ('hist_fit_PSR_ellipse_skymap_sum.Integral() = %0.2e'%(hist_fit_PSR_ellipse_skymap_sum.Integral()))

output_file = ROOT.TFile("output_fitting/J1908_fit_skymap.root","update")
mytree = ROOT.TTree("mytree_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up),"mytree_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up))
PWN_radius_br = array('f',[PWN_radius])
PSR_travel_angle_br = array('f',[PSR_travel_angle])
mytree.Branch("PWN_radius",PWN_radius_br,"PWN_radius/F")
mytree.Branch("PSR_travel_angle",PSR_travel_angle_br,"PSR_travel_angle/F")
mytree.Fill()
mytree.Write()
hist_data_skymap_sum.Write()
hist_bkgd_skymap_sum.Write()
hist_syst_skymap_sum.Write()
hist_flux_skymap_sum.Write()
hist_flux_syst_skymap_sum.Write()
hist_fit_PWN_skymap_sum.Write()
hist_fit_PSR_ellipse_skymap_sum.Write()
hist_fit_CO_north_skymap_sum.Write()
output_file.Close();

other_stars, other_star_coord = CommonPlotFunctions.GetGammaSourceInfo() 
other_star_labels = []
other_star_markers = []
other_star_names = []
star_range = 3.0
source_ra = (MapEdge_left+MapEdge_right)/2.
source_dec = (MapEdge_lower+MapEdge_upper)/2.
for star in range(0,len(other_stars)):
    if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
    other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
    other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
    other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
    other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
    other_star_names += [other_stars[star]]


canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 650, 600)
pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
pad1.SetBottomMargin(0.10)
pad1.SetRightMargin(0.20)
pad1.SetLeftMargin(0.10)
pad1.SetTopMargin(0.10)
pad1.SetBorderMode(0)
pad1.Draw()
pad1.cd()

hist_zscore = CommonPlotFunctions.GetSignificanceMap(hist_data_skymap_sum,hist_bkgd_skymap_sum,hist_syst_skymap_sum,False)
hist_zscore_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore)
hist_zscore_reflect.SetContour(3)
hist_zscore_reflect.SetContourLevel(0,3)
hist_zscore_reflect.SetContourLevel(1,4)
hist_zscore_reflect.SetContourLevel(2,5)

hist_mask_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_mask_skymap_sum)
hist_mask_skymap_reflect.SetContour(1)
hist_mask_skymap_reflect.SetContourLevel(0,0.5)
hist_mask_skymap_reflect.SetLineColor(2)

hist_excess_skymap_sum = ROOT.TH2D("hist_excess_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_excess_skymap_sum.Reset()
hist_excess_skymap_sum.Add(hist_data_skymap_sum)
hist_excess_skymap_sum.Add(hist_bkgd_skymap_sum,-1.)
hist_excess_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_excess_skymap_sum)
x1 = hist_excess_skymap_sum_reflect.GetXaxis().GetXmin()
x2 = hist_excess_skymap_sum_reflect.GetXaxis().GetXmax()
y1 = hist_excess_skymap_sum_reflect.GetYaxis().GetXmin()
y2 = hist_excess_skymap_sum_reflect.GetYaxis().GetXmax()
IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 505, "+")
raLowerAxis.SetLabelSize(hist_excess_skymap_sum_reflect.GetXaxis().GetLabelSize())
raLowerAxis.Draw()
hist_excess_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_excess_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_excess_skymap_sum_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
hist_mask_skymap_reflect.Draw("CONT3 same")
raLowerAxis.Draw()
mycircle_PWN.Draw("same")
mycircle_PSR_head.Draw("same")
mycircle_PSR_tail.Draw("same")
for star in range(0,len(other_star_markers)):
    other_star_markers[star].Draw("same")
    other_star_labels[star].Draw("same")
canvas.SaveAs('output_plots/SkymapDataExcess_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_expo_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_expo_skymap_sum)
hist_expo_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_expo_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_expo_skymap_sum_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
hist_mask_skymap_reflect.Draw("CONT3 same")
raLowerAxis.Draw()
mycircle_PWN.Draw("same")
mycircle_PSR_head.Draw("same")
mycircle_PSR_tail.Draw("same")
canvas.SaveAs('output_plots/SkymapExposure_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_flux_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_flux_skymap_sum)
hist_flux_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_flux_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_flux_skymap_sum_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
hist_mask_skymap_reflect.Draw("CONT3 same")
raLowerAxis.Draw()
mycircle_PWN.Draw("same")
mycircle_PSR_head.Draw("same")
mycircle_PSR_tail.Draw("same")
for star in range(0,len(other_star_markers)):
    other_star_markers[star].Draw("same")
    other_star_labels[star].Draw("same")
canvas.SaveAs('output_plots/SkymapDataFlux_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_fit_all_models_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_fit_all_models_skymap_sum)
hist_fit_all_models_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_fit_all_models_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_fit_all_models_skymap_sum_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
hist_mask_skymap_reflect.Draw("CONT3 same")
raLowerAxis.Draw()
mycircle_PWN.Draw("same")
mycircle_PSR_head.Draw("same")
mycircle_PSR_tail.Draw("same")
canvas.SaveAs('output_plots/SkymapModel_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_fit_error_skymap = ROOT.TH2D("hist_fit_error_skymap","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_error_skymap.Add(hist_flux_skymap_sum)
hist_fit_error_skymap.Add(hist_fit_all_models_skymap_sum,-1.)
hist_fit_error_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_fit_error_skymap)
hist_fit_error_skymap_reflect.GetXaxis().SetLabelOffset(999)
hist_fit_error_skymap_reflect.GetXaxis().SetTickLength(0)
hist_fit_error_skymap_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
hist_mask_skymap_reflect.Draw("CONT3 same")
raLowerAxis.Draw()
mycircle_PWN.Draw("same")
mycircle_PSR_head.Draw("same")
mycircle_PSR_tail.Draw("same")
canvas.SaveAs('output_plots/SkymapFitError_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))


RA_PSR = 286.975
Dec_PSR = 6.03777777778
profile_center_x, profile_center_y, profile_center_z = RA_PSR, Dec_PSR, 0.
profile, profile_err, theta2 = CommonPlotFunctions.FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,2.0)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindExtension(hist_fit_all_models_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindExtension(hist_fit_PWN_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindExtension(hist_fit_CO_north_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindExtension(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x,profile_center_y,2.0)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='g',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('angular distance from source [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsTheta2_PSR'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()


profile_center_x = 286.77
profile_center_y = 7.11
profile, profile_err, theta2 = CommonPlotFunctions.FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,2.0)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindExtension(hist_fit_all_models_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindExtension(hist_fit_PWN_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindExtension(hist_fit_CO_north_skymap_sum,None,profile_center_x,profile_center_y,2.0)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindExtension(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x,profile_center_y,2.0)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='b',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('angular distance from source [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsTheta2_CO_north'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()

PSR_initial_RA = 287.09
PSR_initial_Dec = 6.86
profile_center_x1 = 286.98
profile_center_y1 = 6.04
#profile_center_x2 = 287.2
#profile_center_y2 = 6.6
profile_center_x2 = 287.09
profile_center_y2 = 6.86
profile, profile_err, theta2 = CommonPlotFunctions.FindProjection(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindProjection(hist_fit_all_models_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindProjection(hist_fit_PWN_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindProjection(hist_fit_CO_north_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindProjection(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='g',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('distance along proper motion [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsProjectedX'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()


profile, profile_err, theta2 = CommonPlotFunctions.FindProjection(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindProjection(hist_fit_all_models_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindProjection(hist_fit_PWN_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindProjection(hist_fit_CO_north_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindProjection(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='g',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('distance perpendicular to proper motion [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsProjectedY0'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()

profile, profile_err, theta2 = CommonPlotFunctions.FindProjection(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindProjection(hist_fit_all_models_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindProjection(hist_fit_PWN_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindProjection(hist_fit_CO_north_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindProjection(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='g',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('distance perpendicular to proper motion [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsProjectedY1'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()

profile, profile_err, theta2 = CommonPlotFunctions.FindProjection(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,1.0)
profile_all, profile_err_all, theta2_all = CommonPlotFunctions.FindProjection(hist_fit_all_models_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,1.0)
profile_PWN, profile_err_PWN, theta2_PWN = CommonPlotFunctions.FindProjection(hist_fit_PWN_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,1.0)
profile_CO_north, profile_err_CO_north, theta2_CO_north = CommonPlotFunctions.FindProjection(hist_fit_CO_north_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,1.0)
profile_PSR_ellipse, profile_err_PSR_ellipse, theta2_PSR_ellipse = CommonPlotFunctions.FindProjection(hist_fit_PSR_ellipse_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,1.0)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2,profile,profile_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
axbig.plot(theta2_all,profile_all,color='k',label='all models')
axbig.plot(theta2_PWN,profile_PWN,color='b',label='Halo Region model')
axbig.plot(theta2_CO_north,profile_CO_north,color='g',label='North Region model')
axbig.plot(theta2_PSR_ellipse,profile_PSR_ellipse,color='r',label='PSR Region model')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
axbig.set_xlabel('distance perpendicular to proper motion [degree]')
axbig.legend(loc='best')
plotname = 'ProfileVsProjectedY2'
fig.savefig("output_plots/%s_%s_E%sto%s.png"%(plotname,sys.argv[1],energy_bin_cut_low,energy_bin_cut_up),bbox_inches='tight')
axbig.remove()
