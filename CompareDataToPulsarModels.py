
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

energy_index = CommonPlotFunctions.energy_index_scale
Smoothing = CommonPlotFunctions.Smoothing

def pulsar_model(x,par):

    binx = hist_model_global.GetXaxis().FindBin(x[0])
    biny = hist_model_global.GetYaxis().FindBin(x[1])
    model = hist_model_global.GetBinContent(binx,biny)

    norm = par[0]
    return norm*model

def FitPulsarModel(hist_flux,hist_flux_syst,hist_pwn,hist_co_north,hist_psr,hist_model):

    global hist_model_global
    hist_model_global.Reset()
    hist_model_global.Add(hist_model)

    hist_data = ROOT.TH2D("hist_data","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    hist_data.Add(hist_flux)
    for bx in range(0,hist_data.GetNbinsX()):
        for by in range(0,hist_data.GetNbinsY()):
            syst = hist_flux_syst.GetBinContent(bx+1,by+1)
            hist_data.SetBinError(bx+1,by+1,syst)
    hist_data.Add(hist_pwn,-1.)
    hist_data.Add(hist_co_north,-1.)

    scaling_factor = hist_psr.Integral()/hist_model.Integral()

    npar = 1
    simple_model_2d = ROOT.TF2('simple_model_2d',pulsar_model,MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper,npar)
    simple_model_2d.SetParameter(0,scaling_factor)
    simple_model_2d.SetParLimits(0,0.,100.0*scaling_factor)
    hist_data.Fit('simple_model_2d')
    scale = simple_model_2d.GetParameter(0)

    return scale

def GetPulsarModelChi2(hist_flux,hist_flux_syst,hist_pwn,hist_co_north,hist_psr):

    hist_model = ROOT.TH2D("hist_model","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

    #D_par = [0,1,2,3,4,5,6,7,8,9]
    #V_par = [0,1,2,3,4,5,6,7,8]
    D_par = [0,1,2,3,4,5]
    V_par = [0,1,2,3,4,5]
    chi2_array = np.full((len(D_par),len(V_par)),0.)
    D_axis = [ [0 for x in range(len(D_par))] for y in range(0,2)]
    V_axis = [ [0 for x in range(len(V_par))] for y in range(0,2)]

    for idx_d in range(0,len(D_par)):
        for idx_v in range(0,len(V_par)):

            print ('Fitting model D%sV%s'%(D_par[idx_d],V_par[idx_v]))
            hist_model = ReadPulsarModelMap(hist_model,D_par[idx_d],V_par[idx_v])
            InputFile = ROOT.TFile("output_pulsar_models_3kpc/pulsar_skymap_D%s_V%s.root"%(D_par[idx_d],V_par[idx_v]))
            InfoTree = InputFile.Get("InfoTree")
            InfoTree.GetEntry(0)
            proper_velocity = InfoTree.proper_velocity/1000.
            D0 = InfoTree.D0/1e26
            D_axis[0][idx_d] = '%0.2f'%(D0)
            V_axis[0][idx_v] = '%0.2f'%(proper_velocity)
            D_axis[1][idx_d] = D_par[idx_d]
            V_axis[1][idx_v] = V_par[idx_v]
            InputFile.Close()
            
            if Smoothing:
                hist_model = CommonPlotFunctions.Smooth2DMap(hist_model,CommonPlotFunctions.smooth_size_spectroscopy,False,True)
            print ('hist_flux.Integral() = %s'%(hist_flux.Integral()))
            print ('hist_co_north.Integral() = %s'%(hist_co_north.Integral()))
            scaling_factor = FitPulsarModel(hist_flux,hist_flux_syst,hist_pwn,hist_co_north,hist_psr,hist_model)
            hist_model.Scale(scaling_factor)
            total_chi2 = 0.
            total_bins = 0.
            for binx in range(0,hist_flux.GetNbinsX()):
                for biny in range(0,hist_flux.GetNbinsY()):
                    data_flux = hist_flux.GetBinContent(binx+1,biny+1)
                    data_flux_err = hist_flux.GetBinError(binx+1,biny+1)
                    data_flux_syst = hist_flux_syst.GetBinContent(binx+1,biny+1)
                    pwn_flux = hist_pwn.GetBinContent(binx+1,biny+1)
                    co_north_flux = hist_co_north.GetBinContent(binx+1,biny+1)
                    model_flux = hist_model.GetBinContent(binx+1,biny+1)
                    error_sq = (data_flux_err*data_flux_err+data_flux_syst*data_flux_syst)
                    cell_x = hist_flux.GetXaxis().GetBinCenter(binx+1)
                    cell_y = hist_flux.GetYaxis().GetBinCenter(biny+1)
                    distance_to_pulsar = pow(pow(cell_x-RA_PSR,2)+pow(cell_y-Dec_PSR,2),0.5)
                    mask = 1.
                    if distance_to_pulsar>1.:
                        mask = 0.
                    if error_sq>0.:
                        if mask==0.: continue
                        chi2 = pow(data_flux-pwn_flux-co_north_flux-model_flux,2)/error_sq
                        total_chi2 += chi2
                        total_bins += 1.
            chi2_array[idx_d,idx_v] = pow(total_chi2/total_bins,0.5)

    return chi2_array, D_axis, V_axis

def ReadPulsarModelMap(hist_model_skymap_output,idx_d_best,idx_v_best):

    InputFile = ROOT.TFile("output_pulsar_models_3kpc/pulsar_skymap_D%s_V%s.root"%(idx_d_best,idx_v_best))

    HistName = "hist_pulsar_skymap_631_1585"
    nbins = InputFile.Get(HistName).GetNbinsX()
    MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
    MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
    MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
    hist_model_skymap = ROOT.TH2D("hist_model_skymap","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        HistName = "hist_pulsar_skymap_631_1585"
        if ebin==2:
            HistName = "hist_pulsar_skymap_631_1585"
        elif ebin==3:
            HistName = "hist_pulsar_skymap_1585_3981"
        elif ebin==4:
            HistName = "hist_pulsar_skymap_3981_10000"
        elif ebin==5:
            HistName = "hist_pulsar_skymap_10000_25118"
        hist_model_skymap.Add(InputFile.Get(HistName),pow(energy_bin[ebin]/1e3,energy_index-1))
    print ('hist_model_skymap.Integral() = %s'%(hist_model_skymap.Integral()))

    hist_model_skymap_output.Reset()
    for binx in range(0,hist_model_skymap.GetNbinsX()):
        for biny in range(0,hist_model_skymap.GetNbinsY()):
            cell_x = hist_model_skymap.GetXaxis().GetBinCenter(binx+1)
            cell_y = hist_model_skymap.GetYaxis().GetBinCenter(biny+1)
            cell_z = hist_model_skymap.GetBinContent(binx+1,biny+1)
            binx_out = hist_model_skymap_output.GetXaxis().FindBin(cell_x)
            biny_out = hist_model_skymap_output.GetYaxis().FindBin(cell_y)
            old_content = hist_model_skymap_output.GetBinContent(binx_out,biny_out)
            hist_model_skymap_output.SetBinContent(binx_out,biny_out,cell_z+old_content)

    InputFile.Close()

    return hist_model_skymap_output


fig, ax = plt.subplots()

energy_bin = CommonPlotFunctions.energy_bin

RA_PSR = 286.975
Dec_PSR = 6.03777777778

InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_data_skymap_sum"
nbins = InputFile.Get(HistName).GetNbinsX()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()
 
hist_model_global = ROOT.TH2D("hist_model_global","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_excess_skymap_sum = ROOT.TH2D("hist_excess_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_data_skymap_sum = ROOT.TH2D("hist_data_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_bkgd_skymap_sum = ROOT.TH2D("hist_bkgd_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_syst_skymap_sum = ROOT.TH2D("hist_syst_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_syst_skymap_sum = ROOT.TH2D("hist_flux_syst_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_pwn_skymap_sum = ROOT.TH2D("hist_pwn_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_psr_skymap_sum = ROOT.TH2D("hist_psr_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_co_north_skymap_sum = ROOT.TH2D("hist_co_north_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])

InputFile = ROOT.TFile("output_fitting/J1908_fit_skymap.root")
HistName = "hist_data_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_data_skymap_sum.Add(InputFile.Get(HistName))
hist_excess_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_bkgd_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_bkgd_skymap_sum.Add(InputFile.Get(HistName))
hist_excess_skymap_sum.Add(InputFile.Get(HistName),-1.)
HistName = "hist_syst_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_syst_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_flux_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_flux_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_flux_syst_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_flux_syst_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_fit_PWN_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_pwn_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_fit_PSR_ellipse_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_psr_skymap_sum.Add(InputFile.Get(HistName))
HistName = "hist_fit_CO_north_skymap_sum_E%sto%s"%(energy_bin_cut_low,energy_bin_cut_up)
hist_co_north_skymap_sum.Add(InputFile.Get(HistName))
InputFile.Close()

hist_psr_skymap_sum.Reset()
hist_co_north_skymap_sum.Reset()


hist_model_skymap_sum = ROOT.TH2D("hist_model_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)

find_chi2 = True

if find_chi2:
    chi2_array_pub, D0_axis, Vel_axis = GetPulsarModelChi2(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,hist_pwn_skymap_sum,hist_co_north_skymap_sum,hist_psr_skymap_sum)
    plt.imshow(chi2_array_pub[:,:])
    plt.xticks(Vel_axis[1], Vel_axis[0])
    plt.yticks(D0_axis[1], D0_axis[0])
    plt.xlabel("proper motion 1e3 km/s")
    plt.ylabel("diffusion coefficient 1e26 $\mathrm{cm}^{2}\mathrm{s}^{-1}$")
    plt.colorbar()
    plotname = 'PulsarModelChi2_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    plt.clf()


chi2_min = 1e10
idx_d_best = 2
idx_v_best = 4
if find_chi2:
    for idx_d in range(0,len(D0_axis[1])):
        for idx_v in range(0,len(Vel_axis[1])):
            chi2 = chi2_array_pub[idx_d,idx_v]
            if chi2_min>chi2:
                chi2_min = chi2
                idx_d_best = idx_d
                idx_v_best = idx_v
print ('chi2_min = %s'%(chi2_min))
print ('idx_d_best = %s'%(idx_d_best))
print ('idx_v_best = %s'%(idx_v_best))
hist_model_skymap_sum = ReadPulsarModelMap(hist_model_skymap_sum,idx_d_best,idx_v_best)

if Smoothing:
    hist_model_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_model_skymap_sum,CommonPlotFunctions.smooth_size_spectroscopy,False,True)
scaling_factor = FitPulsarModel(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,hist_pwn_skymap_sum,hist_co_north_skymap_sum,hist_psr_skymap_sum,hist_model_skymap_sum)
hist_model_skymap_sum.Scale(scaling_factor)
print ('scaling_factor = %s'%(scaling_factor))
print ('Done fitting.')

profile_center_x = RA_PSR
profile_center_y = Dec_PSR
profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,2.)
profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindExtension(hist_pwn_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_co_north, profile_co_north_err, theta2_co_north = CommonPlotFunctions.FindExtension(hist_co_north_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindExtension(hist_model_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindExtension(hist_psr_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_total = np.array(profile_pwn)+np.array(profile_co_north)+np.array(profile_model)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='data')
axbig.plot(theta2_pwn,profile_pwn,color='r',label='halo model')
axbig.plot(theta2_co_north,profile_co_north,color='g',label='north model')
axbig.plot(theta2_model,profile_model,color='b',label='pulsar model')
#axbig.plot(theta2_fit,profile_fit,color='g')
axbig.plot(theta2_pwn,profile_total,color='k',label='sum of models')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
axbig.set_xlabel('angular distance from source [degree]')
axbig.legend(loc='best')
plotname = 'ProfilePSRCenter_Sum_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
print ('Done plotting %s.'%(plotname))

profile_center_x = 286.8
profile_center_y = 7.1
profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,2.)
profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindExtension(hist_pwn_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_co_north, profile_co_north_err, theta2_co_north = CommonPlotFunctions.FindExtension(hist_co_north_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindExtension(hist_model_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindExtension(hist_psr_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_total = np.array(profile_pwn)+np.array(profile_co_north)+np.array(profile_model)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='data')
axbig.plot(theta2_pwn,profile_pwn,color='r',label='halo model')
axbig.plot(theta2_co_north,profile_co_north,color='g',label='north model')
axbig.plot(theta2_model,profile_model,color='b',label='pulsar model')
#axbig.plot(theta2_fit,profile_fit,color='g')
axbig.plot(theta2_pwn,profile_total,color='k',label='sum of models')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
axbig.set_xlabel('angular distance from source [degree]')
axbig.legend(loc='best')
plotname = 'ProfileNorthCenter_Sum_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
print ('Done plotting %s.'%(plotname))

travel_angle = 78.02
travel_distance = 1.0
profile_center_x1 = 286.98
profile_center_y1 = 6.04
profile_center_x2 = travel_distance*np.cos(travel_angle*3.14/180.)+profile_center_x1
profile_center_y2 = travel_distance*np.sin(travel_angle*3.14/180.)+profile_center_y1
profile_center_z = 0.

profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindProjection(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindProjection(hist_pwn_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_co_north, profile_co_north_err, theta2_co_north = CommonPlotFunctions.FindProjection(hist_co_north_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindProjection(hist_model_skymap_sum,None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
profile_total = np.array(profile_pwn)+np.array(profile_co_north)+np.array(profile_model)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='data')
axbig.plot(theta2_pwn,profile_pwn,color='r',label='halo model')
axbig.plot(theta2_co_north,profile_co_north,color='g',label='north model')
axbig.plot(theta2_model,profile_model,color='b',label='pulsar model')
axbig.plot(theta2_pwn,profile_total,color='k',label='sum of models')
axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
axbig.set_xlabel('distance along proper motion [degree]')
axbig.legend(loc='best')
plotname = 'ProfileModelProjectedX_Sum_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
print ('Done plotting %s.'%(plotname))


canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 650, 600)
pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
pad1.SetBottomMargin(0.10)
pad1.SetRightMargin(0.20)
pad1.SetLeftMargin(0.10)
pad1.SetTopMargin(0.10)
pad1.SetBorderMode(0)
pad1.Draw()
pad1.cd()

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

hist_zscore = CommonPlotFunctions.GetSignificanceMap(hist_data_skymap_sum,hist_bkgd_skymap_sum,hist_syst_skymap_sum,False)
hist_zscore_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore)
hist_zscore_reflect.SetContour(3)
hist_zscore_reflect.SetContourLevel(0,3)
hist_zscore_reflect.SetContourLevel(1,4)
hist_zscore_reflect.SetContourLevel(2,5)

hist_excess_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_excess_skymap_sum)
CommonPlotFunctions.MatplotlibMap2D(hist_excess_skymap_sum_reflect,hist_zscore_reflect,fig,'RA','Dec','excess','SkymapExcess_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_flux_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_flux_skymap_sum)
x1 = hist_flux_skymap_sum_reflect.GetXaxis().GetXmin()
x2 = hist_flux_skymap_sum_reflect.GetXaxis().GetXmax()
y1 = hist_flux_skymap_sum_reflect.GetYaxis().GetXmin()
y2 = hist_flux_skymap_sum_reflect.GetYaxis().GetXmax()
IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 505, "+")
raLowerAxis.SetLabelSize(hist_flux_skymap_sum_reflect.GetXaxis().GetLabelSize())
raLowerAxis.Draw()
hist_flux_skymap_sum_reflect.GetXaxis().SetLabelOffset(999)
hist_flux_skymap_sum_reflect.GetXaxis().SetTickLength(0)
hist_flux_skymap_sum_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
for star in range(0,len(other_star_markers)):
    other_star_markers[star].Draw("same")
    other_star_labels[star].Draw("same")
raLowerAxis.Draw()
canvas.SaveAs('output_plots/SkymapDataFlux_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_fit_all_models_skymap_sum = ROOT.TH2D("hist_fit_all_models_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_all_models_skymap_sum.Add(hist_model_skymap_sum)
hist_fit_all_models_skymap_sum.Add(hist_pwn_skymap_sum)
hist_fit_all_models_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_fit_all_models_skymap_sum)
y_label = '$E^{%s}$ Flux [$\mathrm{TeV}^{%s}\mathrm{cm}^{-2}\mathrm{s}^{-1}$]'%(energy_index,-1+energy_index)
CommonPlotFunctions.MatplotlibMap2D(hist_fit_all_models_skymap_sum_reflect,hist_zscore_reflect,fig,'RA','Dec',y_label,'SkymapModel_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))

hist_fit_error_skymap = ROOT.TH2D("hist_fit_error_skymap","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_fit_error_skymap.Add(hist_flux_skymap_sum)
hist_fit_error_skymap.Add(hist_fit_all_models_skymap_sum,-1.)
hist_fit_error_skymap.Divide(hist_flux_syst_skymap_sum)
hist_fit_error_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_fit_error_skymap)
hist_fit_error_skymap_reflect.GetXaxis().SetLabelOffset(999)
hist_fit_error_skymap_reflect.GetXaxis().SetTickLength(0)
hist_fit_error_skymap_reflect.Draw("COL4Z")
hist_zscore_reflect.Draw("CONT3 same")
for star in range(0,len(other_star_markers)):
    other_star_markers[star].Draw("same")
    other_star_labels[star].Draw("same")
raLowerAxis.Draw()
canvas.SaveAs('output_plots/SkymapFitError_E%sto%s.png'%(energy_bin_cut_low,energy_bin_cut_up))


