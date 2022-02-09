
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


InputFile = ROOT.TFile("output_fitting/J1908_skymap.root")
HistName = "hist_data_skymap_sum"
nbins = InputFile.Get(HistName).GetNbinsX()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()
 
hist_flux_skymap_sum = ROOT.TH2D("hist_flux_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_syst_skymap_sum = ROOT.TH2D("hist_flux_syst_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_pwn_skymap_sum = ROOT.TH2D("hist_pwn_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_psr_skymap_sum = ROOT.TH2D("hist_psr_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_flux_skymap = []
hist_flux_syst_skymap = []
hist_pwn_skymap = []
hist_psr_skymap = []
for ebin in range(0,len(CommonPlotFunctions.energy_bin)-1):
    hist_flux_skymap += [ROOT.TH2D("hist_flux_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_flux_syst_skymap += [ROOT.TH2D("hist_flux_syst_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_pwn_skymap += [ROOT.TH2D("hist_pwn_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_psr_skymap += [ROOT.TH2D("hist_psr_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])

InputFile = ROOT.TFile("output_fitting/J1908_fit_skymap.root")
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    HistName = "hist_flux_skymap_v2_E%s"%(ebin)
    hist_flux_skymap_sum.Add(InputFile.Get(HistName))
    hist_flux_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_flux_syst_skymap_v2_E%s"%(ebin)
    hist_flux_syst_skymap_sum.Add(InputFile.Get(HistName))
    hist_flux_syst_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_fit_PWN_skymap_E%s"%(ebin)
    hist_pwn_skymap_sum.Add(InputFile.Get(HistName))
    hist_pwn_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_fit_PSR_ellipse_skymap_E%s"%(ebin)
    hist_psr_skymap_sum.Add(InputFile.Get(HistName))
    hist_psr_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

hist_model_skymap_sum = ROOT.TH2D("hist_model_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_model_skymap = []
for ebin in range(0,len(CommonPlotFunctions.energy_bin)-1):
    hist_model_skymap += [ROOT.TH2D("hist_model_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

energy_index = 2
InputFile = ROOT.TFile("output_pulsar_models/pulsar_skymap_D3_V4.root")
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
    hist_model_skymap_sum.Add(InputFile.Get(HistName),pow(CommonPlotFunctions.energy_bin[ebin]/1e3,energy_index-1))
    hist_model_skymap[ebin].Add(InputFile.Get(HistName),pow(CommonPlotFunctions.energy_bin[ebin]/1e3,energy_index-1))
InputFile.Close()

hist_model_skymap_sum = CommonPlotFunctions.Smooth2DMap(hist_model_skymap_sum,CommonPlotFunctions.smooth_size_spectroscopy,False,True)
scaling_factor = hist_psr_skymap_sum.Integral()/hist_model_skymap_sum.Integral()
hist_model_skymap_sum.Scale(scaling_factor)
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_model_skymap[ebin] = CommonPlotFunctions.Smooth2DMap(hist_model_skymap[ebin],CommonPlotFunctions.smooth_size_spectroscopy,False,True)
    scaling_factor = hist_psr_skymap[ebin].Integral()/hist_model_skymap[ebin].Integral()
    hist_model_skymap[ebin].Scale(scaling_factor)
    print ('scaling_factor = %s'%(scaling_factor))

fig, ax = plt.subplots()

RA_PSR = 286.975
Dec_PSR = 6.03777777778
profile_center_x = RA_PSR
profile_center_y = Dec_PSR

profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindExtension(hist_flux_skymap_sum,hist_flux_syst_skymap_sum,profile_center_x,profile_center_y,2.)
profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindExtension(hist_pwn_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindExtension(hist_model_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindExtension(hist_psr_skymap_sum,None,profile_center_x,profile_center_y,2.)
profile_total = np.array(profile_pwn)+np.array(profile_model)
fig.clf()
axbig = fig.add_subplot()
axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='excess')
axbig.plot(theta2_pwn,profile_pwn,color='r')
axbig.plot(theta2_model,profile_model,color='b')
axbig.plot(theta2_fit,profile_fit,color='g')
axbig.plot(theta2_pwn,profile_total,color='k')
axbig.set_ylabel('[counts $\mathrm{deg}^{-2}$]')
axbig.set_xlabel('angular distance from source [degree]')
axbig.legend(loc='best')
plotname = 'ProfileModel_Sum'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

for ebin in range(energy_bin_cut_low,energy_bin_cut_up):

    profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindExtension(hist_flux_skymap[ebin],hist_flux_syst_skymap[ebin],profile_center_x,profile_center_y,2.)
    profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindExtension(hist_pwn_skymap[ebin],None,profile_center_x,profile_center_y,2.)
    profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindExtension(hist_model_skymap[ebin],None,profile_center_x,profile_center_y,2.)
    profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindExtension(hist_psr_skymap[ebin],None,profile_center_x,profile_center_y,2.)
    profile_total = np.array(profile_pwn)+np.array(profile_model)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='excess')
    axbig.plot(theta2_pwn,profile_pwn,color='r')
    axbig.plot(theta2_model,profile_model,color='b')
    axbig.plot(theta2_fit,profile_fit,color='g')
    axbig.plot(theta2_pwn,profile_total,color='k')
    axbig.set_ylabel('[counts $\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileModelRadial_E%s'%(ebin)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()

    travel_angle = 71.94
    travel_distance = 1.0
    profile_center_x1 = 286.98
    profile_center_y1 = 6.04
    profile_center_x2 = travel_distance*np.cos(travel_angle*3.14/180.)+profile_center_x1
    profile_center_y2 = travel_distance*np.sin(travel_angle*3.14/180.)+profile_center_y1
    profile_center_z = 0.
    
    profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindProjection(hist_flux_skymap[ebin],hist_flux_syst_skymap[ebin],profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
    profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindProjection(hist_pwn_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
    profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindProjection(hist_model_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
    profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindProjection(hist_psr_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,False,0.)
    profile_total = np.array(profile_pwn)+np.array(profile_model)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='excess')
    axbig.plot(theta2_pwn,profile_pwn,color='r')
    axbig.plot(theta2_model,profile_model,color='b')
    axbig.plot(theta2_fit,profile_fit,color='g')
    axbig.plot(theta2_pwn,profile_total,color='k')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
    axbig.set_xlabel('distance along proper motion [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileModelProjectedX_E%s'%(ebin)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()
    
    profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindProjection(hist_flux_skymap[ebin],hist_flux_syst_skymap[ebin],profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
    profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindProjection(hist_pwn_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
    profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindProjection(hist_model_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
    profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindProjection(hist_psr_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.)
    profile_total = np.array(profile_pwn)+np.array(profile_model)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='excess')
    axbig.plot(theta2_pwn,profile_pwn,color='r')
    axbig.plot(theta2_fit,profile_fit,color='g')
    axbig.plot(theta2_model,profile_model,color='b')
    axbig.plot(theta2_pwn,profile_total,color='k')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
    axbig.set_xlabel('distance along proper motion [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileModelProjectedY0_E%s'%(ebin)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()
    
    profile_data, profile_data_err, theta2_data = CommonPlotFunctions.FindProjection(hist_flux_skymap[ebin],hist_flux_syst_skymap[ebin],profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
    profile_pwn, profile_pwn_err, theta2_pwn = CommonPlotFunctions.FindProjection(hist_pwn_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
    profile_model, profile_model_err, theta2_model = CommonPlotFunctions.FindProjection(hist_model_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
    profile_fit, profile_fit_err, theta2_fit = CommonPlotFunctions.FindProjection(hist_psr_skymap[ebin],None,profile_center_x1,profile_center_y1,profile_center_x2,profile_center_y2,profile_center_z,0.2,True,0.5)
    profile_total = np.array(profile_pwn)+np.array(profile_model)
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(theta2_data,profile_data,profile_data_err,color='k',marker='s',ls='none',label='excess')
    axbig.plot(theta2_pwn,profile_pwn,color='r')
    axbig.plot(theta2_model,profile_model,color='b')
    axbig.plot(theta2_fit,profile_fit,color='g')
    axbig.plot(theta2_pwn,profile_total,color='k')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}^{%s}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]'%(-1+energy_index))
    axbig.set_xlabel('distance along proper motion [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileModelProjectedY1_E%s'%(ebin)
    fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
    axbig.remove()
    
