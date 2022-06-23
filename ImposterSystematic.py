
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
energy_index_scale = CommonPlotFunctions.energy_index_scale
energy_bin = CommonPlotFunctions.energy_bin
energy_bin_cut_low = 0
energy_bin_cut_up = 4

def FluxBiasCorrection(real_flux,imposter_fluxes, imposter_biases):

    real_flux_correct = real_flux
    imposter_fluxes_correct = imposter_fluxes
    #return real_flux_correct, imposter_fluxes_correct

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for binx in range(0,real_flux[ebin].GetNbinsX()):
            for biny in range(0,real_flux[ebin].GetNbinsY()):
                imposter_bias = 0.
                for imposter in range(0,len(imposter_fluxes)):
                    imposter_bias += imposter_biases[imposter][ebin].GetBinContent(binx+1,biny+1)/float(len(imposter_fluxes))
                real_flux_old = real_flux_correct[ebin].GetBinContent(binx+1,biny+1)
                real_flux_correct[ebin].SetBinContent(binx+1,biny+1,real_flux_old-imposter_bias)
                for imposter in range(0,len(imposter_fluxes)):
                    imposter_flux_old = imposter_fluxes_correct[imposter][ebin].GetBinContent(binx+1,biny+1)
                    imposter_fluxes_correct[imposter][ebin].SetBinContent(binx+1,biny+1,imposter_flux_old-imposter_bias)
    return real_flux_correct, imposter_fluxes_correct

def GetRelSystError(roi_x,roi_y,roi_r,roi_name):

    imposter_data_list = []
    imposter_bkgd_list = []
    for imposter in range(0,5):
        energy_axis, energy_error, imposter_data, imposter_data_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_data_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
        energy_axis, energy_error, imposter_bkgd, imposter_bkgd_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_bkgd_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
        imposter_data_list += [imposter_data]
        imposter_bkgd_list += [imposter_bkgd]

    imposter_s2b_list = []
    imposter_s2b_err_list = []
    for imposter in range(0,5):
        imposter_s2b = []
        imposter_s2b_err = []
        for ebin in range(0,len(energy_axis)):
            s2b = (imposter_data_list[imposter][ebin]-imposter_bkgd_list[imposter][ebin])/imposter_bkgd_list[imposter][ebin]
            s2b_err = pow(imposter_data_list[imposter][ebin],0.5)/imposter_bkgd_list[imposter][ebin]
            imposter_s2b += [s2b]
            imposter_s2b_err += [s2b_err]
        imposter_s2b_list += [imposter_s2b]
        imposter_s2b_err_list += [imposter_s2b_err]

    return imposter_s2b_list


#folder_name = 'output_2x2'
#folder_name = 'output_8x8'
folder_name = 'output_test'
source_name = []
source_name += ['IC443HotSpot']
source_name += ['WComae']
source_name += ['NGC1275']
source_name += ['1ES0229']
source_name += ['MGRO_J1908']
plot_tag = folder_name

imposter_rel_syst_biglist = []
for source in range(0,len(source_name)):

    InputFile = ROOT.TFile("output_fitting/%s_ON_skymap_%s.root"%(source_name[source],folder_name))
    HistName = "hist_data_skymap_0"
    nbins = InputFile.Get(HistName).GetNbinsX()
    MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
    MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
    MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
    InputFile.Close()

    MapPlotSize = 4.
    MapCenter_x = (MapEdge_left+MapEdge_right)/2.
    MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.
    MapPlot_left = MapCenter_x-0.5*MapPlotSize
    MapPlot_right = MapCenter_x+0.5*MapPlotSize
    MapPlot_lower = MapCenter_y-0.5*MapPlotSize
    MapPlot_upper = MapCenter_y+0.5*MapPlotSize

    hist_real_data_skymap = []
    hist_real_bkgd_skymap = []
    hist_real_flux_skymap = []
    hist_real_flux_syst_skymap = []
    for ebin in range(0,len(energy_bin)-1):
        hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_real_flux_syst_skymap += [ROOT.TH2D("hist_real_flux_syst_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    InputFile = ROOT.TFile("output_fitting/%s_ON_skymap_%s.root"%(source_name[source],folder_name))
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        HistName = "hist_data_skymap_%s"%(ebin)
        hist_real_data_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_bkgd_skymap_%s"%(ebin)
        hist_real_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
        HistName = "hist_energy_flux_skymap_%s"%(ebin)
        hist_real_flux_skymap[ebin].Add(InputFile.Get(HistName))
    InputFile.Close()

    hist_imposter_data_skymap = []
    hist_imposter_bkgd_skymap = []
    hist_imposter_bias_skymap = []
    hist_imposter_flux_skymap = []
    for imposter in range(0,5):
        hist_imposter_data_skymap_sublist = []
        hist_imposter_bkgd_skymap_sublist = []
        hist_imposter_bias_skymap_sublist = []
        hist_imposter_flux_skymap_sublist = []
        for ebin in range(0,len(energy_bin)-1):
            hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter%s_data_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
            hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bkgd_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
            hist_imposter_bias_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bias_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
            hist_imposter_flux_skymap_sublist += [ROOT.TH2D("hist_imposter%s_flux_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
        hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
        hist_imposter_bias_skymap += [hist_imposter_bias_skymap_sublist]
        hist_imposter_flux_skymap += [hist_imposter_flux_skymap_sublist]
    for imposter in range(0,5):
        InputFile = ROOT.TFile("output_fitting/%s_Imposter%s_skymap_%s.root"%(source_name[source],imposter+1,folder_name))
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            HistName = "hist_data_skymap_%s"%(ebin)
            hist_imposter_data_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.)
            HistName = "hist_bkgd_skymap_%s"%(ebin)
            hist_imposter_bkgd_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            HistName = "hist_energy_flux_skymap_%s"%(ebin)
            hist_imposter_flux_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        InputFile.Close()
    
    #hist_real_bkgd_skymap, hist_imposter_bkgd_skymap = FluxBiasCorrection(hist_real_bkgd_skymap, hist_imposter_bkgd_skymap, hist_imposter_bias_skymap)
    #hist_real_flux_skymap, hist_imposter_flux_skymap = FluxBiasCorrection(hist_real_flux_skymap, hist_imposter_flux_skymap, hist_imposter_flux_skymap)


    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = 1.5
    region_name = 'Center'
    imposter_rel_syst_biglist += [GetRelSystError(region_x,region_y,region_r,region_name)]

for ebin in range(0,len(energy_bin)-1):
    hist_imposter_rel_syst = ROOT.TH1D("hist_imposter_rel_syst","",20,-0.05,0.05)
    for source in range(0,len(source_name)):
        for imposter in range(0,5):
            rel_syst = imposter_rel_syst_biglist[source][imposter][ebin]
            hist_imposter_rel_syst.Fill(rel_syst)

    x_axis = []
    y_axis = []
    y_error = []
    for binx in range(0,hist_imposter_rel_syst.GetNbinsX()):
        x_axis += [hist_imposter_rel_syst.GetBinCenter(binx+1)]
        y_axis += [hist_imposter_rel_syst.GetBinContent(binx+1)]
        y_error += [pow(hist_imposter_rel_syst.GetBinContent(binx+1),0.5)]
    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(x_axis,y_axis,y_error,color='k',marker='_',ls='none')
    axbig.set_xlabel('relative error')
    axbig.set_ylabel('number of measurements')
    plotname = 'ImposterRelSyst_E%s'%(ebin)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    hist_imposter_rel_syst.Delete()
