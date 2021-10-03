
import os
import sys,ROOT
import array
import math
from prettytable import PrettyTable
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.time import Time
#from scipy import special

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import cycle

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")
ROOT.gROOT.LoadMacro("load_stl.h+") # access vectors and vector of vectors

plt.rcParams["figure.figsize"] = (12,8)
fig, ax = plt.subplots()

energy_bin_cut_low = 0
energy_bin_cut_up = 2

N_bins_for_deconv = 16
#N_bins_for_deconv = 8
#N_bins_for_deconv = 4
#N_bins_for_deconv = 2

#folder_path = 'output_nominal'
folder_path = 'output_16x16'
#folder_path = 'output_8x8'
#folder_path = 'output_4x4'
#folder_path = 'output_2x2'

ErecS_lower_cut = 0
ErecS_upper_cut = 0
Accuracy_source = []
AccuracyErr_source = []
data_exposure = []
validate_data_count = []
validate_bkgd_count = []
validate_rfov_count = []
validate_comb_count = []
data_count = []
bkgd_count = []
bkgd_rank0_count = []
bkgd_rank1_count = []
par8_count = []
par9_count = []
wpar9_count = []
dark_count = []
rank0_count = []
rank1_count = []
rank2_count = []
rank3_count = []
rank4_count = []

method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
method_tag += lowrank_tag

ONOFF_tag = 'OFF'
ONOFF_tag += '_Model0'
sample_list = []
sample_name = []
sample_list += ['1ES0502V6_OFF']
sample_name += ['1ES0502 V6']
sample_list += ['1ES0502V5_OFF']
sample_name += ['1ES0502 V5']
sample_list += ['DracoV6_OFF']
sample_name += ['Draco V6']
sample_list += ['DracoV5_OFF']
sample_name += ['Draco V5']
sample_list += ['PG1553V6_OFF']
sample_name += ['PG1553 V6']
sample_list += ['PG1553V5_OFF']
sample_name += ['PG1553 V5']
sample_list += ['M82V6_OFF']
sample_name += ['M82 V6']
sample_list += ['M82V5_OFF']
sample_name += ['M82 V5']
sample_list += ['1ES0414V5_OFF']
sample_name += ['1ES0410 V5']
sample_list += ['3C273V6_OFF']
sample_name += ['3C273 V6']
sample_list += ['3C273V5_OFF']
sample_name += ['3C273 V5']
sample_list += ['1ES0647V6_OFF']
sample_name += ['1ES0647 V6']
sample_list += ['1ES1011V6_OFF']
sample_name += ['1ES1011 V6']
sample_list += ['OJ287V6_OFF']
sample_name += ['OJ287 V6']
sample_list += ['PKS1424V6_OFF']
sample_name += ['PKS1424 V6']
sample_list += ['PKS1424V5_OFF']
sample_name += ['PKS1424 V5']
sample_list += ['3C264V6_OFF']
sample_name += ['3C264 V6']
sample_list += ['1ES0229V6_OFF']
sample_name += ['1ES0229 V6']
sample_list += ['1ES0229V5_OFF']
sample_name += ['1ES0229 V5']
sample_list += ['Segue1V6_OFF']
sample_name += ['Segue1 V6']
sample_list += ['Segue1V5_OFF']
sample_name += ['Segue1 V5']
sample_list += ['BLLacV6_OFF']
sample_name += ['BLLac V6']
sample_list += ['BLLacV5_OFF']
sample_name += ['BLLac V5']
sample_list += ['H1426V6_OFF']
sample_name += ['H1426 V6']
sample_list += ['RGBJ0710V5_OFF']
sample_name += ['RGBJ0710 V5']

sample_list += ['CrabV5_OFF']
sample_name += ['Crab V5']
sample_list += ['CrabV6_OFF']
sample_name += ['Crab V6']
##sample_list += ['RBS0413V6_OFF']
##sample_name += ['RBS0413 V6']
##sample_list += ['NGC1275V6_OFF']
##sample_name += ['NGC 1275 V6']
##sample_list += ['RBS0413V5_OFF']
##sample_name += ['RBS0413 V5']

#ONOFF_tag = 'ON'
#ONOFF_tag += '_Model0'
#sample_list = []
#sample_name = []
#sample_list += ['NGC1275V6_ON']
#sample_name += ['NGC 1275 V6']
##sample_list += ['GemingaV6_ON']
##sample_name += ['Geminga V6']
##sample_list += ['GemingaV5_ON']
##sample_name += ['Geminga V5']
#sample_list += ['CrabV6_ON']
#sample_name += ['Crab V6']
#sample_list += ['CrabV5_ON']
#sample_name += ['Crab V5']
    
elev_bins = [35,45,55,65,75,85]
#elev_bins = [45,55,65,75,85]
lowrank_tag += 'elev_incl'

theta2_bins = [0,4]
#theta2_bins = [0,6]

number_of_roi = 4

energy_bin = []
energy_bin += [100]
energy_bin += [251]
energy_bin += [631]
energy_bin += [1585]
energy_bin += [3981]
energy_bin += [10000]
energy_bin += [25118]

#energy_mibe_weight = [0.039,0.033,0.032,0.040,0.064,0.073,0.141]
#energy_rfov_weight = [0.061,0.048,0.036,0.050,0.058,0.072,0.055]
energy_mibe_weight = [1,1,1,1,1,1,1]
energy_rfov_weight = [1,1,1,1,1,1,1]

energy_dependent_singularvalue = []

energy_dependent_stat = []
energy_dependent_syst = []
energy_dependent_syst_rank0 = []
energy_dependent_syst_rank1 = []
energy_dependent_syst_init = []
energy_dependent_stat_vali = []
energy_dependent_syst_vbkg = []
energy_dependent_syst_rfov = []
energy_dependent_syst_comb = []
energy_dependent_rank0 = []
energy_dependent_rank1 = []
energy_dependent_rank2 = []
energy_dependent_rank3 = []

root_file_tags = []
mjd_tag = []
mjd_tag += ['']
for elev in range(0,len(elev_bins)-1):
    elev_tag = '_TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

def SystematicErrorMeasurement(list_measurements, list_err_measurements):

    syst_err = 0.
    stat_err = 0.
    sample_weight = 0.
    for entry in range(0,len(list_measurements)):
        if list_err_measurements[entry]==0.: continue
        syst_err += pow(1./list_err_measurements[entry],2)*pow(list_measurements[entry],2)
        stat_err += pow(1./list_err_measurements[entry],2)*pow(list_err_measurements[entry],2)
        sample_weight += pow(1./list_err_measurements[entry],2)
    if sample_weight>0.: 
        syst_err = pow(syst_err/sample_weight,0.5)
        stat_err = pow(stat_err/sample_weight,0.5)


    return syst_err, stat_err

def Make2DPlot(Hist2D,title_x,title_y,name,logz,min_z,max_z):

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
    #lumilab1 = ROOT.TLatex(0.15,0.50,'#tau_{kn} coefficents' )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.3)
    #lumilab1.Draw()
    pad1.cd()
    if logz: pad1.SetLogz()
    Hist2D.GetYaxis().SetTitle(title_y)
    Hist2D.GetXaxis().SetTitle(title_x)
    Hist2D.SetMaximum(max_z)
    Hist2D.SetMinimum(min_z)
    Hist2D.Draw("COL4Z")
    #Hist2D.Draw("TEXT45 same")
    canvas.SaveAs('output_plots/%s.png'%(name))

def MakeMultipleFitPlot(ax,Hists,legends,colors,title_x,title_y):

    func_gauss = []
    for h in range(0,len(Hists)):
        func_gauss += [ROOT.TF1("func_gauss_%s"%(h),"[0]*TMath::Gaus(x,[1],[2]*pow(2,0.5))",-0.2,0.2)]
        func_gauss[h].SetParameters(Hists[h].Integral(),0,0.02)
        Hists[h].Fit("func_gauss_%s"%(h),"N")
        func_gauss[h].SetLineColor(colors[h])
        func_gauss[h].Draw("E same")
        print ("func_gauss[%s].GetNDF() = %s"%(h,func_gauss[h].GetNDF()))
        print ("chi2/NDF = %s"%(func_gauss[h].GetChisquare()/func_gauss[h].GetNDF()))

    hist_xdata = []
    hist_ydata = []
    hist_error = []
    for entry in range(0,len(Hists)):
        xdata = []
        ydata = []
        error = []
        for binx in range(0,Hists[entry].GetNbinsX()):
            xdata += [Hists[entry].GetBinCenter(binx+1)]
            ydata += [Hists[entry].GetBinContent(binx+1)]
            error += [Hists[entry].GetBinError(binx+1)]
        hist_xdata += [xdata]
        hist_ydata += [ydata]
        hist_error += [error]

    func_xdata = []
    func_ydata = []
    for entry in range(0,len(Hists)):
        xdata = []
        ydata = []
        error = []
        for binx in range(0,200):
            xdata += [-0.2+0.4/200.*binx]
            ydata += [func_gauss[entry].Eval(-0.2+0.4/200.*binx)]
        func_xdata += [xdata]
        func_ydata += [ydata]
        hist_error += [error]

    cycol = cycle('brgcmk')
    for entry in range(0,len(Hists)):
        next_color = next(cycol)
        ax.errorbar(hist_xdata[entry], hist_ydata[entry], hist_error[entry], color=next_color, marker='s', ls='none', label='%s'%(legends[entry]))
        ax.plot(func_xdata[entry], func_ydata[entry], color=next_color)

    ax.legend(loc='best')
    ax.set_xlabel(title_x)
    ax.set_ylabel(title_y)
    return(ax)

def MakeMultiplePlot(ax,Hists,legends,colors,title_x,title_y,name,x_min,x_max,logx,logy):
    
    hist_xdata = []
    hist_ydata = []
    for entry in range(0,len(Hists)):
        xdata = []
        ydata = []
        xdata += [Hists[entry].GetBinLowEdge(1)]
        ydata += [0.]
        for binx in range(0,Hists[entry].GetNbinsX()):
            xdata += [Hists[entry].GetBinLowEdge(binx+2)]
            ydata += [Hists[entry].GetBinContent(binx+1)]
        hist_xdata += [xdata]
        hist_ydata += [ydata]

    ax.step(hist_xdata[0], hist_ydata[0], color='black', label='%s'%(legends[0]), alpha=1.0)
    for entry in range(1,len(Hists)):
        ax.step(hist_xdata[entry], hist_ydata[entry], label='%s'%(legends[entry]), alpha=0.5)

    ax.legend(loc='best')
    ax.set_xlabel(title_x)
    ax.set_ylabel(title_y)
    if (x_min!=x_max):
        ax.set_xlim(x_min,x_max)
    return(ax)


def Smooth2DMap(Hist_Old):

    smooth_size = 0.5
    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(2*smooth_size/bin_size) + 1
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
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
    return Hist_Smooth

def GetHistogramsFromFile(FilePath,which_source):
    global data_exposure
    global validate_data_count
    global validate_bkgd_count
    global validate_rfov_count
    global validate_comb_count
    global data_count
    global bkgd_count
    global bkgd_rank0_count
    global bkgd_rank1_count
    global par8_count
    global par9_count
    global wpar9_count
    global dark_count
    global rank0_count
    global rank1_count
    global rank2_count
    global rank3_count
    global rank4_count
    data_gamma_count = ROOT.std.vector("double")(10)
    bkgd_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank0_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank1_gamma_count = ROOT.std.vector("double")(10)
    par8_gamma_count = ROOT.std.vector("double")(10)
    par9_gamma_count = ROOT.std.vector("double")(10)
    wpar9_gamma_count = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)
    rank0_gamma_count = ROOT.std.vector("double")(10)
    rank1_gamma_count = ROOT.std.vector("double")(10)
    rank2_gamma_count = ROOT.std.vector("double")(10)
    rank3_gamma_count = ROOT.std.vector("double")(10)
    rank4_gamma_count = ROOT.std.vector("double")(10)
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    NewInfoTree = InputFile.Get("NewInfoTree")
    InfoTree.GetEntry(0)
    NewInfoTree.SetBranchAddress('data_gamma_count',ROOT.AddressOf(data_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_gamma_count',ROOT.AddressOf(bkgd_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank0_gamma_count',ROOT.AddressOf(bkgd_rank0_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank1_gamma_count',ROOT.AddressOf(bkgd_rank1_gamma_count))
    NewInfoTree.SetBranchAddress('par8_gamma_count',ROOT.AddressOf(par8_gamma_count))
    NewInfoTree.SetBranchAddress('par9_gamma_count',ROOT.AddressOf(par9_gamma_count))
    NewInfoTree.SetBranchAddress('wpar9_gamma_count',ROOT.AddressOf(wpar9_gamma_count))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.SetBranchAddress('rank0_gamma_count',ROOT.AddressOf(rank0_gamma_count))
    NewInfoTree.SetBranchAddress('rank1_gamma_count',ROOT.AddressOf(rank1_gamma_count))
    NewInfoTree.SetBranchAddress('rank2_gamma_count',ROOT.AddressOf(rank2_gamma_count))
    NewInfoTree.SetBranchAddress('rank3_gamma_count',ROOT.AddressOf(rank3_gamma_count))
    NewInfoTree.SetBranchAddress('rank4_gamma_count',ROOT.AddressOf(rank4_gamma_count))
    NewInfoTree.GetEntry(0)
    data_validate_count = NewInfoTree.data_validate_count
    bkgd_validate_count = NewInfoTree.bkgd_validate_count
    rfov_validate_count = NewInfoTree.rfov_validate_count
    #comb_validate_count = NewInfoTree.comb_validate_count
    roi_name = InfoTree.roi_name
    exposure_hours = InfoTree.exposure_hours
    data_exposure[which_source] += exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    for nth_roi in range(0,number_of_roi):
        print ("nth_roi = %s"%(nth_roi))
        validate_data_count[nth_roi][which_source]  += data_validate_count[energy_index][2*nth_roi+1]
        validate_bkgd_count[nth_roi][which_source]  += bkgd_validate_count[energy_index][2*nth_roi+1]
        validate_rfov_count[nth_roi][which_source]  += rfov_validate_count[energy_index][2*nth_roi+1]
        #validate_comb_count[nth_roi][which_source]  += comb_validate_count[energy_index][2*nth_roi+1]
        validate_comb_count[nth_roi][which_source] += (1./energy_mibe_weight[energy_index]*bkgd_validate_count[energy_index][2*nth_roi+1]+1./energy_rfov_weight[energy_index]*rfov_validate_count[energy_index][2*nth_roi+1])/(1./energy_mibe_weight[energy_index]+1./energy_rfov_weight[energy_index])
    data_count[which_source]  += data_gamma_count[energy_index]
    bkgd_count[which_source]  += bkgd_gamma_count[energy_index]
    bkgd_rank0_count[which_source]  += bkgd_rank0_gamma_count[energy_index]
    bkgd_rank1_count[which_source]  += bkgd_rank1_gamma_count[energy_index]
    par8_count[which_source]  += par8_gamma_count[energy_index]
    par9_count[which_source]  += par9_gamma_count[energy_index]
    wpar9_count[which_source]  += wpar9_gamma_count[energy_index]
    dark_count[which_source]  += dark_gamma_count[energy_index]
    rank0_count[which_source]  += rank0_gamma_count[energy_index]
    rank1_count[which_source]  += rank1_gamma_count[energy_index]
    rank2_count[which_source]  += rank2_gamma_count[energy_index]
    rank3_count[which_source]  += rank3_gamma_count[energy_index]
    rank4_count[which_source]  += rank4_gamma_count[energy_index]
    Hist_Bkgd_Optimization[which_source+1].Reset()
    weight = 1./float(len(sample_list))
    #weight = exposure_hours
    HistName = "Hist_OnData_SR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_SR_R2off.Reset()
    Hist_OnData_SR_R2off.Add(InputFile.Get(HistName))
    total_data_count = Hist_OnData_SR_R2off.Integral()
    HistName = "Hist_OnDark_SR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_SR_R2off.Reset()
    Hist_OnDark_SR_R2off.Add(InputFile.Get(HistName))
    total_dark_count = Hist_OnDark_SR_R2off.Integral()
    if total_dark_count>0.:
        Hist_OnDark_SR_R2off.Scale(total_data_count/total_dark_count)
    HistName = "Hist_OnData_CR_R2off_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_CR_R2off.Reset()
    Hist_OnData_CR_R2off.Add(InputFile.Get(HistName))
    total_bkgd_count = Hist_OnData_CR_R2off.Integral()
    if total_bkgd_count>0.:
        Hist_OnData_CR_R2off.Scale(total_data_count/total_bkgd_count)
    HistName = "Hist_OnData_SR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_SR_XYoff.Reset()
    Hist_OnData_SR_XYoff.Add(InputFile.Get(HistName))
    total_data_count = Hist_OnData_SR_XYoff.Integral()
    HistName = "Hist_OnDark_SR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_SR_XYoff.Reset()
    Hist_OnDark_SR_XYoff.Add(InputFile.Get(HistName))
    total_dark_count = Hist_OnDark_SR_XYoff.Integral()
    if total_dark_count>0.:
        Hist_OnDark_SR_XYoff.Scale(total_data_count/total_dark_count)
    HistName = "Hist_OnData_CR_XYoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_CR_XYoff.Reset()
    Hist_OnData_CR_XYoff.Add(InputFile.Get(HistName))
    total_bkgd_count = Hist_OnData_CR_XYoff.Integral()
    if total_bkgd_count>0.:
        Hist_OnData_CR_XYoff.Scale(total_data_count/total_bkgd_count)

    for binx in range (0,Hist_OnBkgd_SystErr_R2off.GetNbinsX()):
        binx_center = Hist_OnBkgd_SystErr_R2off.GetBinCenter(binx+1)
        binx2 = Hist_OnData_SR_R2off.GetXaxis().FindBin(binx_center)
        data_bin_count = Hist_OnData_SR_R2off.GetBinContent(binx2)
        dark_bin_count = Hist_OnDark_SR_R2off.GetBinContent(binx2)
        bkgd_bin_count = Hist_OnData_CR_R2off.GetBinContent(binx2)
        if data_bin_count==0.: continue
        relative_stat = pow(1.*data_bin_count,0.5)/data_bin_count
        relative_bias_dark = (data_bin_count-dark_bin_count)/data_bin_count
        relative_bias_bkgd = (data_bin_count-bkgd_bin_count)/data_bin_count
        old_content_stat = Hist_OnData_StatErr_R2off.GetBinContent(binx+1)
        old_content_weight = Hist_OnData_StatWeight_R2off.GetBinContent(binx+1)
        old_content_dark = Hist_OnDark_InclErr_R2off.GetBinContent(binx+1)
        old_content_bkgd = Hist_OnBkgd_InclErr_R2off.GetBinContent(binx+1)
        old_content_bias = Hist_OnBkgd_Bias_R2off.GetBinContent(binx+1)
        Hist_OnData_StatWeight_R2off.SetBinContent(binx+1,old_content_weight+1./relative_stat)
        Hist_OnData_StatErr_R2off.SetBinContent(binx+1,old_content_stat+1./relative_stat*relative_stat*relative_stat)
        Hist_OnDark_InclErr_R2off.SetBinContent(binx+1,old_content_dark+1./relative_stat*relative_bias_dark*relative_bias_dark)
        Hist_OnBkgd_InclErr_R2off.SetBinContent(binx+1,old_content_bkgd+1./relative_stat*relative_bias_bkgd*relative_bias_bkgd)
        Hist_OnBkgd_Bias_R2off.SetBinContent(binx+1,old_content_bias+1./relative_stat*relative_bias_bkgd)
    for binx in range (0,Hist_OnBkgd_SystErr_XYoff.GetNbinsX()):
        for biny in range (0,Hist_OnBkgd_SystErr_XYoff.GetNbinsY()):
            binx_center = Hist_OnBkgd_SystErr_XYoff.GetXaxis().GetBinCenter(binx+1)
            biny_center = Hist_OnBkgd_SystErr_XYoff.GetYaxis().GetBinCenter(biny+1)
            binx2 = Hist_OnData_SR_XYoff.GetXaxis().FindBin(binx_center)
            biny2 = Hist_OnData_SR_XYoff.GetYaxis().FindBin(biny_center)
            data_bin_count = Hist_OnData_SR_XYoff.GetBinContent(binx2,biny2)
            dark_bin_count = Hist_OnDark_SR_XYoff.GetBinContent(binx2,biny2)
            bkgd_bin_count = Hist_OnData_CR_XYoff.GetBinContent(binx2,biny2)
            if data_bin_count==0.: continue
            relative_stat = pow(2.*data_bin_count,0.5)/data_bin_count
            relative_bias_dark = (data_bin_count-dark_bin_count)/data_bin_count
            relative_bias_bkgd = (data_bin_count-bkgd_bin_count)/data_bin_count
            old_content_stat = Hist_OnData_StatErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_weight = Hist_OnData_StatWeight_XYoff.GetBinContent(binx+1,biny+1)
            old_content_dark = Hist_OnDark_InclErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_bkgd = Hist_OnBkgd_InclErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_bias = Hist_OnBkgd_Bias_XYoff.GetBinContent(binx+1,biny+1)
            Hist_OnData_StatWeight_XYoff.SetBinContent(binx+1,biny+1,old_content_weight+1./relative_stat)
            Hist_OnData_StatErr_XYoff.SetBinContent(binx+1,biny+1,old_content_stat+1./relative_stat*relative_stat*relative_stat)
            Hist_OnDark_InclErr_XYoff.SetBinContent(binx+1,biny+1,old_content_dark+1./relative_stat*relative_bias_dark*relative_bias_dark)
            Hist_OnBkgd_InclErr_XYoff.SetBinContent(binx+1,biny+1,old_content_bkgd+1./relative_stat*relative_bias_bkgd*relative_bias_bkgd)
            Hist_OnBkgd_Bias_XYoff.SetBinContent(binx+1,biny+1,old_content_bias+1./relative_stat*relative_bias_bkgd)

    if energy_index>=energy_bin_cut_low and energy_index<=energy_bin_cut_up:
        if data_gamma_count[energy_index]>0.:
            measured_error_mibe = (data_gamma_count[energy_index]-bkgd_gamma_count[energy_index])/data_gamma_count[energy_index]
            measured_error_init = (data_gamma_count[energy_index]-dark_gamma_count[energy_index])/data_gamma_count[energy_index]
            measured_error_weight = 1.
            #measured_error_weight = 1./(pow(data_gamma_count[energy_index],0.5)/data_gamma_count[energy_index])
            Hist_SystErrDist_MDM.Fill(measured_error_mibe,measured_error_weight)
            Hist_SystErrDist_Init.Fill(measured_error_init,measured_error_weight)
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Optimization[which_source+1].Add(InputFile.Get(HistName))
    Hist_Data_Eigenvalues[which_source+1].Reset()
    HistName = "Hist_Data_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Data_Eigenvalues[0].Add(InputFile.Get(HistName),weight)
    Hist_Data_Eigenvalues[which_source+1].Add(InputFile.Get(HistName))
    Hist_GammaRegion_Contribution[which_source+1].Reset()
    HistName = "Hist_GammaRegion_Contribution_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_GammaRegion_Contribution[0].Add(InputFile.Get(HistName),weight)
    Hist_GammaRegion_Contribution[which_source+1].Add(InputFile.Get(HistName))


R2off_nbins = 9
XYoff_nbins = 12
#XYoff_nbins = 6
#XYoff_nbins = 4
#XYoff_nbins = 3
#XYoff_nbins = 2
#XYoff_nbins = 1
Hist_SystErrDist_MDM = ROOT.TH1D("Hist_SystErrDist_MDM","",20,-0.2,0.2)
Hist_SystErrDist_Init = ROOT.TH1D("Hist_SystErrDist_Init","",20,-0.2,0.2)
Hist_OnData_SR_R2off = ROOT.TH1D("Hist_OnData_SR_R2off","",R2off_nbins,0,9)
Hist_OnDark_SR_R2off = ROOT.TH1D("Hist_OnDark_SR_R2off","",R2off_nbins,0,9)
Hist_OnData_CR_R2off = ROOT.TH1D("Hist_OnData_CR_R2off","",R2off_nbins,0,9)
Hist_OnDark_SystErr_R2off = ROOT.TH1D("Hist_OnDark_SystErr_R2off","",R2off_nbins,0,9)
Hist_OnBkgd_SystErr_R2off = ROOT.TH1D("Hist_OnBkgd_SystErr_R2off","",R2off_nbins,0,9)
Hist_OnDark_InclErr_R2off = ROOT.TH1D("Hist_OnDark_InclErr_R2off","",R2off_nbins,0,9)
Hist_OnBkgd_InclErr_R2off = ROOT.TH1D("Hist_OnBkgd_InclErr_R2off","",R2off_nbins,0,9)
Hist_OnBkgd_Bias_R2off = ROOT.TH1D("Hist_OnBkgd_Bias_R2off","",R2off_nbins,0,9)
Hist_OnData_StatErr_R2off = ROOT.TH1D("Hist_OnData_StatErr_R2off","",R2off_nbins,0,9)
Hist_OnData_StatWeight_R2off = ROOT.TH1D("Hist_OnData_StatWeight_R2off","",R2off_nbins,0,9)
Hist_OnData_SR_XYoff = ROOT.TH2D("Hist_OnData_SR_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnDark_SR_XYoff = ROOT.TH2D("Hist_OnDark_SR_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnData_CR_XYoff = ROOT.TH2D("Hist_OnData_CR_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnDark_SystErr_XYoff = ROOT.TH2D("Hist_OnDark_SystErr_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnBkgd_SystErr_XYoff = ROOT.TH2D("Hist_OnBkgd_SystErr_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnDark_InclErr_XYoff = ROOT.TH2D("Hist_OnDark_InclErr_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnBkgd_InclErr_XYoff = ROOT.TH2D("Hist_OnBkgd_InclErr_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnBkgd_Bias_XYoff = ROOT.TH2D("Hist_OnBkgd_Bias_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnData_StatErr_XYoff = ROOT.TH2D("Hist_OnData_StatErr_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)
Hist_OnData_StatWeight_XYoff = ROOT.TH2D("Hist_OnData_StatWeight_XYoff","",XYoff_nbins,-3,3,XYoff_nbins,-3,3)

Hist_NormStatErr = []
Hist_NormSystErr = []
Hist_NormSystInitErr = []
Hist_ShapeSystErr = []
for elev in range(0,len(elev_bins)-1):
    elev_tag = 'TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    Hist_NormStatErr += [ROOT.TH1D("Hist_NormStatErr_"+elev_tag,"",len(energy_bin)-1,array('d',energy_bin))]
    Hist_NormSystErr += [ROOT.TH1D("Hist_NormSystErr_"+elev_tag,"",len(energy_bin)-1,array('d',energy_bin))]
    Hist_NormSystInitErr += [ROOT.TH1D("Hist_NormSystInitErr_"+elev_tag,"",len(energy_bin)-1,array('d',energy_bin))]
    Hist_ShapeSystErr_ThisElev = []
    for ebin in range(0,len(energy_bin)-1):
        Hist_ShapeSystErr_ThisElev += [ROOT.TH1D("Hist_ShapeSystErr_%s_ErecS%sto%s"%(elev_tag,int(energy_bin[ebin]),int(energy_bin[ebin+1])),"",R2off_nbins,0,9)]
    Hist_ShapeSystErr += [Hist_ShapeSystErr_ThisElev]

optimiz_lower = -5.
optimiz_upper = -1.
Hist_Bkgd_Optimization = []
for e in range(0,len(sample_list)+1):
    Hist_Bkgd_Optimization += [ROOT.TH1D("Hist_Bkgd_Optimization_%s"%(e),"",10,optimiz_lower,optimiz_upper)]
Hist_Data_Eigenvalues = []
for e in range(0,len(sample_list)+1):
    Hist_Data_Eigenvalues += [ROOT.TH1D("Hist_Data_Eigenvalues_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv)]
Hist_GammaRegion_Contribution = []
for e in range(0,len(sample_list)+1):
    Hist_GammaRegion_Contribution += [ROOT.TH1D("Hist_GammaRegion_Contribution_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv)]

for e in range(0,len(energy_bin)-1):
    for elev in range(0,len(root_file_tags)):
        FilePath_List = []
        validate_data_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
        validate_bkgd_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
        validate_rfov_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
        validate_comb_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
        data_count = [0.]*len(sample_list)
        bkgd_count = [0.]*len(sample_list)
        bkgd_rank0_count = [0.]*len(sample_list)
        bkgd_rank1_count = [0.]*len(sample_list)
        par8_count = [0.]*len(sample_list)
        par9_count = [0.]*len(sample_list)
        wpar9_count = [0.]*len(sample_list)
        rank0_count = [0.]*len(sample_list)
        rank1_count = [0.]*len(sample_list)
        rank2_count = [0.]*len(sample_list)
        rank3_count = [0.]*len(sample_list)
        rank4_count = [0.]*len(sample_list)
        data_exposure = [0.]*len(sample_list)
        dark_count = [0.]*len(sample_list)
        for entry in range(0,len(Hist_GammaRegion_Contribution)):
            Hist_GammaRegion_Contribution[entry].Reset()
        for entry in range(0,len(Hist_Bkgd_Optimization)):
            Hist_Bkgd_Optimization[entry].Reset()
        for entry in range(0,len(Hist_Data_Eigenvalues)):
            Hist_Data_Eigenvalues[entry].Reset()
        Hist_OnData_StatWeight_R2off.Reset()
        Hist_OnData_StatErr_R2off.Reset()
        Hist_OnDark_SystErr_R2off.Reset()
        Hist_OnBkgd_SystErr_R2off.Reset()
        Hist_OnDark_InclErr_R2off.Reset()
        Hist_OnBkgd_InclErr_R2off.Reset()
        Hist_OnBkgd_Bias_R2off.Reset()
        Accuracy_source = []
        AccuracyErr_source = []
        for source in range(0,len(sample_list)):
            source_name = sample_list[source]
            FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            print ('Trying to read file %s'%(FilePath_List[len(FilePath_List)-1]))
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1],source)
        for binx in range (0,Hist_OnBkgd_SystErr_R2off.GetNbinsX()):
            old_content_weight = Hist_OnData_StatWeight_R2off.GetBinContent(binx+1)
            if old_content_weight==0.: continue
            old_content_stat = Hist_OnData_StatErr_R2off.GetBinContent(binx+1)
            old_content_dark = Hist_OnDark_InclErr_R2off.GetBinContent(binx+1)
            old_content_bkgd = Hist_OnBkgd_InclErr_R2off.GetBinContent(binx+1)
            old_content_bias = Hist_OnBkgd_Bias_R2off.GetBinContent(binx+1)
            Hist_OnDark_InclErr_R2off.SetBinContent(binx+1,pow(old_content_dark/old_content_weight,0.5))
            Hist_OnBkgd_InclErr_R2off.SetBinContent(binx+1,pow(old_content_bkgd/old_content_weight,0.5))
            Hist_OnBkgd_Bias_R2off.SetBinContent(binx+1,old_content_bias/old_content_weight)
            old_content_dark = max(0.,old_content_dark-old_content_stat)
            old_content_bkgd = max(0.,old_content_bkgd-old_content_stat)
            Hist_OnDark_SystErr_R2off.SetBinContent(binx+1,pow(old_content_dark/old_content_weight,0.5))
            Hist_OnBkgd_SystErr_R2off.SetBinContent(binx+1,pow(old_content_bkgd/old_content_weight,0.5))
            Hist_OnData_StatErr_R2off.SetBinContent(binx+1,pow(old_content_stat/old_content_weight,0.5))
        for binx in range (0,Hist_OnBkgd_SystErr_R2off.GetNbinsX()):
            Hist_ShapeSystErr[elev][e].SetBinContent(binx+1,Hist_OnBkgd_SystErr_R2off.GetBinContent(binx+1))
        AccuracyBkgd_source = []
        AccuracyBkgdSigned_source = []
        AccuracyBkgdErr_source = []
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            if data_count[entry-1]<10.:
                AccuracyBkgd_source += [0.]
                AccuracyBkgdSigned_source += [0.]
                AccuracyBkgdErr_source += [0.]
            else:
                AccuracyBkgd_source += [abs(data_count[entry-1]-bkgd_count[entry-1])/data_count[entry-1]]
                AccuracyBkgdSigned_source += [(data_count[entry-1]-bkgd_count[entry-1])/data_count[entry-1]]
                AccuracyBkgdErr_source += [1./pow(data_count[entry-1],0.5)]
        AccuracyBkgd_mean, AccuracyBkgd_mean_error = SystematicErrorMeasurement(AccuracyBkgd_source,AccuracyBkgdErr_source)
        elev_dependent_syst = pow(max(0.,pow(AccuracyBkgd_mean,2)-pow(AccuracyBkgd_mean_error,2)),0.5)
        Hist_NormSystErr[elev].SetBinContent(e+1,elev_dependent_syst)
        Hist_NormStatErr[elev].SetBinContent(e+1,AccuracyBkgd_mean_error)
        AccuracyInit_source = []
        AccuracyInitErr_source = []
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            if data_count[entry-1]<10.:
                AccuracyInit_source += [0.]
                AccuracyInitErr_source += [0.]
            else:
                AccuracyInit_source += [abs(data_count[entry-1]-dark_count[entry-1])/data_count[entry-1]]
                AccuracyInitErr_source += [1./pow(data_count[entry-1],0.5)]
        AccuracyInit_mean, AccuracyInit_mean_error = SystematicErrorMeasurement(AccuracyInit_source,AccuracyInitErr_source)
        elev_dependent_syst = pow(max(0.,pow(AccuracyInit_mean,2)-pow(AccuracyInit_mean_error,2)),0.5)
        Hist_NormSystInitErr[elev].SetBinContent(e+1,elev_dependent_syst)


    FilePath_List = []
    validate_data_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
    validate_bkgd_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
    validate_rfov_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
    validate_comb_count = [ [0.]*len(sample_list) for i in range(number_of_roi)]
    data_count = [0.]*len(sample_list)
    bkgd_count = [0.]*len(sample_list)
    bkgd_rank0_count = [0.]*len(sample_list)
    bkgd_rank1_count = [0.]*len(sample_list)
    par8_count = [0.]*len(sample_list)
    par9_count = [0.]*len(sample_list)
    wpar9_count = [0.]*len(sample_list)
    rank0_count = [0.]*len(sample_list)
    rank1_count = [0.]*len(sample_list)
    rank2_count = [0.]*len(sample_list)
    rank3_count = [0.]*len(sample_list)
    rank4_count = [0.]*len(sample_list)
    data_exposure = [0.]*len(sample_list)
    dark_count = [0.]*len(sample_list)
    for entry in range(0,len(Hist_GammaRegion_Contribution)):
        Hist_GammaRegion_Contribution[entry].Reset()
    for entry in range(0,len(Hist_Bkgd_Optimization)):
        Hist_Bkgd_Optimization[entry].Reset()
    for entry in range(0,len(Hist_Data_Eigenvalues)):
        Hist_Data_Eigenvalues[entry].Reset()
    Hist_OnData_StatWeight_R2off.Reset()
    Hist_OnData_StatErr_R2off.Reset()
    Hist_OnDark_SystErr_R2off.Reset()
    Hist_OnBkgd_SystErr_R2off.Reset()
    Hist_OnDark_InclErr_R2off.Reset()
    Hist_OnBkgd_InclErr_R2off.Reset()
    Hist_OnBkgd_Bias_R2off.Reset()
    Accuracy_source = []
    AccuracyErr_source = []
    for source in range(0,len(sample_list)):
        for elev in range(0,len(root_file_tags)):
            source_name = sample_list[source]
            FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            print ('Trying to read file %s'%(FilePath_List[len(FilePath_List)-1]))
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1],source)

    for binx in range (0,Hist_OnBkgd_SystErr_R2off.GetNbinsX()):
        old_content_weight = Hist_OnData_StatWeight_R2off.GetBinContent(binx+1)
        if old_content_weight==0.: continue
        old_content_stat = Hist_OnData_StatErr_R2off.GetBinContent(binx+1)
        old_content_dark = Hist_OnDark_InclErr_R2off.GetBinContent(binx+1)
        old_content_bkgd = Hist_OnBkgd_InclErr_R2off.GetBinContent(binx+1)
        old_content_bias = Hist_OnBkgd_Bias_R2off.GetBinContent(binx+1)
        Hist_OnDark_InclErr_R2off.SetBinContent(binx+1,pow(old_content_dark/old_content_weight,0.5))
        Hist_OnBkgd_InclErr_R2off.SetBinContent(binx+1,pow(old_content_bkgd/old_content_weight,0.5))
        Hist_OnBkgd_Bias_R2off.SetBinContent(binx+1,old_content_bias/old_content_weight)
        old_content_dark = max(0.,old_content_dark-old_content_stat)
        old_content_bkgd = max(0.,old_content_bkgd-old_content_stat)
        Hist_OnDark_SystErr_R2off.SetBinContent(binx+1,pow(old_content_dark/old_content_weight,0.5))
        Hist_OnBkgd_SystErr_R2off.SetBinContent(binx+1,pow(old_content_bkgd/old_content_weight,0.5))
        Hist_OnData_StatErr_R2off.SetBinContent(binx+1,pow(old_content_stat/old_content_weight,0.5))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_OnBkgd_InclErr_R2off]
    legends += ['$\sigma_{\mathrm{stat}} \oplus \sigma_{\mathrm{shape}}$']
    colors += [1]
    Hists += [Hist_OnData_StatErr_R2off]
    legends += ['$\sigma_{\mathrm{stat}}$']
    colors += [4]
    Hists += [Hist_OnBkgd_SystErr_R2off]
    legends += ['$\sigma_{\mathrm{shape}}$']
    colors += [2]
    #Hists += [Hist_OnBkgd_Bias_R2off]
    #legends += ['bias']
    #colors += [3]
    ax.cla()
    MakeMultiplePlot(ax,Hists,legends,colors,'$\\theta^{2}$ from camera center','relative uncertainty','Theta2Errors_E%s%s'%(e,folder_path),0.,4.,False,False)
    fig.savefig("output_plots/Theta2Errors_E%s%s.png"%(e,folder_path))
    for binx in range (0,Hist_OnBkgd_SystErr_XYoff.GetNbinsX()):
        for biny in range (0,Hist_OnBkgd_SystErr_XYoff.GetNbinsY()):
            old_content_weight = Hist_OnData_StatWeight_XYoff.GetBinContent(binx+1,biny+1)
            if old_content_weight==0.: continue
            old_content_stat = Hist_OnData_StatErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_dark = Hist_OnDark_InclErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_bkgd = Hist_OnBkgd_InclErr_XYoff.GetBinContent(binx+1,biny+1)
            old_content_bias = Hist_OnBkgd_Bias_XYoff.GetBinContent(binx+1,biny+1)
            Hist_OnDark_InclErr_XYoff.SetBinContent(binx+1,biny+1,pow(old_content_dark/old_content_weight,0.5))
            Hist_OnBkgd_InclErr_XYoff.SetBinContent(binx+1,biny+1,pow(old_content_bkgd/old_content_weight,0.5))
            Hist_OnBkgd_Bias_XYoff.SetBinContent(binx+1,biny+1,old_content_bias/old_content_weight)
            old_content_dark = max(0.,old_content_dark-old_content_stat)
            old_content_bkgd = max(0.,old_content_bkgd-old_content_stat)
            Hist_OnDark_SystErr_XYoff.SetBinContent(binx+1,biny+1,pow(old_content_dark/old_content_weight,0.5))
            Hist_OnBkgd_SystErr_XYoff.SetBinContent(binx+1,biny+1,pow(old_content_bkgd/old_content_weight,0.5))
            Hist_OnData_StatErr_XYoff.SetBinContent(binx+1,biny+1,pow(old_content_stat/old_content_weight,0.5))
    Make2DPlot(Hist_OnData_StatErr_XYoff,'Xoff','Yoff','StatErr_XYoff_Data_E%s%s'%(e,folder_path),False,0.,0.1)
    Make2DPlot(Hist_OnDark_SystErr_XYoff,'Xoff','Yoff','SystErr_XYoff_Dark_E%s%s'%(e,folder_path),False,0.,0.1)
    Make2DPlot(Hist_OnBkgd_SystErr_XYoff,'Xoff','Yoff','SystErr_XYoff_Bkgd_E%s%s'%(e,folder_path),False,0.,0.1)
    Make2DPlot(Hist_OnDark_InclErr_XYoff,'Xoff','Yoff','InclErr_XYoff_Dark_E%s%s'%(e,folder_path),False,0.,0.1)
    Make2DPlot(Hist_OnBkgd_InclErr_XYoff,'Xoff','Yoff','InclErr_XYoff_Bkgd_E%s%s'%(e,folder_path),False,0.,0.1)
    Make2DPlot(Hist_OnBkgd_Bias_XYoff,'Xoff','Yoff','Bias_XYoff_Bkgd_E%s%s'%(e,folder_path),False,-0.05,0.05)
    Hist_OnData_StatErr_XYoff.Reset()
    Hist_OnDark_SystErr_XYoff.Reset()
    Hist_OnBkgd_SystErr_XYoff.Reset()
    Hist_OnDark_InclErr_XYoff.Reset()
    Hist_OnBkgd_InclErr_XYoff.Reset()
    Hist_OnBkgd_Bias_XYoff.Reset()

    for binx in range(1,Hist_Bkgd_Optimization[0].GetNbinsX()+1):
        total_weight = 0.
        y_content = 0.
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            weight = pow(data_count[entry-1],0.5)
            #weight = 1.
            total_weight += weight
            y_content += weight*Hist_Bkgd_Optimization[entry].GetBinContent(binx)
        if total_weight>0.: y_content = y_content/total_weight
        Hist_Bkgd_Optimization[0].SetBinContent(binx,y_content)


    plt.clf()
    ax = fig.add_subplot(111)

    AccuracyInit_source = []
    AccuracyInitSigned_source = []
    AccuracyInitErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyInit_source += [0.]
            AccuracyInitSigned_source += [0.]
            AccuracyInitErr_source += [0.]
        else:
            AccuracyInit_source += [abs(data_count[entry-1]-dark_count[entry-1])/data_count[entry-1]]
            AccuracyInitSigned_source += [(data_count[entry-1]-dark_count[entry-1])/data_count[entry-1]]
            AccuracyInitErr_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyInit_mean, AccuracyInit_mean_error = SystematicErrorMeasurement(AccuracyInit_source,AccuracyInitErr_source)

    AccuracyStat_source = []
    AccuracyStatErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyStat_source += [0.]
            AccuracyStatErr_source += [0.]
        else:
            AccuracyStat_source += [pow(data_count[entry-1],0.5)/data_count[entry-1]]
            AccuracyStatErr_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyStat_mean = 0.
    AccuracyStat_weight = 0.
    AccuracyStat_mean_error = 0.
    for entry in range(0,len(AccuracyStat_source)):
        AccuracyStat_mean += AccuracyStat_source[entry]
        AccuracyStat_weight += 1.
        AccuracyStat_mean_error += pow(AccuracyStatErr_source[entry],2)/len(AccuracyStat_source)
    if AccuracyStat_weight>0.: AccuracyStat_mean = AccuracyStat_mean/AccuracyStat_weight
    AccuracyStat_mean_error = pow(AccuracyStat_mean_error,0.5)

    AccuracyRank0_source = []
    AccuracyRank0Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyRank0_source += [0.]
            AccuracyRank0Err_source += [0.]
        else:
            AccuracyRank0_source += [abs(data_count[entry-1]-rank0_count[entry-1])/data_count[entry-1]]
            AccuracyRank0Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyRank0_mean, AccuracyRank0_mean_error = SystematicErrorMeasurement(AccuracyRank0_source,AccuracyRank0Err_source)

    AccuracyRank1_source = []
    AccuracyRank1Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyRank1_source += [0.]
            AccuracyRank1Err_source += [0.]
        else:
            AccuracyRank1_source += [abs(data_count[entry-1]-rank1_count[entry-1])/data_count[entry-1]]
            AccuracyRank1Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyRank1_mean, AccuracyRank1_mean_error = SystematicErrorMeasurement(AccuracyRank1_source,AccuracyRank1Err_source)

    AccuracyRank2_source = []
    AccuracyRank2Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyRank2_source += [0.]
            AccuracyRank2Err_source += [0.]
        else:
            AccuracyRank2_source += [abs(data_count[entry-1]-rank2_count[entry-1])/data_count[entry-1]]
            AccuracyRank2Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyRank2_mean, AccuracyRank2_mean_error = SystematicErrorMeasurement(AccuracyRank2_source,AccuracyRank2Err_source)

    AccuracyRank3_source = []
    AccuracyRank3Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyRank3_source += [0.]
            AccuracyRank3Err_source += [0.]
        else:
            AccuracyRank3_source += [abs(data_count[entry-1]-rank3_count[entry-1])/data_count[entry-1]]
            AccuracyRank3Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyRank3_mean, AccuracyRank3_mean_error = SystematicErrorMeasurement(AccuracyRank3_source,AccuracyRank3Err_source)

    AccuracyBestPar9_source = []
    AccuracyBestPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyBestPar9_source += [0.]
            AccuracyBestPar9Err_source += [0.]
        else:
            AccuracyBestPar9_source += [pow(Hist_GammaRegion_Contribution[entry].GetBinContent(2+1),0.5)]
            AccuracyBestPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBestPar9_mean, AccuracyBestPar9_mean_error = SystematicErrorMeasurement(AccuracyBestPar9_source,AccuracyBestPar9Err_source)

    AccuracyPar9_source = []
    AccuracyPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyPar9_source += [0.]
            AccuracyPar9Err_source += [0.]
        else:
            AccuracyPar9_source += [abs(data_count[entry-1]-par9_count[entry-1])/data_count[entry-1]]
            AccuracyPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyPar9_mean, AccuracyPar9_mean_error = SystematicErrorMeasurement(AccuracyPar9_source,AccuracyPar9Err_source)

    AccuracyWPar9_source = []
    AccuracyWPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyWPar9_source += [0.]
            AccuracyWPar9Err_source += [0.]
        else:
            AccuracyWPar9_source += [abs(data_count[entry-1]-wpar9_count[entry-1])/data_count[entry-1]]
            AccuracyWPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyWPar9_mean, AccuracyWPar9_mean_error = SystematicErrorMeasurement(AccuracyWPar9_source,AccuracyWPar9Err_source)

    AccuracyBkgd_source = []
    AccuracyBkgdSigned_source = []
    AccuracyBkgdErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyBkgd_source += [0.]
            AccuracyBkgdSigned_source += [0.]
            AccuracyBkgdErr_source += [0.]
        else:
            AccuracyBkgd_source += [abs(data_count[entry-1]-bkgd_count[entry-1])/data_count[entry-1]]
            AccuracyBkgdSigned_source += [(data_count[entry-1]-bkgd_count[entry-1])/data_count[entry-1]]
            AccuracyBkgdErr_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBkgd_mean, AccuracyBkgd_mean_error = SystematicErrorMeasurement(AccuracyBkgd_source,AccuracyBkgdErr_source)

    AccuracyBkgdRank0_source = []
    AccuracyBkgdRank0Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyBkgdRank0_source += [0.]
            AccuracyBkgdRank0Err_source += [0.]
        else:
            AccuracyBkgdRank0_source += [abs(data_count[entry-1]-bkgd_rank0_count[entry-1])/data_count[entry-1]]
            AccuracyBkgdRank0Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBkgdRank0_mean, AccuracyBkgdRank0_mean_error = SystematicErrorMeasurement(AccuracyBkgdRank0_source,AccuracyBkgdRank0Err_source)

    AccuracyBkgdRank1_source = []
    AccuracyBkgdRank1Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyBkgdRank1_source += [0.]
            AccuracyBkgdRank1Err_source += [0.]
        else:
            AccuracyBkgdRank1_source += [abs(data_count[entry-1]-bkgd_rank1_count[entry-1])/data_count[entry-1]]
            AccuracyBkgdRank1Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBkgdRank1_mean, AccuracyBkgdRank1_mean_error = SystematicErrorMeasurement(AccuracyBkgdRank1_source,AccuracyBkgdRank1Err_source)

    AccuracyPar8_source = []
    AccuracyPar8Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<10.:
            AccuracyPar8_source += [0.]
            AccuracyPar8Err_source += [0.]
        else:
            AccuracyPar8_source += [abs(data_count[entry-1]-par8_count[entry-1])/data_count[entry-1]]
            AccuracyPar8Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyPar8_mean, AccuracyPar8_mean_error = SystematicErrorMeasurement(AccuracyPar8_source,AccuracyPar8Err_source)

    #plt.clf()
    #plt.xlabel("Zenith", fontsize=18)
    #plt.ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    #plt.errorbar(Zenith_mean_data,AccuracyInit_source,xerr=Zenith_RMS_data,yerr=AccuracyInitErr_source,fmt='o',color='r')
    #xmin, xmax, ymin, ymax = plt.axis()
    #x = np.linspace(xmin, xmax, 100)
    #y = 0.*x + AccuracyInit_mean
    #plt.plot(x,y,color='#1B2ACC')
    #plt.fill_between(x,y-AccuracyInit_mean_error,y+AccuracyInit_mean_error,edgecolor='#CC4F1B',alpha=0.2,facecolor='#FF9848')
    #plt.title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyInit_mean,AccuracyInit_mean_error))
    #plt.savefig("output_plots/PerformanceInit_Zenith_E%s%s.png"%(e,lowrank_tag))
    #plt.clf()
    #plt.xlabel("NSB", fontsize=18)
    #plt.ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    #plt.errorbar(NSB_mean_data,AccuracyInit_source,xerr=NSB_RMS_data,yerr=AccuracyInitErr_source,fmt='o',color='r')
    #xmin, xmax, ymin, ymax = plt.axis()
    #x = np.linspace(xmin, xmax, 100)
    #y = 0.*x + AccuracyInit_mean
    #plt.plot(x,y,color='#1B2ACC')
    #plt.fill_between(x,y-AccuracyInit_mean_error,y+AccuracyInit_mean_error,edgecolor='#CC4F1B',alpha=0.2,facecolor='#FF9848')
    #plt.title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyInit_mean,AccuracyInit_mean_error))
    #plt.savefig("output_plots/PerformanceInit_NSB_E%s%s.png"%(e,lowrank_tag))


    Accuracy_source = []
    AccuracyErr_source = []
    min_y = 1.0
    min_x = 1.0
    min_bin = 0
    for binx in range(1,Hist_Bkgd_Optimization[0].GetNbinsX()+1):
        if Hist_Bkgd_Optimization[0].GetBinContent(binx)<min_y:
            min_y = Hist_Bkgd_Optimization[0].GetBinContent(binx)
            min_x = Hist_Bkgd_Optimization[0].GetBinCenter(binx)
            min_bin = binx
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]==0.: continue
        Accuracy_source += [Hist_Bkgd_Optimization[entry].GetBinContent(min_bin)]
        AccuracyErr_source += [1./pow(data_count[entry-1],0.5)]
    Accuracy_mean = 0.
    Accuracy_weight = 0.
    Accuracy_mean_error = 0.
    for entry in range(0,len(Accuracy_source)):
        Accuracy_mean += 1./AccuracyErr_source[entry]*Hist_Bkgd_Optimization[0].GetBinContent(min_bin)
        Accuracy_weight += 1./AccuracyErr_source[entry]
        #Accuracy_mean += Hist_Bkgd_Optimization[0].GetBinContent(min_bin)
        #Accuracy_weight += 1.
        Accuracy_mean_error += pow(AccuracyErr_source[entry],2)/len(Accuracy_source)
    if Accuracy_weight>0.: Accuracy_mean = Accuracy_mean/Accuracy_weight
    Accuracy_mean_error = pow(Accuracy_mean_error,0.5)


    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyInit_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyInit_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{data}-N_{bkg})/N_{data}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyInit_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+0.5*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    #ax.legend( (rects1[0]), ('Initial $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("output_plots/PerformanceInit_SourceName_E%s_%s%s.png"%(e,method_tag,folder_path))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyInitSigned_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyInitSigned_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{data}-N_{bkg})/N_{data}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyInit_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+0.5*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    #ax.legend( (rects1[0]), ('Initial $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.25)
    plt.savefig("output_plots/PerformanceInitSigned_SourceName_E%s_%s%s.png"%(e,method_tag,folder_path))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyBestPar9_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyBestPar9_source, width, color='#089FFF', yerr=AccuracyBestPar9Err_source)
    rects2 = ax.bar(ind+width, AccuracyInit_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{data}-N_{bkg})/N_{data}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyBestPar9_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+1.0*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax.legend( (rects1[0], rects2[0]), ('Best 9-par $<\epsilon>=%0.3f$'%(AccuracyBestPar9_mean), 'Initial $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("output_plots/PerformanceBest_SourceName_E%s_%s%s.png"%(e,method_tag,folder_path))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyBkgd_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyBkgd_source, width, color='#089FFF', yerr=AccuracyBkgdErr_source)
    rects2 = ax.bar(ind+width, AccuracyInit_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{data}-N_{bkg})/N_{data}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyBkgd_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+1.0*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax.legend( (rects1[0], rects2[0]), ('MIBE $<\epsilon>=%0.3f$'%(AccuracyBkgd_mean), 'ON/OFF $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.25)
    plt.savefig("output_plots/PerformanceMin_SourceName_E%s_%s%s.png"%(e,method_tag,folder_path))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyBkgdSigned_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyBkgdSigned_source, width, color='#089FFF', yerr=AccuracyBkgdErr_source)
    rects2 = ax.bar(ind+width, AccuracyInitSigned_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{data}-N_{bkg})/N_{data}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyBkgd_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+1.0*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax.legend( (rects1[0], rects2[0]), ('MIBE $<\epsilon>=%0.3f$'%(AccuracyBkgd_mean), 'ON/OFF $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.25)
    plt.savefig("output_plots/PerformanceSigned_SourceName_E%s_%s%s.png"%(e,method_tag,folder_path))

    energy_dependent_stat += [AccuracyBkgd_mean_error]
    energy_dependent_syst += [pow(max(0.,pow(AccuracyBkgd_mean,2)-pow(AccuracyBkgd_mean_error,2)),0.5)]
    energy_dependent_syst_rank0 += [pow(max(0.,pow(AccuracyBkgdRank0_mean,2)-pow(AccuracyBkgdRank0_mean_error,2)),0.5)]
    energy_dependent_syst_rank1 += [pow(max(0.,pow(AccuracyBkgdRank1_mean,2)-pow(AccuracyBkgdRank1_mean_error,2)),0.5)]
    energy_dependent_syst_init += [pow(max(0.,pow(AccuracyInit_mean,2)-pow(AccuracyInit_mean_error,2)),0.5)]
    energy_dependent_rank0 += [AccuracyRank0_mean]
    energy_dependent_rank1 += [AccuracyRank1_mean]
    energy_dependent_rank2 += [AccuracyRank2_mean]
    energy_dependent_rank3 += [AccuracyRank3_mean]


    this_roi_energy_dependent_stat_vali = []
    this_roi_energy_dependent_syst_vbkg = []
    this_roi_energy_dependent_syst_rfov = []
    this_roi_energy_dependent_syst_comb = []
    for nth_roi in range(0,number_of_roi):

        ValidateBkgd_source = []
        ValidateBkgdErr_source = []
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            if validate_data_count[nth_roi][entry-1]<10.:
                ValidateBkgd_source += [0.]
                ValidateBkgdErr_source += [0.]
            else:
                ValidateBkgd_source += [abs(validate_data_count[nth_roi][entry-1]-validate_bkgd_count[nth_roi][entry-1])/validate_data_count[nth_roi][entry-1]]
                ValidateBkgdErr_source += [1./pow(validate_data_count[nth_roi][entry-1],0.5)]
        ValidateBkgd_mean, ValidateBkgd_mean_error = SystematicErrorMeasurement(ValidateBkgd_source,ValidateBkgdErr_source)

        ValidateRFoV_source = []
        ValidateRFoVErr_source = []
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            if validate_data_count[nth_roi][entry-1]<10.:
                ValidateRFoV_source += [0.]
                ValidateRFoVErr_source += [0.]
            else:
                ValidateRFoV_source += [abs(validate_data_count[nth_roi][entry-1]-validate_rfov_count[nth_roi][entry-1])/validate_data_count[nth_roi][entry-1]]
                ValidateRFoVErr_source += [1./pow(validate_data_count[nth_roi][entry-1],0.5)]
        ValidateRFoV_mean, ValidateRFoV_mean_error = SystematicErrorMeasurement(ValidateRFoV_source,ValidateRFoVErr_source)

        ValidateComb_source = []
        ValidateCombErr_source = []
        for entry in range(1,len(Hist_Bkgd_Optimization)):
            if validate_data_count[nth_roi][entry-1]<10.:
                ValidateComb_source += [0.]
                ValidateCombErr_source += [0.]
            else:
                ValidateComb_source += [abs(validate_data_count[nth_roi][entry-1]-validate_comb_count[nth_roi][entry-1])/validate_data_count[nth_roi][entry-1]]
                ValidateCombErr_source += [1./pow(validate_data_count[nth_roi][entry-1],0.5)]
        ValidateComb_mean, ValidateComb_mean_error = SystematicErrorMeasurement(ValidateComb_source,ValidateCombErr_source)

        plt.clf()
        ax = fig.add_subplot(111)
        ind = np.arange(len(ValidateBkgd_source))
        width = 0.35
        rects1 = ax.bar(ind, ValidateBkgd_source, width, color='#089FFF', yerr=ValidateBkgdErr_source)
        rects2 = ax.bar(ind+width, ValidateRFoV_source, width, color='#FF9848', yerr=ValidateBkgdErr_source)
        #rects3 = ax.bar(ind+2*width, ValidateComb_source, width, color='green', yerr=ValidateBkgdErr_source)
        ax.set_xlim(-width,len(ind)+width)
        ax.set_ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
        #combined_mean = 1. / (1./ValidateBkgd_mean + 1./ValidateRFoV_mean)
        combined_mean = ValidateComb_mean
        ax.set_title('$combined <\epsilon>=%0.3f \pm %0.3f$'%(combined_mean,ValidateBkgd_mean_error))
        xTickMarks = sample_name
        ax.set_xticks(ind+width)
        xtickNames = ax.set_xticklabels(xTickMarks)
        plt.setp(xtickNames, rotation=45, fontsize=10)
        ax.legend( (rects1[0], rects2[0]), ('MIBE $<\epsilon>=%0.3f$'%(ValidateBkgd_mean), 'FoV method $<\epsilon>=%0.3f$'%(ValidateRFoV_mean)), loc='best' )
        #ax.legend( (rects1[0], rects2[0], rects3[0]), ('LRR $<\epsilon>=%0.3f$'%(ValidateBkgd_mean), 'FoV method $<\epsilon>=%0.3f$'%(ValidateRFoV_mean), 'combined $<\epsilon>=%0.3f$'%(ValidateComb_mean)), loc='best' )
        plt.subplots_adjust(bottom=0.15)
        plt.savefig("output_plots/PerformanceRFoV_SourceName_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

        this_roi_energy_dependent_stat_vali += [ValidateBkgd_mean_error]
        this_roi_energy_dependent_syst_vbkg += [ValidateBkgd_mean]
        this_roi_energy_dependent_syst_rfov += [ValidateRFoV_mean]
        this_roi_energy_dependent_syst_comb += [ValidateComb_mean]
        #this_roi_energy_dependent_syst_vbkg += [pow(max(0.,pow(ValidateBkgd_mean,2)-pow(ValidateBkgd_mean_error,2)),0.5)]
        #this_roi_energy_dependent_syst_rfov += [pow(max(0.,pow(ValidateRFoV_mean,2)-pow(ValidateBkgd_mean_error,2)),0.5)]
        #this_roi_energy_dependent_syst_comb += [pow(max(0.,pow(ValidateComb_mean,2)-pow(ValidateBkgd_mean_error,2)),0.5)]

    energy_dependent_stat_vali += [this_roi_energy_dependent_stat_vali]
    energy_dependent_syst_vbkg += [this_roi_energy_dependent_syst_vbkg]
    energy_dependent_syst_rfov += [this_roi_energy_dependent_syst_rfov]
    energy_dependent_syst_comb += [this_roi_energy_dependent_syst_comb]


    #RankCounts = []
    #Rank = [1,2,3,4,5]
    #plt.clf()
    #plt.xlabel("Rank (r)", fontsize=18)
    #plt.ylabel("$abs(N_{\gamma bkg}-N^{(r)}_{\gamma bkg})/N_{\gamma bkg}$", fontsize=18)
    #plt.yscale('log')
    #plt.xlim(0,6)
    #for s in range(0,len(sample_list)):
    #    if data_count[s]==0.: continue
    #    RankCounts = []
    #    RankCounts += [abs(data_count[s]-rank0_count[s])/data_count[s]]
    #    RankCounts += [abs(data_count[s]-rank1_count[s])/data_count[s]]
    #    RankCounts += [abs(data_count[s]-rank2_count[s])/data_count[s]]
    #    RankCounts += [abs(data_count[s]-rank3_count[s])/data_count[s]]
    #    RankCounts += [abs(data_count[s]-rank4_count[s])/data_count[s]]
    #    plt.errorbar(Rank,RankCounts,fmt='o')
    #plt.savefig("output_plots/RankCounts_E%s%s.png"%(e,folder_path))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Bkgd_Optimization[0]]
    legends += ['average']
    colors += [1]
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        Hists += [Hist_Bkgd_Optimization[entry]]
        #legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        legends += ['%s'%(sample_name[entry-1])]
        colors += [int(random_gen.Uniform(29.,49.))]
    ax.cla()
    MakeMultiplePlot(ax,Hists,legends,colors,'log10 c','relative error','OptimizationAlpha_E%s%s'%(e,folder_path),0.,0.,False,False)
    fig.savefig("output_plots/OptimizationAlpha_E%s%s.png"%(e,folder_path))


    singularvalue_array = []
    for binx in range(0,Hist_Data_Eigenvalues[0].GetNbinsX()):
        singularvalue_array += [Hist_Data_Eigenvalues[0].GetBinContent(binx+1)]
    energy_dependent_singularvalue += [singularvalue_array]

    #par9_epsilon = []
    #wpar9_epsilon = []
    #for entry in range(0,len(wpar9_count)):
    #    if data_count[entry]==0.: continue
    #    par9_epsilon += [abs(data_count[entry]-par9_count[entry])/data_count[entry]]
    #    wpar9_epsilon += [abs(data_count[entry]-wpar9_count[entry])/data_count[entry]]
    #plt.clf()
    #plt.xlabel("$\epsilon$ (plain Frobenius norm)", fontsize=18)
    #plt.ylabel("$\epsilon$ (weighted Frobenius norm)", fontsize=18)
    #plt.xlim(0.,0.1)
    #plt.ylim(0.,0.1)
    #plt.scatter(par9_epsilon,wpar9_epsilon)
    #line_x = np.arange(0.0, 0.1, 1e-4) # angular size
    #line_y = line_x
    #plt.plot(line_x, line_y, color='r')
    #plt.savefig("output_plots/PlainVsWeightedFrobenius_Count_E%s_%s%s.png"%(e,method_tag,folder_path))

    #par9_epsilon = []
    #wpar9_epsilon = []
    #for entry in range(0,len(wpar9_count)):
    #    par9_epsilon += [par9_chi2[entry]]
    #    wpar9_epsilon += [wpar9_chi2[entry]]
    #plt.clf()
    #plt.xlabel("$\epsilon$ (plain Frobenius norm)", fontsize=18)
    #plt.ylabel("$\epsilon$ (weighted Frobenius norm)", fontsize=18)
    #plt.xlim(0.,0.1)
    #plt.ylim(0.,0.1)
    #plt.scatter(par9_epsilon,wpar9_epsilon)
    #line_x = np.arange(0.0, 0.1, 1e-4) # angular size
    #line_y = line_x
    #plt.plot(line_x, line_y, color='r')
    #plt.savefig("output_plots/PlainVsWeightedFrobenius_Chi2_E%s_%s%s.png"%(e,method_tag,folder_path))

my_table = PrettyTable()
my_table.field_names = ["field","exposure"]
my_table.float_format["exposure"] = ".1"
for entry in range(0,len(data_exposure)):
    my_table.add_row([sample_name[entry],data_exposure[entry]])
print(my_table)

my_table = PrettyTable()
my_table.field_names = ["Syst. err MIBE", "Syst. err init.", "Stat. err"]
my_table.float_format["Syst. err MIBE"] = ".3"
my_table.float_format["Syst. err init."] = ".3"
my_table.float_format["Stat. err"] = ".3"
for entry in range(0,len(energy_dependent_syst)):
    my_table.add_row([energy_dependent_syst[entry],energy_dependent_syst_init[entry],energy_dependent_stat[entry]])
print(my_table)

for nth_roi in range(0,number_of_roi):
    my_table = PrettyTable()
    my_table.field_names = ["Syst. err MIBE(2)", "Syst. err RFoV", "Syst. err MIBE(2)+RFoV", "Stat. err(2)"]
    my_table.float_format["Syst. err MIBE(2)"] = ".3"
    my_table.float_format["Syst. err RFoV"] = ".3"
    my_table.float_format["Syst. err MIBE(2)+RFoV"] = ".3"
    my_table.float_format["Stat. err(2)"] = ".3"
    for entry in range(0,len(energy_dependent_syst)):
        my_table.add_row([energy_dependent_syst_vbkg[entry][nth_roi],energy_dependent_syst_rfov[entry][nth_roi],energy_dependent_syst_comb[entry][nth_roi],energy_dependent_stat_vali[entry][nth_roi]])
    print(my_table)

energy_dependent_stat_vali = np.array(energy_dependent_stat_vali)
energy_dependent_syst_vbkg = np.array(energy_dependent_syst_vbkg)
energy_dependent_syst_rfov = np.array(energy_dependent_syst_rfov)
energy_dependent_syst_comb = np.array(energy_dependent_syst_comb)
energy_array = np.array(energy_bin)
for nth_roi in range(0,number_of_roi):
    plt.clf()
    fig, ax = plt.subplots()
    plt.xlabel("energy [GeV]", fontsize=18)
    plt.ylabel("systematic error", fontsize=18)
    plt.xscale('log')
    plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_syst_vbkg[:,nth_roi], color='b', label='MIBE')
    plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_syst_vbkg[:,nth_roi]-energy_dependent_stat_vali[:,nth_roi], energy_dependent_syst_vbkg[:,nth_roi]+energy_dependent_stat_vali[:,nth_roi], alpha=0.1, color='b')
    plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rfov[:,nth_roi], color='r', label='OFF region')
    plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rfov[:,nth_roi]-energy_dependent_stat_vali[:,nth_roi], energy_dependent_syst_rfov[:,nth_roi]+energy_dependent_stat_vali[:,nth_roi], alpha=0.1, color='r')
    ax.legend(loc='best')
    plt.savefig("output_plots/MIBE_vs_RFoV_SystErr_RoI%s.png"%(nth_roi))

energy_dependent_stat_array = np.array(energy_dependent_stat)
energy_dependent_init_array = np.array(energy_dependent_syst_init)
energy_dependent_syst_array = np.array(energy_dependent_syst)
energy_dependent_syst_rank0_array = np.array(energy_dependent_syst_rank0)
energy_dependent_syst_rank1_array = np.array(energy_dependent_syst_rank1)
plt.clf()
fig, ax = plt.subplots()
plt.xlabel("energy [GeV]", fontsize=18)
plt.ylabel("$\sigma_{\\alpha}$", fontsize=18)
plt.xscale('log')
plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_syst_array, color='b', label='m = d')
plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_syst_array-energy_dependent_stat_array, energy_dependent_syst_array+energy_dependent_stat_array, alpha=0.1, color='b')
plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rank1_array, color='g', label='m = min(2,d)')
plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rank1_array-energy_dependent_stat_array, energy_dependent_syst_rank1_array+energy_dependent_stat_array, alpha=0.1, color='g')
plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rank0_array, color='r', label='m = min(1,d)')
plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_syst_rank0_array-energy_dependent_stat_array, energy_dependent_syst_rank0_array+energy_dependent_stat_array, alpha=0.1, color='r')
#plt.plot(energy_array[0:len(energy_bin)-1], energy_dependent_init_array, color='m', label='initial')
#plt.fill_between(energy_array[0:len(energy_bin)-1], energy_dependent_init_array-energy_dependent_stat_array, energy_dependent_init_array+energy_dependent_stat_array, alpha=0.1, color='m')
ax.legend(loc='best')
plt.savefig("output_plots/MIBE_vs_DiffRanks_SystErr.png")

plt.clf()
fig, ax = plt.subplots()
plt.xlabel("rank $n$", fontsize=18)
plt.ylabel("singular value $\sigma_{n}$", fontsize=18)
plt.yscale('log')
for entry in range(0,len(energy_dependent_singularvalue)):
    plt.plot(energy_dependent_singularvalue[entry],marker='.',label='%s-%s GeV'%(energy_bin[entry],energy_bin[entry+1]))
ax.legend(loc='best')
plt.savefig("output_plots/MatrixSingularValue.png")

plt.clf()
fig, ax = plt.subplots()
plt.xlabel("Energy [GeV]", fontsize=18)
plt.ylabel("relative error", fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank0,marker='.',label='$n \leq 0$')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank1,marker='.',label='$n \leq 1$')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank2,marker='.',label='$n \leq 2$')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank3,marker='.',label='$n \leq 3$')
ax.legend(loc='best')
plt.savefig("output_plots/RankErrors.png")

Hists = []
legends = []
colors = []
Hists += [Hist_SystErrDist_MDM]
legends += ['Refined']
colors += [1]
Hists += [Hist_SystErrDist_Init]
legends += ['Init.']
colors += [2]
ax.cla()
MakeMultipleFitPlot(ax,Hists,legends,colors,'relative error','number of measurements')
fig.savefig("output_plots/SystErrDist.png")

for elev in range(0,len(elev_bins)-1):
    energy_interval = []
    for ebin in range(0,len(energy_bin)-1):
        energy_interval += ["%s-%s"%(energy_bin[ebin],energy_bin[ebin+1])]
        energy_dependent_syst[ebin] = Hist_NormSystErr[elev].GetBinContent(ebin+1)
        energy_dependent_syst_init[ebin] = Hist_NormSystInitErr[elev].GetBinContent(ebin+1)
        energy_dependent_stat[ebin] = Hist_NormStatErr[elev].GetBinContent(ebin+1)
    my_table = PrettyTable()
    my_table.field_names = ["energy","Syst. err init.", "Syst. err refined", "Stat. err"]
    my_table.float_format["Syst. err init."] = ".3"
    my_table.float_format["Syst. err refined"] = ".3"
    my_table.float_format["Stat. err"] = ".3"
    for entry in range(0,len(energy_dependent_syst)):
        my_table.add_row([energy_interval[entry],energy_dependent_syst_init[entry],energy_dependent_syst[entry],energy_dependent_stat[entry]])
    print(my_table)

OutputFile = ROOT.TFile('output_plots/SystErrors.root','recreate')
for elev in range(0,len(elev_bins)-1):
    Hist_NormSystErr[elev].Write()
    for ebin in range(0,len(energy_bin)-1):
        Hist_ShapeSystErr[elev][ebin].Write()
OutputFile.Close()