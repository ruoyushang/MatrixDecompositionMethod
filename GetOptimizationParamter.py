
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

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

energy_bin_cut_low = 3
energy_bin_cut_up = 4

#N_bins_for_deconv = 20
N_bins_for_deconv = 16
#N_bins_for_deconv = 12
#N_bins_for_deconv = 8
ErecS_lower_cut = 0
ErecS_upper_cut = 0
total_exposure_hours = 0.
Zenith_mean_data = []
NSB_mean_data = []
Zenith_RMS_data = []
NSB_RMS_data = []
Zenith_mean_dark = []
NSB_mean_dark = []
Zenith_RMS_dark = []
NSB_RMS_dark = []
Accuracy_source = []
AccuracyErr_source = []
data_exposure = []
validate_data_count = []
validate_bkgd_count = []
validate_rfov_count = []
validate_comb_count = []
data_count = []
bkgd_count = []
par8_count = []
par9_count = []
wpar9_count = []
dark_count = []
bkgd_chi2 = []
par8_chi2 = []
par9_chi2 = []
wpar9_chi2 = []
dark_chi2 = []
rank0_count = []
rank1_count = []
rank2_count = []
rank3_count = []
rank4_count = []

#folder_path = 'output_nominal'
#folder_path = 'output_elev_p5'
folder_path = 'output_nsb_m1'
method_tag = 'tight_mdm_default'
#method_tag = 'tight_mdm_rank3'
#method_tag = 'tight_mdm_rank5'
#method_tag = 'tight_mdm_tikhonov'
#method_tag = 'tight_mdm_cutoff'
#method_tag = 'tight_mdm_cutoff_eigen'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
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

#sample_list += ['CrabV5_OFF']
#sample_name += ['Crab V5']
#sample_list += ['CrabV6_OFF']
#sample_name += ['Crab V6']
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
    
#elev_bins = [45,85]
#elev_bins = [40,70]
#elev_bins = [60,80]
#elev_bins = [40,60]
elev_bins = [50,60,70,80]
#elev_bins = [70,80]
#elev_bins = [60,70]
#elev_bins = [50,60]
#elev_bins = [40,50]
lowrank_tag += 'elev_incl'

theta2_bins = [0,4]

energy_bin = []
energy_bin += [int(pow(10,2.0))]
energy_bin += [int(pow(10,2.33))]
energy_bin += [int(pow(10,2.66))]
energy_bin += [int(pow(10,3.0))]
energy_bin += [int(pow(10,3.33))]
energy_bin += [int(pow(10,3.66))]
energy_bin += [int(pow(10,4.0))]

root_file_tags = []
mjd_tag = []
mjd_tag += ['']
for elev in range(0,len(elev_bins)-1):
    elev_tag = '_TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

def Make2DPlot(Hist2D,title_x,title_y,name,logz):

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
    Hist2D.SetMaximum(0.1)
    Hist2D.SetMinimum(0.)
    Hist2D.Draw("COL4Z")
    Hist2D.Draw("TEXT45 same")
    canvas.SaveAs('output_plots/%s.png'%(name))

def MakeMultiplePlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetLeftMargin(0.2)
    pad3.SetBorderMode(1)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.8)
    pad2.SetBottomMargin(0.2)
    pad2.SetLeftMargin(0.2)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad2.SetGrid()
    pad2.Draw()
    pad3.Draw()

    pad2.cd()
    if logy: pad2.SetLogy()

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

    #gap = 0.1*(max_heigh-min_heigh)
    #    Hists[0].Draw("E")
    #else:
    #    if not logy:
    #        Hists[max_hist].SetMaximum(max_heigh+gap)
    #        Hists[max_hist].SetMinimum(min_heigh-gap)
    #    Hists[max_hist].Draw("E")
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
    #if not logy:
    #    Hists[0].SetMaximum(1e-1)
    #    #Hists[0].SetMinimum(5e-3)
    Hists[0].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(4)
        Hists[h].Draw("E same")
    Hists[0].SetLineWidth(6)
    Hists[0].Draw("E same")

    if Hists[0].GetXaxis().GetLabelOffset()==999:
        x1 = Hists[0].GetXaxis().GetXmin()
        x2 = Hists[0].GetXaxis().GetXmax()
        y1 = Hists[0].GetYaxis().GetXmin()
        y2 = Hists[0].GetYaxis().GetXmax()
        IncValues = ROOT.TF1( "IncValues", "x", 0, 256 )
        raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
        raLowerAxis.SetLabelSize(Hists[0].GetXaxis().GetLabelSize())
        raLowerAxis.Draw()

    pad3.cd()
    legend = ROOT.TLegend(0.2,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.SetNColumns(3)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    min_y = 1.0
    min_bin = 0
    for binx in range(1,Hists[0].GetNbinsX()+1):
        if Hists[0].GetBinContent(binx)<min_y:
            min_y = Hists[0].GetBinContent(binx)
            min_bin = binx

    #lumilab1 = ROOT.TLatex(0.15,0.80,'log10 #alpha = %.1f'%(Hists[0].GetBinCenter(min_bin)) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def GetHistogramsFromFile(FilePath,which_source):
    global total_exposure_hours
    global Zenith_mean_data
    global NSB_mean_data
    global Zenith_RMS_data
    global NSB_RMS_data
    global Zenith_mean_dark
    global NSB_mean_dark
    global Zenith_RMS_dark
    global NSB_RMS_dark
    global data_exposure
    global validate_data_count
    global validate_bkgd_count
    global validate_rfov_count
    global validate_comb_count
    global data_count
    global bkgd_count
    global par8_count
    global par9_count
    global wpar9_count
    global dark_count
    global bkgd_chi2
    global par8_chi2
    global par9_chi2
    global wpar9_chi2
    global dark_chi2
    global rank0_count
    global rank1_count
    global rank2_count
    global rank3_count
    global rank4_count
    data_validate_count = ROOT.std.vector("double")(10)
    bkgd_validate_count = ROOT.std.vector("double")(10)
    rfov_validate_count = ROOT.std.vector("double")(10)
    comb_validate_count = ROOT.std.vector("double")(10)
    data_gamma_count = ROOT.std.vector("double")(10)
    bkgd_gamma_count = ROOT.std.vector("double")(10)
    par8_gamma_count = ROOT.std.vector("double")(10)
    par9_gamma_count = ROOT.std.vector("double")(10)
    wpar9_gamma_count = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)
    bkgd_coeff_chi2 = ROOT.std.vector("double")(10)
    par8_coeff_chi2 = ROOT.std.vector("double")(10)
    par9_coeff_chi2 = ROOT.std.vector("double")(10)
    wpar9_coeff_chi2 = ROOT.std.vector("double")(10)
    dark_coeff_chi2 = ROOT.std.vector("double")(10)
    rank0_gamma_count = ROOT.std.vector("double")(10)
    rank1_gamma_count = ROOT.std.vector("double")(10)
    rank2_gamma_count = ROOT.std.vector("double")(10)
    rank3_gamma_count = ROOT.std.vector("double")(10)
    rank4_gamma_count = ROOT.std.vector("double")(10)
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    NewInfoTree = InputFile.Get("NewInfoTree")
    InfoTree.GetEntry(0)
    NewInfoTree.SetBranchAddress('data_validate_count',ROOT.AddressOf(data_validate_count))
    NewInfoTree.SetBranchAddress('bkgd_validate_count',ROOT.AddressOf(bkgd_validate_count))
    NewInfoTree.SetBranchAddress('rfov_validate_count',ROOT.AddressOf(rfov_validate_count))
    NewInfoTree.SetBranchAddress('comb_validate_count',ROOT.AddressOf(comb_validate_count))
    NewInfoTree.SetBranchAddress('data_gamma_count',ROOT.AddressOf(data_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_gamma_count',ROOT.AddressOf(bkgd_gamma_count))
    NewInfoTree.SetBranchAddress('par8_gamma_count',ROOT.AddressOf(par8_gamma_count))
    NewInfoTree.SetBranchAddress('par9_gamma_count',ROOT.AddressOf(par9_gamma_count))
    NewInfoTree.SetBranchAddress('wpar9_gamma_count',ROOT.AddressOf(wpar9_gamma_count))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_coeff_chi2',ROOT.AddressOf(bkgd_coeff_chi2))
    NewInfoTree.SetBranchAddress('par8_coeff_chi2',ROOT.AddressOf(par8_coeff_chi2))
    NewInfoTree.SetBranchAddress('par9_coeff_chi2',ROOT.AddressOf(par9_coeff_chi2))
    NewInfoTree.SetBranchAddress('wpar9_coeff_chi2',ROOT.AddressOf(wpar9_coeff_chi2))
    NewInfoTree.SetBranchAddress('dark_coeff_chi2',ROOT.AddressOf(dark_coeff_chi2))
    NewInfoTree.SetBranchAddress('rank0_gamma_count',ROOT.AddressOf(rank0_gamma_count))
    NewInfoTree.SetBranchAddress('rank1_gamma_count',ROOT.AddressOf(rank1_gamma_count))
    NewInfoTree.SetBranchAddress('rank2_gamma_count',ROOT.AddressOf(rank2_gamma_count))
    NewInfoTree.SetBranchAddress('rank3_gamma_count',ROOT.AddressOf(rank3_gamma_count))
    NewInfoTree.SetBranchAddress('rank4_gamma_count',ROOT.AddressOf(rank4_gamma_count))
    NewInfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours
    data_exposure[which_source] += exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    if energy_index==0: 
        total_exposure_hours += exposure_hours
    Zenith_mean = NewInfoTree.Zenith_mean_data
    Zenith_RMS = NewInfoTree.Zenith_RMS_data
    Zenith_mean = NewInfoTree.Zenith_mean_data-NewInfoTree.Zenith_mean_dark
    Zenith_RMS = NewInfoTree.Zenith_RMS_dark
    Zenith_mean_data[which_source] += Zenith_mean
    Zenith_RMS_data[which_source] += Zenith_RMS/2.
    Zenith_mean_dark[which_source] += Zenith_mean
    Zenith_RMS_dark[which_source] += Zenith_RMS/2.
    NSB_mean = NewInfoTree.NSB_mean_data
    NSB_RMS = NewInfoTree.NSB_RMS_data
    NSB_mean = NewInfoTree.NSB_mean_data-NewInfoTree.NSB_mean_dark
    NSB_RMS = NewInfoTree.NSB_RMS_dark
    NSB_mean_data[which_source] += NSB_mean
    NSB_RMS_data[which_source] += NSB_RMS/2.
    NSB_mean_dark[which_source] += NSB_mean
    NSB_RMS_dark[which_source] += NSB_RMS/2.
    validate_data_count[which_source]  += data_validate_count[energy_index]
    validate_bkgd_count[which_source]  += bkgd_validate_count[energy_index]
    validate_rfov_count[which_source]  += rfov_validate_count[energy_index]
    validate_comb_count[which_source]  += comb_validate_count[energy_index]
    data_count[which_source]  += data_gamma_count[energy_index]
    bkgd_count[which_source]  += bkgd_gamma_count[energy_index]
    par8_count[which_source]  += par8_gamma_count[energy_index]
    par9_count[which_source]  += par9_gamma_count[energy_index]
    wpar9_count[which_source]  += wpar9_gamma_count[energy_index]
    dark_count[which_source]  += dark_gamma_count[energy_index]
    bkgd_chi2[which_source]  += bkgd_coeff_chi2[energy_index]
    par8_chi2[which_source]  += par8_coeff_chi2[energy_index]
    par9_chi2[which_source]  += par9_coeff_chi2[energy_index]
    wpar9_chi2[which_source]  += wpar9_coeff_chi2[energy_index]
    dark_chi2[which_source]  += dark_coeff_chi2[energy_index]
    rank0_count[which_source]  += rank0_gamma_count[energy_index]
    rank1_count[which_source]  += rank1_gamma_count[energy_index]
    rank2_count[which_source]  += rank2_gamma_count[energy_index]
    rank3_count[which_source]  += rank3_gamma_count[energy_index]
    rank4_count[which_source]  += rank4_gamma_count[energy_index]
    Hist_Bkgd_Optimization[which_source+1].Reset()
    Hist_Bkgd_Optimization_beta[which_source+1].Reset()
    Hist_Bkgd_OptimizationChi2_beta[which_source+1].Reset()
    Hist_Dark_Optimization[which_source+1].Reset()
    weight = 1./float(len(sample_list))
    #weight = exposure_hours
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Optimization[which_source+1].Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Optimization_beta_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization_beta[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Optimization_beta[which_source+1].Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_OptimizationChi2_beta_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_OptimizationChi2_beta[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_OptimizationChi2_beta[which_source+1].Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Dark_Optimization[0].Add(InputFile.Get(HistName),weight)
    Hist_Dark_Optimization[which_source+1].Add(InputFile.Get(HistName))
    Hist_Bkgd_Chi2[which_source+1].Reset()
    HistName = "Hist_Bkgd_Chi2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Chi2[which_source+1].Add(InputFile.Get(HistName))
    Hist_VVV_Eigenvalues[which_source+1].Reset()
    HistName = "Hist_VVV_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_VVV_Eigenvalues[0].Add(InputFile.Get(HistName),weight)
    Hist_VVV_Eigenvalues[which_source+1].Add(InputFile.Get(HistName))
    Hist_Data_Eigenvalues[which_source+1].Reset()
    HistName = "Hist_Data_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Data_Eigenvalues[0].Add(InputFile.Get(HistName),weight)
    Hist_Data_Eigenvalues[which_source+1].Add(InputFile.Get(HistName))
    Hist_U_Proj[which_source+1].Reset()
    HistName = "Hist_U_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_U_Proj[0].Add(InputFile.Get(HistName),weight)
    Hist_U_Proj[which_source+1].Add(InputFile.Get(HistName))
    Hist_V_Proj[which_source+1].Reset()
    HistName = "Hist_V_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_V_Proj[0].Add(InputFile.Get(HistName),weight)
    Hist_V_Proj[which_source+1].Add(InputFile.Get(HistName))
    Hist_GammaRegion_Contribution[which_source+1].Reset()
    HistName = "Hist_GammaRegion_Contribution_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_GammaRegion_Contribution[0].Add(InputFile.Get(HistName),weight)
    Hist_GammaRegion_Contribution[which_source+1].Add(InputFile.Get(HistName))

optimiz_lower = -5.
optimiz_upper = -3.
Hist_Bkgd_Optimization = []
Hist_Bkgd_Optimization_beta = []
Hist_Bkgd_OptimizationChi2_beta = []
Hist_Dark_Optimization = []
for e in range(0,len(sample_list)+1):
    Hist_Bkgd_Optimization += [ROOT.TH1D("Hist_Bkgd_Optimization_%s"%(e),"",10,optimiz_lower,optimiz_upper)]
    Hist_Bkgd_Optimization_beta += [ROOT.TH2D("Hist_Bkgd_Optimization_beta_%s"%(e),"",10,0.,1.,10,0.,1.)]
    Hist_Bkgd_OptimizationChi2_beta += [ROOT.TH2D("Hist_Bkgd_OptimizationChi2_beta_%s"%(e),"",10,0.,1.,10,0.,1.)]
    Hist_Dark_Optimization += [ROOT.TH1D("Hist_Dark_Optimization_%s"%(e),"",int(N_bins_for_deconv/2),0,N_bins_for_deconv/2)]
Hist_Bkgd_Chi2 = []
for e in range(0,len(sample_list)+1):
    Hist_Bkgd_Chi2 += [ROOT.TH1D("Hist_Bkgd_Chi2_%s"%(e),"",10,optimiz_lower,optimiz_upper)]
Hist_VVV_Eigenvalues = []
for e in range(0,len(sample_list)+1):
    Hist_VVV_Eigenvalues += [ROOT.TH1D("Hist_VVV_Eigenvalues_%s"%(e),"",N_bins_for_deconv*N_bins_for_deconv,0,N_bins_for_deconv*N_bins_for_deconv)]
Hist_Data_Eigenvalues = []
for e in range(0,len(sample_list)+1):
    Hist_Data_Eigenvalues += [ROOT.TH1D("Hist_Data_Eigenvalues_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv)]
Hist_U_Proj = []
Hist_V_Proj = []
Hist_GammaRegion_Contribution = []
for e in range(0,len(sample_list)+1):
    Hist_U_Proj += [ROOT.TH2D("Hist_U_Proj_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
    Hist_V_Proj += [ROOT.TH2D("Hist_V_Proj_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
    Hist_GammaRegion_Contribution += [ROOT.TH1D("Hist_GammaRegion_Contribution_%s"%(e),"",N_bins_for_deconv,0,N_bins_for_deconv)]

for e in range(0,len(energy_bin)-1):
    FilePath_List = []
    Zenith_mean_data = [0.]*len(sample_list)
    NSB_mean_data = [0.]*len(sample_list)
    Zenith_RMS_data = [0.]*len(sample_list)
    NSB_RMS_data = [0.]*len(sample_list)
    Zenith_mean_dark = [0.]*len(sample_list)
    NSB_mean_dark = [0.]*len(sample_list)
    Zenith_RMS_dark = [0.]*len(sample_list)
    NSB_RMS_dark = [0.]*len(sample_list)
    validate_data_count = [0.]*len(sample_list)
    validate_bkgd_count = [0.]*len(sample_list)
    validate_rfov_count = [0.]*len(sample_list)
    validate_comb_count = [0.]*len(sample_list)
    data_count = [0.]*len(sample_list)
    bkgd_count = [0.]*len(sample_list)
    par8_count = [0.]*len(sample_list)
    par9_count = [0.]*len(sample_list)
    wpar9_count = [0.]*len(sample_list)
    bkgd_chi2 = [0.]*len(sample_list)
    par8_chi2 = [0.]*len(sample_list)
    par9_chi2 = [0.]*len(sample_list)
    wpar9_chi2 = [0.]*len(sample_list)
    data_exposure = [0.]*len(sample_list)
    dark_count = [0.]*len(sample_list)
    dark_chi2 = [0.]*len(sample_list)
    rank0_count = [0.]*len(sample_list)
    rank1_count = [0.]*len(sample_list)
    rank2_count = [0.]*len(sample_list)
    rank3_count = [0.]*len(sample_list)
    rank4_count = [0.]*len(sample_list)
    for entry in range(0,len(Hist_U_Proj)):
        Hist_U_Proj[entry].Reset()
        Hist_V_Proj[entry].Reset()
        Hist_GammaRegion_Contribution[entry].Reset()
    for entry in range(0,len(Hist_Bkgd_Optimization)):
        Hist_Bkgd_Optimization[entry].Reset()
        Hist_Bkgd_Optimization_beta[entry].Reset()
        Hist_Bkgd_OptimizationChi2_beta[entry].Reset()
        Hist_Dark_Optimization[entry].Reset()
    for entry in range(0,len(Hist_Bkgd_Chi2)):
        Hist_Bkgd_Chi2[entry].Reset()
    for entry in range(0,len(Hist_VVV_Eigenvalues)):
        Hist_VVV_Eigenvalues[entry].Reset()
    for entry in range(0,len(Hist_Data_Eigenvalues)):
        Hist_Data_Eigenvalues[entry].Reset()
    for source in range(0,len(sample_list)):
        source_name = sample_list[source]
        nfiles_used = 0
        for elev in range(0,len(root_file_tags)):
            FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            print ('Trying to read file %s'%(FilePath_List[len(FilePath_List)-1]))
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
            nfiles_used += 1
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1],source)
        if nfiles_used>0:
            Zenith_mean_data[source] = Zenith_mean_data[source]/nfiles_used
            NSB_mean_data[source] = NSB_mean_data[source]/nfiles_used
            Zenith_RMS_data[source] = Zenith_RMS_data[source]/nfiles_used
            NSB_RMS_data[source] = NSB_RMS_data[source]/nfiles_used
            Zenith_mean_dark[source] = Zenith_mean_dark[source]/nfiles_used
            NSB_mean_dark[source] = NSB_mean_dark[source]/nfiles_used
            Zenith_RMS_dark[source] = Zenith_RMS_dark[source]/nfiles_used
            NSB_RMS_dark[source] = NSB_RMS_dark[source]/nfiles_used

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

    for binx in range(1,Hist_Bkgd_Optimization_beta[0].GetNbinsX()+1):
        for biny in range(1,Hist_Bkgd_Optimization_beta[0].GetNbinsY()+1):
            total_weight = 0.
            y_content = 0.
            for entry in range(1,len(Hist_Bkgd_Optimization_beta)):
                weight = pow(data_count[entry-1],0.5)
                #weight = 1.
                total_weight += weight
                y_content += weight*Hist_Bkgd_Optimization_beta[entry].GetBinContent(binx,biny)
            if total_weight>0.: y_content = y_content/total_weight
            Hist_Bkgd_Optimization_beta[0].SetBinContent(binx,biny,y_content)
    for binx in range(1,Hist_Bkgd_OptimizationChi2_beta[0].GetNbinsX()+1):
        for biny in range(1,Hist_Bkgd_OptimizationChi2_beta[0].GetNbinsY()+1):
            total_weight = 0.
            y_content = 0.
            for entry in range(1,len(Hist_Bkgd_OptimizationChi2_beta)):
                weight = pow(data_count[entry-1],0.5)
                #weight = 1.
                total_weight += weight
                y_content += weight*Hist_Bkgd_OptimizationChi2_beta[entry].GetBinContent(binx,biny)
            if total_weight>0.: y_content = y_content/total_weight
            Hist_Bkgd_OptimizationChi2_beta[0].SetBinContent(binx,biny,y_content)

    for binx in range(1,Hist_Dark_Optimization[0].GetNbinsX()+1):
        total_weight = 0.
        y_content = 0.
        for entry in range(1,len(Hist_Dark_Optimization)):
            weight = pow(data_count[entry-1],0.5)
            #weight = 1.
            total_weight += weight
            y_content += weight*Hist_Dark_Optimization[entry].GetBinContent(binx)
        if total_weight>0.: y_content = y_content/total_weight
        Hist_Dark_Optimization[0].SetBinContent(binx,y_content)

    fig, ax = plt.subplots()
    colors = np.random.rand(len(Zenith_mean_data))

    AccuracyInit_source = []
    AccuracyInitErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyInit_source += [0.]
            AccuracyInitErr_source += [0.]
        else:
            AccuracyInit_source += [abs(data_count[entry-1]-dark_count[entry-1])/data_count[entry-1]]
            AccuracyInitErr_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyInit_mean = 0.
    AccuracyInit_weight = 0.
    AccuracyInit_mean_error = 0.
    for entry in range(0,len(AccuracyInit_source)):
        if AccuracyInitErr_source[entry]==0.: continue
        AccuracyInit_mean += 1./AccuracyInitErr_source[entry]*AccuracyInit_source[entry]
        AccuracyInit_weight += 1./AccuracyInitErr_source[entry]
        AccuracyInit_mean_error += pow(AccuracyInitErr_source[entry],2)/len(AccuracyInit_source)
    if AccuracyInit_weight>0.: AccuracyInit_mean = AccuracyInit_mean/AccuracyInit_weight
    AccuracyInit_mean_error = pow(AccuracyInit_mean_error,0.5)

    AccuracyStat_source = []
    AccuracyStatErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
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

    AccuracyRank2_source = []
    AccuracyRank2Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyRank2_source += [0.]
            AccuracyRank2Err_source += [0.]
        else:
            AccuracyRank2_source += [abs(data_count[entry-1]-rank2_count[entry-1])/data_count[entry-1]]
            AccuracyRank2Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyRank2_mean = 0.
    AccuracyRank2_weight = 0.
    AccuracyRank2_mean_error = 0.
    for entry in range(0,len(AccuracyRank2_source)):
        if AccuracyRank2Err_source[entry]==0.: continue
        AccuracyRank2_mean += 1./AccuracyRank2Err_source[entry]*AccuracyRank2_source[entry]
        AccuracyRank2_weight += 1./AccuracyRank2Err_source[entry]
        AccuracyRank2_mean_error += pow(AccuracyRank2Err_source[entry],2)/len(AccuracyRank2_source)
    if AccuracyRank2_weight>0.: AccuracyRank2_mean = AccuracyRank2_mean/AccuracyRank2_weight
    AccuracyRank2_mean_error = pow(AccuracyRank2_mean_error,0.5)

    AccuracyBestPar9_source = []
    AccuracyBestPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyBestPar9_source += [0.]
            AccuracyBestPar9Err_source += [0.]
        else:
            AccuracyBestPar9_source += [pow(Hist_GammaRegion_Contribution[entry].GetBinContent(2+1),0.5)]
            AccuracyBestPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBestPar9_mean = 0.
    AccuracyBestPar9_weight = 0.
    AccuracyBestPar9_mean_error = 0.
    for entry in range(0,len(AccuracyBestPar9_source)):
        if AccuracyBestPar9Err_source[entry]==0.: continue
        AccuracyBestPar9_mean += 1./AccuracyBestPar9Err_source[entry]*AccuracyBestPar9_source[entry]
        AccuracyBestPar9_weight += 1./AccuracyBestPar9Err_source[entry]
        AccuracyBestPar9_mean_error += pow(AccuracyBestPar9Err_source[entry],2)/len(AccuracyBestPar9_source)
    if AccuracyBestPar9_weight>0.: AccuracyBestPar9_mean = AccuracyBestPar9_mean/AccuracyBestPar9_weight
    AccuracyBestPar9_mean_error = pow(AccuracyBestPar9_mean_error,0.5)

    AccuracyPar9_source = []
    AccuracyPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyPar9_source += [0.]
            AccuracyPar9Err_source += [0.]
        else:
            AccuracyPar9_source += [abs(data_count[entry-1]-par9_count[entry-1])/data_count[entry-1]]
            AccuracyPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyPar9_mean = 0.
    AccuracyPar9_weight = 0.
    AccuracyPar9_mean_error = 0.
    for entry in range(0,len(AccuracyPar9_source)):
        if AccuracyPar9Err_source[entry]==0.: continue
        AccuracyPar9_mean += 1./AccuracyPar9Err_source[entry]*AccuracyPar9_source[entry]
        AccuracyPar9_weight += 1./AccuracyPar9Err_source[entry]
        AccuracyPar9_mean_error += pow(AccuracyPar9Err_source[entry],2)/len(AccuracyPar9_source)
    if AccuracyPar9_weight>0.: AccuracyPar9_mean = AccuracyPar9_mean/AccuracyPar9_weight
    AccuracyPar9_mean_error = pow(AccuracyPar9_mean_error,0.5)

    AccuracyWPar9_source = []
    AccuracyWPar9Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyWPar9_source += [0.]
            AccuracyWPar9Err_source += [0.]
        else:
            AccuracyWPar9_source += [abs(data_count[entry-1]-wpar9_count[entry-1])/data_count[entry-1]]
            AccuracyWPar9Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyWPar9_mean = 0.
    AccuracyWPar9_weight = 0.
    AccuracyWPar9_mean_error = 0.
    for entry in range(0,len(AccuracyWPar9_source)):
        if AccuracyWPar9Err_source[entry]==0.: continue
        AccuracyWPar9_mean += 1./AccuracyWPar9Err_source[entry]*AccuracyWPar9_source[entry]
        AccuracyWPar9_weight += 1./AccuracyWPar9Err_source[entry]
        AccuracyWPar9_mean_error += pow(AccuracyWPar9Err_source[entry],2)/len(AccuracyWPar9_source)
    if AccuracyWPar9_weight>0.: AccuracyWPar9_mean = AccuracyWPar9_mean/AccuracyWPar9_weight
    AccuracyWPar9_mean_error = pow(AccuracyWPar9_mean_error,0.5)

    AccuracyBkgd_source = []
    AccuracyBkgdErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyBkgd_source += [0.]
            AccuracyBkgdErr_source += [0.]
        else:
            AccuracyBkgd_source += [abs(data_count[entry-1]-bkgd_count[entry-1])/data_count[entry-1]]
            AccuracyBkgdErr_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyBkgd_mean = 0.
    AccuracyBkgd_weight = 0.
    AccuracyBkgd_mean_error = 0.
    for entry in range(0,len(AccuracyBkgd_source)):
        if AccuracyBkgdErr_source[entry]==0.: continue
        AccuracyBkgd_mean += 1./AccuracyBkgdErr_source[entry]*AccuracyBkgd_source[entry]
        AccuracyBkgd_weight += 1./AccuracyBkgdErr_source[entry]
        AccuracyBkgd_mean_error += pow(AccuracyBkgdErr_source[entry],2)/len(AccuracyBkgd_source)
    if AccuracyBkgd_weight>0.: AccuracyBkgd_mean = AccuracyBkgd_mean/AccuracyBkgd_weight
    AccuracyBkgd_mean_error = pow(AccuracyBkgd_mean_error,0.5)

    ValidateBkgd_source = []
    ValidateBkgdErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if validate_data_count[entry-1]<100.:
            ValidateBkgd_source += [0.]
            ValidateBkgdErr_source += [0.]
        else:
            ValidateBkgd_source += [abs(validate_data_count[entry-1]-validate_bkgd_count[entry-1])/validate_data_count[entry-1]]
            ValidateBkgdErr_source += [1./pow(validate_data_count[entry-1],0.5)]
    ValidateBkgd_mean = 0.
    ValidateBkgd_weight = 0.
    ValidateBkgd_mean_error = 0.
    for entry in range(0,len(ValidateBkgd_source)):
        if ValidateBkgdErr_source[entry]==0.: continue
        ValidateBkgd_mean += 1./ValidateBkgdErr_source[entry]*ValidateBkgd_source[entry]
        ValidateBkgd_weight += 1./ValidateBkgdErr_source[entry]
        ValidateBkgd_mean_error += pow(ValidateBkgdErr_source[entry],2)/len(ValidateBkgd_source)
    if ValidateBkgd_weight>0.: ValidateBkgd_mean = ValidateBkgd_mean/ValidateBkgd_weight
    ValidateBkgd_mean_error = pow(ValidateBkgd_mean_error,0.5)

    ValidateRFoV_source = []
    ValidateRFoVErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if validate_data_count[entry-1]<100.:
            ValidateRFoV_source += [0.]
            ValidateRFoVErr_source += [0.]
        else:
            ValidateRFoV_source += [abs(validate_data_count[entry-1]-validate_rfov_count[entry-1])/validate_data_count[entry-1]]
            ValidateRFoVErr_source += [1./pow(validate_data_count[entry-1],0.5)]
    ValidateRFoV_mean = 0.
    ValidateRFoV_weight = 0.
    ValidateRFoV_mean_error = 0.
    for entry in range(0,len(ValidateRFoV_source)):
        if ValidateRFoVErr_source[entry]==0.: continue
        ValidateRFoV_mean += 1./ValidateRFoVErr_source[entry]*ValidateRFoV_source[entry]
        ValidateRFoV_weight += 1./ValidateRFoVErr_source[entry]
        ValidateRFoV_mean_error += pow(ValidateRFoVErr_source[entry],2)/len(ValidateRFoV_source)
    if ValidateRFoV_weight>0.: ValidateRFoV_mean = ValidateRFoV_mean/ValidateRFoV_weight
    ValidateRFoV_mean_error = pow(ValidateRFoV_mean_error,0.5)

    ValidateComb_source = []
    ValidateCombErr_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if validate_data_count[entry-1]<100.:
            ValidateComb_source += [0.]
            ValidateCombErr_source += [0.]
        else:
            ValidateComb_source += [abs(validate_data_count[entry-1]-validate_comb_count[entry-1])/validate_data_count[entry-1]]
            ValidateCombErr_source += [1./pow(validate_data_count[entry-1],0.5)]
    ValidateComb_mean = 0.
    ValidateComb_weight = 0.
    ValidateComb_mean_error = 0.
    for entry in range(0,len(ValidateComb_source)):
        if ValidateCombErr_source[entry]==0.: continue
        ValidateComb_mean += 1./ValidateCombErr_source[entry]*ValidateComb_source[entry]
        ValidateComb_weight += 1./ValidateCombErr_source[entry]
        ValidateComb_mean_error += pow(ValidateCombErr_source[entry],2)/len(ValidateComb_source)
    if ValidateComb_weight>0.: ValidateComb_mean = ValidateComb_mean/ValidateComb_weight
    ValidateComb_mean_error = pow(ValidateComb_mean_error,0.5)

    AccuracyPar8_source = []
    AccuracyPar8Err_source = []
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        if data_count[entry-1]<100.:
            AccuracyPar8_source += [0.]
            AccuracyPar8Err_source += [0.]
        else:
            AccuracyPar8_source += [abs(data_count[entry-1]-par8_count[entry-1])/data_count[entry-1]]
            AccuracyPar8Err_source += [1./pow(data_count[entry-1],0.5)]
    AccuracyPar8_mean = 0.
    AccuracyPar8_weight = 0.
    AccuracyPar8_mean_error = 0.
    for entry in range(0,len(AccuracyPar8_source)):
        if AccuracyPar8Err_source[entry]==0.: continue
        AccuracyPar8_mean += 1./AccuracyPar8Err_source[entry]*AccuracyPar8_source[entry]
        AccuracyPar8_weight += 1./AccuracyPar8Err_source[entry]
        AccuracyPar8_mean_error += pow(AccuracyPar8Err_source[entry],2)/len(AccuracyPar8_source)
    if AccuracyPar8_weight>0.: AccuracyPar8_mean = AccuracyPar8_mean/AccuracyPar8_weight
    AccuracyPar8_mean_error = pow(AccuracyPar8_mean_error,0.5)

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
    ax.set_ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyInit_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+0.5*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    #ax.legend( (rects1[0]), ('Initial $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("output_plots/PerformanceInit_SourceName_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyBestPar9_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyBestPar9_source, width, color='#089FFF', yerr=AccuracyBestPar9Err_source)
    rects2 = ax.bar(ind+width, AccuracyInit_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyBestPar9_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+1.0*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax.legend( (rects1[0], rects2[0]), ('Best 9-par $<\epsilon>=%0.3f$'%(AccuracyBestPar9_mean), 'Initial $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("output_plots/PerformanceBest_SourceName_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    plt.clf()
    ax = fig.add_subplot(111)
    ind = np.arange(len(AccuracyBkgd_source))
    width = 0.35
    rects1 = ax.bar(ind, AccuracyBkgd_source, width, color='#089FFF', yerr=AccuracyBkgdErr_source)
    rects2 = ax.bar(ind+width, AccuracyInit_source, width, color='#FF9848', yerr=AccuracyInitErr_source)
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    ax.set_title('$<\epsilon>=%0.3f \pm %0.3f$'%(AccuracyBkgd_mean,Accuracy_mean_error))
    xTickMarks = sample_name
    ax.set_xticks(ind+1.0*width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax.legend( (rects1[0], rects2[0]), ('MIBE $<\epsilon>=%0.3f$'%(AccuracyBkgd_mean), 'ON/OFF $<\epsilon>=%0.3f$'%(AccuracyInit_mean)), loc='best' )
    plt.subplots_adjust(bottom=0.25)
    plt.savefig("output_plots/PerformanceMin_SourceName_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    #plt.clf()
    #ax = fig.add_subplot(111)
    #ind = np.arange(len(ValidateBkgd_source))
    #width = 0.35
    #rects1 = ax.bar(ind, ValidateBkgd_source, width, color='#089FFF', yerr=ValidateBkgdErr_source)
    #rects2 = ax.bar(ind+width, ValidateRFoV_source, width, color='#FF9848', yerr=ValidateBkgdErr_source)
    ##rects3 = ax.bar(ind+2*width, ValidateComb_source, width, color='green', yerr=ValidateBkgdErr_source)
    #ax.set_xlim(-width,len(ind)+width)
    #ax.set_ylabel("$abs(N_{\gamma bkg}-N_{model})/N_{\gamma bkg}$", fontsize=18)
    #combined_mean = 1. / (1./ValidateBkgd_mean + 1./ValidateRFoV_mean)
    #ax.set_title('$combined <\epsilon>=%0.3f \pm %0.3f$'%(combined_mean,ValidateBkgd_mean_error))
    #xTickMarks = sample_name
    #ax.set_xticks(ind+width)
    #xtickNames = ax.set_xticklabels(xTickMarks)
    #plt.setp(xtickNames, rotation=45, fontsize=10)
    #ax.legend( (rects1[0], rects2[0]), ('LRR $<\epsilon>=%0.3f$'%(ValidateBkgd_mean), 'FoV method $<\epsilon>=%0.3f$'%(ValidateRFoV_mean)), loc='best' )
    ##ax.legend( (rects1[0], rects2[0], rects3[0]), ('LRR $<\epsilon>=%0.3f$'%(ValidateBkgd_mean), 'FoV method $<\epsilon>=%0.3f$'%(ValidateRFoV_mean), 'combined $<\epsilon>=%0.3f$'%(ValidateComb_mean)), loc='best' )
    #plt.subplots_adjust(bottom=0.15)
    #plt.savefig("output_plots/PerformanceRFoV_SourceName_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    RankCounts = []
    Rank = [1,2,3,4,5]
    plt.clf()
    plt.xlabel("Rank (r)", fontsize=18)
    plt.ylabel("$abs(N_{\gamma bkg}-N^{(r)}_{\gamma bkg})/N_{\gamma bkg}$", fontsize=18)
    plt.yscale('log')
    plt.xlim(0,6)
    for s in range(0,len(sample_list)):
        if data_count[s]==0.: continue
        RankCounts = []
        RankCounts += [abs(data_count[s]-rank0_count[s])/data_count[s]]
        RankCounts += [abs(data_count[s]-rank1_count[s])/data_count[s]]
        RankCounts += [abs(data_count[s]-rank2_count[s])/data_count[s]]
        RankCounts += [abs(data_count[s]-rank3_count[s])/data_count[s]]
        RankCounts += [abs(data_count[s]-rank4_count[s])/data_count[s]]
        plt.errorbar(Rank,RankCounts,fmt='o')
    plt.savefig("output_plots/RankCounts_E%s%s.png"%(e,lowrank_tag))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Bkgd_Optimization[0]]
    legends += ['average']
    colors += [1]
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        Hists += [Hist_Bkgd_Optimization[entry]]
        legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        colors += [int(random_gen.Uniform(29.,49.))]
    MakeMultiplePlot(Hists,legends,colors,'log10 #beta','abs(N_{#gamma bkg}-N_{model})/N_{#gamma bkg}','OptimizationAlpha_E%s%s'%(e,lowrank_tag),1e-3,0.1,False,False)

    Make2DPlot(Hist_Bkgd_Optimization_beta[0],'C_{1}','C_{2}','OptimizationBeta_E%s%s'%(e,lowrank_tag),False)
    Make2DPlot(Hist_Bkgd_OptimizationChi2_beta[0],'C_{1}','C_{2}','OptimizationChi2Beta_E%s%s'%(e,lowrank_tag),False)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Dark_Optimization[0]]
    legends += ['average']
    colors += [1]
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_Dark_Optimization)):
        Hists += [Hist_Dark_Optimization[entry]]
        legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        colors += [int(random_gen.Uniform(29.,49.))]
    MakeMultiplePlot(Hists,legends,colors,'normalization bins','abs(N_{#gamma bkg}-N_{model})/N_{#gamma bkg}','OptimizationNormalization_E%s%s'%(e,lowrank_tag),1e-3,0.1,False,False)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Bkgd_Chi2[0]]
    legends += ['average']
    colors += [1]
    #Hist_Bkgd_Chi2[0].GetXaxis().SetLabelOffset(999)
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_Bkgd_Chi2)):
        Hists += [Hist_Bkgd_Chi2[entry]]
        legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        #colors += [entry+1]
        colors += [int(random_gen.Uniform(29.,49.))]
    MakeMultiplePlot(Hists,legends,colors,'log10 #beta','#chi^{2} = #sum #deltaM_{ij}^{2} in CR','Chi2_E%s%s'%(e,lowrank_tag),pow(10.,-6.5),pow(10.,-4.5),False,False)

    Hists = []
    legends = []
    colors = []
    Hist_VVV_Eigenvalues[0].SetMinimum(1e-3)
    Hists += [Hist_VVV_Eigenvalues[0]]
    legends += ['average']
    colors += [1]
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_VVV_Eigenvalues)):
        Hist_VVV_Eigenvalues[entry].SetMinimum(1e-3)
        Hists += [Hist_VVV_Eigenvalues[entry]]
        legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        #colors += [entry+1]
        colors += [int(random_gen.Uniform(29.,49.))]
    MakeMultiplePlot(Hists,legends,colors,'entry','singular value','SingularValue_E%s%s'%(e,lowrank_tag),0,0,False,True)

    Hists = []
    legends = []
    colors = []
    Hist_Data_Eigenvalues[0].SetMinimum(1e-3)
    Hists += [Hist_Data_Eigenvalues[0]]
    legends += ['average']
    colors += [1]
    random_gen = ROOT.TRandom3()
    for entry in range(1,len(Hist_Data_Eigenvalues)):
        Hist_Data_Eigenvalues[entry].SetMinimum(1e-3)
        Hists += [Hist_Data_Eigenvalues[entry]]
        legends += ['exposure %0.1f'%(data_exposure[entry-1])]
        #colors += [entry+1]
        colors += [int(random_gen.Uniform(29.,49.))]
    MakeMultiplePlot(Hists,legends,colors,'Rank (r)','singular value','SingularValue_Mon_E%s%s'%(e,lowrank_tag),0,0,False,True)

    print ('Energy %s'%(energy_bin[e]))
    for entry in range(0,len(Hist_U_Proj)):
        source_name = 'Average'
        if entry>0:
            source_name = sample_name[entry-1]
        print (source_name)
        tab_U_proj = [['1',0.,0.,0.,0.],['2',0.,0.,0.,0.],['3',0.,0.,0.,0.]]
        tab_V_proj = [['1',0.,0.,0.,0.],['2',0.,0.,0.,0.],['3',0.,0.,0.,0.]]
        for row in range(0,Hist_U_Proj[entry].GetNbinsX()):
            for col in range(0,Hist_U_Proj[entry].GetNbinsY()):
                content_u = Hist_U_Proj[entry].GetBinContent(row+1,col+1)
                content_v = Hist_V_Proj[entry].GetBinContent(row+1,col+1)
                if row>=3: continue
                if col<3:
                    tab_U_proj[row][col+1] = pow(content_u,0.5)
                    tab_V_proj[row][col+1] = pow(content_v,0.5)
                else:
                    tab_U_proj[row][3+1] += content_u
                    tab_V_proj[row][3+1] += content_v
        for row in range(0,Hist_U_Proj[entry].GetNbinsX()):
            if row>=3: continue
            tab_U_proj[row][3+1] = pow(tab_U_proj[row][3+1],0.5)
            tab_V_proj[row][3+1] = pow(tab_V_proj[row][3+1],0.5)
        tab_GammaRegion_Contribution = ['$\epsilon$ in gamma region',0.,0.,0.,0.]
        for row in range(0,Hist_GammaRegion_Contribution[entry].GetNbinsX()):
            if row<3:
                tab_GammaRegion_Contribution[row+1] = pow(Hist_GammaRegion_Contribution[entry].GetBinContent(row+1),0.5)
        tab_GammaRegion_Contribution[3+1] = pow(Hist_GammaRegion_Contribution[entry].GetBinContent(Hist_GammaRegion_Contribution[entry].GetNbinsX()),0.5)
        my_table = PrettyTable()
        my_table.field_names = ["Rank (r)", "$E^{U}_{r1}$", "$E^{U}_{r2}$", "$E^{U}_{r3}$", "$\sum_{4}^{16} E^{U}_{ri}$"]
        for row in range(0,len(tab_U_proj)):
            my_table.add_row(tab_U_proj[row])
        print(my_table)
        my_table = PrettyTable()
        my_table.field_names = ["Rank (r)", "$E^{V}_{r1}$", "$E^{V}_{r2}$", "$E^{V}_{r3}$", "$\sum_{4}^{16} E^{V}_{ri}$"]
        for row in range(0,len(tab_V_proj)):
            my_table.add_row(tab_V_proj[row])
        print(my_table)
        sum_proj = 1./6.*(pow(tab_U_proj[0][4],2)+pow(tab_U_proj[1][4],2)+pow(tab_U_proj[2][4],2))
        sum_proj += 1./6.*(pow(tab_V_proj[0][4],2)+pow(tab_V_proj[1][4],2)+pow(tab_V_proj[2][4],2))
        print ('eta = %s, eta^{0.5} = %s'%(sum_proj,pow(sum_proj,0.5)))
        my_table = PrettyTable()
        my_table.field_names = ["k<=3", "n<=1", "n<=2", "n<=3", "n<=16"]
        my_table.add_row(tab_GammaRegion_Contribution)
        print(my_table)
    eta_array = []
    best_array = []
    par9_array = []
    par9_sig_array = []
    my_table = PrettyTable()
    my_table.field_names = ["source","data","sqrt data epsilon","best rank3 epsilon","best 9-par epsilon","initial epsilon","8(9)-par epsilon", "weighted 9-par epsilon", "low-rank epsilon","eta"]
    my_table.float_format["sqrt data epsilon"] = ".2e"
    my_table.float_format["best rank3 epsilon"] = ".2e"
    my_table.float_format["best 9-par epsilon"] = ".2e"
    my_table.float_format["initial epsilon"] = ".2e"
    my_table.float_format["9-par epsilon"] = ".2e"
    my_table.float_format["weighted 9-par epsilon"] = ".2e"
    my_table.float_format["low-rank epsilon"] = ".2e"
    my_table.float_format["eta"] = ".2e"
    for entry in range(0,len(Hist_U_Proj)):
        tab_U_proj = [['1',0.,0.,0.,0.],['2',0.,0.,0.,0.],['3',0.,0.,0.,0.]]
        tab_V_proj = [['1',0.,0.,0.,0.],['2',0.,0.,0.,0.],['3',0.,0.,0.,0.]]
        for row in range(0,Hist_U_Proj[entry].GetNbinsX()):
            for col in range(0,Hist_U_Proj[entry].GetNbinsY()):
                content_u = Hist_U_Proj[entry].GetBinContent(row+1,col+1)
                content_v = Hist_V_Proj[entry].GetBinContent(row+1,col+1)
                if row>=3: continue
                if col<3:
                    tab_U_proj[row][col+1] = pow(content_u,0.5)
                    tab_V_proj[row][col+1] = pow(content_v,0.5)
                else:
                    tab_U_proj[row][3+1] += content_u
                    tab_V_proj[row][3+1] += content_v
        for row in range(0,Hist_U_Proj[entry].GetNbinsX()):
            if row>=3: continue
            tab_U_proj[row][3+1] = pow(tab_U_proj[row][3+1],0.5)
            tab_V_proj[row][3+1] = pow(tab_V_proj[row][3+1],0.5)
        sum_proj = 1./6.*(pow(tab_U_proj[0][4],2)+pow(tab_U_proj[1][4],2)+pow(tab_U_proj[2][4],2))
        sum_proj += 1./6.*(pow(tab_V_proj[0][4],2)+pow(tab_V_proj[1][4],2)+pow(tab_V_proj[2][4],2))
        if entry>0 and data_count[entry-1]>0.:
            eta_array += [sum_proj]
            best_array += [abs(1.-rank2_count[entry-1]/data_count[entry-1])]
            par9_array += [pow(Hist_GammaRegion_Contribution[entry].GetBinContent(2+1),0.5)]
            par9_sig_array += [pow(Hist_GammaRegion_Contribution[entry].GetBinContent(2+1),0.5)/(pow(data_count[entry-1],0.5)/data_count[entry-1])]
        source_name = 'Average'
        if entry>0 and data_count[entry-1]>0.:
            source_name = sample_name[entry-1]
            my_table.add_row([source_name,data_count[entry-1],pow(data_count[entry-1],0.5)/data_count[entry-1],abs(1.-rank2_count[entry-1]/data_count[entry-1]),pow(Hist_GammaRegion_Contribution[entry].GetBinContent(2+1),0.5),abs(1.-dark_count[entry-1]/data_count[entry-1]),abs(1.-par9_count[entry-1]/data_count[entry-1]),abs(1.-wpar9_count[entry-1]/data_count[entry-1]),abs(1.-bkgd_count[entry-1]/data_count[entry-1]),sum_proj])
        else:
            my_table.add_row([source_name,1.,AccuracyStat_mean,AccuracyRank2_mean,AccuracyBestPar9_mean,AccuracyInit_mean,AccuracyPar9_mean,AccuracyWPar9_mean,Accuracy_mean,sum_proj])
    print(my_table)

    plt.clf()
    plt.xlabel("$eta$", fontsize=18)
    plt.ylabel("$abs(N_{\gamma bkg}-N_{best-9-par})/N_{\gamma bkg}$", fontsize=18)
    plt.scatter(eta_array,par9_array)
    plt.savefig("output_plots/EtaVsBestPar9_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    plt.clf()
    plt.xlabel("$eta$", fontsize=18)
    plt.ylabel("$abs(N_{\gamma bkg}-N_{best-9-par})/ sqrt(N_{\gamma bkg})$", fontsize=18)
    plt.scatter(eta_array,par9_sig_array)
    plt.savefig("output_plots/EtaVsBestPar9Sig_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    par9_epsilon = []
    wpar9_epsilon = []
    for entry in range(0,len(wpar9_count)):
        if data_count[entry]==0.: continue
        par9_epsilon += [abs(data_count[entry]-par9_count[entry])/data_count[entry]]
        wpar9_epsilon += [abs(data_count[entry]-wpar9_count[entry])/data_count[entry]]
    plt.clf()
    plt.xlabel("$\epsilon$ (plain Frobenius norm)", fontsize=18)
    plt.ylabel("$\epsilon$ (weighted Frobenius norm)", fontsize=18)
    plt.xlim(0.,0.1)
    plt.ylim(0.,0.1)
    plt.scatter(par9_epsilon,wpar9_epsilon)
    line_x = np.arange(0.0, 0.1, 1e-4) # angular size
    line_y = line_x
    plt.plot(line_x, line_y, color='r')
    plt.savefig("output_plots/PlainVsWeightedFrobenius_Count_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    par9_epsilon = []
    wpar9_epsilon = []
    for entry in range(0,len(wpar9_count)):
        par9_epsilon += [par9_chi2[entry]]
        wpar9_epsilon += [wpar9_chi2[entry]]
    plt.clf()
    plt.xlabel("$\epsilon$ (plain Frobenius norm)", fontsize=18)
    plt.ylabel("$\epsilon$ (weighted Frobenius norm)", fontsize=18)
    plt.xlim(0.,0.1)
    plt.ylim(0.,0.1)
    plt.scatter(par9_epsilon,wpar9_epsilon)
    line_x = np.arange(0.0, 0.1, 1e-4) # angular size
    line_y = line_x
    plt.plot(line_x, line_y, color='r')
    plt.savefig("output_plots/PlainVsWeightedFrobenius_Chi2_E%s_%s%s.png"%(e,method_tag,lowrank_tag))

    #Hists = []
    #legends = []
    #colors = []
    #Hist_Bkgd_Optimization[0].GetXaxis().SetLabelOffset(999)
    #Hists += [Hist_Bkgd_Optimization[0]]
    #legends += ['average']
    #colors += [1]
    #for entry in range(1,len(Hist_Bkgd_Optimization)):
    #    Hists += [Hist_Bkgd_Optimization[entry]]
    #    legends += ['source %s'%(entry)]
    #    colors += [entry+1]
    #MakeMultiplePlot(Hists,legends,colors,'number of entries included','abs(N_{#gamma bkg}-N_{model})/N_{#gamma bkg}','OptimizationParameter_Entry_E%s'%(e),0,0,False,False)

    #Hists = []
    #legends = []
    #colors = []
    #Hist_Bkgd_Chi2[0].GetXaxis().SetLabelOffset(999)
    #Hists += [Hist_Bkgd_Chi2[0]]
    #legends += ['average']
    #colors += [1]
    #for entry in range(1,len(Hist_Bkgd_Chi2)):
    #    Hists += [Hist_Bkgd_Chi2[entry]]
    #    legends += ['source %s'%(entry)]
    #    colors += [entry+1]
    #MakeMultiplePlot(Hists,legends,colors,'number of entries included','#chi^{2} in CR','Chi2_Entry_E%s'%(e),0,0,False,False)

