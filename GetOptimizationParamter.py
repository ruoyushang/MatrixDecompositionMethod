
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

N_bins_for_deconv = 16
ErecS_lower_cut = 0
ErecS_upper_cut = 0
total_exposure_hours = 0.
Zenith_mean_data = []
NSB_mean_data = []
Accuracy_source = []
color_scatter = []
area_scatter = []
legend_scatter = []

folder_path = 'output_test'
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
sample_list = []
sample_list += ['H1426V6']
sample_list += ['OJ287V6']
sample_list += ['1ES0229V6']
sample_list += ['1ES0229V5']
sample_list += ['PKS1424V6']
sample_list += ['3C264V6']
sample_list += ['1ES1011V6']
sample_list += ['Segue1V6']
sample_list += ['Segue1V5']
sample_list += ['BLLacV6']
sample_list += ['BLLacV5']
sample_list += ['CrabV5']
sample_list += ['CrabV6']
#sample_list += ['M82V5']

#ONOFF_tag = 'ON'
#sample_list = []
#sample_list += ['WComaeV6']
#sample_list += ['WComaeV5']
#sample_list += ['WComaeV4']
    
elev_bins = [45,85]
theta2_bins = [0,4]

energy_bin = []
energy_bin += [int(pow(10,2.0))]
energy_bin += [int(pow(10,2.33))]
energy_bin += [int(pow(10,2.66))]
energy_bin += [int(pow(10,3.0))]
energy_bin += [int(pow(10,3.33))]
energy_bin += [int(pow(10,3.66))]
energy_bin += [int(pow(10,4.0))]

signal_tag = '_S0'
root_file_tags = []
mjd_tag = []
mjd_tag += ['']
for elev in range(0,len(elev_bins)-1):
    elev_tag = '_TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+signal_tag+mjd_tag[d]+'_'+ONOFF_tag]

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
        Hists[h].SetLineWidth(2)
        Hists[h].Draw("E same")
    Hists[0].SetLineWidth(4)
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
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
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

    min_y = 1.0
    min_bin = 0
    for binx in range(1,Hists[0].GetNbinsX()+1):
        if Hists[0].GetBinContent(binx)<min_y:
            min_y = Hists[0].GetBinContent(binx)
            min_bin = binx

    lumilab1 = ROOT.TLatex(0.15,0.80,'log10 #alpha = %.1f'%(Hists[0].GetBinCenter(min_bin)) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def GetHistogramsFromFile(FilePath,which_source):
    global total_exposure_hours
    global Zenith_mean_data
    global NSB_mean_data
    global area_scatter
    global legend_scatter
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    NewInfoTree = InputFile.Get("NewInfoTree")
    InfoTree.GetEntry(0)
    NewInfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    if energy_index==0: 
        total_exposure_hours += exposure_hours
    Zenith_mean = NewInfoTree.Zenith_mean_data
    Zenith_mean_data += [Zenith_mean]
    area_scatter += [pow(exposure_hours/10.,2)]
    legend_scatter += ['%0.1f hrs'%(exposure_hours)]
    NSB_mean = NewInfoTree.NSB_mean_data
    NSB_mean_data += [NSB_mean]
    Hist_Bkgd_Optimization[which_source+1].Reset()
    weight = 1./float(len(sample_list))
    #weight = exposure_hours
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Optimization[which_source+1].Add(InputFile.Get(HistName))
    Hist_Bkgd_Chi2[which_source+1].Reset()
    HistName = "Hist_Bkgd_Chi2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2[0].Add(InputFile.Get(HistName),weight)
    Hist_Bkgd_Chi2[which_source+1].Add(InputFile.Get(HistName))
    Hist_VVV_Eigenvalues[which_source+1].Reset()
    HistName = "Hist_VVV_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_VVV_Eigenvalues[0].Add(InputFile.Get(HistName),weight)
    Hist_VVV_Eigenvalues[which_source+1].Add(InputFile.Get(HistName))

optimiz_lower = -10.
optimiz_upper = 0.
Hist_Bkgd_Optimization = []
for e in range(0,len(sample_list)+1):
    Hist_Bkgd_Optimization += [ROOT.TH1D("Hist_Bkgd_Optimization_%s"%(e),"",N_bins_for_deconv*N_bins_for_deconv,optimiz_lower,optimiz_upper)]
Hist_Bkgd_Chi2 = []
for e in range(0,len(sample_list)+1):
    Hist_Bkgd_Chi2 += [ROOT.TH1D("Hist_Bkgd_Chi2_%s"%(e),"",N_bins_for_deconv*N_bins_for_deconv,optimiz_lower,optimiz_upper)]
Hist_VVV_Eigenvalues = []
for e in range(0,len(sample_list)+1):
    Hist_VVV_Eigenvalues += [ROOT.TH1D("Hist_VVV_Eigenvaluesi_%s"%(e),"",N_bins_for_deconv*N_bins_for_deconv,0,N_bins_for_deconv*N_bins_for_deconv)]

for e in range(0,len(energy_bin)-1):
    FilePath_List = []
    Zenith_mean_data = []
    NSB_mean_data = []
    Accuracy_source = []
    color_scatter = []
    area_scatter = []
    legend_scatter = []
    for entry in range(0,len(Hist_Bkgd_Optimization)):
        Hist_Bkgd_Optimization[entry].Reset()
    for entry in range(0,len(Hist_Bkgd_Chi2)):
        Hist_Bkgd_Chi2[entry].Reset()
    for entry in range(0,len(Hist_VVV_Eigenvalues)):
        Hist_VVV_Eigenvalues[entry].Reset()
    for source in range(0,len(sample_list)):
        source_name = sample_list[source]
        for elev in range(0,len(root_file_tags)):
            FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1],source)

    #Hist_Bkgd_Optimization[0].Scale(1./total_exposure_hours)
    min_y = 1.0
    min_bin = 0
    for binx in range(1,Hist_Bkgd_Optimization[0].GetNbinsX()+1):
        if Hist_Bkgd_Optimization[0].GetBinContent(binx)<min_y:
            min_y = Hist_Bkgd_Optimization[0].GetBinContent(binx)
            min_bin = binx
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        Accuracy_source += [Hist_Bkgd_Optimization[entry].GetBinContent(min_bin)]

    fig, ax = plt.subplots()
    colors = np.random.rand(len(Zenith_mean_data))
    plt.clf()
    plt.xlabel("Zenith")
    plt.ylabel("$abs(N_{\gamma}-N_{model})/N_{\gamma}$")
    plt.scatter(Zenith_mean_data,Accuracy_source,s=area_scatter,c=colors,alpha=0.5)
    plt.savefig("output_plots/Performance_Zenith_E%s%s.png"%(e,lowrank_tag))

    plt.clf()
    plt.xlabel("NSB")
    plt.ylabel("$abs(N_{\gamma}-N_{model})/N_{\gamma}$")
    plt.scatter(NSB_mean_data,Accuracy_source,s=area_scatter,c=colors,alpha=0.5)
    plt.savefig("output_plots/Performance_NSB_E%s%s.png"%(e,lowrank_tag))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Bkgd_Optimization[0]]
    legends += ['average']
    colors += [1]
    #Hist_Bkgd_Optimization[0].GetXaxis().SetLabelOffset(999)
    for entry in range(1,len(Hist_Bkgd_Optimization)):
        Hists += [Hist_Bkgd_Optimization[entry]]
        legends += ['source %s'%(entry)]
        colors += [entry+1]
    MakeMultiplePlot(Hists,legends,colors,'log10 #alpha','abs(N_{#gamma}-N_{model})/N_{#gamma}','OptimizationParameter_E%s%s'%(e,lowrank_tag),1e-3,0.1,False,False)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Bkgd_Chi2[0]]
    legends += ['average']
    colors += [1]
    #Hist_Bkgd_Chi2[0].GetXaxis().SetLabelOffset(999)
    for entry in range(1,len(Hist_Bkgd_Chi2)):
        Hists += [Hist_Bkgd_Chi2[entry]]
        legends += ['source %s'%(entry)]
        colors += [entry+1]
    MakeMultiplePlot(Hists,legends,colors,'log10 #alpha','#chi^{2} in CR','Chi2_E%s%s'%(e,lowrank_tag),1e-7,1e-4,False,True)

    Hists = []
    legends = []
    colors = []
    Hist_VVV_Eigenvalues[0].SetMinimum(1e-3)
    Hists += [Hist_VVV_Eigenvalues[0]]
    legends += ['average']
    colors += [1]
    for entry in range(1,len(Hist_VVV_Eigenvalues)):
        Hist_VVV_Eigenvalues[entry].SetMinimum(1e-3)
        Hists += [Hist_VVV_Eigenvalues[entry]]
        legends += ['source %s'%(entry)]
        colors += [entry+1]
    MakeMultiplePlot(Hists,legends,colors,'entry','singular value','SingularValue_E%s%s'%(e,lowrank_tag),0,0,False,True)

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
    #MakeMultiplePlot(Hists,legends,colors,'number of entries included','abs(N_{#gamma}-N_{model})/N_{#gamma}','OptimizationParameter_Entry_E%s'%(e),0,0,False,False)

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

