
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
from prettytable import PrettyTable
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

target_energy_index = 0
N_bins_for_deconv = 12
gamma_hadron_dim_ratio_w = 1.
gamma_hadron_dim_ratio_l = 1.
MSCW_blind_cut = 0.5
MSCL_blind_cut = 0.6
MSCW_plot_lower = -0.6
MSCL_plot_lower = -0.6
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
ErecS_lower_cut = 0
ErecS_upper_cut = 0
total_exposure_hours = 0.
Skymap_nbins = 120
Skymap_size = 3.

folder_path = 'output_root'
#folder_path = 'output_test'
#method_tag = 'loose_mdm_default'
method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

ONOFF_tag = 'OFF'
ONOFF_tag += '_Model0'
sample_list = []
sample_name = []
sample_list += ['1ES0647V6_OFF']
sample_name += ['1ES0647 V6']
sample_list += ['1ES1011V6_OFF']
sample_name += ['1ES1011 V6']
sample_list += ['OJ287V6_OFF']
sample_name += ['OJ287 V6']
sample_list += ['PKS1424V6_OFF']
sample_name += ['PKS1424 V6']
sample_list += ['3C264V6_OFF']
sample_name += ['3C264 V6']
sample_list += ['1ES0229V6_OFF']
sample_name += ['1ES0229 V6']
#sample_list += ['1ES0229V5_OFF']
#sample_name += ['1ES0229 V5']
sample_list += ['Segue1V6_OFF']
sample_name += ['Segue1 V6']
#sample_list += ['Segue1V5_OFF']
#sample_name += ['Segue1 V5']
#sample_list += ['CrabV5_OFF']
#sample_name += ['Crab V5']
#sample_list += ['CrabV6_OFF']
#sample_name += ['Crab V6']
#sample_list += ['BLLacV6_OFF']
#sample_name += ['BLLac V6']
#sample_list += ['BLLacV5_OFF']
#sample_name += ['BLLac V5']
#sample_list += ['PG1553V5_OFF']
#sample_name += ['PG1553 V5']
sample_list += ['H1426V6_OFF']
sample_name += ['H1426 V6']
#sample_list += ['CasAV6_OFF']
#sample_name += ['CasA V6']
#sample_list += ['RBS0413V6_OFF']
#sample_name += ['RBS0413 V6']
#sample_list += ['NGC1275V6_OFF']
#sample_name += ['NGC 1275 V6']
sample_list += ['M82V6_OFF']
sample_name += ['M82 V6']

    
elev_bins = [45,85]
theta2_bins = [0,9]

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

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy,normalized):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.1)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

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

    gap = 0.1*(max_heigh-min_heigh)
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        if normalized: Hists[0].DrawNormalized("E")
        else: Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        if normalized: Hists[max_hist].DrawNormalized("E")
        else: Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(2)
        if normalized: Hists[h].DrawNormalized("E same")
        else: Hists[h].Draw("E same")

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

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))


def GetHistogramsFromFile(FilePath):
    global total_exposure_hours
    global Hist_OnData_PerElev_MSCW
    global Hist_OnData_PerElev_MSCL
    global Hist_OnData_PerElev_Xoff
    global Hist_OnData_PerElev_Yoff
    global Hist_CRData_PerElev_Yoff
    global Hist_OnData_PerYear_MSCW
    global Hist_OnData_PerYear_MSCL
    global Hist_OnData_PerYear_Xoff
    global Hist_OnData_PerYear_Yoff
    InputFile = ROOT.TFile(FilePath)
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy_index + elev
        HistName = "Hist_OnData_ThisElev_MSCW_V%s_ErecS%sto%s"%(elev,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerElev_MSCW[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisElev_MSCL_V%s_ErecS%sto%s"%(elev,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerElev_MSCL[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisElev_Xoff_V%s_ErecS%sto%s"%(elev,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerElev_Xoff[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisElev_Yoff_V%s_ErecS%sto%s"%(elev,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerElev_Yoff[index].Add(InputFile.Get(HistName))
        HistName = "Hist_CRData_ThisElev_Yoff_V%s_ErecS%sto%s"%(elev,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_CRData_PerElev_Yoff[index].Add(InputFile.Get(HistName))
    for year in range(0,len(MJD_bins)-1):
        index = (len(MJD_bins)-1)*energy_index + year
        HistName = "Hist_OnData_ThisYear_MSCW_V%s_ErecS%sto%s"%(year,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerYear_MSCW[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisYear_MSCL_V%s_ErecS%sto%s"%(year,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerYear_MSCL[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisYear_Xoff_V%s_ErecS%sto%s"%(year,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerYear_Xoff[index].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_ThisYear_Yoff_V%s_ErecS%sto%s"%(year,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_PerYear_Yoff[index].Add(InputFile.Get(HistName))

elev_bins = [45,50,55,60,65,70,75,80,85]
MJD_bins = [53613,55074,56535,57996,59457]

Hist_OnData_PerElev_MSCW = []
Hist_OnData_PerElev_MSCL = []
Hist_OnData_PerElev_Xoff = []
Hist_OnData_PerElev_Yoff = []
Hist_CRData_PerElev_Yoff = []
for energy in range(0,len(energy_bin)-1):
    for elev in range(0,len(elev_bins)-1):
        Hist_OnData_PerElev_MSCW += [ROOT.TH1D("Hist_OnData_PerElev_MSCW_E%s_Z%s"%(energy,elev),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_OnData_PerElev_MSCL += [ROOT.TH1D("Hist_OnData_PerElev_MSCL_E%s_Z%s"%(energy,elev),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
        Hist_OnData_PerElev_Xoff += [ROOT.TH1D("Hist_OnData_PerElev_Xoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        Hist_OnData_PerElev_Yoff += [ROOT.TH1D("Hist_OnData_PerElev_Yoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        Hist_CRData_PerElev_Yoff += [ROOT.TH1D("Hist_OnData_PerElev_Yoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]

Hist_OnData_DiffElev_Yoff = []
Hist_OnData_DiffRefElev_Yoff = []
Hist_CRData_DiffRefElev_Yoff = []
Hist_OnData_PerYear_MSCW = []
Hist_OnData_PerYear_MSCL = []
Hist_OnData_PerYear_Xoff = []
Hist_OnData_PerYear_Yoff = []
for energy in range(0,len(energy_bin)-1):
    for year in range(0,len(MJD_bins)-1):
        Hist_OnData_PerYear_MSCW += [ROOT.TH1D("Hist_OnData_PerYear_MSCW_E%s_Y%s"%(energy,year),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_OnData_PerYear_MSCL += [ROOT.TH1D("Hist_OnData_PerYear_MSCL_E%s_Y%s"%(energy,year),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
        Hist_OnData_PerYear_Xoff += [ROOT.TH1D("Hist_OnData_PerYear_Xoff_E%s_Y%s"%(energy,year),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        Hist_OnData_PerYear_Yoff += [ROOT.TH1D("Hist_OnData_PerYear_Yoff_E%s_Y%s"%(energy,year),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]

FilePath_List = []
#for source in range(0,len(sample_list)):
#    source_name = sample_list[source]
#    for the_file in range(0,len(root_file_tags)):
#        FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[the_file])+".root"
#        FilePath_List += [FilePath]
#        print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
#        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
#        for e in range(0,len(energy_bin)-1):
#            ErecS_lower_cut = energy_bin[e]
#            ErecS_upper_cut = energy_bin[e+1]
#            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1])

FilePath = "output_root/ObservingEffect/ObservingEffect.root"
FilePath_List += [FilePath]
print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
for energy in range(0,len(energy_bin)-1):
    ErecS_lower_cut = energy_bin[energy]
    ErecS_upper_cut = energy_bin[energy+1]
    GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1])

for energy in range(0,len(energy_bin)-1):
    Hists = []
    legends = []
    colors = []
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hists += [Hist_OnData_PerElev_MSCW[index]]
        legends += ['%s-%s'%(elev_bins[elev],elev_bins[elev+1])]
        colors += [elev+1]
    MakeComparisonPlot(Hists,legends,colors,'MSCW','Normalized counts','PerElev_MSCW_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hists += [Hist_OnData_PerElev_MSCL[index]]
        legends += ['%s-%s'%(elev_bins[elev],elev_bins[elev+1])]
        colors += [elev+1]
    MakeComparisonPlot(Hists,legends,colors,'MSCL','Normalized counts','PerElev_MSCL_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hists += [Hist_OnData_PerElev_Xoff[index]]
        legends += ['%s-%s'%(elev_bins[elev],elev_bins[elev+1])]
        colors += [elev+1]
    MakeComparisonPlot(Hists,legends,colors,'X off','Normalized counts','PerElev_Xoff_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hists += [Hist_OnData_PerElev_Yoff[index]]
        legends += ['%s-%s'%(elev_bins[elev],elev_bins[elev+1])]
        colors += [elev+1]
    MakeComparisonPlot(Hists,legends,colors,'Y off','Normalized counts','PerElev_Yoff_E%s'%(energy),0,0,False,False,True)
    
    Hists = []
    legends = []
    colors = []
    for year in range(0,len(MJD_bins)-1):
        index = (len(MJD_bins)-1)*energy + year
        Hists += [Hist_OnData_PerYear_MSCW[index]]
        legends += ['%s-%s'%(MJD_bins[year],MJD_bins[year+1])]
        colors += [year+1]
    MakeComparisonPlot(Hists,legends,colors,'MSCW','Normalized counts','PerYear_MSCW_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for year in range(0,len(MJD_bins)-1):
        index = (len(MJD_bins)-1)*energy + year
        Hists += [Hist_OnData_PerYear_MSCL[index]]
        legends += ['%s-%s'%(MJD_bins[year],MJD_bins[year+1])]
        colors += [year+1]
    MakeComparisonPlot(Hists,legends,colors,'MSCL','Normalized counts','PerYear_MSCL_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for year in range(0,len(MJD_bins)-1):
        index = (len(MJD_bins)-1)*energy + year
        Hists += [Hist_OnData_PerYear_Xoff[index]]
        legends += ['%s-%s'%(MJD_bins[year],MJD_bins[year+1])]
        colors += [year+1]
    MakeComparisonPlot(Hists,legends,colors,'X off','Normalized counts','PerYear_Xoff_E%s'%(energy),0,0,False,False,True)
    Hists = []
    legends = []
    colors = []
    for year in range(0,len(MJD_bins)-1):
        index = (len(MJD_bins)-1)*energy + year
        Hists += [Hist_OnData_PerYear_Yoff[index]]
        legends += ['%s-%s'%(MJD_bins[year],MJD_bins[year+1])]
        colors += [year+1]
    MakeComparisonPlot(Hists,legends,colors,'Y off','Normalized counts','PerYear_Yoff_E%s'%(energy),0,0,False,False,True)
    
    
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        norm = Hist_OnData_PerElev_Yoff[index].Integral()
        if norm>0.: Hist_OnData_PerElev_Yoff[index].Scale(1./norm)
        norm = Hist_CRData_PerElev_Yoff[index].Integral()
        if norm>0.: Hist_CRData_PerElev_Yoff[index].Scale(1./norm)
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hist_OnData_DiffElev_Yoff += [ROOT.TH1D("Hist_OnData_DiffElev_Yoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        ref_bin = elev - 1 
        if ref_bin<0: continue
        print "ref_bin = %s"%(ref_bin)
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hist_OnData_DiffElev_Yoff[index].Add(Hist_OnData_PerElev_Yoff[ref_index])
        Hist_OnData_DiffElev_Yoff[index].Add(Hist_OnData_PerElev_Yoff[index],-1.)
        Hist_OnData_DiffElev_Yoff[index].Divide(Hist_OnData_PerElev_Yoff[index])
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hist_OnData_DiffRefElev_Yoff += [ROOT.TH1D("Hist_OnData_DiffRefElev_Yoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        ref_bin = len(elev_bins)-3 
        if ref_bin<0: continue
        print "ref_bin = %s"%(ref_bin)
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hist_OnData_DiffRefElev_Yoff[index].Add(Hist_OnData_PerElev_Yoff[ref_index])
        Hist_OnData_DiffRefElev_Yoff[index].Add(Hist_OnData_PerElev_Yoff[index],-1.)
        Hist_OnData_DiffRefElev_Yoff[index].Divide(Hist_OnData_PerElev_Yoff[index])
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hist_CRData_DiffRefElev_Yoff += [ROOT.TH1D("Hist_CRData_DiffRefElev_Yoff_E%s_Z%s"%(energy,elev),"",Skymap_nbins/6,-Skymap_size,Skymap_size)]
        ref_bin = len(elev_bins)-3 
        if ref_bin<0: continue
        print "ref_bin = %s"%(ref_bin)
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hist_CRData_DiffRefElev_Yoff[index].Add(Hist_CRData_PerElev_Yoff[ref_index])
        Hist_CRData_DiffRefElev_Yoff[index].Add(Hist_CRData_PerElev_Yoff[index],-1.)
        Hist_CRData_DiffRefElev_Yoff[index].Divide(Hist_CRData_PerElev_Yoff[index])
    
    Hists = []
    legends = []
    colors = []
    for elev in range(3,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        ref_bin = elev - 1 
        if ref_bin<0: continue
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hists += [Hist_OnData_DiffElev_Yoff[index]]
        legends += ['%s-%s'%(elev_bins[ref_bin],elev_bins[elev])]
        colors += [elev-2]
    MakeComparisonPlot(Hists,legends,colors,'Y off','Correction','DiffElev_Yoff_E%s'%(energy),-0.5,0.5,False,False,False)
    Hists = []
    legends = []
    colors = []
    for elev in range(3,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        ref_bin = len(elev_bins)-3 
        if ref_bin<0: continue
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hists += [Hist_OnData_DiffRefElev_Yoff[index]]
        legends += ['%s-%s'%(elev_bins[ref_bin],elev_bins[elev])]
        colors += [elev-2]
    MakeComparisonPlot(Hists,legends,colors,'Y off','Correction','DiffRefElev_Yoff_E%s'%(energy),-0.5,0.5,False,False,False)
    Hists = []
    legends = []
    colors = []
    for elev in range(3,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        ref_bin = len(elev_bins)-3 
        if ref_bin<0: continue
        ref_index = (len(elev_bins)-1)*energy + ref_bin
        Hists += [Hist_CRData_DiffRefElev_Yoff[index]]
        legends += ['%s-%s'%(elev_bins[ref_bin],elev_bins[elev])]
        colors += [elev-2]
    MakeComparisonPlot(Hists,legends,colors,'Y off','Correction','CRDiffRefElev_Yoff_E%s'%(energy),-0.5,0.5,False,False,False)

OutputFilePath = "output_root/ObservingEffect/YoffCorrection.root"
OutputFile = ROOT.TFile(OutputFilePath,"recreate")
for energy in range(0,len(energy_bin)-1):
    for elev in range(0,len(elev_bins)-1):
        index = (len(elev_bins)-1)*energy + elev
        Hist_OnData_DiffElev_Yoff[index].Write()
OutputFile.Close()
