
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

folder_path = 'output_test'
method_tag = []
legends = []
colors = []
#method_tag += ['loose_mdm_default']
#legends += ['Gen. Tik.']
#method_tag += ['loose_mdm_rank3']
#legends += ['9-Par']
#method_tag += ['loose_mdm_rank5']
#legends += ['21-Par']
#method_tag += ['loose_mdm_cutoff']
#legends += ['Cutoff signular']
#method_tag += ['loose_mdm_tikhonov']
#legends += ['Tikhonov']
method_tag += ['tight_mdm_default']
legends += ['Gen. Tik.']
colors += [1]
#method_tag += ['tight_mdm_cutoff']
#legends += ['Cutoff signular']
#method_tag += ['tight_mdm_tikhonov']
#legends += ['Tikhonov']
#colors += [1]
#method_tag += ['tight_mdm_cutoff_eigen']
#legends += ['Cutoff Eigen']
#method_tag += ['tight_mdm_rank3']
#legends += ['9-Par']
#colors += [4]
#method_tag += ['tight_mdm_rank5']
#legends += ['21-Par']
#colors += [2]

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
for method in range(0,len(method_tag)):
    method_tag[method] += lowrank_tag

ONOFF_tag = 'OFF'
sample_list = []
sample_list += ['OJ287V6']
sample_list += ['1ES0229V6']
sample_list += ['1ES0229V5']
sample_list += ['H1426V6']
sample_list += ['PKS1424V6']
sample_list += ['3C264V6']
sample_list += ['1ES1011V6']
sample_list += ['Segue1V6']
sample_list += ['Segue1V5']
sample_list += ['BLLacV6']
sample_list += ['BLLacV5']
sample_list += ['CrabV6']
sample_list += ['CrabV5']

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
            root_file_tags += [elev_tag+theta2_tag+signal_tag+mjd_tag[d]+'_'+ONOFF_tag]

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
    Hists[0][0].GetXaxis().SetTitleOffset(0.8)
    Hists[0][0].GetXaxis().SetTitleSize(0.06)
    Hists[0][0].GetXaxis().SetLabelSize(0.06)
    Hists[0][0].GetYaxis().SetLabelSize(0.06)
    Hists[0][0].GetYaxis().SetTitleOffset(1.2)
    Hists[0][0].GetYaxis().SetTitleSize(0.06)
    Hists[0][0].GetXaxis().SetTitle(title_x)
    Hists[0][0].GetYaxis().SetTitle(title_y)
    Hists[0][0].SetMaximum(y_max)
    Hists[0][0].SetMinimum(y_min)
    Hists[0][0].Draw("E")

    for method in range(0,len(Hists)):
        for source in range(0,len(Hists[method])):
            print '%s'%(Hists[method][source].GetBinContent(1))
            Hists[method][source].SetLineColor(colors[method])
            if source!=0:
                if colors[method]==1:
                    Hists[method][source].SetLineColor(15)
                elif colors[method]==2:
                    Hists[method][source].SetLineColor(45)
                elif colors[method]==4:
                    Hists[method][source].SetLineColor(38)
            if source==0:
                Hists[method][source].SetLineWidth(4)
            else:
                Hists[method][source].SetLineWidth(2)
            #if source!=0:
            #    continue
            Hists[method][source].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.2,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for method in range(0,len(Hists)):
        legend.AddEntry(Hists[method][0],'%s'%(legends[method]),"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def GetHistogramsFromFile(FilePath,which_method,which_source):
    global total_exposure_hours
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    if energy_index==0: 
        total_exposure_hours += exposure_hours
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization[which_method][0].Add(InputFile.Get(HistName),1./float(len(sample_list)))
    Hist_Bkgd_Optimization[which_method][which_source+1].Add(InputFile.Get(HistName))

optimiz_lower = -10.
optimiz_upper = 0.
Hist_Bkgd_Optimization = []
for e in range(0,len(method_tag)):
    Hist_Bkgd_Optimization_ThisMethod = []
    for s in range(0,len(sample_list)+1):
        Hist_Bkgd_Optimization_ThisMethod += [ROOT.TH1D("Hist_Bkgd_Optimization_%s_%s"%(e,s),"",N_bins_for_deconv*N_bins_for_deconv,optimiz_lower,optimiz_upper)]
    Hist_Bkgd_Optimization += [Hist_Bkgd_Optimization_ThisMethod]

for e in range(0,len(energy_bin)-1):
    FilePath_List = []
    for entry in range(0,len(method_tag)):
        Hist_Bkgd_Optimization[entry][0].Reset()
        for source in range(0,len(sample_list)):
            Hist_Bkgd_Optimization[entry][source+1].Reset()
            source_name = sample_list[source]
            for elev in range(0,len(root_file_tags)):
                FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s%s"%(method_tag[entry],root_file_tags[elev])+".root"
                FilePath_List += [FilePath]
                if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
                print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
                ErecS_lower_cut = energy_bin[e]
                ErecS_upper_cut = energy_bin[e+1]
                GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1],entry,source)

    MakeMultiplePlot(Hist_Bkgd_Optimization,legends,colors,'log10 #alpha','abs(N_{#gamma}-N_{model})/N_{#gamma}','PerformanceEvaluation%s_E%s'%(lowrank_tag,e),0,0.1,False,False)
    #MakeMultiplePlot(Hist_Bkgd_Optimization,legends,colors,'log10 #alpha','abs(N_{#gamma}-N_{model})/N_{#gamma}','PerformanceEvaluation%s_E%s'%(lowrank_tag,e),1e-3,0.1,False,True)

