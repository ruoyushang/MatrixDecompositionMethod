
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

N_bins_for_deconv = 16
ErecS_lower_cut = 0
ErecS_upper_cut = 0
total_exposure_hours = 0.

folder_path = 'output_test'
#method_tag = 'loose_mdm_default'
method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

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
sample_list += ['CrabV5']
sample_list += ['CrabV6']
    
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

def GetHistogramsFromFile(FilePath):
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
    Hist2D_Coeff_Data.Reset()
    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Coeff_Data.Add(InputFile.Get(HistName))
    for binx in range(0,Hist2D_Coeff_Data.GetNbinsX()):
        for biny in range(0,Hist2D_Coeff_Data.GetNbinsY()):
            old_content = Hist2D_Regularization[energy_index].GetBinContent(binx+1,biny+1)
            new_content = Hist2D_Coeff_Data.GetBinContent(binx+1,biny+1)
            Hist2D_Regularization[energy_index].SetBinContent(binx+1,biny+1,old_content+new_content*new_content)
            old_content = Hist2D_Regularization_rev[energy_index].GetBinContent(binx+1,biny+1)
            new_content = Hist2D_Coeff_Data.GetBinContent(Hist2D_Coeff_Data.GetNbinsX()-binx,Hist2D_Coeff_Data.GetNbinsY()-biny)
            Hist2D_Regularization_rev[energy_index].SetBinContent(binx+1,biny+1,old_content+new_content*new_content)

Hist2D_Coeff_Data = ROOT.TH2D("Hist2D_Coeff_Data","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Regularization = []
Hist2D_Regularization_rev = []
for e in range(0,len(energy_bin)-1):
    ErecS_lower_cut = energy_bin[e]
    ErecS_upper_cut = energy_bin[e+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    Hist2D_Regularization += [ROOT.TH2D("Hist2D_Regularization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
    Hist2D_Regularization_rev += [ROOT.TH2D("Hist2D_Regularization_rev_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]

FilePath_List = []
for source in range(0,len(sample_list)):
    source_name = sample_list[source]
    for elev in range(0,len(root_file_tags)):
        FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
        FilePath_List += [FilePath]
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
        print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
        for e in range(0,len(energy_bin)-1):
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1])

OutputFilePath = "%s/Regularization%s.root"%(folder_path,lowrank_tag)
OutputFile = ROOT.TFile(OutputFilePath,"recreate")
print "total_exposure_hours = %s"%(total_exposure_hours)
for e in range(0,len(energy_bin)-1):
    #Hist2D_Regularization[e].Scale(1./(total_exposure_hours*total_exposure_hours))
    Hist2D_Regularization[e].Write()
OutputFile.Close()


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
for e in range(0,len(energy_bin)-1):
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sum_{all fields} t_{kn}^{2} coefficents' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    pad1.SetLogz()
    #Hist2D_Coeff_Data_Sum.SetMaximum(1e-1);
    #Hist2D_Coeff_Data_Sum.SetMinimum(1.0);
    Hist2D_Regularization[e].GetYaxis().SetTitle('n')
    Hist2D_Regularization[e].GetXaxis().SetTitle('k')
    Hist2D_Regularization[e].Draw("COL4Z")
    #Hist2D_Regularization_rev[e].GetYaxis().SetTitle('n')
    #Hist2D_Regularization_rev[e].GetXaxis().SetTitle('k')
    #Hist2D_Regularization_rev[e].Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_t2_data_%s_%s.png'%(e,method_tag))
