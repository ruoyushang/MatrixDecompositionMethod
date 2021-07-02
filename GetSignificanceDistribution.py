
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


folder_path = 'output_root'
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
sample_list += ['1ES1011_OFF']
sample_list += ['H1426_OFF']
sample_list += ['OJ287_OFF']
sample_list += ['1ES0229_OFF']
sample_list += ['PKS1424_OFF']
sample_list += ['3C264_OFF']
sample_list += ['Segue1_OFF']
sample_list += ['Crab_OFF']
sample_list += ['BLLac_OFF']
sample_list += ['NGC1275_OFF']
sample_list += ['PG1553_OFF']
sample_list += ['CasA_OFF']
sample_list += ['RBS0413_OFF']
sample_list += ['RGBJ0710_OFF']
sample_list += ['1ES0647_OFF']

#ONOFF_tag = 'ON'
#sample_list = []
#sample_list += ['WComaeV6']
#sample_list += ['WComaeV5']
#sample_list += ['WComaeV4']
    
elev_bins = [55,85]
theta2_bins = [0,4]

energy_bin = []
energy_bin += [int(pow(10,2.0))]
energy_bin += [int(pow(10,2.33))]
energy_bin += [int(pow(10,2.66))]
energy_bin += [int(pow(10,3.0))]
energy_bin += [int(pow(10,3.33))]
energy_bin += [int(pow(10,3.66))]
energy_bin += [int(pow(10,4.0))]

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
        if legends[h]=='reference':
            Hists[h].Draw("hist same")
        else:
            Hists[h].Draw("E same")
    Hists[0].SetLineWidth(4)
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
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

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))


Hist_Sig_MDM = ROOT.TH1D("Hist_Sig_MDM","",65,-5,8)
Hist_Sig_Raw = ROOT.TH1D("Hist_Sig_Raw","",65,-5,8)
Hist_Sig_Dark = ROOT.TH1D("Hist_Sig_Dark","",65,-5,8)
Hist_Sig_Ref = ROOT.TH1D("Hist_Sig_Ref","",65,-5,8)

InputFile = ROOT.TFile("output_plots/SigDist.root")
for entry in range(0,len(sample_list)):
    HistName = 'Hist_Sig_SigDist_MDM_%s'%(sample_list[entry])
    print 'Get %s'%(HistName)
    Hist_Sig_MDM.Add(InputFile.Get(HistName))
    HistName = 'Hist_Sig_SigDist_Raw_%s'%(sample_list[entry])
    print 'Get %s'%(HistName)
    Hist_Sig_Raw.Add(InputFile.Get(HistName))
    HistName = 'Hist_Sig_SigDist_Dark_%s'%(sample_list[entry])
    print 'Get %s'%(HistName)
    Hist_Sig_Dark.Add(InputFile.Get(HistName))
    HistName = 'Hist_Model_SigDist_MDM_%s'%(sample_list[entry])
    print 'Get %s'%(HistName)
    Hist_Sig_Ref.Add(InputFile.Get(HistName))



Hists = []
legends = []
colors = []
Hists += [Hist_Sig_MDM]
legends += ['all corrections']
colors += [1]
Hists += [Hist_Sig_Dark]
legends += ['accptance corrected']
colors += [4]
Hists += [Hist_Sig_Raw]
legends += ['no correction']
colors += [3]
Hists += [Hist_Sig_Ref]
legends += ['reference']
colors += [2]
MakeMultiplePlot(Hists,legends,colors,'significance (#sigma)','counts','AllDarkFieldsSignificanceDistribution',0.,0.,False,True)

bin_8sigma = Hist_Sig_MDM.FindBin(8.0)+1
bin_5sigma = Hist_Sig_MDM.FindBin(5.0)
bin_4sigma = Hist_Sig_MDM.FindBin(4.0)
bin_3sigma = Hist_Sig_MDM.FindBin(3.0)
bin_2sigma = Hist_Sig_MDM.FindBin(2.0)
bin_1sigma = Hist_Sig_MDM.FindBin(1.0)

prob_ref_total = Hist_Sig_Ref.Integral()
prob_mdm_total = Hist_Sig_MDM.Integral()
prob_dark_total = Hist_Sig_Dark.Integral()
prob_raw_total = Hist_Sig_Raw.Integral()
prob_ref_1sigma = Hist_Sig_Ref.Integral(bin_1sigma,bin_8sigma)/prob_ref_total
prob_ref_2sigma = Hist_Sig_Ref.Integral(bin_2sigma,bin_8sigma)/prob_ref_total
prob_ref_3sigma = Hist_Sig_Ref.Integral(bin_3sigma,bin_8sigma)/prob_ref_total
prob_ref_4sigma = Hist_Sig_Ref.Integral(bin_4sigma,bin_8sigma)/prob_ref_total
prob_ref_5sigma = Hist_Sig_Ref.Integral(bin_5sigma,bin_8sigma)/prob_ref_total
prob_mdm_1sigma = Hist_Sig_MDM.Integral(bin_1sigma,bin_8sigma)/prob_mdm_total
prob_mdm_2sigma = Hist_Sig_MDM.Integral(bin_2sigma,bin_8sigma)/prob_mdm_total
prob_mdm_3sigma = Hist_Sig_MDM.Integral(bin_3sigma,bin_8sigma)/prob_mdm_total
prob_mdm_4sigma = Hist_Sig_MDM.Integral(bin_4sigma,bin_8sigma)/prob_mdm_total
prob_mdm_5sigma = Hist_Sig_MDM.Integral(bin_5sigma,bin_8sigma)/prob_mdm_total
prob_mdm_1sigma_err = pow(Hist_Sig_MDM.Integral(bin_1sigma,bin_8sigma),0.5)/prob_mdm_total
prob_mdm_2sigma_err = pow(Hist_Sig_MDM.Integral(bin_2sigma,bin_8sigma),0.5)/prob_mdm_total
prob_mdm_3sigma_err = pow(Hist_Sig_MDM.Integral(bin_3sigma,bin_8sigma),0.5)/prob_mdm_total
prob_mdm_4sigma_err = pow(Hist_Sig_MDM.Integral(bin_4sigma,bin_8sigma),0.5)/prob_mdm_total
prob_mdm_5sigma_err = pow(Hist_Sig_MDM.Integral(bin_5sigma,bin_8sigma),0.5)/prob_mdm_total
prob_dark_1sigma = Hist_Sig_Dark.Integral(bin_1sigma,bin_8sigma)/prob_dark_total
prob_dark_2sigma = Hist_Sig_Dark.Integral(bin_2sigma,bin_8sigma)/prob_dark_total
prob_dark_3sigma = Hist_Sig_Dark.Integral(bin_3sigma,bin_8sigma)/prob_dark_total
prob_dark_4sigma = Hist_Sig_Dark.Integral(bin_4sigma,bin_8sigma)/prob_dark_total
prob_dark_5sigma = Hist_Sig_Dark.Integral(bin_5sigma,bin_8sigma)/prob_dark_total
prob_dark_1sigma_err = pow(Hist_Sig_Dark.Integral(bin_1sigma,bin_8sigma),0.5)/prob_dark_total
prob_dark_2sigma_err = pow(Hist_Sig_Dark.Integral(bin_2sigma,bin_8sigma),0.5)/prob_dark_total
prob_dark_3sigma_err = pow(Hist_Sig_Dark.Integral(bin_3sigma,bin_8sigma),0.5)/prob_dark_total
prob_dark_4sigma_err = pow(Hist_Sig_Dark.Integral(bin_4sigma,bin_8sigma),0.5)/prob_dark_total
prob_dark_5sigma_err = pow(Hist_Sig_Dark.Integral(bin_5sigma,bin_8sigma),0.5)/prob_dark_total
prob_raw_1sigma = Hist_Sig_Raw.Integral(bin_1sigma,bin_8sigma)/prob_raw_total
prob_raw_2sigma = Hist_Sig_Raw.Integral(bin_2sigma,bin_8sigma)/prob_raw_total
prob_raw_3sigma = Hist_Sig_Raw.Integral(bin_3sigma,bin_8sigma)/prob_raw_total
prob_raw_4sigma = Hist_Sig_Raw.Integral(bin_4sigma,bin_8sigma)/prob_raw_total
prob_raw_5sigma = Hist_Sig_Raw.Integral(bin_5sigma,bin_8sigma)/prob_raw_total
prob_raw_1sigma_err = pow(Hist_Sig_Raw.Integral(bin_1sigma,bin_8sigma),0.5)/prob_raw_total
prob_raw_2sigma_err = pow(Hist_Sig_Raw.Integral(bin_2sigma,bin_8sigma),0.5)/prob_raw_total
prob_raw_3sigma_err = pow(Hist_Sig_Raw.Integral(bin_3sigma,bin_8sigma),0.5)/prob_raw_total
prob_raw_4sigma_err = pow(Hist_Sig_Raw.Integral(bin_4sigma,bin_8sigma),0.5)/prob_raw_total
prob_raw_5sigma_err = pow(Hist_Sig_Raw.Integral(bin_5sigma,bin_8sigma),0.5)/prob_raw_total

my_table = PrettyTable()
my_table.field_names = ["deviation","normal","no correction","after accep. correction","full correction"]
table_row = []
table_row += ["> 1 sigma"]
table_row += ["%.2e"%(prob_ref_1sigma)]
table_row += ["%.2e +/- %.2e"%(prob_raw_1sigma,prob_raw_1sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_dark_1sigma,prob_dark_1sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_mdm_1sigma,prob_mdm_1sigma_err)]
my_table.add_row(table_row)
table_row = []
table_row += ["> 2 sigma"]
table_row += ["%.2e"%(prob_ref_2sigma)]
table_row += ["%.2e +/- %.2e"%(prob_raw_2sigma,prob_raw_2sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_dark_2sigma,prob_dark_2sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_mdm_2sigma,prob_mdm_2sigma_err)]
my_table.add_row(table_row)
table_row = []
table_row += ["> 3 sigma"]
table_row += ["%.2e"%(prob_ref_3sigma)]
table_row += ["%.2e +/- %.2e"%(prob_raw_3sigma,prob_raw_3sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_dark_3sigma,prob_dark_3sigma_err)]
table_row += ["%.2e +/- %.2e"%(prob_mdm_3sigma,prob_mdm_3sigma_err)]
my_table.add_row(table_row)
print(my_table)

#print 'prob_ref_2sigma = %s'%(prob_ref_2sigma)
#print 'prob_mdm_2sigma = %s +/- %s'%(prob_mdm_2sigma,prob_mdm_2sigma_err)
#print 'prob_dark_2sigma = %s +/- %s'%(prob_dark_2sigma,prob_dark_2sigma_err)
#print 'prob_raw_2sigma = %s +/- %s'%(prob_raw_2sigma,prob_raw_2sigma_err)
#print 'prob_ref_3sigma = %s'%(prob_ref_3sigma)
#print 'prob_mdm_3sigma = %s +/- %s'%(prob_mdm_3sigma,prob_mdm_3sigma_err)
#print 'prob_dark_3sigma = %s +/- %s'%(prob_dark_3sigma,prob_dark_3sigma_err)
#print 'prob_raw_3sigma = %s +/- %s'%(prob_raw_3sigma,prob_raw_3sigma_err)
#print 'prob_ref_4sigma = %s'%(prob_ref_4sigma)
#print 'prob_mdm_4sigma = %s +/- %s'%(prob_mdm_4sigma,prob_mdm_4sigma_err)
#print 'prob_dark_4sigma = %s +/- %s'%(prob_dark_4sigma,prob_dark_4sigma_err)
#print 'prob_raw_4sigma = %s +/- %s'%(prob_raw_4sigma,prob_raw_4sigma_err)
