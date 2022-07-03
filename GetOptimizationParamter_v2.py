
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

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")
np.set_printoptions(precision=2)




sample_list = []
sample_list += ['3C273V6_OFF']
sample_list += ['3C273V5_OFF']
sample_list += ['1ES0502V6_OFF']
sample_list += ['1ES0502V5_OFF']
sample_list += ['DracoV6_OFF']
sample_list += ['DracoV5_OFF']
sample_list += ['1ES0647V6_OFF']
sample_list += ['1ES1011V6_OFF']
sample_list += ['PKS1424V6_OFF']
sample_list += ['PKS1424V5_OFF']
sample_list += ['H1426V6_OFF']
sample_list += ['NGC1275V6_OFF']
sample_list += ['OJ287V6_OFF']
sample_list += ['UrsaMajorIIV6_OFF']
sample_list += ['BLLacV6_OFF']
sample_list += ['BLLacV5_OFF']
sample_list += ['RGBJ0710V5_OFF']
sample_list += ['1ES0229V6_OFF']
sample_list += ['1ES0229V5_OFF']
sample_list += ['1ES0414V5_OFF']
sample_list += ['PG1553V6_OFF']
sample_list += ['PG1553V5_OFF']
sample_list += ['Segue1V6_OFF']
sample_list += ['Segue1V5_OFF']
sample_list += ['M82V6_OFF']
sample_list += ['M82V5_OFF']
sample_list += ['SNR_G150p3Plus04p5_V6_OFF']

folder_path = CommonPlotFunctions.folder_path
elev_range = CommonPlotFunctions.elev_range
theta2_bins = [0,4]
skymap_zoomin_scale = CommonPlotFunctions.skymap_zoomin_scale
smooth_size_spectroscopy = CommonPlotFunctions.smooth_size_spectroscopy
Skymap_size = CommonPlotFunctions.Skymap_size
Skymap_nbins = CommonPlotFunctions.Skymap_nbins
energy_bin = CommonPlotFunctions.energy_bin
calibration_radius = CommonPlotFunctions.calibration_radius
elev_range = CommonPlotFunctions.elev_range
energy_index_scale = CommonPlotFunctions.energy_index_scale
Smoothing = CommonPlotFunctions.Smoothing
Skymap_normalization_nbins = CommonPlotFunctions.Skymap_normalization_nbins

method_tag = 'tight_mdm_default'
lowrank_tag = '_svd'
method_tag += lowrank_tag

ONOFF_tag = 'OFF'
ONOFF_tag += '_Model0'

root_file_tags = []
for elev in range(0,len(elev_range)-1):
    elev_tag = '_TelElev%sto%s'%(elev_range[elev],elev_range[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        root_file_tags += [method_tag+elev_tag+theta2_tag+'_'+ONOFF_tag]

def FindSourceIndex(source_name):
    for source in range(0,len(sample_list)):
        if source_name==sample_list[source]:
            return source
    return 0

def GetGammaCounts(file_path,ebin):

    InputFile = ROOT.TFile(file_path)
    NewInfoTree = InputFile.Get("NewInfoTree")
    data_gamma_count = ROOT.std.vector("double")(10)
    bkgd_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank0_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank1_gamma_count = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)
    NewInfoTree.SetBranchAddress('data_gamma_count',ROOT.AddressOf(data_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_gamma_count',ROOT.AddressOf(bkgd_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank0_gamma_count',ROOT.AddressOf(bkgd_rank0_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank1_gamma_count',ROOT.AddressOf(bkgd_rank1_gamma_count))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.GetEntry(0)

    return data_gamma_count[ebin], bkgd_gamma_count[ebin], dark_gamma_count[ebin]

def GetOptimizationHistogram(file_path,ebin,hist_optimization):

    InputFile = ROOT.TFile(file_path)
    ErecS_lower_cut = energy_bin[ebin]
    ErecS_upper_cut = energy_bin[ebin+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_optimization.Add(InputFile.Get(HistName))

def MakeMultipleFitPlot(ax,Hists,legends,colors,title_x,title_y):

    func_gauss = []
    for h in range(0,len(Hists)):
        func_gauss += [ROOT.TF1("func_gauss_%s"%(h),"[0]*TMath::Gaus(x,[1],[2]*pow(2,0.5))",-1.0,1.0)]
        func_gauss[h].SetParameters(Hists[h].Integral(),0,0.02)
        Hists[h].Fit("func_gauss_%s"%(h),"N")
        func_gauss[h].SetLineColor(colors[h])
        func_gauss[h].Draw("E same")
        #print ("func_gauss[%s].GetNDF() = %s"%(h,func_gauss[h].GetNDF()))
        #print ("chi2/NDF = %s"%(func_gauss[h].GetChisquare()/func_gauss[h].GetNDF()))

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

def MakeMultiplePlot(ax,Hists,colors,title_x,title_y,name,y_min,y_max,logx,logy):
    
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

    for entry in range(1,len(Hists)):
        ax.plot(hist_xdata[entry], hist_ydata[entry], alpha=0.5)
    ax.errorbar(hist_xdata[0], hist_ydata[0], hist_error[0], color='k', marker='s', ls='none')

    ax.set_xlabel(title_x)
    ax.set_ylabel(title_y)
    return(ax)

def LoopOverFiles():

    global FilePath_Folder
    global Hist_Bkgd_Optimization
    global data_count
    global bkgd_count
    global dark_count
    global n_measurements


    n_measurements = 0
    for source in range(0,len(sample_list)):
        source_idx = FindSourceIndex(sample_list[source])
        for elev in range(0,len(root_file_tags)):
            file_exists = True
            n_groups = 0
            while file_exists:
                SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(n_groups,0,0)+".root"
                print ('Read file: %s'%(SourceFilePath))
                if os.path.exists(SourceFilePath):
                    n_groups += 1
                    print ('file exists.')
                else:
                    file_exists = False
                    print ('file does not exist.')
            for x_idx in range(0,Skymap_normalization_nbins):
                for y_idx in range(0,Skymap_normalization_nbins):
                    for g_idx in range(0,n_groups):
                        SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+"_G%d_X%d_Y%d"%(g_idx,x_idx,y_idx)+".root"
                        FilePath_Folder += [SourceFilePath]
                        print ('Get %s...'%(FilePath_Folder[len(FilePath_Folder)-1]))
                        if not os.path.isfile(FilePath_Folder[len(FilePath_Folder)-1]): 
                            print ('Found no file!!')
                            continue
                        else:
                            print ('Found a file.')
                        Hist_Bkgd_Optimization_E = []
                        data_count_E = []
                        bkgd_count_E = []
                        dark_count_E = []
                        for eb in range(0,len(energy_bin)-1):
                            Hist_Bkgd_Optimization_E += [ROOT.TH1D("Hist_Bkgd_Optimization_M%s_E%s"%(n_measurements,eb),"",20,optimiz_lower,optimiz_upper)]
                            Hist_Bkgd_Optimization_E[eb].Reset()
                            GetOptimizationHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Bkgd_Optimization_E[eb])
                            data, bkgd, dark = GetGammaCounts(FilePath_Folder[len(FilePath_Folder)-1],eb)
                            data_count_E += [data]
                            bkgd_count_E += [bkgd]
                            dark_count_E += [dark]
                        Hist_Bkgd_Optimization += [Hist_Bkgd_Optimization_E]
                        data_count += [data_count_E]
                        bkgd_count += [bkgd_count_E]
                        dark_count += [dark_count_E]
                        n_measurements += 1

optimiz_lower = -7.
optimiz_upper = -0.

n_measurements = 0
FilePath_Folder = []
Hist_Bkgd_Optimization = []
data_count = []
bkgd_count = []
dark_count = []
LoopOverFiles()

Hist_Bkgd_Optimization_Mean = ROOT.TH1D("Hist_Bkgd_Optimization_Mean","",20,optimiz_lower,optimiz_upper)
for eb in range(0,len(energy_bin)-1):
    Hist_Bkgd_Optimization_Mean.Reset()
    for binx in range(0,Hist_Bkgd_Optimization_Mean.GetNbinsX()):
        mean = 0.
        for entry in range(0,n_measurements-1):
            error = Hist_Bkgd_Optimization[entry][eb].GetBinContent(binx+1)
            mean += error
        mean = mean/float(n_measurements-1)
        rms = 0.
        for entry in range(0,n_measurements-1):
            error = Hist_Bkgd_Optimization[entry][eb].GetBinContent(binx+1)
            rms += pow(error-mean,2)
        rms = rms/float(n_measurements-1)
        rms = pow(rms,0.5)
        Hist_Bkgd_Optimization_Mean.SetBinContent(binx+1,mean)
        Hist_Bkgd_Optimization_Mean.SetBinError(binx+1,rms)
    Hists = []
    colors = []
    Hists += [Hist_Bkgd_Optimization_Mean]
    colors += [1]
    #random_gen = ROOT.TRandom3()
    #for entry in range(0,n_measurements-1):
    #    Hists += [Hist_Bkgd_Optimization[entry][eb]]
    #    colors += [int(random_gen.Uniform(29.,49.))]
    ax.cla()
    MakeMultiplePlot(ax,Hists,colors,'log10 c','relative error','OptimizationAlpha_E%s%s'%(eb,folder_path),1e-3,0.1,False,False)
    fig.savefig("output_plots/OptimizationAlpha_E%s%s.png"%(eb,folder_path))

Hist_SystErrDist_MDM = []
Hist_SystErrDist_Init = []
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E0","",10,-0.2,0.2)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E0","",10,-0.2,0.2)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E1","",10,-0.2,0.2)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E1","",10,-0.2,0.2)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E3","",10,-0.5,0.5)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E3","",10,-0.5,0.5)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E4","",10,-1.,1.)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E4","",10,-1.,1.)]
for eb in range(0,len(energy_bin)-1):
    Hist_SystErrDist_MDM[eb].Reset()
    Hist_SystErrDist_Init[eb].Reset()
    for entry in range(0,n_measurements-1):
        if data_count[entry][eb]==0.: continue
        Hist_SystErrDist_MDM[eb].Fill(1.-bkgd_count[entry][eb]/data_count[entry][eb])
        Hist_SystErrDist_Init[eb].Fill(1.-dark_count[entry][eb]/data_count[entry][eb])
    Hists = []
    legends = []
    colors = []
    Hists += [Hist_SystErrDist_MDM[eb]]
    legends += ['MIBE']
    colors += [1]
    Hists += [Hist_SystErrDist_Init[eb]]
    legends += ['Init']
    colors += [2]
    ax.cla()
    MakeMultipleFitPlot(ax,Hists,legends,colors,'relative error','number of measurements')
    fig.savefig("output_plots/SystErrDist_E%s.png"%(eb))

