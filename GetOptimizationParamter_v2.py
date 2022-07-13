
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
np.set_printoptions(precision=4)


analysis_type = 'inclusive'
#analysis_type = '1ES0229 OFF'
#analysis_type = '1ES0229 Imposter'
#analysis_type = 'UrsaMajorII OFF'
#analysis_type = 'UrsaMajorII Imposter'
#analysis_type = 'MGRO_J1908 Imposter'
#analysis_type = 'MGRO_J2019 Imposter'
#imposter_ID = 2

sample_list = []

if analysis_type=='inclusive':
    ONOFF_tag = 'OFF'
    sample_list += ['1ES0502V6_OFF']
    sample_list += ['1ES0502V5_OFF']
    sample_list += ['DracoV6_OFF']
    sample_list += ['DracoV5_OFF']
    sample_list += ['BLLacV6_OFF']
    sample_list += ['BLLacV5_OFF']
    sample_list += ['1ES0414V5_OFF']
    sample_list += ['M82V6_OFF']
    sample_list += ['M82V5_OFF']
    sample_list += ['1ES0647V6_OFF']
    sample_list += ['1ES1011V6_OFF']
    sample_list += ['NGC1275V6_OFF']
    sample_list += ['OJ287V6_OFF']
    sample_list += ['RGBJ0710V5_OFF']
    sample_list += ['Segue1V6_OFF']
    sample_list += ['Segue1V5_OFF']
    sample_list += ['3C264V6_OFF']
    sample_list += ['3C273V6_OFF']
    sample_list += ['3C273V5_OFF']
    sample_list += ['PG1553V6_OFF']
    sample_list += ['PG1553V5_OFF']
    #sample_list += ['H1426V6_OFF']
    sample_list += ['UrsaMajorIIV6_OFF']
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
    sample_list += ['PKS1424V6_OFF']
    sample_list += ['PKS1424V5_OFF']
if analysis_type=='1ES0229 OFF':
    ONOFF_tag = 'OFF'
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
if analysis_type=='1ES0229 Imposter':
    ONOFF_tag = 'ON'
    sample_list += ['1ES0229V6_Imposter%s'%(imposter_ID)]
    sample_list += ['1ES0229V5_Imposter%s'%(imposter_ID)]
if analysis_type=='UrsaMajorII OFF':
    ONOFF_tag = 'OFF'
    sample_list += ['UrsaMajorIIV6_OFF']
if analysis_type=='UrsaMajorII Imposter':
    ONOFF_tag = 'ON'
    sample_list += ['UrsaMajorIIV6_Imposter%s'%(imposter_ID)]
if analysis_type=='MGRO_J1908 Imposter':
    ONOFF_tag = 'ON'
    for imposter_ID in range(0,5):
        sample_list += ['MGRO_J1908_V6_Imposter%s'%(imposter_ID)]
        sample_list += ['MGRO_J1908_V5_Imposter%s'%(imposter_ID)]
if analysis_type=='MGRO_J2019 Imposter':
    ONOFF_tag = 'ON'
    for imposter_ID in range(0,5):
        sample_list += ['MGRO_J2019_V6_Imposter%s'%(imposter_ID)]
        sample_list += ['MGRO_J2019_V5_Imposter%s'%(imposter_ID)]


#ONOFF_tag = 'ON'
#sample_list += ['UrsaMajorIIV6_Imposter1']
#sample_list += ['UrsaMajorIIV6_Imposter2']
#sample_list += ['UrsaMajorIIV6_Imposter3']
#sample_list += ['UrsaMajorIIV6_Imposter4']
#sample_list += ['UrsaMajorIIV6_Imposter5']
#
#ONOFF_tag = 'ON'
#sample_list += ['PG1553V6_Imposter1']
#sample_list += ['PG1553V6_Imposter2']
#sample_list += ['PG1553V6_Imposter3']
#sample_list += ['PG1553V6_Imposter4']
#sample_list += ['PG1553V6_Imposter5']
#sample_list += ['PG1553V5_Imposter1']
#sample_list += ['PG1553V5_Imposter2']
#sample_list += ['PG1553V5_Imposter3']
#sample_list += ['PG1553V5_Imposter4']
#sample_list += ['PG1553V5_Imposter5']
#
#ONOFF_tag = 'ON'
#sample_list += ['3C273V6_Imposter1']
#sample_list += ['3C273V6_Imposter2']
#sample_list += ['3C273V6_Imposter3']
#sample_list += ['3C273V6_Imposter4']
#sample_list += ['3C273V6_Imposter5']
#sample_list += ['3C273V5_Imposter1']
#sample_list += ['3C273V5_Imposter2']
#sample_list += ['3C273V5_Imposter3']
#sample_list += ['3C273V5_Imposter4']
#sample_list += ['3C273V5_Imposter5']
#

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

def GetMatrixCoefficients(hist_mtx,gamma_count):

    mtx_CDE = []
    for row in range(0,stable_rank):
        for col in range(0,stable_rank):
            mtx_CDE += [0.]
    if gamma_count<100.: return mtx_CDE
    for row in range(0,stable_rank):
        for col in range(0,stable_rank):
            idx = row*stable_rank+col
            content = hist_mtx.GetBinContent(row+1,col+1)
            mtx_CDE[idx] = content
    return mtx_CDE

def GetGammaCounts(file_path,ebin):

    dark_stable_rank = ROOT.std.vector("int")(10)
    data_gamma_count = ROOT.std.vector("double")(10)
    bkgd_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank0_gamma_count = ROOT.std.vector("double")(10)
    bkgd_rank1_gamma_count = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)

    InputFile = ROOT.TFile(file_path)
    NewInfoTree = InputFile.Get("NewInfoTree")
    NewInfoTree.SetBranchAddress('dark_stable_rank',ROOT.AddressOf(dark_stable_rank))
    NewInfoTree.SetBranchAddress('data_gamma_count',ROOT.AddressOf(data_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_gamma_count',ROOT.AddressOf(bkgd_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank0_gamma_count',ROOT.AddressOf(bkgd_rank0_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_rank1_gamma_count',ROOT.AddressOf(bkgd_rank1_gamma_count))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.GetEntry(0)

    print ('Open %s, E%s, data_gamma_count = %s'%(file_path,ebin,data_gamma_count[ebin]))

    return dark_stable_rank[ebin], data_gamma_count[ebin], bkgd_gamma_count[ebin], dark_gamma_count[ebin]

def GetCoefficientHistogram(file_path,ebin,hist_data,hist_bkgd):

    InputFile = ROOT.TFile(file_path)
    ErecS_lower_cut = energy_bin[ebin]
    ErecS_upper_cut = energy_bin[ebin+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_data.Add(InputFile.Get(HistName))
    HistName = "Hist_Coeff_Bkgd_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_bkgd.Add(InputFile.Get(HistName))

def GetOptimizationHistogram(file_path,ebin,hist_optimization):

    InputFile = ROOT.TFile(file_path)
    ErecS_lower_cut = energy_bin[ebin]
    ErecS_upper_cut = energy_bin[ebin+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_optimization.Add(InputFile.Get(HistName))

def gaussian_func(x,a,b,c):
    return a*np.exp(-0.5*pow((x-b)/c,2))

def MakeMultipleFitPlot(ax,Hists,legends,colors,title_x,title_y):

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
            error += [max(1.,Hists[entry].GetBinError(binx+1))]
        hist_xdata += [xdata]
        hist_ydata += [ydata]
        hist_error += [error]

    cycol = cycle('brgcmk')
    for entry in range(0,len(Hists)):
        next_color = next(cycol)
        #start = (50., 0., 0.1)
        #popt, pcov = curve_fit(gaussian_func, np.array(hist_xdata[entry]), np.array(hist_ydata[entry]), p0=start, sigma=np.array(hist_error[entry]))
        #xaxis = np.linspace(hist_xdata[entry][0],hist_xdata[entry][len(hist_xdata[entry])-1],50)
        #ax.plot(xaxis, gaussian_func(np.array(xaxis), *popt),color=next_color)
        #ax.errorbar(hist_xdata[entry], hist_ydata[entry], hist_error[entry], color=next_color, marker='s', ls='none', label='%s, $\mu = %0.3f$, $\sigma = %0.3f$'%(legends[entry],popt[1],popt[2]))
        ax.errorbar(hist_xdata[entry], hist_ydata[entry], hist_error[entry], color=next_color, marker='s', ls='none', label='%s'%(legends[entry]))

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

def PrincipalComponentAnalysis(list_var, ebin):

    n_variables = len(list_var)
    n_samples = len(mtx_CDE_data)
    mtx_var_data = np.zeros((n_samples,n_variables))
    mtx_var_bkgd = np.zeros((n_samples,n_variables))
    chi2 = np.zeros(n_variables)
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            row = list_var[var][0]-1
            col = list_var[var][1]-1
            idx = row*stable_rank+col
            mtx_var_data[sample][var] = mtx_CDE_data[sample][ebin][idx]
            mtx_var_bkgd[sample][var] = mtx_CDE_bkgd[sample][ebin][idx]
            chi2[var] += pow(mtx_var_data[sample][var]-mtx_var_bkgd[sample][var],2)

    #mtx_var_centered = mtx_var - np.mean(mtx_var , axis = 0)
    mtx_var_square = np.square(mtx_var_data)
    mtx_var_mean = np.mean(mtx_var_square , axis = 0)
    mtx_var_rms = np.sqrt(mtx_var_mean)
    for var in range(0,n_variables):
        chi2[var] = chi2[var]/mtx_var_rms[var]

    mtx_var_data_norm = np.zeros((n_samples,n_variables))
    mtx_var_bkgd_norm = np.zeros((n_samples,n_variables))
    n_good_samples = 0
    is_good_sample = []
    for sample in range(0,n_samples):
        is_good_sample += [True]
        chi2 = 0.
        for var in range(0,n_variables):
            mtx_var_data_norm[sample][var] = mtx_var_data[sample][var]/mtx_var_rms[var]
            mtx_var_bkgd_norm[sample][var] = mtx_var_bkgd[sample][var]/mtx_var_rms[var]
            chi2 += pow(mtx_var_data_norm[sample][var],2)
            if mtx_var_data_norm[sample][var]==0.:
                is_good_sample[sample] = False
        #if data_count[sample][ebin]<500.:
        #    is_good_sample[sample] = False
        if chi2>25.:
            is_good_sample[sample] = False
        if is_good_sample[sample]:
            n_good_samples += 1

    mtx_var_data_good = np.zeros((n_good_samples,n_variables))
    mtx_var_bkgd_good = np.zeros((n_good_samples,n_variables))
    good_sample = 0
    for sample in range(0,n_samples):
        if is_good_sample[sample]:
            for var in range(0,n_variables):
                mtx_var_data_good[good_sample][var] = mtx_var_data_norm[sample][var]
                mtx_var_bkgd_good[good_sample][var] = mtx_var_bkgd_norm[sample][var]
                if mtx_var_data_good[good_sample][var]==0.:
                    mtx_var_data_good[good_sample][0] = 0
                    mtx_var_data_good[good_sample][1] = 0
                    mtx_var_bkgd_good[good_sample][0] = 0
                    mtx_var_bkgd_good[good_sample][1] = 0
            good_sample += 1

    mtx_cov = np.cov(mtx_var_data_good, rowvar = False)
    eigen_values , eigen_vectors = np.linalg.eigh(mtx_cov)
    print ('E%s'%(ebin))
    if eigen_values[0]/eigen_values[n_variables-1]<1.0:
        print ('eigen_values = \n {0}'.format(eigen_values))
        print ('eigen_values ratio = \n {0}'.format(eigen_values[0]/eigen_values[n_variables-1]))
        print ('mtx_var_rms = \n {0}'.format(mtx_var_rms))
        print ('chi2 = \n {0}'.format(chi2))
        print ('primary eigen_vectors = ')
        for var in range(0,n_variables):
            print ('({1},{2}) {0}'.format(eigen_vectors[n_variables-1][var],list_var[var][0],list_var[var][1]))
        #print ('sencondary eigen_vectors = ')
        #for var in range(0,n_variables):
        #    print ('({1},{2}) {0}'.format(eigen_vectors[n_variables-2][var],list_var[var][0],list_var[var][1]))
    #return eigen_values[n_variables-1]/eigen_values[n_variables-2], chi2

def MakeCorrelationPlot(list_var,ebin):

    par1_row = list_var[0][0]
    par1_col = list_var[0][1]
    par2_row = list_var[1][0]
    par2_col = list_var[1][1]

    n_variables = len(list_var)
    n_samples = len(mtx_CDE_data)
    mtx_var_data = np.zeros((n_samples,n_variables))
    mtx_var_bkgd = np.zeros((n_samples,n_variables))
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            row = list_var[var][0]-1
            col = list_var[var][1]-1
            idx = row*stable_rank+col
            mtx_var_data[sample][var] = mtx_CDE_data[sample][ebin][idx]
            mtx_var_bkgd[sample][var] = mtx_CDE_bkgd[sample][ebin][idx]

    #mtx_var_centered = mtx_var - np.mean(mtx_var , axis = 0)
    mtx_var_square = np.square(mtx_var_data)
    mtx_var_mean = np.mean(mtx_var_square , axis = 0)
    mtx_var_rms = np.sqrt(mtx_var_mean)

    mtx_var_data_norm = np.zeros((n_samples,n_variables))
    mtx_var_bkgd_norm = np.zeros((n_samples,n_variables))
    n_good_samples = 0
    is_good_sample = []
    for sample in range(0,n_samples):
        is_good_sample += [True]
        chi2 = 0.
        for var in range(0,n_variables):
            mtx_var_data_norm[sample][var] = mtx_var_data[sample][var]/mtx_var_rms[var]
            mtx_var_bkgd_norm[sample][var] = mtx_var_bkgd[sample][var]/mtx_var_rms[var]
            chi2 += pow(mtx_var_data_norm[sample][var],2)
            if mtx_var_data_norm[sample][var]==0.:
                is_good_sample[sample] = False
        #if data_count[sample][ebin]<500.:
        #    is_good_sample[sample] = False
        if chi2>25.:
            is_good_sample[sample] = False
        if is_good_sample[sample]:
            n_good_samples += 1

    mtx_var_data_good = np.zeros((n_good_samples,n_variables))
    mtx_var_bkgd_good = np.zeros((n_good_samples,n_variables))
    good_sample = 0
    for sample in range(0,n_samples):
        if is_good_sample[sample]:
            for var in range(0,n_variables):
                mtx_var_data_good[good_sample][var] = mtx_var_data_norm[sample][var]
                mtx_var_bkgd_good[good_sample][var] = mtx_var_bkgd_norm[sample][var]
                if mtx_var_data_good[good_sample][var]==0.:
                    mtx_var_data_good[good_sample][0] = 0
                    mtx_var_data_good[good_sample][1] = 0
                    mtx_var_bkgd_good[good_sample][0] = 0
                    mtx_var_bkgd_good[good_sample][1] = 0
            good_sample += 1

    plt.clf()
    x_var = mtx_var_data_good.transpose()[0]
    y_var = mtx_var_data_good.transpose()[1]
    plt.xlabel("$t_{%s,%s}$ (arbitrary unit)"%(par1_row,par1_col), fontsize=16)
    plt.ylabel("$t_{%s,%s}$ (arbitrary unit)"%(par2_row,par2_col), fontsize=16)
    plt.scatter(x_var,y_var,color='b',alpha=0.2)

    x_var = mtx_var_bkgd_good.transpose()[0]
    y_var = mtx_var_bkgd_good.transpose()[1]
    plt.scatter(x_var,y_var,color='r',alpha=0.2)

    plt.savefig("output_plots/par_correlation_%s%s_%s%s_E%s.png"%(par1_row,par1_col,par2_row,par2_col,ebin))


def LoopOverFiles():

    global FilePath_Folder
    global Hist_Bkgd_Optimization
    global data_rank
    global data_count
    global bkgd_count
    global dark_count
    global n_measurements
    global Hist_Coeff_Data
    global Hist_Coeff_Bkgd
    global mtx_CDE_data
    global mtx_CDE_bkgd

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
                        rank_E = []
                        data_count_E = []
                        bkgd_count_E = []
                        dark_count_E = []
                        Hist_Coeff_Data_E = []
                        Hist_Coeff_Bkgd_E = []
                        mtx_CDE_data_E = []
                        mtx_CDE_bkgd_E = []
                        for eb in range(0,len(energy_bin)-1):
                            Hist_Bkgd_Optimization_E += [ROOT.TH1D("Hist_Bkgd_Optimization_M%s_E%s"%(n_measurements,eb),"",optimiz_nbins,optimiz_lower[eb],optimiz_upper[eb])]
                            Hist_Bkgd_Optimization_E[eb].Reset()
                            GetOptimizationHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Bkgd_Optimization_E[eb])
                            rank, data, bkgd, dark = GetGammaCounts(FilePath_Folder[len(FilePath_Folder)-1],eb)
                            rank_E += [rank]
                            data_count_E += [data]
                            bkgd_count_E += [bkgd]
                            dark_count_E += [dark]
                            Hist_Coeff_Data_E += [ROOT.TH2D("Hist_Coeff_Data_M%s_E%s"%(n_measurements,eb),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
                            Hist_Coeff_Bkgd_E += [ROOT.TH2D("Hist_Coeff_Bkgd_M%s_E%s"%(n_measurements,eb),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
                            GetCoefficientHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Coeff_Data_E[eb],Hist_Coeff_Bkgd_E[eb])
                            mtx_data = GetMatrixCoefficients(Hist_Coeff_Data_E[eb],data)
                            mtx_bkgd = GetMatrixCoefficients(Hist_Coeff_Bkgd_E[eb],data)
                            mtx_CDE_data_E += [mtx_data]
                            mtx_CDE_bkgd_E += [mtx_bkgd]
                        Hist_Bkgd_Optimization += [Hist_Bkgd_Optimization_E]
                        data_rank += [rank_E]
                        data_count += [data_count_E]
                        bkgd_count += [bkgd_count_E]
                        dark_count += [dark_count_E]
                        Hist_Coeff_Data += [Hist_Coeff_Data_E]
                        Hist_Coeff_Bkgd += [Hist_Coeff_Bkgd_E]
                        mtx_CDE_data += [mtx_CDE_data_E]
                        mtx_CDE_bkgd += [mtx_CDE_bkgd_E]
                        n_measurements += 1

optimiz_lower = [-4.,-5.,-4.,-4.]
optimiz_upper = [-2.,-3.,-2.,-2.]
optimiz_nbins = 20
N_bins_for_deconv = 8
stable_rank = 3

n_measurements = 0
plot_rank = 2
FilePath_Folder = []
Hist_Bkgd_Optimization = []
data_rank = []
data_count = []
bkgd_count = []
dark_count = []
Hist_Coeff_Data = []
Hist_Coeff_Bkgd = []
mtx_CDE_data = []
mtx_CDE_bkgd = []

LoopOverFiles()

Hist_Bkgd_Optimization_Mean = []
for eb in range(0,len(energy_bin)-1):
    Hist_Bkgd_Optimization_Mean += [ROOT.TH1D("Hist_Bkgd_Optimization_Mean_E%s"%(eb),"",optimiz_nbins,optimiz_lower[eb],optimiz_upper[eb])]
for eb in range(0,len(energy_bin)-1):
    Hist_Bkgd_Optimization_Mean[eb].Reset()
    for binx in range(0,Hist_Bkgd_Optimization_Mean[eb].GetNbinsX()):
        mean = 0.
        n_entries = 0.
        for entry in range(0,n_measurements):
            if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
            if eb==0 and data_count[entry][eb]<1000.: continue
            if eb==1 and data_count[entry][eb]<500.: continue
            error = Hist_Bkgd_Optimization[entry][eb].GetBinContent(binx+1)
            mean += error
            n_entries += 1.
        mean = mean/n_entries
        rms = 0.
        n_entries = 0.
        for entry in range(0,n_measurements):
            if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
            if eb==0 and data_count[entry][eb]<1000.: continue
            if eb==1 and data_count[entry][eb]<500.: continue
            error = Hist_Bkgd_Optimization[entry][eb].GetBinContent(binx+1)
            rms += pow(error-mean,2)
            #rms += pow(error,2)
            n_entries += 1.
        rms = rms/n_entries
        rms = pow(rms,0.5)
        #Hist_Bkgd_Optimization_Mean[eb].SetBinContent(binx+1,mean)
        #Hist_Bkgd_Optimization_Mean[eb].SetBinError(binx+1,rms)
        Hist_Bkgd_Optimization_Mean[eb].SetBinContent(binx+1,rms)
    Hists = []
    colors = []
    Hists += [Hist_Bkgd_Optimization_Mean[eb]]
    colors += [1]
    #random_gen = ROOT.TRandom3()
    #for entry in range(0,n_measurements):
    #    Hists += [Hist_Bkgd_Optimization[entry][eb]]
    #    colors += [int(random_gen.Uniform(29.,49.))]
    ax.cla()
    MakeMultiplePlot(ax,Hists,colors,'log10 c','relative error','OptimizationAlpha_E%s_%s'%(eb,folder_path),1e-3,0.1,False,False)
    fig.savefig("output_plots/OptimizationAlpha_E%s_%s.png"%(eb,folder_path))

Hist_SystErrDist_MDM = []
Hist_SystErrDist_Init = []
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E0","",21,-0.2,0.2)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E0","",21,-0.2,0.2)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E1","",21,-0.2,0.2)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E1","",21,-0.2,0.2)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E3","",21,-0.4,0.4)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E3","",21,-0.4,0.4)]
Hist_SystErrDist_MDM += [ROOT.TH1D("Hist_SystErrDist_MDM_E4","",21,-0.8,0.8)]
Hist_SystErrDist_Init += [ROOT.TH1D("Hist_SystErrDist_Init_E4","",21,-0.8,0.8)]
for eb in range(0,len(energy_bin)-1):
    Hist_SystErrDist_MDM[eb].Reset()
    Hist_SystErrDist_Init[eb].Reset()
    MDM_mean = 0.
    Init_mean = 0.
    MDM_rms = 0.
    Init_rms = 0.
    n_entries = 0.
    for entry in range(0,n_measurements):
        if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
        if data_count[entry][eb]==0.: continue
        if eb==0 and data_count[entry][eb]<1000.: continue
        if eb==1 and data_count[entry][eb]<500.: continue
        Hist_SystErrDist_MDM[eb].Fill(1.-bkgd_count[entry][eb]/data_count[entry][eb])
        Hist_SystErrDist_Init[eb].Fill(1.-dark_count[entry][eb]/data_count[entry][eb])
        MDM_mean += 1.-bkgd_count[entry][eb]/data_count[entry][eb]
        Init_mean += 1.-dark_count[entry][eb]/data_count[entry][eb]
        MDM_rms += pow(1.-bkgd_count[entry][eb]/data_count[entry][eb],2)
        Init_rms += pow(1.-dark_count[entry][eb]/data_count[entry][eb],2)
        n_entries += 1.
    MDM_mean = MDM_mean/n_entries
    Init_mean = Init_mean/n_entries
    MDM_rms = pow(MDM_rms/n_entries,0.5)
    Init_rms = pow(Init_rms/n_entries,0.5)
    Hists = []
    legends = []
    colors = []
    Hists += [Hist_SystErrDist_MDM[eb]]
    legends += ['%0.2f-%0.2f TeV, mean = %0.3f, sigma = %0.3f'%(energy_bin[eb]/1000.,energy_bin[eb+1]/1000.,MDM_mean,MDM_rms)]
    #legends += ['matrix method, %0.2f-%0.2f TeV, mean = %0.3f, sigma = %0.3f'%(energy_bin[eb]/1000.,energy_bin[eb+1]/1000.,MDM_mean,MDM_rms)]
    colors += [1]
    Hists += [Hist_SystErrDist_Init[eb]]
    legends += ['initial matching, mean = %0.3f, sigma = %0.3f'%(Init_mean,Init_rms)]
    colors += [2]
    ax.cla()
    MakeMultipleFitPlot(ax,Hists,legends,colors,'relative error $\epsilon$','number of measurements')
    fig.savefig("output_plots/SystErrDist_E%s_%s.png"%(eb,folder_path))

list_var_pair = []
good_var_pair = []
good_eigenvalue = []
for eb in range(2,3):
    for row1 in range(0,stable_rank):
        for col1 in range(0,stable_rank):
            for row2 in range(0,stable_rank):
                for col2 in range(0,stable_rank):
                    idx1 = row1*stable_rank+col1
                    idx2 = row2*stable_rank+col2
                    list_var_pair = [[row1+1,col1+1]]
                    list_var_pair += [[row2+1,col2+1]]
                    if idx1<idx2:
                        print('=======================================================')
                        PrincipalComponentAnalysis(list_var_pair,eb)
                        #if math.isnan(max_eigenvalue): continue
                        #chi2_ratio = min(chi2[0]/chi2[1],chi2[1]/chi2[0])
                        MakeCorrelationPlot(list_var_pair,eb)
                        #good_var_pair += [list_var_pair]
                        #good_eigenvalue += [max_eigenvalue]
