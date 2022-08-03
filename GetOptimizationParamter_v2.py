
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

N_bins_for_deconv = CommonPlotFunctions.N_bins_for_deconv
gamma_hadron_dim_ratio_w = CommonPlotFunctions.gamma_hadron_dim_ratio_w
gamma_hadron_dim_ratio_l = CommonPlotFunctions.gamma_hadron_dim_ratio_l
MSCW_blind_cut = CommonPlotFunctions.MSCW_blind_cut
MSCL_blind_cut = CommonPlotFunctions.MSCL_blind_cut

n_measures_per_entry = 4

analysis_type = 'inclusive'
#analysis_type = 'SelectOFF'
#analysis_type='UrsaMajorII'
#analysis_type='1ES0229'
#analysis_type='PG1553'

#observing_condition = 'all'
#observing_condition = 'north'
observing_condition = 'south'
#observing_condition = 'eastwest'
#observing_condition = 'sza'
#observing_condition = 'lza'
#observing_condition = 'hnsb'
#observing_condition = 'lnsb'

sample_list = []
ONOFF_tag_sample = 'OFF'
imposter_list = []
ONOFF_tag_imposter = 'OFF'

if analysis_type=='inclusive':
    ONOFF_tag_sample = 'OFF'
    sample_list += ['1ES0414V5_OFF']
    sample_list += ['1ES0647V6_OFF']
    sample_list += ['OJ287V6_OFF']
    sample_list += ['RGBJ0710V5_OFF']
    sample_list += ['Segue1V6_OFF']
    sample_list += ['Segue1V5_OFF']
    sample_list += ['3C264V6_OFF']
    sample_list += ['3C273V6_OFF']
    sample_list += ['3C273V5_OFF']
    sample_list += ['PG1553V6_OFF']
    sample_list += ['PG1553V5_OFF']
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
    sample_list += ['PKS1424V6_OFF']
    sample_list += ['PKS1424V5_OFF']
    sample_list += ['UrsaMajorIIV6_OFF']
    sample_list += ['H1426V6_OFF']
    sample_list += ['NGC1275V6_OFF']
    sample_list += ['DracoV6_OFF']
    sample_list += ['DracoV5_OFF']
    sample_list += ['BLLacV6_OFF']
    sample_list += ['BLLacV5_OFF']
    sample_list += ['1ES0502V6_OFF']
    sample_list += ['1ES0502V5_OFF']
    sample_list += ['M82V6_OFF']
    sample_list += ['M82V5_OFF']
    sample_list += ['1ES1011V6_OFF']
if analysis_type=='SelectOFF':
    ONOFF_tag_sample = 'OFF'
    sample_list += ['UrsaMajorIIV6_OFF']
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
    sample_list += ['PG1553V6_OFF']
    sample_list += ['PG1553V5_OFF']
    sample_list += ['PKS1424V6_OFF']
    sample_list += ['PKS1424V5_OFF']
    sample_list += ['3C273V6_OFF']
    sample_list += ['3C273V5_OFF']
    sample_list += ['Segue1V6_OFF']
    sample_list += ['Segue1V5_OFF']
    ONOFF_tag_imposter = 'ON'
    for imposter_ID in range(0,5):
        imposter_list += ['UrsaMajorIIV6_Imposter%s'%(imposter_ID)]
        imposter_list += ['1ES0229V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['1ES0229V5_Imposter%s'%(imposter_ID)]
        imposter_list += ['PG1553V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['PG1553V5_Imposter%s'%(imposter_ID)]
        imposter_list += ['PKS1424V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['PKS1424V5_Imposter%s'%(imposter_ID)]
        imposter_list += ['3C273V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['3C273V5_Imposter%s'%(imposter_ID)]
        imposter_list += ['Segue1V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['Segue1V5_Imposter%s'%(imposter_ID)]
if analysis_type=='UrsaMajorII':
    ONOFF_tag_sample = 'OFF'
    sample_list += ['UrsaMajorIIV6_OFF']
    ONOFF_tag_imposter = 'ON'
    for imposter_ID in range(0,5):
        imposter_list += ['UrsaMajorIIV6_Imposter%s'%(imposter_ID)]
if analysis_type=='1ES0229':
    ONOFF_tag_sample = 'OFF'
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
    ONOFF_tag_imposter = 'ON'
    for imposter_ID in range(0,5):
        imposter_list += ['1ES0229V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['1ES0229V5_Imposter%s'%(imposter_ID)]
if analysis_type=='PG1553':
    ONOFF_tag_sample = 'OFF'
    sample_list += ['PG1553V6_OFF']
    sample_list += ['PG1553V5_OFF']
    ONOFF_tag_imposter = 'ON'
    for imposter_ID in range(0,5):
        imposter_list += ['PG1553V6_Imposter%s'%(imposter_ID)]
        imposter_list += ['PG1553V5_Imposter%s'%(imposter_ID)]



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

ONOFF_tag_sample += '_Model0'
ONOFF_tag_imposter += '_Model0'

sample_file_tags = []
imposter_file_tags = []
for elev in range(0,len(elev_range)-1):
    elev_tag = '_TelElev%sto%s'%(elev_range[elev],elev_range[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        sample_file_tags += [method_tag+elev_tag+theta2_tag+'_'+ONOFF_tag_sample]
        imposter_file_tags += [method_tag+elev_tag+theta2_tag+'_'+ONOFF_tag_imposter]

def FindSourceIndex(source_name,the_list):
    for source in range(0,len(the_list)):
        if source_name==the_list[source]:
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

def GetObservingCondition(file_path):

    InputFile = ROOT.TFile(file_path)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    mean_elev = InfoTree.mean_elev
    mean_azim = InfoTree.mean_azim
    mean_nsb = InfoTree.mean_nsb

    print ('Open %s, mean_azim = %0.2f'%(file_path,mean_azim))

    return mean_elev, mean_azim, mean_nsb

def GetGammaCounts(file_path,ebin):

    dark_stable_rank = ROOT.std.vector("int")(10)
    data_gamma_count = ROOT.std.vector("double")(10)
    bkgd_gamma_count = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)
    rank0_gamma_count = ROOT.std.vector("double")(10)
    rank1_gamma_count = ROOT.std.vector("double")(10)
    rank2_gamma_count = ROOT.std.vector("double")(10)
    rank3_gamma_count = ROOT.std.vector("double")(10)

    InputFile = ROOT.TFile(file_path)
    NewInfoTree = InputFile.Get("NewInfoTree")
    NewInfoTree.SetBranchAddress('dark_stable_rank',ROOT.AddressOf(dark_stable_rank))
    NewInfoTree.SetBranchAddress('data_gamma_count',ROOT.AddressOf(data_gamma_count))
    NewInfoTree.SetBranchAddress('bkgd_gamma_count',ROOT.AddressOf(bkgd_gamma_count))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.SetBranchAddress('rank0_gamma_count',ROOT.AddressOf(rank0_gamma_count))
    NewInfoTree.SetBranchAddress('rank1_gamma_count',ROOT.AddressOf(rank1_gamma_count))
    NewInfoTree.SetBranchAddress('rank2_gamma_count',ROOT.AddressOf(rank2_gamma_count))
    NewInfoTree.SetBranchAddress('rank3_gamma_count',ROOT.AddressOf(rank3_gamma_count))
    NewInfoTree.GetEntry(0)

    print ('Open %s, E%s, data_gamma_count = %s'%(file_path,ebin,data_gamma_count[ebin]))

    return dark_stable_rank[ebin], data_gamma_count[ebin], bkgd_gamma_count[ebin], dark_gamma_count[ebin], rank0_gamma_count[ebin], rank1_gamma_count[ebin], rank2_gamma_count[ebin]

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

def GetShowerShapeHistogram(file_path,ebin,hist_data,hist_dark,hist_bkgd):

    InputFile = ROOT.TFile(file_path)
    ErecS_lower_cut = energy_bin[ebin]
    ErecS_upper_cut = energy_bin[ebin+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_data.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_dark.Add(InputFile.Get(HistName))
    HistName = "Hist_OnBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_bkgd.Add(InputFile.Get(HistName))

def GetSingularValueHistogram(file_path,ebin,hist_optimization):

    InputFile = ROOT.TFile(file_path)
    ErecS_lower_cut = energy_bin[ebin]
    ErecS_upper_cut = energy_bin[ebin+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    HistName = "Hist_Data_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    hist_optimization.Add(InputFile.Get(HistName))

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
        if ebin==0 or ebin==1:
            if plot_rank!=0 and data_rank[sample][ebin]!=plot_rank:
                is_good_sample[sample] = False
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

    plt.savefig("output_plots/par_correlation_%s%s_%s%s_E%s_%s_%s.png"%(par1_row,par1_col,par2_row,par2_col,ebin,analysis_type,observing_condition))


def LoopOverFiles():

    global FilePath_Folder
    global Hist_Data_Eigenvalues
    global Hist_Bkgd_Optimization
    global Hist_Data_ShowerShape
    global data_rank
    global data_count
    global bkgd_count
    global dark_count
    global rank0_count
    global rank1_count
    global rank2_count
    global imposter_data_rank
    global imposter_data_count
    global imposter_bkgd_count
    global imposter_dark_count
    global n_measurements
    global n_imposter_measurements
    global Hist_Coeff_Data
    global Hist_Coeff_Bkgd
    global mtx_CDE_data
    global mtx_CDE_bkgd

    n_measurements = 0
    for source in range(0,len(sample_list)):
        source_idx = FindSourceIndex(sample_list[source],sample_list)
        for elev in range(0,len(sample_file_tags)):
            file_exists = True
            n_groups = 0
            while file_exists:
                SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(sample_file_tags[elev])+"_G%d_X%d_Y%d"%(n_groups,0,0)+".root"
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
                        SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(sample_file_tags[elev])+"_G%d_X%d_Y%d"%(g_idx,x_idx,y_idx)+".root"
                        FilePath_Folder += [SourceFilePath]
                        print ('Get %s...'%(FilePath_Folder[len(FilePath_Folder)-1]))
                        if not os.path.isfile(FilePath_Folder[len(FilePath_Folder)-1]): 
                            print ('Found no file!!')
                            continue
                        else:
                            print ('Found a file.')
                        mean_elev, mean_azim, mean_nsb = GetObservingCondition(FilePath_Folder[len(FilePath_Folder)-1])
                        if observing_condition=='north':
                            if abs(mean_azim-0.)>45.: continue
                        if observing_condition=='south':
                            if abs(mean_azim-180.)>45.: continue
                        if observing_condition=='eastwest':
                            if abs(mean_azim-0.)<45. or abs(mean_azim-180.)<45.: continue
                        if observing_condition=='sza':
                            if mean_elev<65.: continue
                        if observing_condition=='lza':
                            if mean_elev>65.: continue
                        if observing_condition=='hnsb':
                            if mean_nsb<5.: continue
                        if observing_condition=='lnsb':
                            if mean_nsb>5.: continue
                        Hist_Data_Eigenvalues_E = []
                        Hist_Bkgd_Optimization_E = []
                        rank_E = []
                        data_count_E = []
                        bkgd_count_E = []
                        dark_count_E = []
                        rank0_count_E = []
                        rank1_count_E = []
                        rank2_count_E = []
                        Hist_Coeff_Data_E = []
                        Hist_Coeff_Bkgd_E = []
                        mtx_CDE_data_E = []
                        mtx_CDE_bkgd_E = []
                        for eb in range(0,len(energy_bin)-1):
                            GetShowerShapeHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Data_ShowerShape[eb],Hist_Dark_ShowerShape[eb],Hist_Bkgd_ShowerShape[eb])
                            Hist_Data_Eigenvalues_E += [ROOT.TH1D("Hist_Data_Eigenvalues_M%s_E%s"%(n_measurements,eb),"",N_bins_for_deconv,0,N_bins_for_deconv)]
                            Hist_Data_Eigenvalues_E[eb].Reset()
                            GetSingularValueHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Data_Eigenvalues_E[eb])
                            Hist_Bkgd_Optimization_E += [ROOT.TH2D("Hist_Bkgd_Optimization_M%s_E%s"%(n_measurements,eb),"",optimiz_nbins,optimiz_alpha_lower[eb],optimiz_alpha_upper[eb],optimiz_nbins,optimiz_beta_lower[eb],optimiz_beta_upper[eb])]
                            Hist_Bkgd_Optimization_E[eb].Reset()
                            GetOptimizationHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Bkgd_Optimization_E[eb])
                            rank, data, bkgd, dark, rank0, rank1, rank2 = GetGammaCounts(FilePath_Folder[len(FilePath_Folder)-1],eb)
                            rank_E += [rank]
                            data_count_E += [data]
                            bkgd_count_E += [bkgd]
                            dark_count_E += [dark]
                            rank0_count_E += [rank0]
                            rank1_count_E += [rank1]
                            rank2_count_E += [rank2]
                            Hist_Coeff_Data_E += [ROOT.TH2D("Hist_Coeff_Data_M%s_E%s"%(n_measurements,eb),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
                            Hist_Coeff_Bkgd_E += [ROOT.TH2D("Hist_Coeff_Bkgd_M%s_E%s"%(n_measurements,eb),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
                            GetCoefficientHistogram(FilePath_Folder[len(FilePath_Folder)-1],eb,Hist_Coeff_Data_E[eb],Hist_Coeff_Bkgd_E[eb])
                            mtx_data = GetMatrixCoefficients(Hist_Coeff_Data_E[eb],data)
                            mtx_bkgd = GetMatrixCoefficients(Hist_Coeff_Bkgd_E[eb],data)
                            mtx_CDE_data_E += [mtx_data]
                            mtx_CDE_bkgd_E += [mtx_bkgd]
                        Hist_Data_Eigenvalues += [Hist_Data_Eigenvalues_E]
                        Hist_Bkgd_Optimization += [Hist_Bkgd_Optimization_E]
                        data_rank += [rank_E]
                        data_count += [data_count_E]
                        bkgd_count += [bkgd_count_E]
                        dark_count += [dark_count_E]
                        rank0_count += [rank0_count_E]
                        rank1_count += [rank1_count_E]
                        rank2_count += [rank2_count_E]
                        Hist_Coeff_Data += [Hist_Coeff_Data_E]
                        Hist_Coeff_Bkgd += [Hist_Coeff_Bkgd_E]
                        mtx_CDE_data += [mtx_CDE_data_E]
                        mtx_CDE_bkgd += [mtx_CDE_bkgd_E]
                        n_measurements += 1

    n_imposter_measurements = 0
    for source in range(0,len(imposter_list)):
        source_idx = FindSourceIndex(imposter_list[source],imposter_list)
        for elev in range(0,len(imposter_file_tags)):
            file_exists = True
            n_groups = 0
            while file_exists:
                SourceFilePath = "%s/Netflix_"%(folder_path)+imposter_list[source_idx]+"_%s"%(imposter_file_tags[elev])+"_G%d_X%d_Y%d"%(n_groups,0,0)+".root"
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
                        SourceFilePath = "%s/Netflix_"%(folder_path)+imposter_list[source_idx]+"_%s"%(imposter_file_tags[elev])+"_G%d_X%d_Y%d"%(g_idx,x_idx,y_idx)+".root"
                        FilePath_Folder += [SourceFilePath]
                        print ('Get %s...'%(FilePath_Folder[len(FilePath_Folder)-1]))
                        if not os.path.isfile(FilePath_Folder[len(FilePath_Folder)-1]): 
                            print ('Found no file!!')
                            continue
                        else:
                            print ('Found a file.')
                        rank_E = []
                        data_count_E = []
                        bkgd_count_E = []
                        dark_count_E = []
                        for eb in range(0,len(energy_bin)-1):
                            rank, data, bkgd, dark, rank0, rank1, rank2 = GetGammaCounts(FilePath_Folder[len(FilePath_Folder)-1],eb)
                            rank_E += [rank]
                            data_count_E += [data]
                            bkgd_count_E += [bkgd]
                            dark_count_E += [dark]
                        imposter_data_rank += [rank_E]
                        imposter_data_count += [data_count_E]
                        imposter_bkgd_count += [bkgd_count_E]
                        imposter_dark_count += [dark_count_E]
                        n_imposter_measurements += 1


optimiz_alpha_lower = [-1.5,-1.5,-1.5,-1.5]
optimiz_alpha_upper = [1.5,1.5,1.5,1.5]
optimiz_beta_lower = [-1.5,-1.5,-1.5,-1.5]
optimiz_beta_upper = [1.5,1.5,1.5,1.5]
optimiz_nbins = 10
stable_rank = 2

n_measurements = 0
n_imposter_measurements = 0
plot_rank = 0
FilePath_Folder = []
Hist_Data_Eigenvalues = []
Hist_Bkgd_Optimization = []
data_rank = []
data_count = []
bkgd_count = []
dark_count = []
rank0_count = []
rank1_count = []
rank2_count = []
imposter_data_rank = []
imposter_data_count = []
imposter_bkgd_count = []
imposter_dark_count = []
Hist_Coeff_Data = []
Hist_Coeff_Bkgd = []
mtx_CDE_data = []
mtx_CDE_bkgd = []

MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))+MSCL_blind_cut
MSCW_plot_lower = -0.5*gamma_hadron_dim_ratio_w*(MSCW_blind_cut-(-1.*MSCW_blind_cut))-MSCW_blind_cut
MSCL_plot_lower = -0.5*gamma_hadron_dim_ratio_l*(MSCL_blind_cut-(-1.*MSCL_blind_cut))-MSCL_blind_cut
Hist_Data_ShowerShape = []
Hist_Dark_ShowerShape = []
Hist_Bkgd_ShowerShape = []
Hist_Diff_ShowerShape = []
for eb in range(0,len(energy_bin)-1):
    Hist_Data_ShowerShape += [ROOT.TH2D("Hist_Data_ShowerShape_E%s"%(eb),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_Dark_ShowerShape += [ROOT.TH2D("Hist_Dark_ShowerShape_E%s"%(eb),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_Bkgd_ShowerShape += [ROOT.TH2D("Hist_Bkgd_ShowerShape_E%s"%(eb),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_Diff_ShowerShape += [ROOT.TH2D("Hist_Diff_ShowerShape_E%s"%(eb),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]

LoopOverFiles()

for eb in range(0,len(energy_bin)-1):
    Hist_Diff_ShowerShape[eb].Reset()
    Hist_Diff_ShowerShape[eb].Add(Hist_Data_ShowerShape[eb])
    Hist_Diff_ShowerShape[eb].Add(Hist_Dark_ShowerShape[eb],-1.)
    CommonPlotFunctions.MatplotlibHist2D(Hist_Diff_ShowerShape[eb],fig,'scaled length','scaled width','count difference','DiffMatrixDark_E%s_%s_%s_%s'%(eb,analysis_type,observing_condition,folder_path))
    Hist_Diff_ShowerShape[eb].Reset()
    Hist_Diff_ShowerShape[eb].Add(Hist_Data_ShowerShape[eb])
    Hist_Diff_ShowerShape[eb].Add(Hist_Bkgd_ShowerShape[eb],-1.)
    CommonPlotFunctions.MatplotlibHist2D(Hist_Diff_ShowerShape[eb],fig,'scaled length','scaled width','count difference','DiffMatrixBkgd_E%s_%s_%s_%s'%(eb,analysis_type,observing_condition,folder_path))

Hist_Data_Eigenvalues_Mean = []
for eb in range(0,len(energy_bin)-1):
    Hist_Data_Eigenvalues_Mean += [ROOT.TH1D("Hist_Data_Eigenvalues_%s"%(eb),"",N_bins_for_deconv,0,N_bins_for_deconv)]
for eb in range(0,len(energy_bin)-1):
    Hist_Data_Eigenvalues_Mean[eb].Reset()
    for binx in range(0,Hist_Data_Eigenvalues_Mean[eb].GetNbinsX()):
        singularvalue_avg = 0.
        n_entries = 0.
        for entry in range(0,n_measurements):
            if eb==0 or eb==1:
                if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
            if eb==0 and data_count[entry][eb]<2000.: continue
            if eb==1 and data_count[entry][eb]<500.: continue
            if data_count[entry][eb]<10.: continue
            singularvalue = Hist_Data_Eigenvalues[entry][eb].GetBinContent(binx+1)
            singularvalue_avg += singularvalue
            n_entries += 1.
        singularvalue_avg = singularvalue_avg/n_entries
        Hist_Data_Eigenvalues_Mean[eb].SetBinContent(binx+1,singularvalue_avg)
energy_dependent_singularvalue = []
for eb in range(0,len(energy_bin)-1):
    singularvalue_array = []
    for binx in range(0,Hist_Data_Eigenvalues_Mean[eb].GetNbinsX()):
        singularvalue_array += [Hist_Data_Eigenvalues_Mean[eb].GetBinContent(binx+1)]
    energy_dependent_singularvalue += [singularvalue_array]
plt.clf()
fig, ax = plt.subplots()
plt.xlabel("rank $n$", fontsize=18)
plt.ylabel("singular value $\sigma_{n}$", fontsize=18)
plt.yscale('log')
for entry in range(0,len(energy_dependent_singularvalue)):
    plt.plot(energy_dependent_singularvalue[entry],marker='.',label='%s-%s GeV'%(energy_bin[entry],energy_bin[entry+1]))
ax.legend(loc='best')
plt.savefig("output_plots/MatrixSingularValue.png")

energy_array = np.array(energy_bin)
energy_dependent_rank0 = []
energy_dependent_rank1 = []
energy_dependent_rank2 = []
energy_dependent_stat = []
for eb in range(0,len(energy_bin)-1):
    epsilon_rank0_avg = 0.
    epsilon_rank1_avg = 0.
    epsilon_rank2_avg = 0.
    epsilon_stat_avg = 0.
    n_entries = 0.
    for entry in range(0,n_measurements):
        if eb==0 or eb==1:
            if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
        if eb==0 and data_count[entry][eb]<2000.: continue
        if eb==1 and data_count[entry][eb]<500.: continue
        if data_count[entry][eb]<10.: continue
        epsilon_rank0 = (data_count[entry][eb]-rank0_count[entry][eb])/data_count[entry][eb]
        epsilon_rank1 = (data_count[entry][eb]-rank1_count[entry][eb])/data_count[entry][eb]
        epsilon_rank2 = (data_count[entry][eb]-rank2_count[entry][eb])/data_count[entry][eb]
        epsilon_rank0_avg += epsilon_rank0*epsilon_rank0
        epsilon_rank1_avg += epsilon_rank1*epsilon_rank1
        epsilon_rank2_avg += epsilon_rank2*epsilon_rank2
        epsilon_stat_avg += pow(pow(data_count[entry][eb],0.5)/data_count[entry][eb],2)
        n_entries += 1.
    epsilon_rank0_avg = pow(epsilon_rank0_avg/n_entries,0.5)
    epsilon_rank1_avg = pow(epsilon_rank1_avg/n_entries,0.5)
    epsilon_rank2_avg = pow(epsilon_rank2_avg/n_entries,0.5)
    epsilon_stat_avg = pow(epsilon_stat_avg/n_entries,0.5)
    energy_dependent_rank0 += [epsilon_rank0_avg]
    energy_dependent_rank1 += [epsilon_rank1_avg]
    energy_dependent_rank2 += [epsilon_rank2_avg]
    energy_dependent_stat += [epsilon_stat_avg]
plt.clf()
fig, ax = plt.subplots()
plt.xlabel("Energy [GeV]", fontsize=18)
plt.ylabel("RMS of $\\epsilon$", fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank0,marker='.',label='$n \leq 0$')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank1,marker='.',label='$n \leq 1$')
#plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_rank2,marker='.',label='$n \leq 2$')
plt.plot(energy_array[0:len(energy_bin)-1],energy_dependent_stat,marker='.',label='stat. unc.')
ax.legend(loc='best')
plt.savefig("output_plots/RankErrors.png")


Hist_Bkgd_Optimization_Mean = []
for eb in range(0,len(energy_bin)-1):
    Hist_Bkgd_Optimization_Mean += [ROOT.TH2D("Hist_Bkgd_Optimization_Mean_E%s"%(eb),"",optimiz_nbins,optimiz_alpha_lower[eb],optimiz_alpha_upper[eb],optimiz_nbins,optimiz_beta_lower[eb],optimiz_beta_upper[eb])]
for eb in range(0,len(energy_bin)-1):
    Hist_Bkgd_Optimization_Mean[eb].Reset()
    for binx in range(0,Hist_Bkgd_Optimization_Mean[eb].GetNbinsX()):
        for biny in range(0,Hist_Bkgd_Optimization_Mean[eb].GetNbinsY()):
            rms = 0.
            n_entries = 0.
            for entry in range(0,n_measurements):
                if eb==0 or eb==1:
                    if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
                if eb==0 and data_count[entry][eb]<2000.: continue
                if eb==1 and data_count[entry][eb]<500.: continue
                if data_count[entry][eb]<10.: continue
                error = Hist_Bkgd_Optimization[entry][eb].GetBinContent(binx+1,biny+1)
                rms += pow(error,2)
                n_entries += 1.
            rms = rms/n_entries
            rms = pow(rms,0.5)
            Hist_Bkgd_Optimization_Mean[eb].SetBinContent(binx+1,biny+1,rms)
    CommonPlotFunctions.MatplotlibHist2D(Hist_Bkgd_Optimization_Mean[eb],fig,'log10 $\\beta_{11}$','log10 $\\beta_{22}$','RMS of $\\epsilon$','Optimization_E%s_%s_%s'%(eb,analysis_type,observing_condition))

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
Hist_Imposter_SystErrDist_MDM = []
Hist_Imposter_SystErrDist_Init = []
Hist_Imposter_SystErrDist_MDM += [ROOT.TH1D("Hist_Imposter_SystErrDist_MDM_E0","",21,-0.2,0.2)]
Hist_Imposter_SystErrDist_Init += [ROOT.TH1D("Hist_Imposter_SystErrDist_Init_E0","",21,-0.2,0.2)]
Hist_Imposter_SystErrDist_MDM += [ROOT.TH1D("Hist_Imposter_SystErrDist_MDM_E1","",21,-0.2,0.2)]
Hist_Imposter_SystErrDist_Init += [ROOT.TH1D("Hist_Imposter_SystErrDist_Init_E1","",21,-0.2,0.2)]
Hist_Imposter_SystErrDist_MDM += [ROOT.TH1D("Hist_Imposter_SystErrDist_MDM_E3","",21,-0.4,0.4)]
Hist_Imposter_SystErrDist_Init += [ROOT.TH1D("Hist_Imposter_SystErrDist_Init_E3","",21,-0.4,0.4)]
Hist_Imposter_SystErrDist_MDM += [ROOT.TH1D("Hist_Imposter_SystErrDist_MDM_E4","",21,-0.8,0.8)]
Hist_Imposter_SystErrDist_Init += [ROOT.TH1D("Hist_Imposter_SystErrDist_Init_E4","",21,-0.8,0.8)]

for eb in range(0,len(energy_bin)-1):
    Hist_SystErrDist_MDM[eb].Reset()
    Hist_SystErrDist_Init[eb].Reset()
    MDM_mean = 0.
    Init_mean = 0.
    MDM_rms = 0.
    Init_rms = 0.
    n_entries = 0.
    n_measures = 0.
    data_count_n_measures = 0.
    bkgd_count_n_measures = 0.
    dark_count_n_measures = 0.
    for entry in range(0,n_measurements):
        if eb==0 or eb==1:
            if plot_rank!=0 and data_rank[entry][eb]!=plot_rank: continue
        if data_count[entry][eb]==0.: continue
        if eb==0 and data_count[entry][eb]<2000.: continue
        if eb==1 and data_count[entry][eb]<500.: continue
        if data_count[entry][eb]<10.: continue
        n_measures += 1.
        data_count_n_measures += data_count[entry][eb]
        bkgd_count_n_measures += bkgd_count[entry][eb]
        dark_count_n_measures += dark_count[entry][eb]
        if n_measures==n_measures_per_entry:
            Hist_SystErrDist_MDM[eb].Fill(1.-bkgd_count_n_measures/data_count_n_measures)
            Hist_SystErrDist_Init[eb].Fill(1.-dark_count_n_measures/data_count_n_measures)
            MDM_mean += 1.-bkgd_count_n_measures/data_count_n_measures
            Init_mean += 1.-dark_count_n_measures/data_count_n_measures
            MDM_rms += pow(1.-bkgd_count_n_measures/data_count_n_measures,2)
            Init_rms += pow(1.-dark_count_n_measures/data_count_n_measures,2)
            n_entries += 1.
            n_measures = 0.
            data_count_n_measures = 0.
            bkgd_count_n_measures = 0.
            dark_count_n_measures = 0.
    MDM_mean = MDM_mean/n_entries
    Init_mean = Init_mean/n_entries
    MDM_rms = pow(MDM_rms/n_entries,0.5)
    Init_rms = pow(Init_rms/n_entries,0.5)

    Hist_Imposter_SystErrDist_MDM[eb].Reset()
    Hist_Imposter_SystErrDist_Init[eb].Reset()
    MDM_mean_imposter = 0.
    Init_mean_imposter = 0.
    MDM_rms_imposter = 0.
    Init_rms_imposter = 0.
    n_entries_imposter = 0.
    n_measures = 0.
    data_count_n_measures = 0.
    bkgd_count_n_measures = 0.
    dark_count_n_measures = 0.
    for entry in range(0,n_imposter_measurements):
        if eb==0 or eb==1:
            if plot_rank!=0 and imposter_data_rank[entry][eb]!=plot_rank: continue
        if imposter_data_count[entry][eb]==0.: continue
        if eb==0 and imposter_data_count[entry][eb]<2000.: continue
        if eb==1 and imposter_data_count[entry][eb]<500.: continue
        if data_count[entry][eb]<10.: continue
        n_measures += 1.
        data_count_n_measures += imposter_data_count[entry][eb]
        bkgd_count_n_measures += imposter_bkgd_count[entry][eb]
        dark_count_n_measures += imposter_dark_count[entry][eb]
        if n_measures==n_measures_per_entry:
            Hist_Imposter_SystErrDist_MDM[eb].Fill(1.-bkgd_count_n_measures/data_count_n_measures,1./5.)
            Hist_Imposter_SystErrDist_Init[eb].Fill(1.-dark_count_n_measures/data_count_n_measures,1./5.)
            MDM_mean_imposter += 1.-bkgd_count_n_measures/data_count_n_measures
            Init_mean_imposter += 1.-dark_count_n_measures/data_count_n_measures
            MDM_rms_imposter += pow(1.-bkgd_count_n_measures/data_count_n_measures,2)
            Init_rms_imposter += pow(1.-dark_count_n_measures/data_count_n_measures,2)
            n_entries_imposter += 1.
            n_measures = 0.
            data_count_n_measures = 0.
            bkgd_count_n_measures = 0.
            dark_count_n_measures = 0.
    if n_entries_imposter>0.:
        MDM_mean_imposter = MDM_mean_imposter/n_entries_imposter
        Init_mean_imposter = Init_mean_imposter/n_entries_imposter
        MDM_rms_imposter = pow(MDM_rms_imposter/n_entries_imposter,0.5)
        Init_rms_imposter = pow(Init_rms_imposter/n_entries_imposter,0.5)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_SystErrDist_MDM[eb]]
    legends += ['%0.2f-%0.2f TeV, ON data, RMS = $%0.3f \pm %0.3f$'%(energy_bin[eb]/1000.,energy_bin[eb+1]/1000.,MDM_rms,MDM_rms/pow(n_entries,0.5))]
    colors += [1]
    if analysis_type=='inclusive':
        Hists += [Hist_SystErrDist_Init[eb]]
        legends += ['initial matching, RMS = $%0.3f \pm %0.3f$'%(Init_rms,Init_rms/pow(n_entries,0.5))]
        colors += [2]
    else:
        Hists += [Hist_Imposter_SystErrDist_MDM[eb]]
        legends += ['mimic data (entries scaled by 1/5), RMS = $%0.3f \pm %0.3f$'%(MDM_rms_imposter,MDM_rms_imposter/pow(n_entries_imposter,0.5))]
        colors += [2]
    fig.clf()
    ax = fig.add_subplot()
    MakeMultipleFitPlot(ax,Hists,legends,colors,'relative error $\epsilon$','number of entries')
    fig.savefig("output_plots/SystErrDist_E%s_M%s_%s_%s_%s.png"%(eb,n_measures_per_entry,analysis_type,observing_condition,folder_path))

list_var_pair = []
good_var_pair = []
good_eigenvalue = []
for eb in range(1,2):
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
