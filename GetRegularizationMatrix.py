
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

target_energy_index = 1
#N_bins_for_deconv = 16
N_bins_for_deconv = 8
#N_bins_for_deconv = 4

folder_path = 'output_loose'

gamma_hadron_dim_ratio_w = 1.
gamma_hadron_dim_ratio_l = 1.
MSCW_blind_cut = 0.6
MSCL_blind_cut = 0.6
MSCW_plot_lower = -0.6
MSCL_plot_lower = -0.6
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
ErecS_lower_cut = 0
ErecS_upper_cut = 0
total_exposure_hours = 0.

method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

ONOFF_tag = 'OFF'
ONOFF_tag += '_Model0'
sample_list = []
sample_name = []
sample_list += ['UrsaMajorIIV6_OFF']
sample_name += ['UrsaMajorII V6']
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
sample_list += ['NGC1275V6_OFF']
sample_name += ['NGC 1275 V6']
    
elev_bins = [45,55,65,75,85]
#elev_bins = [65,75,85]
#elev_bins = [45,55,65]
theta2_bins = [0,4]

#energy_bin_ref = 1
stable_rank = 3
#stable_rank = 2

energy_bin = []
energy_bin += [100]
energy_bin += [251]
energy_bin += [631]
energy_bin += [1585]
energy_bin += [3981]
energy_bin += [10000]
energy_bin += [25118]

root_file_tags = []
mjd_tag = []
mjd_tag += ['']
for elev in range(0,len(elev_bins)-1):
    elev_tag = '_TelElev%sto%s'%(elev_bins[elev],elev_bins[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

def MakeOneHistPlot(Hist,title_x,title_y,name,logy):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad3.cd()
    pad1.cd()
    if logy: pad1.SetLogy()

    Hist.GetXaxis().SetTitleOffset(0.8)
    Hist.GetXaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetLabelSize(0.06)
    Hist.GetYaxis().SetTitleOffset(1.2)
    Hist.GetYaxis().SetTitleSize(0.06)
    Hist.GetXaxis().SetTitle(title_x)
    Hist.GetYaxis().SetTitle(title_y)
    Hist.Draw("E")

    pad3.cd()

    c_both.SaveAs('output_plots/%s.png'%(name))

def PrincipalComponentAnalysis(list_var, output_type):

    n_variables = len(list_var)
    n_samples = len(mtx_CDE_all_sources)
    mtx_var = np.zeros((n_samples,n_variables))
    mtx_var_bkg = np.zeros((n_samples,n_variables))
    chi2 = np.zeros(n_variables)
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            row = list_var[var][0]-1
            col = list_var[var][1]-1
            idx = row*3+col
            mtx_var[sample][var] = mtx_CDE_all_sources[sample][idx]
            mtx_var_bkg[sample][var] = mtx_CDE_bkgd_all_sources[sample][idx]
            chi2[var] += pow(mtx_var[sample][var]-mtx_var_bkg[sample][var],2)

    #mtx_var_centered = mtx_var - np.mean(mtx_var , axis = 0)
    mtx_var_square = np.square(mtx_var)
    mtx_var_mean = np.mean(mtx_var_square , axis = 0)
    mtx_var_rms = np.sqrt(mtx_var_mean)
    for var in range(0,n_variables):
        chi2[var] = chi2[var]/mtx_var_rms[var]
    if output_type==1: return mtx_var_rms

    mtx_var_norm = np.zeros((n_samples,n_variables))
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            if mtx_var_rms[var]==0.:
                mtx_var_norm[sample][var] = 0.
            else:
                mtx_var_norm[sample][var] = mtx_var[sample][var]/mtx_var_rms[var]

    mtx_cov = np.cov(mtx_var_norm, rowvar = False)
    eigen_values , eigen_vectors = np.linalg.eigh(mtx_cov)
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
    return eigen_values[n_variables-1]/eigen_values[n_variables-2], chi2

def MakeCorrelationPlot(list_var):

    par1_row = list_var[0][0]
    par1_col = list_var[0][1]
    par2_row = list_var[1][0]
    par2_col = list_var[1][1]

    n_variables = len(list_var)
    n_samples = len(mtx_CDE_all_sources)
    mtx_var = np.zeros((n_samples,n_variables))
    mtx_var_bkgd = np.zeros((n_samples,n_variables))
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            row = list_var[var][0]-1
            col = list_var[var][1]-1
            idx = row*3+col
            mtx_var[sample][var] = mtx_CDE_all_sources[sample][idx]
            mtx_var_bkgd[sample][var] = mtx_CDE_bkgd_all_sources[sample][idx]

    #mtx_var_centered = mtx_var - np.mean(mtx_var , axis = 0)
    mtx_var_square = np.square(mtx_var)
    mtx_var_mean = np.mean(mtx_var_square , axis = 0)
    mtx_var_rms = np.sqrt(mtx_var_mean)

    mtx_var_norm = np.zeros((n_samples,n_variables))
    mtx_var_bkgd_norm = np.zeros((n_samples,n_variables))
    for sample in range(0,n_samples):
        for var in range(0,n_variables):
            mtx_var_norm[sample][var] = mtx_var[sample][var]/mtx_var_rms[var]
            mtx_var_bkgd_norm[sample][var] = mtx_var_bkgd[sample][var]/mtx_var_rms[var]


    plt.clf()
    x_var = mtx_var_norm.transpose()[0]
    y_var = mtx_var_norm.transpose()[1]
    plt.xlabel("$t_{%s,%s}$ (arbitrary unit)"%(par1_row,par1_col), fontsize=16)
    plt.ylabel("$t_{%s,%s}$ (arbitrary unit)"%(par2_row,par2_col), fontsize=16)
    plt.scatter(x_var,y_var)

    x_var = mtx_var_bkgd_norm.transpose()[0]
    y_var = mtx_var_bkgd_norm.transpose()[1]
    plt.scatter(x_var,y_var,color='r')

    plt.savefig("output_plots/parameter_correlation_%s%s_%s%s_data.png"%(par1_row,par1_col,par2_row,par2_col))


def GetHistogramsFromFile(FilePath):
    global total_exposure_hours
    global mtx_CDE_all_sources
    global mtx_CDE_bkgd_all_sources
    global sigma_rank0_all_sources
    global sigma_rank1_all_sources
    global sigma_rank2_all_sources
    dark_stable_rank = ROOT.std.vector("int")(10)
    dark_sigma_rank0 = ROOT.std.vector("double")(10)
    dark_sigma_rank1 = ROOT.std.vector("double")(10)
    dark_sigma_rank2 = ROOT.std.vector("double")(10)
    dark_gamma_count = ROOT.std.vector("double")(10)
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    NewInfoTree = InputFile.Get("NewInfoTree")
    NewInfoTree.SetBranchAddress('dark_stable_rank',ROOT.AddressOf(dark_stable_rank))
    NewInfoTree.SetBranchAddress('dark_sigma_rank0',ROOT.AddressOf(dark_sigma_rank0))
    NewInfoTree.SetBranchAddress('dark_sigma_rank1',ROOT.AddressOf(dark_sigma_rank1))
    NewInfoTree.SetBranchAddress('dark_sigma_rank2',ROOT.AddressOf(dark_sigma_rank2))
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    if energy_index==0: 
        total_exposure_hours += exposure_hours
    Hist2D_Coeff_Data.Reset()
    Hist2D_Coeff_Bkgd.Reset()
    Hist2D_MSCLW_Data.Reset()
    Hist2D_MSCLW_Best.Reset()
    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Get %s'%(HistName))
    Hist2D_Coeff_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Coeff_Bkgd_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Get %s'%(HistName))
    Hist2D_Coeff_Bkgd.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Get %s'%(HistName))
    Hist2D_MSCLW_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_OnBest_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print ('Get %s'%(HistName))
    Hist2D_MSCLW_Best.Add(InputFile.Get(HistName))
    if Hist2D_MSCLW_Data.Integral()>0.:
        Hist2D_MSCLW_Best.Add(Hist2D_MSCLW_Data,-1.)
        Hist2D_MSCLW_Best.Divide(Hist2D_MSCLW_Data)
        for binx in range(0,Hist2D_MSCLW_Best.GetNbinsX()):
            for biny in range(0,Hist2D_MSCLW_Best.GetNbinsY()):
                old_content = Hist2D_MSCLW_Best.GetBinContent(binx+1,biny+1)
                Hist2D_MSCLW_Best.SetBinContent(binx+1,biny+1,pow(old_content,2))
        Hist2D_MSCLW_Best_Sum[energy_index].Add(Hist2D_MSCLW_Best,1./len(sample_list))
    for binx in range(0,Hist2D_Coeff_Data.GetNbinsX()):
        for biny in range(0,Hist2D_Coeff_Data.GetNbinsY()):
            old_content = Hist2D_Regularization[energy_index].GetBinContent(binx+1,biny+1)
            new_content = pow(Hist2D_Coeff_Data.GetBinContent(binx+1,biny+1),2)
            Hist2D_Regularization[energy_index].SetBinContent(binx+1,biny+1,old_content+new_content)
    #if dark_gamma_count[energy_index]>100. and dark_stable_rank[energy_index]==stable_rank and energy_index==energy_bin_ref: 
    if dark_gamma_count[energy_index]>1000. and dark_stable_rank[energy_index]==stable_rank: 
        mtx_CDE = []
        mtx_CDE_bkgd = []
        for row in range(0,3):
            for col in range(0,3):
                mtx_CDE += [0.]
                mtx_CDE_bkgd += [0.]
        for row in range(0,3):
            for col in range(0,3):
                idx = row*3+col
                content = Hist2D_Coeff_Data.GetBinContent(row+1,col+1)
                Hist_CDE[idx].Fill(content)
                mtx_CDE[idx] = content
                content = Hist2D_Coeff_Bkgd.GetBinContent(row+1,col+1)
                mtx_CDE_bkgd[idx] = content
        mtx_CDE_all_sources += [mtx_CDE]
        mtx_CDE_bkgd_all_sources += [mtx_CDE_bkgd]
    if target_energy_index==energy_index:
        sigma_rank0_all_sources += [dark_sigma_rank0[energy_index]]
        sigma_rank1_all_sources += [dark_sigma_rank1[energy_index]]
        sigma_rank2_all_sources += [dark_sigma_rank2[energy_index]]

Hist2D_MSCLW_Data = ROOT.TH2D("Hist2D_MSCLW_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_MSCLW_Best = ROOT.TH2D("Hist2D_MSCLW_Best","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Coeff_Data = ROOT.TH2D("Hist2D_Coeff_Data","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Bkgd = ROOT.TH2D("Hist2D_Coeff_Bkgd","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Regularization = []
Hist2D_MSCLW_Best_Sum = []
for e in range(0,len(energy_bin)-1):
    ErecS_lower_cut = energy_bin[e]
    ErecS_upper_cut = energy_bin[e+1]
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    Hist2D_Regularization += [ROOT.TH2D("Hist2D_Regularization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int),"",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)]
    Hist2D_MSCLW_Best_Sum += [ROOT.TH2D("Hist2D_MSCLW_Best_Sum_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]

Hist_CDE = []
mtx_CDE_all_sources = []
mtx_CDE_bkgd_all_sources = []
sigma_rank0_all_sources = []
sigma_rank1_all_sources = []
sigma_rank2_all_sources = []
for row in range(0,3):
    for col in range(0,3):
        idx = row*3+col
        Hist_CDE += [ROOT.TH1D("Hist_CDE_%s"%(idx),"",20,-0.1,0.1)]

FilePath_List = []
for source in range(0,len(sample_list)):
    source_name = sample_list[source]
    for elev in range(0,len(root_file_tags)):
        FilePath = "%s/Netflix_"%(folder_path)+sample_list[source]+"_%s"%(root_file_tags[elev])+".root"
        FilePath_List += [FilePath]
        print ('Reading file %s'%(FilePath_List[len(FilePath_List)-1]))
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):
            print ('cannot find file %s'%(FilePath_List[len(FilePath_List)-1]))
            continue
        for e in range(0,len(energy_bin)-1):
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1])


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
    lumilab1 = ROOT.TLatex(0.15,0.50,'#tau_{kn} coefficents' )
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
    canvas.SaveAs('output_plots/Coeff_t2_data_%s_%s.png'%(e,method_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'best 9-par model residual' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    pad1.SetLogz()
    Hist2D_MSCLW_Best_Sum[e].GetYaxis().SetTitle('MSCL')
    Hist2D_MSCLW_Best_Sum[e].GetXaxis().SetTitle('MSCW')
    Hist2D_MSCLW_Best_Sum[e].Draw("COL4Z")
    canvas.SaveAs('output_plots/Best_9Par_Residual_%s_%s.png'%(e,method_tag))

for row in range(0,3):
    for col in range(0,3):
        idx = row*3+col
        MakeOneHistPlot(Hist_CDE[idx],'perturbation','number of runs','CDE_%s_%s'%(row,col),False)

#my_table = PrettyTable()
#my_table.field_names = ["source","t11","t12","t13","t21","t22","t23","t31","t32","t33","sigma1","sigma2","sigma3"]
#my_table.float_format["t11"] = ".2e"
#my_table.float_format["t12"] = ".2e"
#my_table.float_format["t13"] = ".2e"
#my_table.float_format["t21"] = ".2e"
#my_table.float_format["t22"] = ".2e"
#my_table.float_format["t23"] = ".2e"
#my_table.float_format["t31"] = ".2e"
#my_table.float_format["t32"] = ".2e"
#my_table.float_format["t33"] = ".2e"
#my_table.float_format["sigma1"] = ".2e"
#my_table.float_format["sigma2"] = ".2e"
#my_table.float_format["sigma3"] = ".2e"
#for entry in range(0,len(sample_list)):
#    table_row = []
#    source_name = sample_name[entry]
#    table_row += [source_name]
#    for row in range(0,3):
#        for col in range(0,3):
#            idx = row*3+col
#            table_row += [mtx_CDE_all_sources[entry][idx]]
#    table_row += [sigma_rank0_all_sources[entry]]
#    table_row += [sigma_rank1_all_sources[entry]]
#    table_row += [sigma_rank2_all_sources[entry]]
#    my_table.add_row(table_row)
#print(my_table)


list_var_pair = []
good_var_pair = []
good_eigenvalue = []
for row1 in range(0,stable_rank):
    for col1 in range(0,stable_rank):
        for row2 in range(0,stable_rank):
            for col2 in range(0,stable_rank):
                idx1 = row1*stable_rank+col1
                idx2 = row2*stable_rank+col2
                list_var_pair = [[row1+1,col1+1]]
                list_var_pair += [[row2+1,col2+1]]
                if idx1<idx2:
                #if row1==1 and col1==1:
                    print('=======================================================')
                    max_eigenvalue, chi2 = PrincipalComponentAnalysis(list_var_pair,0)
                    if math.isnan(max_eigenvalue): continue
                    chi2_ratio = min(chi2[0]/chi2[1],chi2[1]/chi2[0])
                    #if chi2_ratio>0.8: continue
                    #if max_eigenvalue<5.0: continue
                    #if max_eigenvalue<2.0: continue
                    MakeCorrelationPlot(list_var_pair)
                    good_var_pair += [list_var_pair]
                    good_eigenvalue += [max_eigenvalue]
for entry in range(0,len(good_eigenvalue)):
    print('{0}, {1}'.format(good_var_pair[entry],good_eigenvalue[entry]))

list_var = []
list_rms = []
for row in range(0,stable_rank):
    for col in range(0,stable_rank):
        idx = row*3+col
        tmp_list_var = [[row+1,col+1]]
        tmp_list_rms = PrincipalComponentAnalysis(tmp_list_var,1)
        if math.isnan(tmp_list_rms[0]): continue
        list_var += [[row+1,col+1]]
        list_rms += [tmp_list_rms[0]]
for entry in range(0,len(list_var)):
    print('RMS: \n {0}, {1}'.format(list_var[entry],list_rms[entry]))

fig, ax = plt.subplots()

plt.clf()
ax = fig.add_subplot(111)
ind = np.arange(len(list_rms))
width = 0.35
rects1 = ax.bar(ind, list_rms, width, color='#089FFF')
ax.set_xlim(-width,len(ind)+width)
ax.set_ylabel("RMS", fontsize=18)
xTickMarks = list_var
ax.set_xticks(ind+0.5*width)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=45, fontsize=10)
plt.subplots_adjust(bottom=0.15)
plt.savefig("output_plots/VariableRMS.png")

plt.clf()
ax = fig.add_subplot(111)
ind = np.arange(len(good_eigenvalue))
width = 0.35
rects1 = ax.bar(ind, good_eigenvalue, width, color='#089FFF')
ax.set_xlim(-width,len(ind)+width)
#ax.set_ylim(1.5,2.2)
ax.set_ylabel("eigenvalue", fontsize=18)
xTickMarks = good_var_pair
ax.set_xticks(ind+0.5*width)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=45, fontsize=10)
plt.subplots_adjust(bottom=0.15)
plt.savefig("output_plots/VarPairEigenvalue.png")
