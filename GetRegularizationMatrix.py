
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

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

N_bins_for_deconv = 16
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

folder_path = 'output_root'
#folder_path = 'output_test'
#method_tag = 'loose_mdm_default'
method_tag = 'tight_mdm_default'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

ONOFF_tag = 'OFF'
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
sample_list += ['1ES0229V5_OFF']
sample_name += ['1ES0229 V5']
sample_list += ['Segue1V6_OFF']
sample_name += ['Segue1 V6']
sample_list += ['Segue1V5_OFF']
sample_name += ['Segue1 V5']
sample_list += ['CrabV5_OFF']
sample_name += ['Crab V5']
sample_list += ['CrabV6_OFF']
sample_name += ['Crab V6']
sample_list += ['BLLacV6_OFF']
sample_name += ['BLLac V6']
sample_list += ['BLLacV5_OFF']
sample_name += ['BLLac V5']
sample_list += ['PG1553V5_OFF']
sample_name += ['PG1553 V5']
sample_list += ['H1426V6_OFF']
sample_name += ['H1426 V6']
sample_list += ['CasAV6_OFF']
sample_name += ['CasA V6']
sample_list += ['RBS0413V6_OFF']
sample_name += ['RBS0413 V6']

    
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

def GetHistogramsFromFile(FilePath):
    global total_exposure_hours
    global mtx_CDE_all_sources
    dark_gamma_count = ROOT.std.vector("double")(10)
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    NewInfoTree = InputFile.Get("NewInfoTree")
    NewInfoTree.SetBranchAddress('dark_gamma_count',ROOT.AddressOf(dark_gamma_count))
    NewInfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours
    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)
    energy_index = energy_bin.index(ErecS_lower_cut)
    if energy_index==0: 
        total_exposure_hours += exposure_hours
    Hist2D_Coeff_Data.Reset()
    Hist2D_MSCLW_Data.Reset()
    Hist2D_MSCLW_Best.Reset()
    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Get %s'%(HistName)
    Hist2D_Coeff_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Get %s'%(HistName)
    Hist2D_MSCLW_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_OnBest_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Get %s'%(HistName)
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
    if energy_index==1: 
        mtx_CDE = []
        for row in range(0,3):
            for col in range(0,3):
                mtx_CDE += [0.]
        for row in range(0,3):
            for col in range(0,3):
                idx = row*3+col
                content = Hist2D_Coeff_Data.GetBinContent(row+1,col+1)
                Hist_CDE[idx].Fill(content)
                mtx_CDE[idx] = content
        mtx_CDE_all_sources += [mtx_CDE]

Hist2D_MSCLW_Data = ROOT.TH2D("Hist2D_MSCLW_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_MSCLW_Best = ROOT.TH2D("Hist2D_MSCLW_Best","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Coeff_Data = ROOT.TH2D("Hist2D_Coeff_Data","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
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
        print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
        for e in range(0,len(energy_bin)-1):
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            GetHistogramsFromFile(FilePath_List[len(FilePath_List)-1])

OutputFilePath = "%s/Regularization%s.root"%(folder_path,lowrank_tag)
OutputFile = ROOT.TFile(OutputFilePath,"recreate")
print "total_exposure_hours = %s"%(total_exposure_hours)
for e in range(0,len(energy_bin)-1):
    Hist2D_Regularization[e].Write()
    Hist2D_MSCLW_Best_Sum[e].Write()
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

my_table = PrettyTable()
my_table.field_names = ["source","t11","D12","D13","C21","t22","D23","C31","C32","t33"]
my_table.float_format["t11"] = ".3f"
my_table.float_format["D12"] = ".3f"
my_table.float_format["D13"] = ".3f"
my_table.float_format["C21"] = ".3f"
my_table.float_format["t22"] = ".3f"
my_table.float_format["D23"] = ".3f"
my_table.float_format["C31"] = ".3f"
my_table.float_format["C32"] = ".3f"
my_table.float_format["t33"] = ".3f"
for entry in range(0,len(sample_list)):
    table_row = []
    source_name = sample_name[entry]
    table_row += [source_name]
    for row in range(0,3):
        for col in range(0,3):
            idx = row*3+col
            table_row += [mtx_CDE_all_sources[entry][idx]]
    my_table.add_row(table_row)
print(my_table)

