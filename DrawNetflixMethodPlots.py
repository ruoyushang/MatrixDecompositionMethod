
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from scipy import special

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder_path = 'output_root'
PercentCrab = '_Crab0'

energy_fine_bin_cut_low = 3
energy_fine_bin_cut_up = 20

N_bins_for_deconv = 12
gamma_hadron_dim_ratio = 1.
MSCW_blind_cut = 0.5
MSCL_blind_cut = 0.5
MSCW_lower_cut = -1.
MSCL_lower_cut = -1.
MSCW_plot_lower = -1.
MSCL_plot_lower = -1.
MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
ErecS_lower_cut = 0
ErecS_upper_cut = 0

n_bad_matches = 0
exposure_hours = 0.
source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
n_control_samples = 5

Syst_MDM = 0.02

elev_range = []
#elev_range += [[25,45]]
elev_range += [[45,85]]

energy_list = []
energy_list += [int(pow(10,2.3))]
energy_list += [int(pow(10,4.0))]

energy_fine_bin = []
energy_fine_bin += [pow(10,2.0)]
energy_fine_bin += [pow(10,2.1)]
energy_fine_bin += [pow(10,2.2)]
energy_fine_bin += [pow(10,2.3)]
energy_fine_bin += [pow(10,2.4)]
energy_fine_bin += [pow(10,2.5)]
energy_fine_bin += [pow(10,2.6)]
energy_fine_bin += [pow(10,2.7)]
energy_fine_bin += [pow(10,2.8)]
energy_fine_bin += [pow(10,2.9)]
energy_fine_bin += [pow(10,3.0)]
energy_fine_bin += [pow(10,3.1)]
energy_fine_bin += [pow(10,3.2)]
energy_fine_bin += [pow(10,3.3)]
energy_fine_bin += [pow(10,3.4)]
energy_fine_bin += [pow(10,3.5)]
energy_fine_bin += [pow(10,3.6)]
energy_fine_bin += [pow(10,3.7)]
energy_fine_bin += [pow(10,3.8)]
energy_fine_bin += [pow(10,3.9)]
energy_fine_bin += [pow(10,4.0)]

def ResetStackedShowerHistograms():

    Hist2D_OnData_Sum.Reset()
    Hist2D_OnBkgd_Sum.Reset()
    Hist2D_OnDark_Sum.Reset()
    Hist_OnData_MSCL_Sum.Reset()
    Hist_OnBkgd_MSCL_Sum.Reset()
    Hist_OnDark_MSCL_Sum.Reset()
    Hist_OnData_MSCW_Sum.Reset()
    Hist_OnBkgd_MSCW_Sum.Reset()
    Hist_OnDark_MSCW_Sum.Reset()

    for nth_sample in range(0,n_control_samples-1):

        Hist2D_OffData_Sum[nth_sample].Reset()
        Hist2D_OffBkgd_Sum[nth_sample].Reset()
        Hist2D_OffDark_Sum[nth_sample].Reset()
        Hist_OffData_MSCL_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCL_Sum[nth_sample].Reset()
        Hist_OffDark_MSCL_Sum[nth_sample].Reset()
        Hist_OffData_MSCW_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCW_Sum[nth_sample].Reset()
        Hist_OffDark_MSCW_Sum[nth_sample].Reset()

def GetSourceInfo(file_list):

    global MSCW_blind_cut
    global MSCL_blind_cut
    global n_bad_matches
    global exposure_hours
    global source_ra
    global source_dec
    global source_l
    global source_b

    n_bad_matches = 0
    exposure_hours = 0.
    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        n_bad_matches += InfoTree.n_bad_matches
        exposure_hours += InfoTree.exposure_hours
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        source_ra = InfoTree.mean_tele_point_ra
        source_dec = InfoTree.mean_tele_point_dec
        source_l = InfoTree.mean_tele_point_l
        source_b = InfoTree.mean_tele_point_b
        HistName = "Hist_EffArea"
        Hist_EffArea.Reset()
        Hist_EffArea.Add(InputFile.Get(HistName))
        Hist_EffArea.Scale(exposure_hours*3600.)
        Hist_EffArea_Sum.Add(Hist_EffArea)
        InputFile.Close()

def GetShowerHistogramsFromFile(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut

    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnData.Reset()
    Hist2D_OnData.Add(InputFile.Get(HistName))

    HistName = "Hist_OnDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnDark.Reset()
    Hist2D_OnDark.Add(InputFile.Get(HistName))

    HistName = "Hist_OnBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnBkgd.Reset()
    Hist2D_OnBkgd.Add(InputFile.Get(HistName))

    Hist_OnData_MSCL.Reset()
    Hist_OnData_MSCL.Add(Hist2D_OnData.ProjectionX("Hist1D_OnData_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnData_MSCW.Reset()
    Hist_OnData_MSCW.Add(Hist2D_OnData.ProjectionY("Hist1D_OnData_MSCW",bin_lower_x,bin_upper_x))

    Hist_OnDark_MSCL.Reset()
    Hist_OnDark_MSCL.Add(Hist2D_OnDark.ProjectionX("Hist1D_OnDark_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnDark_MSCW.Reset()
    Hist_OnDark_MSCW.Add(Hist2D_OnDark.ProjectionY("Hist1D_OnDark_MSCW",bin_lower_x,bin_upper_x))

    Hist_OnBkgd_MSCL.Reset()
    Hist_OnBkgd_MSCL.Add(Hist2D_OnBkgd.ProjectionX("Hist1D_OnBkgd_MSCL",bin_lower_y,bin_upper_y))
    Hist_OnBkgd_MSCW.Reset()
    Hist_OnBkgd_MSCW.Add(Hist2D_OnBkgd.ProjectionY("Hist1D_OnBkgd_MSCW",bin_lower_x,bin_upper_x))

    for nth_sample in range(0,n_control_samples-1):

        HistName = "Hist_OffData_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffData[nth_sample].Reset()
        Hist2D_OffData[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OffDark_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffDark[nth_sample].Reset()
        Hist2D_OffDark[nth_sample].Add(InputFile.Get(HistName))

        HistName = "Hist_OffBkgd_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffBkgd[nth_sample].Reset()
        Hist2D_OffBkgd[nth_sample].Add(InputFile.Get(HistName))

        Hist_OffData_MSCL[nth_sample].Reset()
        Hist_OffData_MSCL[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionX("Hist1D_OffData_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffDark_MSCL[nth_sample].Reset()
        Hist_OffDark_MSCL[nth_sample].Add(Hist2D_OffDark[nth_sample].ProjectionX("Hist1D_OffDark_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffBkgd_MSCL[nth_sample].Reset()
        Hist_OffBkgd_MSCL[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionX("Hist1D_OffBkgd_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffData_MSCW[nth_sample].Reset()
        Hist_OffData_MSCW[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionY("Hist1D_OffData_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

        Hist_OffDark_MSCW[nth_sample].Reset()
        Hist_OffDark_MSCW[nth_sample].Add(Hist2D_OffDark[nth_sample].ProjectionY("Hist1D_OffDark_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

        Hist_OffBkgd_MSCW[nth_sample].Reset()
        Hist_OffBkgd_MSCW[nth_sample].Add(Hist2D_OffBkgd[nth_sample].ProjectionY("Hist1D_OffBkgd_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))

def StackShowerHistograms():

    Hist2D_OnData_Sum.Add(Hist2D_OnData)
    Hist2D_OnDark_Sum.Add(Hist2D_OnDark)
    Hist2D_OnBkgd_Sum.Add(Hist2D_OnBkgd)

    Hist_OnData_MSCL_Sum.Add(Hist_OnData_MSCL)
    Hist_OnData_MSCW_Sum.Add(Hist_OnData_MSCW)
    Hist_OnDark_MSCL_Sum.Add(Hist_OnDark_MSCL)
    Hist_OnDark_MSCW_Sum.Add(Hist_OnDark_MSCW)
    Hist_OnBkgd_MSCL_Sum.Add(Hist_OnBkgd_MSCL)
    Hist_OnBkgd_MSCW_Sum.Add(Hist_OnBkgd_MSCW)

    for nth_sample in range(0,n_control_samples-1):

        Hist2D_OffData_Sum[nth_sample].Add(Hist2D_OffData[nth_sample])
        Hist2D_OffDark_Sum[nth_sample].Add(Hist2D_OffDark[nth_sample])
        Hist2D_OffBkgd_Sum[nth_sample].Add(Hist2D_OffBkgd[nth_sample])

        Hist_OffData_MSCL_Sum[nth_sample].Add(Hist_OffData_MSCL[nth_sample])
        Hist_OffData_MSCW_Sum[nth_sample].Add(Hist_OffData_MSCW[nth_sample])
        Hist_OffDark_MSCL_Sum[nth_sample].Add(Hist_OffDark_MSCL[nth_sample])
        Hist_OffDark_MSCW_Sum[nth_sample].Add(Hist_OffDark_MSCW[nth_sample])
        Hist_OffBkgd_MSCL_Sum[nth_sample].Add(Hist_OffBkgd_MSCL[nth_sample])
        Hist_OffBkgd_MSCW_Sum[nth_sample].Add(Hist_OffBkgd_MSCW[nth_sample])

def CalculateSignificance(s,b,err):
    if (b*b+(s+b)*err*err)==0.: return 0.
    if (s+b)*(b+err*err)==0.: return 0.
    if ((s+b)*(b+err*err)/(b*b+(s+b)*err*err))<=0.: return 0.
    first_term = (s+b)*math.log((s+b)*(b+err*err)/(b*b+(s+b)*err*err))
    if err>0. and b>0:
        second_term = b*b/(err*err)*math.log(1.+err*err*s/(b*(b+err*err)))
    else: 
        second_term = 0.
    result = 0.
    if first_term>second_term: result = pow(2*(first_term-second_term),0.5)
    else: result = pow(2*(-first_term+second_term),0.5)
    if s>0: return result
    else: return -1.*result

def IntegralAndError(Hist,bin1,bin2):
    
    integral = 0
    error = 0
    for b in range(1,Hist.GetNbinsX()+1):
        if b<bin1: continue
        if b>bin2: continue
        integral += Hist.GetBinContent(b)
        error += pow(Hist.GetBinError(b),2)
    error = pow(error,0.5)
    error = max(error,1.)
    if math.isnan(integral) or math.isnan(error):
        integral = 0
        error = 1.
    return integral, error

def MakeChi2Plot(Hists,legends,colors,title,name,doSum,range_lower,range_upper,syst):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,0.8)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
    pad2.SetBottomMargin(0.39)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad1.SetGrid()
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    norm_bin_low_target = Hists[0].FindBin(range_lower)
    norm_bin_up_target = Hists[0].FindBin(range_upper)-1

    max_heigh = 0
    max_hist = 0
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
    Hists[max_hist].Draw("E")

    fill_color = [0,0,0,30,38,46]
    if doSum:
        stack = ROOT.THStack("stack", "")
        Hist_Sum = Hists[1].Clone()
        Hist_Sum.Reset()
        Hist_Sys = Hists[1].Clone()
        Hist_Sys.Reset()
        for h in range(0,len(Hists)):
            if colors[h]==1: continue
            if colors[h]==2: continue
            if not legends[h]=='syst.':
                set_histStyle( Hists[h] , fill_color[colors[h]])
                stack.Add( Hists[h] )
                Hist_Sum.Add( Hists[h] )
            else:
                Hist_Sys.Add( Hists[h] )
        stack.Draw("hist same")
        Hist_Err = Hist_Sum.Clone()
        for binx in range(0,Hist_Err.GetNbinsX()):
            old_err = Hist_Err.GetBinError(binx+1)
            new_err = Hist_Sys.GetBinContent(binx+1)
            Hist_Err.SetBinError(binx+1,pow(old_err*old_err+new_err*new_err,0.5))
        Hist_Err.SetFillColor(1)
        Hist_Err.SetFillStyle(3004)
        Hist_Err.SetMarkerSize(0)
        Hist_Err.Draw("e2 same")

    Hist_Bkgd = Hists[1].Clone()
    Hist_Bkgd.Reset()
    Hist_Syst = Hists[1].Clone()
    Hist_Syst.Reset()
    for h in range(0,len(Hists)):
        if colors[h]==0 or colors[h]==4:
            Hist_Syst.Add(Hists[h].Clone())
        if colors[h]==4:
            Hist_Bkgd.Add(Hists[h].Clone())

    for h in range(0,len(Hists)):
        if colors[h]==0: continue
        Hists[h].SetLineWidth(3)
        Hists[h].Draw("E same")
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    #legend.SetTextSize(0.2)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hists[0],legends[0],"pl")
    for h in range(1,len(Hists)):
        if doSum:
            legend.AddEntry(Hists[h],legends[h],"f")
        else:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(ErecS_lower_cut,exposure_hours) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()

    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hists[0],norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkgd,norm_bin_low_target,norm_bin_up_target)
    syst_err = Syst_MDM
    if not 'MDM' in name: syst_err = 0.
    err_bkg = pow(err_bkg*err_bkg+(syst_err*predict_bkg)*(syst_err*predict_bkg),0.5)
    Sig = 1.*CalculateSignificance(data_SR-predict_bkg,predict_bkg,err_bkg)
    lumilab2 = ROOT.TLatex(0.15,0.60,'Excess = %0.1f#pm%0.1f (%0.1f#sigma)'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5),Sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.40,'Bkg = %0.1f#pm%0.1f'%(predict_bkg,err_bkg) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0: 
        sbratio = (data_SR-predict_bkg)/(predict_bkg)
    if not data_SR-predict_bkg==0 and not predict_bkg==0:
        sbratio_err = sbratio*pow(pow(pow(err_SR*err_SR+err_bkg*err_bkg,0.5)/(data_SR-predict_bkg),2)+pow(err_bkg/predict_bkg,2),0.5)
    else: sbratio_err = 0
    lumilab4 = ROOT.TLatex(0.15,0.20,'S/B = %0.3f#pm%0.3f'%(sbratio,sbratio_err) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()

    pad2.cd()
    Hist_Band = Hist_Syst.Clone()
    Hist_Band.Divide(Hist_Bkgd)
    Hist_Band.SetFillColor(1)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    Hist_Band.GetYaxis().SetTitleSize(0.10)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(1.2)
    Hist_Band.SetMinimum(0.8)
    Hist_Band.Draw("e2")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    Hist_Ratio = []
    for h in range(1,len(Hists)):
        Hist_Bkg = Hists[h].Clone()
        Hist_Ratio += [Hists[0].Clone()]
        Hist_Ratio[h-1].Divide(Hist_Bkg)
        Hist_Ratio[h-1].SetLineWidth(2)
        Hist_Ratio[h-1].SetLineColor(colors[h])
        Hist_Ratio[h-1].Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def PlotsStackedHistograms(tag):

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_OnDark_MSCW_Sum]
    legends += ['init. bkg.']
    colors += [2]
    Hists += [Hist_OnBkgd_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_SystErr_MSCW]
    legends += ['syst. error']
    colors += [0]
    plotname = 'Stack_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,title,plotname,False,MSCW_lower_cut,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_OnDark_MSCL_Sum]
    legends += ['init. bkg.']
    colors += [2]
    Hists += [Hist_OnBkgd_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_SystErr_MSCL]
    legends += ['syst. error']
    colors += [0]
    plotname = 'Stack_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,title,plotname,False,MSCL_lower_cut,MSCL_blind_cut,-1)

    #for nth_sample in range(0,n_control_samples-1):

    #    Hists = []
    #    legends = []
    #    colors = []
    #    Hists += [Hist_OffData_MSCW_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    Hists += [Hist_OffDark_MSCW_Sum[nth_sample]]
    #    legends += ['init. bkg.']
    #    colors += [2]
    #    Hists += [Hist_OffBkgd_MSCW_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    plotname = 'Stack_MSCW_MDM_Control%s_%s'%(nth_sample,tag)
    #    title = 'MSCW'
    #    MakeChi2Plot(Hists,legends,colors,title,plotname,False,MSCW_lower_cut,MSCW_blind_cut,-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    Hists += [Hist_OffData_MSCL_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    Hists += [Hist_OffDark_MSCL_Sum[nth_sample]]
    #    legends += ['init. bkg.']
    #    colors += [2]
    #    Hists += [Hist_OffBkgd_MSCL_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    plotname = 'Stack_MSCL_MDM_Control%s_%s'%(nth_sample,tag)
    #    title = 'MSCL'
    #    MakeChi2Plot(Hists,legends,colors,title,plotname,False,MSCL_lower_cut,MSCL_blind_cut,-1)

def CalculateSystError():

    global Syst_MDM
    Syst_MDM = 0.
    Hist_SystErr_MSCL.Reset()
    Hist_SystErr_MSCW.Reset()

    for nth_sample in range(0,n_control_samples-1):

        norm_bin_low_target = Hist_OffData_MSCL_Sum[nth_sample].FindBin(MSCL_lower_cut)
        norm_bin_up_target = Hist_OffData_MSCL_Sum[nth_sample].FindBin(MSCL_blind_cut)-1
        Total_Bkgd = Hist_OffData_MSCL_Sum[nth_sample].Integral(norm_bin_low_target,norm_bin_up_target)
        Hist_Diff_MSCL = Hist_OffData_MSCL_Sum[nth_sample].Clone()
        Hist_Diff_MSCL.Add(Hist_OffBkgd_MSCL_Sum[nth_sample],-1.)
        Syst_MDM += pow(Hist_Diff_MSCL.Integral(norm_bin_low_target,norm_bin_up_target)/Total_Bkgd,2)
        for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
            syst_err = max(0.,pow(Hist_Diff_MSCL.GetBinContent(binx+1),2)-pow(Hist_Diff_MSCL.GetBinError(binx+1),2))
            old_syst = Hist_SystErr_MSCL.GetBinError(binx+1)
            new_syst = old_syst+syst_err
            Hist_SystErr_MSCL.SetBinError(binx+1,new_syst)

        Hist_Diff_MSCW = Hist_OffData_MSCW_Sum[nth_sample].Clone()
        Hist_Diff_MSCW.Add(Hist_OffBkgd_MSCW_Sum[nth_sample],-1.)
        for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
            syst_err = max(0.,pow(Hist_Diff_MSCW.GetBinContent(binx+1),2)-pow(Hist_Diff_MSCW.GetBinError(binx+1),2))
            old_syst = Hist_SystErr_MSCW.GetBinError(binx+1)
            new_syst = old_syst+syst_err
            Hist_SystErr_MSCW.SetBinError(binx+1,new_syst)

    Syst_MDM = pow(Syst_MDM/n_control_samples,0.5)
    print "Syst_MDM = %s"%(Syst_MDM) 

    for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
            old_syst = Hist_SystErr_MSCL.GetBinError(binx+1)
            new_syst = pow(old_syst/n_control_samples,0.5)
            Hist_SystErr_MSCL.SetBinError(binx+1,new_syst)

    for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
            old_syst = Hist_SystErr_MSCW.GetBinError(binx+1)
            new_syst = pow(old_syst/n_control_samples,0.5)
            Hist_SystErr_MSCW.SetBinError(binx+1,new_syst)

def SingleSourceAnalysis(source_list,doMap):

    global ErecS_lower_cut
    global ErecS_upper_cut

    FilePath_List = []
    ResetStackedShowerHistograms()
    for source in range(0,len(source_list)):
        source_name = source_list[source]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+".root";
            FilePath_List += [FilePath]
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
            InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
            InfoTree = InputFile.Get("InfoTree")
            InfoTree.GetEntry(0)
            MSCW_blind_cut = InfoTree.MSCW_cut_blind
            MSCL_blind_cut = InfoTree.MSCL_cut_blind
            bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_lower_cut)
            bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
            bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_lower_cut)
            bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                if ErecS_upper_cut<=energy_fine_bin[energy_fine_bin_cut_low]: continue
                if ErecS_lower_cut>=energy_fine_bin[energy_fine_bin_cut_up]: continue
                GetShowerHistogramsFromFile(FilePath_List[len(FilePath_List)-1])
                StackShowerHistograms()
    GetSourceInfo(FilePath_List)

    CalculateSystError()
    PlotsStackedHistograms('%s%s'%(source_list[0],PercentCrab))

Hist2D_OnData_Sum = ROOT.TH2D("Hist2D_OnData_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd_Sum = ROOT.TH2D("Hist2D_OnBkgd_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark_Sum = ROOT.TH2D("Hist2D_OnDark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData = ROOT.TH2D("Hist2D_OnData","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd = ROOT.TH2D("Hist2D_OnBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark = ROOT.TH2D("Hist2D_OnDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCL_Sum = ROOT.TH1D("Hist_OnData_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL_Sum = ROOT.TH1D("Hist_OnDark_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCL = ROOT.TH1D("Hist_OnData_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL = ROOT.TH1D("Hist_OnBkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL = ROOT.TH1D("Hist_OnDark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCW_Sum = ROOT.TH1D("Hist_OnData_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnDark_MSCW_Sum = ROOT.TH1D("Hist_OnDark_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCW = ROOT.TH1D("Hist_OnData_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW = ROOT.TH1D("Hist_OnBkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnDark_MSCW = ROOT.TH1D("Hist_OnDark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist_SystErr_MSCL = ROOT.TH1D("Hist_SystErr_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_SystErr_MSCW = ROOT.TH1D("Hist_SystErr_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist_EffArea = ROOT.TH1D("Hist_EffArea","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_EffArea_Sum = ROOT.TH1D("Hist_EffArea_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

Hist2D_OffData = []
Hist2D_OffData_Sum = []
Hist2D_OffDark = []
Hist2D_OffDark_Sum = []
Hist2D_OffBkgd = []
Hist2D_OffBkgd_Sum = []
Hist_OffData_MSCL = []
Hist_OffData_MSCL_Sum = []
Hist_OffDark_MSCL = []
Hist_OffDark_MSCL_Sum = []
Hist_OffBkgd_MSCL = []
Hist_OffBkgd_MSCL_Sum = []
Hist_OffData_MSCW = []
Hist_OffData_MSCW_Sum = []
Hist_OffDark_MSCW = []
Hist_OffDark_MSCW_Sum = []
Hist_OffBkgd_MSCW = []
Hist_OffBkgd_MSCW_Sum = []
for nth_sample in range(0,n_control_samples-1):
    Hist2D_OffData += [ROOT.TH2D("Hist2D_OffData_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffData_Sum += [ROOT.TH2D("Hist2D_OffData_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffDark += [ROOT.TH2D("Hist2D_OffDark_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffDark_Sum += [ROOT.TH2D("Hist2D_OffDark_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffBkgd += [ROOT.TH2D("Hist2D_OffBkgd_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffBkgd_Sum += [ROOT.TH2D("Hist2D_OffBkgd_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCL += [ROOT.TH1D("Hist_OnData_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCL_Sum += [ROOT.TH1D("Hist_OnData_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCW += [ROOT.TH1D("Hist_OnData_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCW_Sum += [ROOT.TH1D("Hist_OnData_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffDark_MSCL += [ROOT.TH1D("Hist_OnDark_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffDark_MSCL_Sum += [ROOT.TH1D("Hist_OnDark_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffDark_MSCW += [ROOT.TH1D("Hist_OnDark_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffDark_MSCW_Sum += [ROOT.TH1D("Hist_OnDark_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCL += [ROOT.TH1D("Hist_OnBkgd_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCL_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCW += [ROOT.TH1D("Hist_OnBkgd_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCW_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]


sample_list = []
sky_coord = []
sample_list += ['H1426']
sky_coord += ['14 28 32.609 +42 40 21.05']
#sample_list += ['M82']
#sky_coord += ['09 55 52.7 +69 40 46']
#sample_list += ['PKS1424']
#sky_coord += ['14 27 00 +23 47 00']

SingleSourceAnalysis(sample_list,False)
