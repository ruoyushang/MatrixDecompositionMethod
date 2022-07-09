
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
energy_index_scale = CommonPlotFunctions.energy_index_scale
energy_bin = CommonPlotFunctions.energy_bin
energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])

def GetFluxCalibration(map_x,map_y,energy):


    binx = hist_real_data_skymap[0].GetXaxis().FindBin(map_x)
    biny = hist_real_data_skymap[0].GetYaxis().FindBin(map_y)
    data_count = hist_real_data_skymap[0].GetBinContent(binx,biny)
    if data_count==0.: return 0.
    
    #return 1.
    # 9 bins
    flux_calibration = [3.2538602333424753e-09, 1.280723442573218e-10, 2.7261479794937783e-12, 4.6984480245008585e-14]

    return flux_calibration[energy]

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn):

    return CommonPlotFunctions.GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn)

def MakeSignificanceMap():

    hist_total_err_skymap = []
    hist_syst_err_skymap = []
    hist_avg_bkgd_skymap = []
    hist_zscore_skymap = []
    hist_excess_skymap = []
    for ebin in range(0,len(energy_bin)-1):
        hist_total_err_skymap += [ROOT.TH2D("hist_total_err_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_syst_err_skymap += [ROOT.TH2D("hist_syst_err_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_avg_bkgd_skymap += [ROOT.TH2D("hist_avg_bkgd_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_zscore_skymap += [ROOT.TH2D("hist_zscore_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_excess_skymap += [ROOT.TH2D("hist_excess_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for binx in range(0,hist_imposter_data_skymap[0][ebin].GetNbinsX()):
            for biny in range(0,hist_imposter_data_skymap[0][ebin].GetNbinsY()):
                sum_square = 0.
                sum_bkgd = 0.
                for imposter in range(0,5):
                    imposter_data = hist_imposter_data_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    imposter_bkgd = hist_imposter_bkgd_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    sum_square += pow(imposter_data-imposter_bkgd,2)
                    sum_bkgd += imposter_bkgd
                total_error = pow(sum_square/n_imposters,0.5)
                avg_bkgd = sum_bkgd/5.
                #print ('ebin %s, binx %s, biny %s, avg_bkgd %s '%(ebin,binx,biny,avg_bkgd))
                syst_error = pow(max(0.,total_error*total_error-avg_bkgd),0.5)
                if avg_bkgd==0.: continue
                hist_total_err_skymap[ebin].SetBinContent(binx+1,biny+1,total_error)
                hist_syst_err_skymap[ebin].SetBinContent(binx+1,biny+1,syst_error)
                hist_avg_bkgd_skymap[ebin].SetBinContent(binx+1,biny+1,avg_bkgd)

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        skymap = GetSignificanceMap(hist_real_data_skymap[ebin],hist_real_bkgd_skymap[ebin],hist_syst_err_skymap[ebin],False)
        hist_zscore_skymap[ebin].Add(skymap)
        hist_excess_skymap[ebin].Add(hist_real_data_skymap[ebin])
        hist_excess_skymap[ebin].Add(hist_real_bkgd_skymap[ebin],-1.)

    hist_zscore_skymap_sum = ROOT.TH2D("hist_zscore_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    hist_real_data_skymap_sum = ROOT.TH2D("hist_real_data_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    hist_real_bkgd_skymap_sum = ROOT.TH2D("hist_real_bkgd_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    hist_real_syst_skymap_sum = ROOT.TH2D("hist_real_syst_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_real_data_skymap_sum.Add(hist_real_data_skymap[ebin])
        hist_real_bkgd_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
        for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
            for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
                syst_err = hist_syst_err_skymap[ebin].GetBinContent(binx+1,biny+1)
                old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,old_syst_err+syst_err*syst_err)
    for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
        for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
            old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
            hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,pow(old_syst_err,0.5))
    hist_zscore_skymap_sum.Add(GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum,hist_real_syst_skymap_sum,False))

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        #hist_total_err_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_total_err_skymap[ebin])
        #CommonPlotFunctions.MatplotlibMap2D(hist_total_err_skymap_reflect,None,fig,'RA','Dec','total error','SkymapTotalErr_E%s'%(ebin))
        #hist_syst_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_syst_err_skymap[ebin])
        #CommonPlotFunctions.MatplotlibMap2D(hist_syst_skymap_reflect,None,fig,'RA','Dec','systematic error','SkymapSyst_E%s'%(ebin))
        #hist_avg_bkgd_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_avg_bkgd_skymap[ebin])
        #CommonPlotFunctions.MatplotlibMap2D(hist_avg_bkgd_skymap_reflect,None,fig,'RA','Dec','avg. bkg','SkymapBkgd_E%s'%(ebin))
        hist_zscore_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_reflect,None,fig,'RA','Dec','Z score','SkymapZscore_E%s_%s'%(ebin,plot_tag))
        hist_excess_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_excess_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_excess_skymap_reflect,None,fig,'RA','Dec','excess count','SkymapExcess_E%s_%s'%(ebin,plot_tag))

    hist_real_zscore_skymap_sum.Reset()
    hist_real_zscore_skymap_sum.Add(hist_zscore_skymap_sum)
    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_sum_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Z score','SkymapZscore_Sum_%s'%(plot_tag))

    CommonPlotFunctions.SaveAsFITS(hist_zscore_skymap_sum_reflect,'skymap_zscore_sum_%s'%(plot_tag))

def diffusion_func(x,A,d):
    return A*1.22/(pow(3.14,1.5)*d*(x+0.06*d))*np.exp(-x*x/(d*d))

def power_law_func(x,a,b):
    return a*pow(x*1./1000.,b)

def flux_crab_func(x):
    # Crab https://arxiv.org/pdf/1508.06442.pdf
    return 37.5*pow(10,-12)*pow(x*1./1000.,-2.467-0.16*log(x/1000.))

def flux_hawc_j1908_func(x):
    # MGRO J1908  TeV^{-1}cm^{-2}s^{-1}
    return 0.95*pow(10,-13)*pow(x*1./10000.,-2.46-0.11*log(x/10000.))

def GetHAWCFluxJ1908(energy_index):
    energies = [pow(10.,12.00),pow(10.,12.20),pow(10.,12.47),pow(10.,12.74),pow(10.,13.01),pow(10.,13.29),pow(10.,13.55),pow(10.,13.78),pow(10.,14.00),pow(10.,14.21)]
    fluxes = [pow(10.,-10.60),pow(10.,-10.62),pow(10.,-10.64),pow(10.,-10.74),pow(10.,-10.86),pow(10.,-11.03),pow(10.,-11.18),pow(10.,-11.36),pow(10.,-11.37),pow(10.,-11.55)]
    flux_errs = [pow(10.,-10.60),pow(10.,-10.62),pow(10.,-10.64),pow(10.,-10.74),pow(10.,-10.86),pow(10.,-11.03),pow(10.,-11.18),pow(10.,-11.36),pow(10.,-11.37),pow(10.,-11.55)]
    flux_errs_up = [pow(10.,-10.57),pow(10.,-10.58),pow(10.,-10.61),pow(10.,-10.70),pow(10.,-10.82),pow(10.,-10.99),pow(10.,-11.14),pow(10.,-11.29),pow(10.,-11.28),pow(10.,-11.38)]
    flux_errs_low = [pow(10.,-10.64),pow(10.,-10.65),pow(10.,-10.68),pow(10.,-10.77),pow(10.,-10.89),pow(10.,-11.07),pow(10.,-11.22),pow(10.,-11.43),pow(10.,-11.48),pow(10.,-11.81)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e9
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxJ1908(energy_index):
    energies = [pow(10.,10.6),pow(10.,10.86),pow(10.,11.12),pow(10.,11.37)]
    fluxes = [pow(10.,-11.31),pow(10.,-11.29),pow(10.,-11.12),pow(10.,-10.88)]
    flux_errs = [pow(10.,-11.31),pow(10.,-11.29),pow(10.,-11.12),pow(10.,-10.88)]
    flux_errs_up = [pow(10.,-11.17),pow(10.,-11.12),pow(10.,-10.96),pow(10.,-10.76)]
    flux_errs_low = [pow(10.,-11.53),pow(10.,-11.57),pow(10.,-11.38),pow(10.,-11.05)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e9
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiUpperLimitFluxGeminga(energy_index):
    energies = [pow(10.,2.09),pow(10.,2.35),pow(10.,2.61),pow(10.,2.87)]
    fluxes = [pow(10.,-7.35),pow(10.,-7.23),pow(10.,-7.34),pow(10.,-7.18)]
    fluxes_err = [pow(10.,-7.35),pow(10.,-7.23),pow(10.,-7.34),pow(10.,-7.18)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        fluxes_err[entry] = fluxes[entry]*0.3

    return energies, fluxes, fluxes_err

def GetFermiHAWCFluxGeminga(energy_index):
    energies = [pow(10.,3.90),pow(10.,4.60)]
    fluxes = [pow(10.,-8.13),pow(10.,-8.36)]
    flux_errs = [0.,0.]
    flux_errs_up = [pow(10.,-8.06),pow(10.,-8.30)]
    flux_errs_low = [pow(10.,-8.24),pow(10.,-8.47)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxGeminga(energy_index):
    energies = [pow(10.,1.03),pow(10.,1.30),pow(10.,1.56),pow(10.,1.82)]
    fluxes = [pow(10.,-7.27),pow(10.,-7.29),pow(10.,-7.41),pow(10.,-7.35)]
    flux_errs = [pow(10.,-7.27),pow(10.,-7.29),pow(10.,-7.41),pow(10.,-7.35)]
    flux_errs_up = [pow(10.,-7.14),pow(10.,-7.16),pow(10.,-7.29),pow(10.,-7.23)]
    flux_errs_low = [pow(10.,-7.46),pow(10.,-7.49),pow(10.,-7.58),pow(10.,-7.50)]

    GeV_to_TeV = 1e-3
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]
        fluxes[entry] = fluxes[entry]*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*GeV_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def flux_ic443_func(x):
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    #return 0.838*pow(10,-12)*pow(x*1./1000.,-2.99)
    # IC 443 https://arxiv.org/pdf/1512.01911.pdf
    return 9.92*pow(10,-12)*pow(x*1./550.,-2.8)
    #return 9.92*pow(10,-13)*pow(x*1./1000.,-2.8)
def flux_ic443_hawc_func(x):
    # 3HWC J0617+224 https://arxiv.org/pdf/2007.08582.pdf
    return 4.5*pow(10,-15)*pow(x*1./7000.,-3.05)

def ConvertGalacticToRaDec(l, b):
    return CommonPlotFunctions.ConvertGalacticToRaDec(l, b)

def ConvertRaDecToGalactic(ra, dec):
    return CommonPlotFunctions.ConvertRaDecToGalactic(ra, dec)

def ConvertRaDecToGalacticMap(hist_map_radec,hist_map_galactic):

    hist_map_galactic.Reset()
    for bx1 in range(1,hist_map_galactic.GetNbinsX()+1):
        for by1 in range(1,hist_map_galactic.GetNbinsY()+1):
            locationx1 = hist_map_galactic.GetXaxis().GetBinCenter(bx1)
            locationy1 = hist_map_galactic.GetYaxis().GetBinCenter(by1)
            locationx2, locationy2 = ConvertGalacticToRaDec(locationx1,locationy1)
            bx2 = hist_map_radec.GetXaxis().FindBin(locationx2)
            by2 = hist_map_radec.GetYaxis().FindBin(locationy2)
            content = hist_map_radec.GetBinContent(bx2,by2)
            hist_map_galactic.SetBinContent(bx1,by1,content)
    #return hist_map_galactic

def SumFluxMap():

    hist_real_flux_skymap_sum.Reset()
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_real_flux_skymap_sum.Add(hist_real_flux_skymap[ebin])

    for imposter in range(0,5):
        hist_imposter_flux_skymap_sum[imposter].Reset()
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            hist_imposter_flux_skymap_sum[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)
    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_%s'%(plot_tag))

    MapCenter_l, MapCenter_b = ConvertRaDecToGalactic(MapCenter_x, MapCenter_y)
    hist_real_flux_skymap_sum_galactic = ROOT.TH2D("hist_real_flux_skymap_sum_galactic","",nbins,MapCenter_l-0.5*MapPlotSize,MapCenter_l+0.5*MapPlotSize,nbins,MapCenter_b-0.5*MapPlotSize,MapCenter_b+0.5*MapPlotSize)
    ConvertRaDecToGalacticMap(hist_real_flux_skymap_sum,hist_real_flux_skymap_sum_galactic)
    hist_real_flux_skymap_galactic_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum_galactic)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_galactic_reflect,None,fig,'gal. l','gal. b','$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFluxGal_%s'%(plot_tag))

    CommonPlotFunctions.SaveAsFITS(hist_real_flux_skymap_reflect,'skymap_flux_sum_%s'%(plot_tag))

def MakeSpectrum(roi_x,roi_y,roi_r,roi_name):

    energy_axis, energy_error, real_flux, real_flux_stat_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_real_flux_skymap,hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
    imposter_flux_list = []
    imposter_flux_err_list = []
    for imposter in range(0,5):
        energy_axis, energy_error, imposter_flux, imposter_flux_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_flux_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
        imposter_flux_list += [imposter_flux]
        imposter_flux_err_list += [imposter_flux_err]

    real_flux_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_err = 0.
        for imposter in range(0,5):
            syst_err += pow(imposter_flux_list[imposter][ebin],2)
        syst_err = pow(syst_err/n_imposters,0.5)
        real_flux_syst_err += [syst_err]

    imposter_data_list = []
    imposter_bkgd_list = []
    for imposter in range(0,5):
        energy_axis, energy_error, imposter_data, imposter_data_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_data_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
        energy_axis, energy_error, imposter_bkgd, imposter_bkgd_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_bkgd_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r)
        imposter_data_list += [imposter_data]
        imposter_bkgd_list += [imposter_bkgd]
    imposter_s2b_list = []
    imposter_s2b_err_list = []
    for imposter in range(0,5):
        imposter_s2b = []
        imposter_s2b_err = []
        for ebin in range(0,len(energy_axis)):
            s2b = (imposter_data_list[imposter][ebin]-imposter_bkgd_list[imposter][ebin])/imposter_bkgd_list[imposter][ebin]
            s2b_err = pow(imposter_data_list[imposter][ebin],0.5)/imposter_bkgd_list[imposter][ebin]
            imposter_s2b += [s2b]
            imposter_s2b_err += [s2b_err]
        imposter_s2b_list += [imposter_s2b]
        imposter_s2b_err_list += [imposter_s2b_err]

    real_s2b_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_err = 0.
        for imposter in range(0,5):
            syst_err += pow(imposter_s2b_list[imposter][ebin],2)
        syst_err = pow(syst_err/n_imposters,0.5)
        real_s2b_syst_err += [syst_err]

    if 'Crab' in source_name:
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_crab = np.vectorize(flux_crab_func)
        ydata_crab = pow(xdata/1e3,energy_index_scale)*vectorize_f_crab(xdata)
        xdata_array = []
        for binx in range(0,len(energy_axis)):
            xdata_array += [energy_axis[binx]]
        ydata_veritas = pow(np.array(xdata_array)/1e3,energy_index_scale)*vectorize_f_crab(xdata_array)
        calibration_new = []
        for binx in range(0,len(energy_axis)):
            if real_flux[binx]>0.:
                calibration_new += [ydata_veritas[binx]/real_flux[binx]]
            else:
                calibration_new += [0.]
        print ('=======================================================================')
        print ('new flux_calibration = %s'%(calibration_new))
        print ('=======================================================================')
    if source_name=='MGRO_J1908':
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
        ydata_hawc = pow(xdata_ref/1e3,energy_index_scale)*vectorize_f_hawc(xdata_ref)
        # HAWC systematic uncertainty, The Astrophysical Journal 881, 134. Fig 13
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1908(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1908(energy_index_scale)
    if source_name=='IC443HotSpot':
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_ic443_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)
        vectorize_f_hawc = np.vectorize(flux_ic443_hawc_func)
        ydata_hawc = pow(xdata/1e3,energy_index_scale)*vectorize_f_hawc(xdata)
    if source_name=='Geminga':
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetFermiHAWCFluxGeminga(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxGeminga(energy_index_scale)
        Fermi_UL_energies, Fermi_UL_fluxes, Fermi_UL_err = GetFermiUpperLimitFluxGeminga(energy_index_scale)

    energy_axis = np.array(energy_axis)
    energy_error = np.array(energy_error)
    real_flux = np.array(real_flux)
    real_flux_stat_err = np.array(real_flux_stat_err)
    real_flux_syst_err = np.array(real_flux_syst_err)
    real_s2b_syst_err = np.array(real_s2b_syst_err)

    fig.clf()
    axbig = fig.add_subplot()
    cycol = cycle('krgbcmy')
    for imposter in range(0,5):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_flux_list[imposter],imposter_flux_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=0.-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    plotname = 'ImposterSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    cycol = cycle('krgbcmy')
    for imposter in range(0,5):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_s2b_list[imposter],imposter_s2b_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_s2b_syst_err, bottom=0.-real_s2b_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    for ebin in range(0,len(energy_axis)):
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.*real_s2b_syst_err[len(real_s2b_syst_err)-1], '%0.1f %%'%(real_s2b_syst_err[ebin]*100.))
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('relative error')
    axbig.set_xscale('log')
    plotname = 'ImposterRelError_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    axbig = fig.add_subplot()
    if source_name=='MGRO_J1908':
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='g',marker='s',ls='none',label='Fermi',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC',zorder=3)
        axbig.plot(xdata_ref, ydata_hawc,'r-',label='1909.08609 (HAWC)',zorder=2)
        axbig.fill_between(xdata_ref, ydata_hawc-0.15*ydata_hawc, ydata_hawc+0.15*ydata_hawc, alpha=0.2, color='r',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS',zorder=6)
    elif source_name=='IC443HotSpot':
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (0905.3291)',zorder=1)
        axbig.plot(xdata, ydata_hawc,'g-',label='HAWC (2007.08582)',zorder=2)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=3)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=4)
    elif 'Crab' in source_name:
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    elif 'Geminga' in source_name:
        axbig.plot(HAWC_energies, HAWC_fluxes,'g-',label='HAWC')
        axbig.fill_between(HAWC_energies, np.array(HAWC_fluxes)-np.array(HAWC_flux_errs), np.array(HAWC_fluxes)+np.array(HAWC_flux_errs), alpha=0.2, color='g')
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='r',marker='_',ls='none',label='Fermi')
        uplims = np.array([1,1,1,1], dtype=bool)
        axbig.errorbar(Fermi_UL_energies,Fermi_UL_fluxes,Fermi_UL_err,color='r',marker='_',ls='none',uplims=uplims)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    else:
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'RealSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

def MakeDiffusionSpectrum(energy_axis,energy_error,r_axis,y_axis,y_axis_stat_err,y_axis_syst_err):

    real_flux = np.zeros(len(energy_axis))
    real_flux_stat_err = np.zeros(len(energy_axis))
    real_flux_syst_err = np.zeros(len(energy_axis))
    for eb in range(0,len(energy_axis)):
        flux = 0.
        flux_stat_err = 0.
        flux_syst_err = 0.
        for rb in range(0,len(r_axis[eb])):
            radius = r_axis[eb][rb]
            delta_r = r_axis[eb][1]-r_axis[eb][0]
            brightness = y_axis[eb][rb]
            brightness_stat_err = y_axis_stat_err[eb][rb]
            brightness_syst_err = y_axis_syst_err[eb][rb]
            flux += brightness*2.*3.14*radius*delta_r
            flux_stat_err += pow(brightness_stat_err*2.*3.14*radius*delta_r,2)
            flux_syst_err += brightness_syst_err*2.*3.14*radius*delta_r
        real_flux[eb] = flux
        real_flux_stat_err[eb] = pow(flux_stat_err,0.5)
        real_flux_syst_err[eb] = flux_syst_err

    if source_name=='MGRO_J1908':
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
        ydata_hawc = pow(xdata_ref/1e3,energy_index_scale)*vectorize_f_hawc(xdata_ref)
        # HAWC systematic uncertainty, The Astrophysical Journal 881, 134. Fig 13
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1908(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1908(energy_index_scale)
    if source_name=='IC443HotSpot':
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_ic443_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)
        vectorize_f_hawc = np.vectorize(flux_ic443_hawc_func)
        ydata_hawc = pow(xdata/1e3,energy_index_scale)*vectorize_f_hawc(xdata)
    if source_name=='Geminga':
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetFermiHAWCFluxGeminga(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxGeminga(energy_index_scale)
        Fermi_UL_energies, Fermi_UL_fluxes, Fermi_UL_err = GetFermiUpperLimitFluxGeminga(energy_index_scale)

    energy_axis = np.array(energy_axis)
    energy_error = np.array(energy_error)
    real_flux = np.array(real_flux)
    real_flux_stat_err = np.array(real_flux_stat_err)
    real_flux_syst_err = np.array(real_flux_syst_err)
    fig.clf()
    axbig = fig.add_subplot()
    if source_name=='MGRO_J1908':
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='g',marker='s',ls='none',label='Fermi',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC',zorder=3)
        axbig.plot(xdata_ref, ydata_hawc,'r-',label='1909.08609 (HAWC)',zorder=2)
        axbig.fill_between(xdata_ref, ydata_hawc-0.15*ydata_hawc, ydata_hawc+0.15*ydata_hawc, alpha=0.2, color='r',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS',zorder=6)
    elif source_name=='IC443HotSpot':
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (0905.3291)',zorder=1)
        axbig.plot(xdata, ydata_hawc,'g-',label='HAWC (2007.08582)',zorder=2)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=3)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=4)
    elif 'Crab' in source_name:
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    elif 'Geminga' in source_name:
        axbig.plot(HAWC_energies, HAWC_fluxes,'g-',label='HAWC')
        axbig.fill_between(HAWC_energies, np.array(HAWC_fluxes)-np.array(HAWC_flux_errs), np.array(HAWC_fluxes)+np.array(HAWC_flux_errs), alpha=0.2, color='g')
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='r',marker='_',ls='none',label='Fermi')
        uplims = np.array([1,1,1,1], dtype=bool)
        axbig.errorbar(Fermi_UL_energies,Fermi_UL_fluxes,Fermi_UL_err,color='r',marker='_',ls='none',uplims=uplims)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    else:
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'ProfileSpectrum'
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

def MakeExtensionProfile(roi_x,roi_y,roi_r,fit_profile,roi_name):

    r_axis_allE = []
    y_axis_allE = []
    y_axis_stat_err_allE = []
    y_axis_syst_err_allE = []
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(hist_real_flux_skymap[ebin],hist_real_flux_syst_skymap[ebin],roi_x,roi_y,2.0)
        imposter_profile_list = []
        imposter_profile_err_list = []
        for imposter in range(0,5):
            imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(hist_imposter_flux_skymap[imposter][ebin],hist_real_flux_syst_skymap[ebin],roi_x,roi_y,2.0)
            imposter_profile_list += [imposter_profile]
            imposter_profile_err_list += [imposter_profile_stat_err]

        real_profile_syst_err = []
        for ubin in range(0,len(theta2)):
            syst_err = 0.
            for imposter in range(0,5):
                syst_err += pow(imposter_profile_list[imposter][ubin],2)
            syst_err = pow(syst_err/n_imposters,0.5)
            real_profile_syst_err += [syst_err]

        theta2 = np.array(theta2)
        theta2_err = np.array(theta2_err)
        real_profile = np.array(real_profile)
        real_profile_stat_err = np.array(real_profile_stat_err)
        real_profile_syst_err = np.array(real_profile_syst_err)

        if fit_profile:
            start = (real_profile[0], 0.5)
            popt, pcov = curve_fit(diffusion_func,theta2,real_profile,p0=start,sigma=real_profile_stat_err)
            profile_fit = diffusion_func(theta2, *popt)
            residual = real_profile - profile_fit
            chisq = np.sum((residual/real_profile_stat_err)**2)
            dof = len(theta2)-2
            print ('diffusion radius = %0.2f deg (chi2/dof = %0.2f)'%(popt[1],chisq/dof))
            r_axis_allE += [theta2]
            y_axis_allE += [profile_fit]
            y_axis_stat_err_allE += [real_profile_stat_err]
            y_axis_syst_err_allE += [real_profile_syst_err]

        fig.clf()
        axbig = fig.add_subplot()
        #axbig.bar(theta2, 2.*real_profile_syst_err, bottom=real_profile-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
        axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
        axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[ebin],energy_bin[ebin+1]))
        if fit_profile:
            axbig.plot(theta2,diffusion_func(theta2,*popt),color='r')
        axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
        axbig.set_xlabel('angular distance from source [degree]')
        axbig.legend(loc='best')
        plotname = 'ProfileVsTheta2_%s'%(roi_name)
        fig.savefig("output_plots/%s_E%s_%s.png"%(plotname,ebin,plot_tag),bbox_inches='tight')
        axbig.remove()

    real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(hist_real_flux_skymap_sum,hist_real_flux_syst_skymap[0],roi_x,roi_y,2.0)
    imposter_profile_list = []
    imposter_profile_err_list = []
    for imposter in range(0,5):
        imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(hist_imposter_flux_skymap_sum[imposter],hist_real_flux_syst_skymap[0],roi_x,roi_y,2.0)
        imposter_profile_list += [imposter_profile]
        imposter_profile_err_list += [imposter_profile_stat_err]

    real_profile_syst_err = []
    for ubin in range(0,len(theta2)):
        syst_err = 0.
        for imposter in range(0,5):
            syst_err += pow(imposter_profile_list[imposter][ubin],2)
        syst_err = pow(syst_err/n_imposters,0.5)
        real_profile_syst_err += [syst_err]

    theta2 = np.array(theta2)
    theta2_err = np.array(theta2_err)
    real_profile = np.array(real_profile)
    real_profile_stat_err = np.array(real_profile_stat_err)
    real_profile_syst_err = np.array(real_profile_syst_err)

    if fit_profile:
        start = (real_profile[0], 0.5)
        popt, pcov = curve_fit(diffusion_func,theta2,real_profile,p0=start,sigma=real_profile_stat_err)
        profile_fit = diffusion_func(theta2, *popt)
        residual = real_profile - profile_fit
        chisq = np.sum((residual/real_profile_stat_err)**2)
        dof = len(theta2)-2
        print ('diffusion radius = %0.2f deg (chi2/dof = %0.2f)'%(popt[1],chisq/dof))

    fig.clf()
    axbig = fig.add_subplot()
    #axbig.bar(theta2, 2.*real_profile_syst_err, bottom=real_profile-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
    if fit_profile:
        axbig.plot(theta2,diffusion_func(theta2,*popt),color='r')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileVsTheta2_%s'%(roi_name)
    fig.savefig("output_plots/%s_sum_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    if fit_profile:
        energy_axis = []
        energy_axis_err = []
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            energy_axis += [0.5*(energy_bin[ebin]+energy_bin[ebin+1])]
            energy_axis_err += [0.5*(energy_bin[ebin+1]-energy_bin[ebin])]
        MakeDiffusionSpectrum(energy_axis,energy_axis_err,r_axis_allE,y_axis_allE,y_axis_stat_err_allE,y_axis_syst_err_allE)

def MakeFluxMap(flux_map, data_map, bkgd_map, expo_map):

    skymap_bin_size_x = data_map[0].GetXaxis().GetBinCenter(2)-data_map[0].GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = data_map[0].GetYaxis().GetBinCenter(2)-data_map[0].GetYaxis().GetBinCenter(1)
    calibration_radius = CommonPlotFunctions.calibration_radius
    for ebin in range(0,len(energy_bin)-1):
        flux_map[ebin].Reset()
        expo_content_max = expo_map[ebin].GetMaximum()
        for binx in range(0,bkgd_map[ebin].GetNbinsX()):
            for biny in range(0,bkgd_map[ebin].GetNbinsY()):
                data_content = data_map[ebin].GetBinContent(binx+1,biny+1)
                bkgd_content = bkgd_map[ebin].GetBinContent(binx+1,biny+1)
                expo_content = expo_map[ebin].GetBinContent(binx+1,biny+1)
                expo_content = max(expo_content,0.5*expo_content_max)
                map_x = data_map[ebin].GetXaxis().GetBinCenter(binx+1)
                map_y = data_map[ebin].GetYaxis().GetBinCenter(biny+1)
                flux_calibration = GetFluxCalibration(map_x,map_y,ebin)
                #correction = flux_calibration*(skymap_bin_size_x*skymap_bin_size_y)/(3.14*calibration_radius*calibration_radius)
                correction = flux_calibration*map_bin_area/calibration_bin_area
                stat_data_err = data_map[ebin].GetBinError(binx+1,biny+1)
                stat_bkgd_err = bkgd_map[ebin].GetBinError(binx+1,biny+1)
                flux_stat_err = pow(stat_data_err*stat_data_err+stat_bkgd_err*stat_bkgd_err,0.5)/expo_content*correction*pow(energy_bin[ebin]/1e3,energy_index_scale)
                flux_content = (data_content-bkgd_content)/expo_content*correction*pow(energy_bin[ebin]/1e3,energy_index_scale)
                flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content)
                flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err)



folder_name = CommonPlotFunctions.folder_path
source_name = sys.argv[1]
source_name = source_name.split('_ON')[0]
plot_tag = source_name+'_'+folder_name

InputFile = ROOT.TFile("output_fitting/%s_ON_skymap_%s.root"%(source_name,folder_name))
HistName = "hist_data_skymap_0"
nbins = InputFile.Get(HistName).GetNbinsX()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()

MapPlotSize = 4.
MapCenter_x = (MapEdge_left+MapEdge_right)/2.
MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.
MapPlot_left = MapCenter_x-0.5*MapPlotSize
MapPlot_right = MapCenter_x+0.5*MapPlotSize
MapPlot_lower = MapCenter_y-0.5*MapPlotSize
MapPlot_upper = MapCenter_y+0.5*MapPlotSize

calibration_bin_area = pow(MapPlotSize/9.,2)
map_bin_area = pow(MapPlotSize/float(nbins),2)

hist_real_zscore_skymap_sum = ROOT.TH2D("hist_real_zscore_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_sum = ROOT.TH2D("hist_real_flux_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_expo_skymap = []
hist_real_bias_skymap = []
hist_real_flux_skymap = []
hist_real_flux_syst_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_expo_skymap += [ROOT.TH2D("hist_real_expo_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_bias_skymap += [ROOT.TH2D("hist_real_bias_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_syst_skymap += [ROOT.TH2D("hist_real_flux_syst_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
InputFile = ROOT.TFile("output_fitting/%s_ON_skymap_%s.root"%(source_name,folder_name))
for ebin in range(0,len(energy_bin)-1):
    HistName = "hist_data_skymap_%s"%(ebin)
    hist_real_data_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    hist_real_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
    #HistName = "hist_expo_skymap_%s"%(ebin)
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    hist_real_expo_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

hist_imposter_data_skymap = []
hist_imposter_bkgd_skymap = []
hist_imposter_expo_skymap = []
hist_imposter_bias_skymap = []
hist_imposter_flux_skymap = []
hist_imposter_flux_skymap_sum = []
for imposter in range(0,5):
    hist_imposter_flux_skymap_sum += [ROOT.TH2D("hist_imposter%s_flux_skymap_sum"%(imposter),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sublist = []
    hist_imposter_bkgd_skymap_sublist = []
    hist_imposter_expo_skymap_sublist = []
    hist_imposter_bias_skymap_sublist = []
    hist_imposter_flux_skymap_sublist = []
    for ebin in range(0,len(energy_bin)-1):
        hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter%s_data_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bkgd_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_expo_skymap_sublist += [ROOT.TH2D("hist_imposter%s_expo_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bias_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bias_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_flux_skymap_sublist += [ROOT.TH2D("hist_imposter%s_flux_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
    hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
    hist_imposter_expo_skymap += [hist_imposter_expo_skymap_sublist]
    hist_imposter_bias_skymap += [hist_imposter_bias_skymap_sublist]
    hist_imposter_flux_skymap += [hist_imposter_flux_skymap_sublist]
for imposter in range(0,5):
    InputFile = ROOT.TFile("output_fitting/%s_Imposter%s_skymap_%s.root"%(source_name,imposter+1,folder_name))
    for ebin in range(0,len(energy_bin)-1):
        HistName = "hist_data_skymap_%s"%(ebin)
        hist_imposter_data_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        #HistName = "hist_expo_skymap_%s"%(ebin)
        HistName = "hist_bkgd_skymap_%s"%(ebin)
        hist_imposter_expo_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        HistName = "hist_bkgd_skymap_%s"%(ebin)
        hist_imposter_bkgd_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        data_norm = hist_imposter_data_skymap[imposter][ebin].Integral()
        bkgd_norm = hist_imposter_bkgd_skymap[imposter][ebin].Integral()
        hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.*data_norm/bkgd_norm)
        #hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.)
    InputFile.Close()

correct_bias = False
n_imposters = 5.
if correct_bias:
    for ebin in range(0,len(energy_bin)-1):
        for imposter in range(0,5):
            hist_real_bias_skymap[ebin].Add(hist_imposter_bias_skymap[imposter][ebin],1./5.)
    for ebin in range(0,len(energy_bin)-1):
        hist_real_bkgd_skymap[ebin].Add(hist_real_bias_skymap[ebin])
        for imposter in range(0,5):
            hist_imposter_bkgd_skymap[imposter][ebin].Add(hist_real_bias_skymap[ebin])
    n_imposters = 4.
else:
    plot_tag += '_NC'


MakeFluxMap(hist_real_flux_skymap, hist_real_data_skymap, hist_real_bkgd_skymap, hist_real_expo_skymap)
for imposter in range(0,5):
    MakeFluxMap(hist_imposter_flux_skymap[imposter], hist_imposter_data_skymap[imposter], hist_imposter_bkgd_skymap[imposter], hist_imposter_expo_skymap[imposter])


MakeSignificanceMap()
SumFluxMap()

if 'Crab' in source_name:
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = CommonPlotFunctions.calibration_radius
    region_name = 'Center'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,True,region_name)
elif source_name=='MGRO_J1908':
    #3HWC J1908+063, 287.05, 6.39
    region_x = 287.05
    region_y = 6.39
    region_r = 1.5
    region_name = '3HWC'
    #PSR J1907+0602
    #region_x = 286.975
    #region_y = 6.03777777778
    #region_r = 1.2
    #region_name = 'PSR'
    #PSR J1907+0602
    #North hot spot
    #region_x = 286.8
    #region_y = 7.1
    #region_r = 0.3
    #SS 433 e1
    #region_x = 288.404
    #region_y = 4.930
    #region_r = 0.3
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,True,region_name)

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)
    Hist_fermi5 = ROOT.TH2D("Hist_fermi5","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    Hist_fermi5 = CommonPlotFunctions.GetSkyViewMap("MWL_maps/skv1826930706371_j1908_fermi5.txt", Hist_fermi5, True)
    # 3-300 GeV Band 5, Atwood et al. 2009
    Hist_fermi5_reflect = CommonPlotFunctions.reflectXaxis(Hist_fermi5)
    CommonPlotFunctions.MatplotlibMap2D(Hist_fermi5_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','cnts$/s/cm^{2}/sr$','SkymapFermi5.png')
    Hist_hawc = ROOT.TH2D("Hist_hawc","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    Hist_hawc = CommonPlotFunctions.GetHawcSkymap(Hist_hawc, True)
    Hist_hawc_reflect = CommonPlotFunctions.reflectXaxis(Hist_hawc)
    CommonPlotFunctions.MatplotlibMap2D(Hist_hawc_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Significance','SkymapHAWC.png')

    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_25_50_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_55_85_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    Hist_mc_intensity.Scale(FITS_correction)
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapMolecularColumn.png')

elif source_name=='IC443HotSpot':
    region_x = 94.213
    region_y = 22.503
    region_r = 0.5
    region_name = 'Center'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,False,region_name)
elif source_name=='WComae':
    # 1ES 1218+304
    region_x = 185.360
    region_y = 30.191
    region_r = 0.2
    region_name = '1ES1218'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,False,region_name)
    ## 1ES 1215+303
    #region_x = 184.616
    #region_y = 30.130
    #region_r = 0.2
    #region_name = '1ES1215'
    #MakeSpectrum(region_x,region_y,region_r,region_name)
    #MakeExtensionProfile(region_x,region_y,region_r,False,region_name)
    ## W Comae
    #region_x = 185.382
    #region_y = 28.233
    #region_r = 0.2
    #region_name = 'WComae'
    #MakeSpectrum(region_x,region_y,region_r,region_name)
    #MakeExtensionProfile(region_x,region_y,region_r,False,region_name)
else:
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = 1.5
    region_name = 'Center'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,False,region_name)
