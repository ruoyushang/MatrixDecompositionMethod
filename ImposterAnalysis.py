
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

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn):

    return CommonPlotFunctions.GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isZoomIn)

def FluxBiasCorrection(real_flux,imposter_fluxes, imposter_biases):

    real_flux_correct = real_flux
    imposter_fluxes_correct = imposter_fluxes
    #return real_flux_correct, imposter_fluxes_correct

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for binx in range(0,real_flux[ebin].GetNbinsX()):
            for biny in range(0,real_flux[ebin].GetNbinsY()):
                imposter_bias = 0.
                for imposter in range(0,len(imposter_fluxes)):
                    imposter_bias += imposter_biases[imposter][ebin].GetBinContent(binx+1,biny+1)/float(len(imposter_fluxes))
                real_flux_old = real_flux_correct[ebin].GetBinContent(binx+1,biny+1)
                real_flux_correct[ebin].SetBinContent(binx+1,biny+1,real_flux_old-imposter_bias)
                for imposter in range(0,len(imposter_fluxes)):
                    imposter_flux_old = imposter_fluxes_correct[imposter][ebin].GetBinContent(binx+1,biny+1)
                    imposter_fluxes_correct[imposter][ebin].SetBinContent(binx+1,biny+1,imposter_flux_old-imposter_bias)
    return real_flux_correct, imposter_fluxes_correct

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
                total_error = pow(sum_square/5.,0.5)
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
        CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_reflect,None,fig,'RA','Dec','Z score','SkymapZscore_E%s'%(ebin))
        hist_excess_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_excess_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_excess_skymap_reflect,None,fig,'RA','Dec','excess count','SkymapExcess_E%s'%(ebin))
    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_sum_reflect,None,fig,'RA','Dec','Z score','SkymapZscore_Sum')

def power_law_func(x,a,b):
    return a*pow(x*1./1000.,b)

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

def flux_ic443_func(x):
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    #return 0.838*pow(10,-12)*pow(x*1./1000.,-2.99)
    # IC 443 https://arxiv.org/pdf/1512.01911.pdf
    return 9.92*pow(10,-12)*pow(x*1./550.,-2.8)
    #return 9.92*pow(10,-13)*pow(x*1./1000.,-2.8)
def flux_ic443_hawc_func(x):
    # 3HWC J0617+224 https://arxiv.org/pdf/2007.08582.pdf
    return 4.5*pow(10,-15)*pow(x*1./7000.,-3.05)

def MakeFluxMap():

    hist_real_flux_skymap_sum = ROOT.TH2D("hist_real_flux_skymap_sum","",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)
    for ebin in range(0,len(energy_bin)-1):
        hist_real_flux_skymap_sum.Add(hist_real_flux_skymap[ebin])

    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,None,fig,'RA','Dec','$E^{2}$ Flux [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux')


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
        syst_err = pow(syst_err/5.,0.5)
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
        syst_err = pow(syst_err/5.,0.5)
        real_s2b_syst_err += [syst_err]

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

def MakeExtensionProfile(roi_x,roi_y,roi_r,roi_name):

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
            syst_err = pow(syst_err/5.,0.5)
            real_profile_syst_err += [syst_err]

        theta2 = np.array(theta2)
        theta2_err = np.array(theta2_err)
        real_profile = np.array(real_profile)
        real_profile_stat_err = np.array(real_profile_stat_err)
        real_profile_syst_err = np.array(real_profile_syst_err)

        fig.clf()
        axbig = fig.add_subplot()
        axbig.bar(theta2, 2.*real_profile_syst_err, bottom=real_profile-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
        axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[ebin],energy_bin[ebin+1]))
        axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
        axbig.set_xlabel('angular distance from source [degree]')
        axbig.legend(loc='best')
        plotname = 'ProfileVsTheta2_%s'%(roi_name)
        fig.savefig("output_plots/%s_E%s_%s.png"%(plotname,ebin,plot_tag),bbox_inches='tight')
        axbig.remove()



#folder_name = 'output_2x2'
#folder_name = 'output_8x8'
folder_name = 'output_test'
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


hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_flux_skymap = []
hist_real_flux_syst_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_syst_skymap += [ROOT.TH2D("hist_real_flux_syst_skymap_E%s"%(ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
InputFile = ROOT.TFile("output_fitting/%s_ON_skymap_%s.root"%(source_name,folder_name))
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    HistName = "hist_data_skymap_%s"%(ebin)
    hist_real_data_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    hist_real_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_energy_flux_skymap_%s"%(ebin)
    hist_real_flux_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

hist_imposter_data_skymap = []
hist_imposter_bkgd_skymap = []
hist_imposter_bias_skymap = []
hist_imposter_flux_skymap = []
for imposter in range(0,5):
    hist_imposter_data_skymap_sublist = []
    hist_imposter_bkgd_skymap_sublist = []
    hist_imposter_bias_skymap_sublist = []
    hist_imposter_flux_skymap_sublist = []
    for ebin in range(0,len(energy_bin)-1):
        hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter%s_data_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bkgd_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bias_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bias_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
        hist_imposter_flux_skymap_sublist += [ROOT.TH2D("hist_imposter%s_flux_skymap_E%s"%(imposter,ebin),"",nbins,MapEdge_left,MapEdge_right,nbins,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
    hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
    hist_imposter_bias_skymap += [hist_imposter_bias_skymap_sublist]
    hist_imposter_flux_skymap += [hist_imposter_flux_skymap_sublist]
for imposter in range(0,5):
    InputFile = ROOT.TFile("output_fitting/%s_Imposter%s_skymap_%s.root"%(source_name,imposter+1,folder_name))
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        HistName = "hist_data_skymap_%s"%(ebin)
        hist_imposter_data_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.)
        HistName = "hist_bkgd_skymap_%s"%(ebin)
        hist_imposter_bkgd_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName))
        HistName = "hist_energy_flux_skymap_%s"%(ebin)
        hist_imposter_flux_skymap[imposter][ebin].Add(InputFile.Get(HistName))
    InputFile.Close()

hist_real_bkgd_skymap, hist_imposter_bkgd_skymap = FluxBiasCorrection(hist_real_bkgd_skymap, hist_imposter_bkgd_skymap, hist_imposter_bias_skymap)
hist_real_flux_skymap, hist_imposter_flux_skymap = FluxBiasCorrection(hist_real_flux_skymap, hist_imposter_flux_skymap, hist_imposter_flux_skymap)





MakeSignificanceMap()
MakeFluxMap()

#if source_name=='MGRO_J1908':
#    #3HWC J1908+063, 287.05, 6.39
#    region_x = 287.05
#    region_y = 6.39
#    region_r = 1.5
#    #PSR J1907+0602
#    #region_x = 286.975
#    #region_y = 6.03777777778
#    #region_r = 1.5
#    #North hot spot
#    #region_x = 286.8
#    #region_y = 7.1
#    #region_r = 0.3
#    #SS 433 e1
#    #region_x = 288.404
#    #region_y = 4.930
#    #region_r = 0.3
#if source_name=='1ES0229':
#    region_x = 38.222
#    region_y = 20.273
#    region_r = 1.5
#if source_name=='NGC1275':
#    # NGC 1275
#    region_x = 49.950
#    region_y = 41.512
#    region_r = 1.5
#    # IC 310
#    #region_x = 49.179
#    #region_y = 41.325
#    #region_r = 0.2
if source_name=='IC443HotSpot':
    region_x = 94.213
    region_y = 22.503
    region_r = 0.5
    region_name = 'Center'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,region_name)
elif source_name=='WComae':
    # 1ES 1218+304
    region_x = 185.360
    region_y = 30.191
    region_r = 0.2
    region_name = '1ES1218'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,region_name)
    # 1ES 1215+303
    region_x = 184.616
    region_y = 30.130
    region_r = 0.2
    region_name = '1ES1215'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,region_name)
    # W Comae
    region_x = 185.382
    region_y = 28.233
    region_r = 0.2
    region_name = 'WComae'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,region_name)
else:
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = 1.5
    region_name = 'Center'
    MakeSpectrum(region_x,region_y,region_r,region_name)
    MakeExtensionProfile(region_x,region_y,region_r,region_name)
