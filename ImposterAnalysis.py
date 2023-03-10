
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

import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

import CommonPlotFunctions

fig, ax = plt.subplots()
#figsize_x = 8.6
#figsize_y = 4.2
figsize_x = 8.6
figsize_y = 6.4
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

UseEffectiveArea = CommonPlotFunctions.UseEffectiveArea
energy_index_scale = CommonPlotFunctions.energy_index_scale
energy_bin = CommonPlotFunctions.energy_bin
energy_bin_cut_low = int(sys.argv[2])
energy_bin_cut_up = int(sys.argv[3])
doImposter = int(sys.argv[4])

plot_energy_break_bin = 3

#n_imposters = 3
#n_imposters = 4
n_imposters = 5
#n_imposters = 6

correct_bias = True
#correct_bias = False
#use_rfov = True
use_rfov = False

# Our function to fit is going to be a sum of two-dimensional Gaussians
def gaussian(x, y, x0, y0, sigma, A):
    #return A * np.exp( -((x-x0)/(2.*sigma))**2 -((y-y0)/(2.*sigma))**2)
    return A * np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma*sigma))/(2*np.pi*sigma*sigma)
# https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
# This is the callable that is passed to curve_fit. M is a (2,N) array
# where N is the total number of data points in Z, which will be ravelled
# to one dimension.
def _gaussian(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(len(args)//4):
       arr += gaussian(x, y, *args[i*4:i*4+4])
    return arr

def fit_2d_model(hist_map_data, hist_map_bkgd, src_x, src_y):

    print ('+++++++++++++++++++++++++++++++++++++++++++++')
    nbins_x = hist_map_data.GetNbinsX()
    nbins_y = hist_map_data.GetNbinsY()
    lon_min = MapEdge_left
    lon_max = MapEdge_right
    lat_min = MapEdge_lower
    lat_max = MapEdge_upper
    x_axis = np.linspace(lon_min,lon_max,nbins_x)
    y_axis = np.linspace(lat_min,lat_max,nbins_y)
    X_grid, Y_grid = np.meshgrid(x_axis, y_axis)
    # We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
    XY_stack = np.vstack((X_grid.ravel(), Y_grid.ravel()))

    image_data = np.zeros((hist_map_data.GetNbinsX(),hist_map_data.GetNbinsY()))
    image_error = np.zeros((hist_map_data.GetNbinsX(),hist_map_data.GetNbinsY()))
    for binx in range (0,hist_map_data.GetNbinsX()):
        for biny in range (0,hist_map_data.GetNbinsY()):
            image_data[biny,binx] = hist_map_data.GetBinContent(binx+1,biny+1) - hist_map_bkgd.GetBinContent(binx+1,biny+1)
            image_error[biny,binx] = max(1.,pow(hist_map_data.GetBinContent(binx+1,biny+1),0.5))

    #print ('set initial avlues and bounds')
    initial_prms = []
    bound_upper_prms = []
    bound_lower_prms = []
    lon = src_x
    lat = src_y
    sigma = 0.03807
    initial_prms += [(lon,lat,sigma,10.)]
    bound_lower_prms += [(lon-0.5,lat-0.5,sigma+0.0,0.)]
    bound_upper_prms += [(lon+0.5,lat+0.5,sigma+2.0,1e10)]
    # Flatten the initial guess parameter list.
    p0 = [p for prms in initial_prms for p in prms]
    p0_lower = [p for prms in bound_lower_prms for p in prms]
    p0_upper = [p for prms in bound_upper_prms for p in prms]
    print ('p0 = %s'%(p0))

    popt, pcov = curve_fit(_gaussian, XY_stack, image_data.ravel(), p0, sigma=image_error.ravel(), absolute_sigma=True, bounds=(p0_lower,p0_upper))
    fit_src_x = popt[0*4+0]
    fit_src_x_err = pow(pcov[0*4+0][0*4+0],0.5)
    print ('fit_src_x = %0.3f +/- %0.3f'%(fit_src_x,fit_src_x_err))
    fit_src_y = popt[0*4+1]
    fit_src_y_err = pow(pcov[0*4+1][0*4+1],0.5)
    print ('fit_src_y = %0.3f +/- %0.3f'%(fit_src_y,fit_src_y_err))
    fit_src_sigma = popt[0*4+2]
    fit_src_sigma_err = pow(pcov[0*4+2][0*4+2],0.5)
    print ('fit_src_sigma = %0.3f +/- %0.3f'%(fit_src_sigma,fit_src_sigma_err))
    fit_src_A = popt[0*4+3]
    print ('fit_src_A = %0.1e'%(fit_src_A))

def MakeDiagnisticPlots():

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)

    hist_real_elev_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_elev_skymap)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_elev_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','elevation [deg]','SkymapElev_%s.png'%(plot_tag),fill_gaps=True,rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    hist_real_azim_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_azim_skymap)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_azim_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','azim [deg]','SkymapAzim_%s.png'%(plot_tag),fill_gaps=True,rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    hist_real_nsb_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_nsb_skymap)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_nsb_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','NSB','SkymapNSB_%s.png'%(plot_tag),fill_gaps=True,rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    hist_real_mjd_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_mjd_skymap)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_mjd_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','MJD','SkymapMJD_%s.png'%(plot_tag),fill_gaps=True,rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    hist_real_expo_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_expo_skymap[0])
    CommonPlotFunctions.MatplotlibMap2D(hist_real_expo_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Hours','SkymapExpo_E%s_%s'%(0,plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    for ebin in range(0,len(energy_bin)-1):
        hist_real_norm_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_norm_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_real_norm_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Count','SkymapNorm_E%s_%s'%(ebin,plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    if doImposter==1:
        for imposter in range(0,n_imposters):
            print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print ('imposter set %s'%(imposter))
            MakeFluxMap(hist_imposter_flux_skymap[imposter], hist_imposter_data_skymap[imposter], hist_imposter_bkgd_skymap[imposter], hist_imposter_norm_skymap[imposter], hist_imposter_Aeff_skymap[imposter], hist_imposter_expo_skymap[imposter])

def GetFluxCalibration(energy):

    #return 1.

    # elevation = 70, energy threshold = 200 GeV
    #flux_calibration = [2.5579605915713342e-11, 9.417832458337564e-12, 4.186357340959233e-12, 2.1646972910242393e-12, 1.142544538194112e-12, 5.433497919384933e-13]

    # elevation = 60, energy threshold = 100 GeV
    flux_calibration = [1.7506329704556129e-10, 1.4669366489694362e-11, 5.575025320745524e-12, 2.6520194512426122e-12, 1.3856458120508363e-12, 6.633980234341466e-13]


    return flux_calibration[energy]

def GetEffectiveAreaCorrection(energy):

    return 1.
    flux_calibration = [0.7567235002074865, 0.3387015370792934, 0.33462622803554154, 0.43292054363330446, 0.495372487346452, 0.20459414829433548]

    return flux_calibration[energy]

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isON):

    return CommonPlotFunctions.GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,isON)

def FindDetectedPWN(hist_on_data_skymap,hist_on_bkgd_skymap,hist_mimic_data_skymap,hist_mimic_bkgd_skymap):

    Skymap_size = CommonPlotFunctions.Skymap_size_x
    target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = CommonPlotFunctions.ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
    if CommonPlotFunctions.doGalacticCoord:
        for src in range(0,len(target_psr_name)):
            target_psr_ra[src], target_psr_dec[src] = ConvertRaDecToGalactic(target_psr_ra[src], target_psr_dec[src])

    target_psr_detect = []
    for psr in range(0,len(target_psr_name)):
        target_psr_detect += [False]
        dist_to_map_center = pow(pow(target_psr_ra[psr]-MapCenter_x,2)+pow(target_psr_dec[psr]-MapCenter_y,2),0.5)
        if dist_to_map_center>Skymap_size: continue

        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            cnt_on_data = 0.
            cnt_on_bkgd = 0.
            for binx in range(0,hist_on_data_skymap[ebin].GetNbinsX()):
                for biny in range(0,hist_on_data_skymap[ebin].GetNbinsY()):
                    cell_x = hist_on_data_skymap[ebin].GetXaxis().GetBinCenter(binx+1)
                    cell_y = hist_on_data_skymap[ebin].GetYaxis().GetBinCenter(biny+1)
                    dist_deg_to_psr = pow(pow(target_psr_ra[psr]-cell_x,2)+pow(target_psr_dec[psr]-cell_y,2),0.5)
                    roi_radius_pc = 60.
                    roi_radius_deg = roi_radius_pc/(1000.*target_psr_dist[psr])*180./3.14
                    if dist_deg_to_psr>roi_radius_deg: continue
                    cnt_on_data += hist_on_data_skymap[ebin].GetBinContent(binx+1,biny+1)
                    cnt_on_bkgd += hist_on_bkgd_skymap[ebin].GetBinContent(binx+1,biny+1)
            if cnt_on_data==0.: continue
            imp_syst_err = 0.
            for imp in range(0,n_imposters):
                cnt_imp_data = 0.
                cnt_imp_bkgd = 0.
                for binx in range(0,hist_mimic_data_skymap[imp][ebin].GetNbinsX()):
                    for biny in range(0,hist_mimic_data_skymap[imp][ebin].GetNbinsY()):
                        cell_x = hist_mimic_data_skymap[imp][ebin].GetXaxis().GetBinCenter(binx+1)
                        cell_y = hist_mimic_data_skymap[imp][ebin].GetYaxis().GetBinCenter(biny+1)
                        dist_deg_to_psr = pow(pow(target_psr_ra[psr]-cell_x,2)+pow(target_psr_dec[psr]-cell_y,2),0.5)
                        roi_radius_pc = 60.
                        roi_radius_deg = roi_radius_pc/(1000.*target_psr_dist[psr])*180./3.14
                        if dist_deg_to_psr>roi_radius_deg: continue
                        cnt_imp_data += hist_mimic_data_skymap[imp][ebin].GetBinContent(binx+1,biny+1)
                        cnt_imp_bkgd += hist_mimic_bkgd_skymap[imp][ebin].GetBinContent(binx+1,biny+1)
                imp_syst_err += pow(cnt_imp_data-cnt_imp_bkgd,2)
            imp_syst_err = pow(imp_syst_err/float(n_imposters-1),0.5)
            if imp_syst_err>0.:
                zscore = (cnt_on_data-cnt_on_bkgd)/imp_syst_err
                if zscore>5.:
                    target_psr_detect[psr] = True

    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = CommonPlotFunctions.ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
    target_snr_name, target_snr_ra, target_snr_dec = CommonPlotFunctions.ReadSNRTargetListFromCSVFile()
    target_tev_name = []
    target_tev_ra = []
    target_tev_dec = []
    inputFile = open('TeVCat_RaDec_w_Names.txt')
    for line in inputFile:
        target_tev_name += [line.split(',')[0]]
        target_tev_ra += [float(line.split(',')[1])]
        target_tev_dec += [float(line.split(',')[2])]
    short_list_psr_name = []
    short_list_psr_ra = []
    short_list_psr_dec = []
    for psr in range(0,len(target_psr_name)):
        if target_psr_detect[psr]:
            print ('PSR %s is detected. Dist = %0.2f kpc. Age = %0.2f kyr. RA = %0.2f Dec = %0.2f'%(target_psr_name[psr],target_psr_dist[psr],target_psr_age[psr]/1000.,target_psr_ra[psr],target_psr_dec[psr]))
            short_list_psr_name += [target_psr_name[psr]]
            short_list_psr_ra += [target_psr_ra[psr]]
            short_list_psr_dec += [target_psr_dec[psr]]
            for snr in range(0,len(target_snr_name)):
                dist_psr_snr = pow(pow(target_psr_ra[psr]-target_snr_ra[snr],2)+pow(target_psr_dec[psr]-target_snr_dec[snr],2),0.5)
                if dist_psr_snr<0.5:
                    print ('SNR %s might be associated.'%(target_snr_name[snr]))
            for tev in range(0,len(target_tev_name)):
                dist_psr_tev = pow(pow(target_psr_ra[psr]-target_tev_ra[tev],2)+pow(target_psr_dec[psr]-target_tev_dec[tev],2),0.5)
                if dist_psr_tev<0.5:
                    print ('TeV source %s might be associated.'%(target_tev_name[tev]))
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    return short_list_psr_name, short_list_psr_ra, short_list_psr_dec

def GetMapNormalDistribution(hist_mean):

    zscores = []
    hist_mean_skymap = ROOT.TH2D("hist_mean_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    hist_noise_skymap = ROOT.TH2D("hist_noise_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    hist_mean_skymap.Add(hist_mean)

    n_trials = 10
    for trial in range(0,n_trials):
        hist_noise_skymap.Reset()
        for binx in range(0,hist_noise_skymap.GetNbinsX()):
            for biny in range(0,hist_noise_skymap.GetNbinsY()):
                error = pow(hist_mean_skymap.GetBinContent(binx+1,biny+1),0.5)
                if error<1.: continue
                random_number = np.random.normal(loc=0, scale=error)
                hist_noise_skymap.SetBinContent(binx+1,biny+1,random_number)
        if CommonPlotFunctions.Smoothing:
            hist_mean_skymap = CommonPlotFunctions.Smooth2DMap(hist_mean_skymap,CommonPlotFunctions.smooth_size_spectroscopy,False,False)
            hist_noise_skymap = CommonPlotFunctions.Smooth2DMap(hist_noise_skymap,CommonPlotFunctions.smooth_size_spectroscopy,False,False)
        for binx in range(0,hist_noise_skymap.GetNbinsX()):
            for biny in range(0,hist_noise_skymap.GetNbinsY()):
                error = pow(hist_mean_skymap.GetBinContent(binx+1,biny+1),0.5)
                if error<1.: continue
                content = hist_noise_skymap.GetBinContent(binx+1,biny+1)/error
                zscores += [content]

    return zscores

def GetMapChi2Distribution(hist_z,hist_cnt):

    count_max = hist_cnt.GetMaximum()
    zscores = []
    for binx in range(0,hist_z.GetNbinsX()):
        for biny in range(0,hist_z.GetNbinsY()):
            content = hist_z.GetBinContent(binx+1,biny+1)
            count = hist_cnt.GetBinContent(binx+1,biny+1)
            if count/count_max<0.3: continue
            zscores += [content]
    return zscores

def MakeSignificanceMap(hist_on_data_skymap,hist_on_bkgd_skymap,hist_mimic_data_skymap,hist_mimic_bkgd_skymap):

    hist_total_err_skymap = []
    hist_syst_err_skymap = []
    hist_avg_bkgd_skymap = []
    hist_zscore_skymap = []
    hist_excess_skymap = []
    for ebin in range(0,len(energy_bin)-1):
        hist_total_err_skymap += [ROOT.TH2D("hist_total_err_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_syst_err_skymap += [ROOT.TH2D("hist_syst_err_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_avg_bkgd_skymap += [ROOT.TH2D("hist_avg_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_zscore_skymap += [ROOT.TH2D("hist_zscore_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_excess_skymap += [ROOT.TH2D("hist_excess_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for binx in range(0,hist_mimic_data_skymap[0][ebin].GetNbinsX()):
            for biny in range(0,hist_mimic_data_skymap[0][ebin].GetNbinsY()):
                sum_square = 0.
                sum_bkgd = 0.
                for imposter in range(0,n_imposters):
                    imposter_data = hist_mimic_data_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    imposter_bkgd = hist_mimic_bkgd_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    sum_square += pow(imposter_data-imposter_bkgd,2)
                    sum_bkgd += imposter_bkgd
                total_error = pow(sum_square/float(n_imposters-1),0.5)
                avg_bkgd = sum_bkgd/float(n_imposters)
                #print ('ebin %s, binx %s, biny %s, avg_bkgd %s '%(ebin,binx,biny,avg_bkgd))
                syst_error = total_error
                if avg_bkgd==0.: continue
                hist_total_err_skymap[ebin].SetBinContent(binx+1,biny+1,total_error)
                hist_syst_err_skymap[ebin].SetBinContent(binx+1,biny+1,syst_error)
                hist_avg_bkgd_skymap[ebin].SetBinContent(binx+1,biny+1,avg_bkgd)

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        skymap = GetSignificanceMap(hist_on_data_skymap[ebin],hist_on_bkgd_skymap[ebin],hist_syst_err_skymap[ebin],True)
        hist_zscore_skymap[ebin].Add(skymap)
        hist_excess_skymap[ebin].Add(hist_on_data_skymap[ebin])
        hist_excess_skymap[ebin].Add(hist_on_bkgd_skymap[ebin],-1.)

    hist_zscore_skymap_sum = ROOT.TH2D("hist_zscore_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)

    hist_zscore_skymap_sum.Reset()
    hist_real_data_skymap_sum.Reset()
    hist_real_bkgd_skymap_sum.Reset()
    hist_real_syst_skymap_sum.Reset()
    for ebin in range(max(0,energy_bin_cut_low),plot_energy_break_bin):
        hist_real_data_skymap_sum.Add(hist_on_data_skymap[ebin])
        hist_real_bkgd_skymap_sum.Add(hist_on_bkgd_skymap[ebin])
        for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
            for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
                syst_err = hist_syst_err_skymap[ebin].GetBinContent(binx+1,biny+1)
                old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,old_syst_err+syst_err*syst_err)
    for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
        for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
            old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
            hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,pow(old_syst_err,0.5))
    hist_zscore_skymap_sum.Add(GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum,hist_real_syst_skymap_sum,True))

    hist_real_zscore_skymap_le.Reset()
    hist_real_zscore_skymap_le.Add(hist_zscore_skymap_sum)

    hist_zscore_skymap_sum.Reset()
    hist_real_data_skymap_sum.Reset()
    hist_real_bkgd_skymap_sum.Reset()
    hist_real_syst_skymap_sum.Reset()
    for ebin in range(max(plot_energy_break_bin,energy_bin_cut_low),6):
        hist_real_data_skymap_sum.Add(hist_on_data_skymap[ebin])
        hist_real_bkgd_skymap_sum.Add(hist_on_bkgd_skymap[ebin])
        for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
            for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
                syst_err = hist_syst_err_skymap[ebin].GetBinContent(binx+1,biny+1)
                old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,old_syst_err+syst_err*syst_err)
    for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
        for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
            old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
            hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,pow(old_syst_err,0.5))
    hist_zscore_skymap_sum.Add(GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum,hist_real_syst_skymap_sum,True))

    hist_real_zscore_skymap_he.Reset()
    hist_real_zscore_skymap_he.Add(hist_zscore_skymap_sum)

    hist_zscore_skymap_sum.Reset()
    hist_real_data_skymap_sum.Reset()
    hist_real_bkgd_skymap_sum.Reset()
    hist_real_syst_skymap_sum.Reset()
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_real_data_skymap_sum.Add(hist_on_data_skymap[ebin])
        hist_real_bkgd_skymap_sum.Add(hist_on_bkgd_skymap[ebin])
        for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
            for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
                syst_err = hist_syst_err_skymap[ebin].GetBinContent(binx+1,biny+1)
                old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,old_syst_err+syst_err*syst_err)
    for binx in range(0,hist_real_syst_skymap_sum.GetNbinsX()):
        for biny in range(0,hist_real_syst_skymap_sum.GetNbinsY()):
            old_syst_err = hist_real_syst_skymap_sum.GetBinContent(binx+1,biny+1)
            hist_real_syst_skymap_sum.SetBinContent(binx+1,biny+1,pow(old_syst_err,0.5))
    hist_zscore_skymap_sum.Add(GetSignificanceMap(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum,hist_real_syst_skymap_sum,True))

    hist_real_zscore_skymap_sum.Reset()
    hist_real_zscore_skymap_sum.Add(hist_zscore_skymap_sum)

    hist_mimic_data_skymap_sum = []
    hist_mimic_bkgd_skymap_sum = []
    for imposter in range(0,n_imposters):
        hist_mimic_data_skymap_sum += [ROOT.TH2D("hist_mimic_data_skymap_sum_%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_mimic_bkgd_skymap_sum += [ROOT.TH2D("hist_mimic_bkgd_skymap_sum_%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            hist_mimic_data_skymap_sum[imposter].Add(hist_mimic_data_skymap[imposter][ebin])
            hist_mimic_bkgd_skymap_sum[imposter].Add(hist_mimic_bkgd_skymap[imposter][ebin])
    data_cnt, data_cnt_err, bkgd_cnt, bkgd_cnt_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v1(hist_real_data_skymap_sum,hist_real_bkgd_skymap_sum,region_x,region_y,2.0*region_r[2])
    mimic_data_cnt_list = []
    mimic_bkgd_cnt_list = []
    for imposter in range(0,n_imposters):
        mimic_data_cnt, mimic_data_cnt_err, mimic_bkgd_cnt, mimic_bkgd_cnt_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v1(hist_mimic_data_skymap_sum[imposter],hist_mimic_bkgd_skymap_sum[imposter],region_x,region_y,2.0*region_r[2])
        mimic_data_cnt_list += [mimic_data_cnt]
        mimic_bkgd_cnt_list += [mimic_bkgd_cnt]
    syst_cnt = []
    for u in range(0,len(data_cnt)):
        mimic_square = 0.
        for imposter in range(0,n_imposters):
            mimic_square += pow(mimic_data_cnt_list[imposter][u]-mimic_bkgd_cnt_list[imposter][u],2)
        mimic_square = pow(mimic_square/float(n_imposters-1),0.5)
        syst_cnt += [mimic_square]

    theta2 = np.array(theta2)
    theta2_err = np.array(theta2_err)
    syst_cnt = np.array(syst_cnt)
    data_cnt = np.array(data_cnt)
    bkgd_cnt = np.array(bkgd_cnt)
    data_cnt_err = np.array(data_cnt_err)
    bkgd_cnt_err = np.array(bkgd_cnt_err)
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.bar(theta2, 2.*syst_cnt, bottom=bkgd_cnt-syst_cnt, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.errorbar(theta2,data_cnt,data_cnt_err,color='k',marker='s',ls='none',label='Data')
    axbig.plot(theta2,bkgd_cnt,color='r',ls='solid',label='Background')
    axbig.set_ylabel('event count per degree square')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'CountVsTheta2_%s'%(region_name)
    fig.savefig("output_plots/%s_sum_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()


    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap_sum)
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_zscore_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Z score','SkymapZscore_E%s_%s'%(ebin,plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
        hist_excess_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_excess_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_excess_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','excess count','SkymapExcess_E%s_%s'%(ebin,plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
        hist_on_bkgd_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_on_bkgd_skymap[ebin])
        CommonPlotFunctions.MatplotlibMap2D(hist_on_bkgd_skymap_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','excess count','SkymapBkgdCnt_E%s_%s'%(ebin,plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_sum_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Z score','SkymapZscore_Sum_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    CommonPlotFunctions.SaveAsFITS(hist_zscore_skymap_sum_reflect,'skymap_zscore_sum_%s'%(plot_tag))
    hist_zscore_skymap_sum_reflect.Delete()

    hist_imposter_zscore_skymap_sum = []
    for imposter in range(0,n_imposters):
        hist_imposter_zscore_skymap_sum += [ROOT.TH2D("hist_imposter_zscore_skymap_sum_%s"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_mimic_data_skymap_sum = ROOT.TH2D("hist_mimic_data_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    hist_mimic_bkgd_skymap_sum = ROOT.TH2D("hist_mimic_bkgd_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    for imposter in range(0,n_imposters):
        hist_mimic_data_skymap_sum.Reset()
        hist_mimic_bkgd_skymap_sum.Reset()
        hist_zscore_skymap_sum.Reset()
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            hist_syst_err_skymap[ebin].Reset()
            for binx in range(0,hist_mimic_data_skymap[0][ebin].GetNbinsX()):
                for biny in range(0,hist_mimic_data_skymap[0][ebin].GetNbinsY()):
                    sum_square = 0.
                    sum_bkgd = 0.
                    for imposter2 in range(0,n_imposters):
                        if imposter==imposter2: continue
                        imposter_data = hist_mimic_data_skymap[imposter2][ebin].GetBinContent(binx+1,biny+1)
                        imposter_bkgd = hist_mimic_bkgd_skymap[imposter2][ebin].GetBinContent(binx+1,biny+1)
                        sum_square += pow(imposter_data-imposter_bkgd,2)
                        sum_bkgd += imposter_bkgd
                    total_error = pow(sum_square/float(n_imposters-1-1),0.5)
                    avg_bkgd = sum_bkgd/float(n_imposters-1)
                    syst_error = total_error
                    if avg_bkgd==0.: continue
                    hist_syst_err_skymap[ebin].SetBinContent(binx+1,biny+1,syst_error)
        hist_mimic_syst_skymap_sum.Reset()
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            for binx in range(0,hist_mimic_syst_skymap_sum.GetNbinsX()):
                for biny in range(0,hist_mimic_syst_skymap_sum.GetNbinsY()):
                    syst_err = hist_syst_err_skymap[ebin].GetBinContent(binx+1,biny+1)
                    old_syst_err = hist_mimic_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                    hist_mimic_syst_skymap_sum.SetBinContent(binx+1,biny+1,old_syst_err+syst_err*syst_err)
        for binx in range(0,hist_mimic_syst_skymap_sum.GetNbinsX()):
            for biny in range(0,hist_mimic_syst_skymap_sum.GetNbinsY()):
                old_syst_err = hist_mimic_syst_skymap_sum.GetBinContent(binx+1,biny+1)
                hist_mimic_syst_skymap_sum.SetBinContent(binx+1,biny+1,pow(old_syst_err,0.5))
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            hist_mimic_data_skymap_sum.Add(hist_mimic_data_skymap[imposter][ebin])
            hist_mimic_bkgd_skymap_sum.Add(hist_mimic_bkgd_skymap[imposter][ebin])
        hist_zscore_skymap_sum.Add(GetSignificanceMap(hist_mimic_data_skymap_sum,hist_mimic_bkgd_skymap_sum,hist_mimic_syst_skymap_sum,False))
        hist_imposter_zscore_skymap_sum[imposter].Reset()
        hist_imposter_zscore_skymap_sum[imposter].Add(hist_zscore_skymap_sum)
        hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_zscore_skymap_sum)
        CommonPlotFunctions.MatplotlibMap2D(hist_zscore_skymap_sum_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Z score','SkymapZscore_Imposter%s_Sum_%s'%(imposter,plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
        hist_zscore_skymap_sum_reflect.Delete()


    fig.clf()
    fig.set_figheight(figsize_x)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    #zscores = GetMapChi2Distribution(hist_real_zscore_skymap_sum,hist_real_bkgd_skymap_sum)
    #z_max_range = max(4.,max(zscores))
    z_max_range = 6.
    z_min_range = -3.
    zscore_normal = GetMapNormalDistribution(hist_real_bkgd_skymap_sum)
    axbig.hist(zscore_normal, bins=40, range=[z_min_range, z_max_range], density=True, color='gray', alpha=0.5, label='Normal')
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        zscores = GetMapChi2Distribution(hist_imposter_zscore_skymap_sum[imposter],hist_imposter_bkgd_skymap_sum[imposter])
        axbig.hist(zscores, bins=40, range=[z_min_range, z_max_range], density=True, histtype='step', color=next_color, alpha=1.0, label='Mimic #%s'%(imposter+1))
    zscores = GetMapChi2Distribution(hist_real_zscore_skymap_sum,hist_real_bkgd_skymap_sum)
    axbig.hist(zscores, bins=40, range=[z_min_range, z_max_range], density=True, histtype='step', color='black', alpha=1.0, label='ON data')
    axbig.set_xlabel('Z scores')
    axbig.set_ylabel('Entries')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'Chi2Distribution'
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

def diffusion_func(x,A,d):
    return A*1.22/(pow(3.14,1.5)*d*(x+0.06*d))*np.exp(-x*x/(d*d))

def gauss_func(x,A,d):
    return A*np.exp(-x*x/(2.*d*d))

def power_law_func(x,a,b):
    return a*pow(x*1./1000.,b)

def flux_crab_func(x):
    # TeV^{-1}cm^{-2}s^{-1}
    # Crab https://arxiv.org/pdf/1508.06442.pdf
    return 37.5*pow(10,-12)*pow(x*1./1000.,-2.467-0.16*log(x/1000.))

def flux_1es0229_func(x):
    # TeV^{-1}cm^{-2}s^{-1}
    # https://articles.adsabs.harvard.edu/pdf/2013ICRC...33.1003C
    return 5.54*pow(10,-13)*pow(x*1./1000.,-2.59)

def flux_mgro_j2019_func(x):
    # TeV^{-1}cm^{-2}s^{-1}
    # https://arxiv.org/pdf/1404.1841.pdf
    return 8.1*pow(10,-14)*pow(x*1./5000.,-1.75)

def flux_hawc_j1908_func(x):
    # MGRO J1908  TeV^{-1}cm^{-2}s^{-1}
    return 0.95*pow(10,-13)*pow(x*1./10000.,-2.46-0.11*log(x/10000.))

def GetHessFlux1ES0229(energy_index):
    # https://arxiv.org/pdf/0709.4584.pdf
    energies = [-0.221,-0.055,0.110,0.277,0.444,0.611,0.812,1.062]
    fluxes = [-11.487,-12.029,-12.632, -12.786, -13.132, -13.606, -14.418, -14.846]
    flux_errs = [-11.487,-12.029,-12.632, -12.786, -13.132, -13.606, -14.418, -14.846]
    flux_errs_up = [-11.341, -11.914, -12.470, -12.675, -13.008, -13.444, -14.153, -14.602]
    flux_errs_low = [-11.726, -12.196, -12.897, -12.948, -13.299, -13.884, -15.226, -15.401]

    for entry in range(0,len(energies)):
        energies[entry] = pow(10.,energies[entry])
        fluxes[entry] = pow(10.,fluxes[entry])*pow(energies[entry],energy_index)
        flux_errs_up[entry] = pow(10.,flux_errs_up[entry])*pow(energies[entry],energy_index)
        flux_errs_low[entry] = pow(10.,flux_errs_low[entry])*pow(energies[entry],energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])
        energies[entry] = energies[entry]*1e3

    return energies, fluxes, flux_errs

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

def GetHawcDiffusionFluxJ1908(energy_index):
    energies = [pow(10.,0.075),pow(10.,0.259),pow(10.,0.492),pow(10.,0.740),pow(10.,0.997),pow(10.,1.270),pow(10.,1.532),pow(10.,1.776),pow(10.,2.011),pow(10.,2.243)]
    fluxes = [pow(10.,-10.710),pow(10.,-10.696),pow(10.,-10.699),pow(10.,-10.801),pow(10.,-10.924),pow(10.,-11.137),pow(10.,-11.329),pow(10.,-11.560),pow(10.,-11.669),pow(10.,-11.860)]
    flux_errs = [pow(10.,-10.60),pow(10.,-10.62),pow(10.,-10.64),pow(10.,-10.74),pow(10.,-10.86),pow(10.,-11.03),pow(10.,-11.18),pow(10.,-11.36),pow(10.,-11.37),pow(10.,-11.55)]
    flux_errs_up = [pow(10.,-10.670),pow(10.,-10.656),pow(10.,-10.660),pow(10.,-10.764),pow(10.,-10.884),pow(10.,-11.101),pow(10.,-11.289),pow(10.,-11.491),pow(10.,-11.585),pow(10.,-11.712)]
    flux_errs_low = [pow(10.,-10.750),pow(10.,-10.739),pow(10.,-10.739),pow(10.,-10.844),pow(10.,-10.967),pow(10.,-11.180),pow(10.,-11.372),pow(10.,-11.636),pow(10.,-11.777),pow(10.,-12.095)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHessFluxJ1908(energy_index):
    energies = [pow(10.,-0.332),pow(10.,0.022),pow(10.,0.396),pow(10.,0.769),pow(10.,1.124),pow(10.,1.478)]
    fluxes = [pow(10.,-10.981),pow(10.,-10.967),pow(10.,-11.057),pow(10.,-11.169),pow(10.,-11.188),pow(10.,-11.386)]
    flux_errs = [pow(10.,-0.332),pow(10.,0.022),pow(10.,0.396),pow(10.,0.769),pow(10.,1.124),pow(10.,1.478)]
    flux_errs_up = [pow(10.,-10.895),pow(10.,-10.916),pow(10.,-11.003),pow(10.,-11.101),pow(10.,-11.101),pow(10.,-11.264)]
    flux_errs_low = [pow(10.,-11.086),pow(10.,-11.010),pow(10.,-11.126),pow(10.,-11.264),pow(10.,-11.292),pow(10.,-11.556)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetVeritasFluxJ1908(energy_index):
    energies = [pow(10.,-0.270),pow(10.,-0.126),pow(10.,0.022),pow(10.,0.175),pow(10.,0.323),pow(10.,0.467),pow(10.,0.618),pow(10.,0.776),pow(10.,0.922),pow(10.,1.070),pow(10.,1.219)]
    fluxes = [pow(10.,-11.061),pow(10.,-11.028),pow(10.,-11.036),pow(10.,-11.097),pow(10.,-11.448),pow(10.,-11.166),pow(10.,-11.213),pow(10.,-11.068),pow(10.,-11.209),pow(10.,-11.231),pow(10.,-11.318)]
    flux_errs = [pow(10.,-0.270),pow(10.,-0.126),pow(10.,0.022),pow(10.,0.175),pow(10.,0.323),pow(10.,0.467),pow(10.,0.618),pow(10.,0.776),pow(10.,0.922),pow(10.,1.070),pow(10.,1.219)]
    flux_errs_up = [pow(10.,-10.952),pow(10.,-10.960),pow(10.,-10.974),pow(10.,-11.028),pow(10.,-11.303),pow(10.,-11.083),pow(10.,-11.112),pow(10.,-10.974),pow(10.,-11.083),pow(10.,-11.097),pow(10.,-11.141)]
    flux_errs_low = [pow(10.,-11.245),pow(10.,-11.155),pow(10.,-11.137),pow(10.,-11.195),pow(10.,-11.661),pow(10.,-11.282),pow(10.,-11.339),pow(10.,-11.162),pow(10.,-11.372),pow(10.,-11.408),pow(10.,-11.538)]

    scale_factor = 1./2.03
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*scale_factor/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*scale_factor/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetLHAASOFluxJ1908(energy_index):
    energies = [pow(10.,1.102),pow(10.,1.302),pow(10.,1.498),pow(10.,1.700),pow(10.,1.900),pow(10.,2.099),pow(10.,2.299),pow(10.,2.498),pow(10.,2.697)]
    fluxes = [pow(10.,-11.033),pow(10.,-10.988),pow(10.,-11.201),pow(10.,-11.324),pow(10.,-11.553),pow(10.,-11.860),pow(10.,-11.921),pow(10.,-12.346),pow(10.,-12.653)]
    flux_errs = [pow(10.,1.102),pow(10.,1.302),pow(10.,1.498),pow(10.,1.700),pow(10.,1.900),pow(10.,2.099),pow(10.,2.299),pow(10.,2.498),pow(10.,2.697)]
    flux_errs_up = [pow(10.,-10.966),pow(10.,-10.949),pow(10.,-11.167),pow(10.,-11.296),pow(10.,-11.513),pow(10.,-11.798),pow(10.,-11.854),pow(10.,-12.173),pow(10.,-12.391)]
    flux_errs_low = [pow(10.,-11.094),pow(10.,-11.027),pow(10.,-11.240),pow(10.,-11.368),pow(10.,-11.597),pow(10.,-11.944),pow(10.,-12.022),pow(10.,-12.536),pow(10.,-13.128)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetLHAASOFluxJ2226(energy_index):
    energies = [pow(10.,1.302),pow(10.,1.5),pow(10.,1.700),pow(10.,1.901),pow(10.,2.102),pow(10.,2.299),pow(10.,2.5),pow(10.,2.70)]
    fluxes = [pow(10.,-11.391),pow(10.,-11.486),pow(10.,-11.525),pow(10.,-11.849),pow(10.,-11.916),pow(10.,-12.463),pow(10.,-12.888),pow(10.,-13.156)]
    flux_errs = [pow(10.,1.302),pow(10.,1.5),pow(10.,1.700),pow(10.,1.901),pow(10.,2.102),pow(10.,2.299),pow(10.,2.5),pow(10.,2.70)]
    flux_errs_up = [pow(10.,-11.312),pow(10.,-11.424),pow(10.,-11.480),pow(10.,-11.793),pow(10.,-11.860),pow(10.,-12.318),pow(10.,-12.608),pow(10.,-12.687)]
    flux_errs_low = [pow(10.,-11.491),pow(10.,-11.558),pow(10.,-11.569),pow(10.,-11.916),pow(10.,-12),pow(10.,-12.648),pow(10.,-13.391),pow(10.,-13.653)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHAWCFluxJ2019(energy_index):
    #https://arxiv.org/pdf/1909.08609.pdf
    energies = [pow(10.,12.230),pow(10.,12.412),pow(10.,12.613),pow(10.,12.807),pow(10.,13.034),pow(10.,13.286),pow(10.,13.538),pow(10.,13.771),pow(10.,14.009)]
    fluxes = [pow(10.,-11.659),pow(10.,-11.389),pow(10.,-11.420),pow(10.,-11.345),pow(10.,-11.324),pow(10.,-11.358),pow(10.,-11.642),pow(10.,-11.752),pow(10.,-11.827)]
    flux_errs = [pow(10.,-11.659),pow(10.,-11.389),pow(10.,-11.420),pow(10.,-11.345),pow(10.,-11.324),pow(10.,-11.358),pow(10.,-11.642),pow(10.,-11.752),pow(10.,-11.827)]
    flux_errs_up = [pow(10.,-11.540),pow(10.,-11.338),pow(10.,-11.372),pow(10.,-11.314),pow(10.,-11.300),pow(10.,-11.335),pow(10.,-11.598),pow(10.,-11.690),pow(10.,-11.741)]
    flux_errs_low = [pow(10.,-11.820),pow(10.,-11.447),pow(10.,-11.482),pow(10.,-11.379),pow(10.,-11.348),pow(10.,-11.382),pow(10.,-11.694),pow(10.,-11.827),pow(10.,-11.936)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e9
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetVeritasFluxJ2019(energy_index):

    #https://arxiv.org/pdf/1805.05989.pdf
    # integration radius = 0.23 deg
    energies = [pow(10.,2.75),pow(10.,3.),pow(10.,3.25),pow(10.,3.5),pow(10.,3.75),pow(10.,4.0),pow(10.,4.25),pow(10.,4.5)]
    fluxes = [pow(10.,-12.291),pow(10.,-12.175),pow(10.,-11.980),pow(10.,-11.832),pow(10.,-11.868),pow(10.,-12.097),pow(10.,-11.982),pow(10.,-12.270)]
    flux_errs = [pow(10.,-12.291),pow(10.,-12.175),pow(10.,-11.980),pow(10.,-11.832),pow(10.,-11.868),pow(10.,-12.097),pow(10.,-11.982),pow(10.,-12.270)]
    flux_errs_up = [pow(10.,-11.982),pow(10.,-12.046),pow(10.,-11.895),pow(10.,-11.763),pow(10.,-11.771),pow(10.,-11.936),pow(10.,-11.812),pow(10.,-12.029)]
    flux_errs_low = [pow(10.,-13.299),pow(10.,-12.362),pow(10.,-12.087),pow(10.,-11.919),pow(10.,-11.995),pow(10.,-12.355),pow(10.,-12.270),pow(10.,-12.861)]

    for entry in range(0,len(energies)):
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxJ1856(energy_index):

    #https://arxiv.org/pdf/2209.11855.pdf, table 1
    energies = [535.,1705.,5432.,17303.,55116.,175560.,559204.]
    fluxes = [11.1, 10.1, 7, 7.55, 8.1, 6.36, 6.27]
    flux_errs = [2.53, 1.08, 1.28, 1.79, 2.55, 3.91, 6.38]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e3
        fluxes[entry] = fluxes[entry]*1e-12
        flux_errs[entry] = flux_errs[entry]*1e-12

    return energies, fluxes, flux_errs

def GetFermiFluxJ2019(energy_index):
    #https://arxiv.org/pdf/2209.11855.pdf, table 1
    energies = [535.,1705.,5432.,17303.,55116.,175560.]
    fluxes = [3.25, 2.02, 1.77, 1.24, 0.229, 3.38]
    flux_errs = [1.3, 0.511, 0.417, 0.452, 0.41, 1.97]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]/1e3
        fluxes[entry] = fluxes[entry]*1e-12
        flux_errs[entry] = flux_errs[entry]*1e-12

    return energies, fluxes, flux_errs

def GetTibetASFluxJ2227(energy_index):
    energies = [pow(10.,0.824),pow(10.,1.032),pow(10.,1.207),pow(10.,1.441),pow(10.,1.659),pow(10.,1.861),pow(10.,2.062)]
    fluxes = [pow(10.,-11.218),pow(10.,-11.206),pow(10.,-11.502),pow(10.,-11.549),pow(10.,-11.904),pow(10.,-11.968),pow(10.,-12.526)]
    flux_errs = [pow(10.,-11.218),pow(10.,-11.206),pow(10.,-11.502),pow(10.,-11.549),pow(10.,-11.904),pow(10.,-11.968),pow(10.,-12.526)]
    flux_errs_up = [pow(10.,-10.979),pow(10.,-11.049),pow(10.,-11.343),pow(10.,-11.444),pow(10.,-11.764),pow(10.,-11.834),pow(10.,-12.238)]
    flux_errs_low = [pow(10.,-11.761),pow(10.,-11.441),pow(10.,-11.764),pow(10.,-11.694),pow(10.,-12.084),pow(10.,-12.116),pow(10.,-12.933)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxJ2227(energy_index):
    energies = [pow(10.,-2.299),pow(10.,-1.851),pow(10.,-1.409),pow(10.,-0.961),pow(10.,-0.513)]
    fluxes = [pow(10.,-12.343),pow(10.,-12.142),pow(10.,-12.200),pow(10.,-11.959),pow(10.,-11.950)]
    flux_errs = [pow(10.,-12.343),pow(10.,-12.142),pow(10.,-12.200),pow(10.,-11.959),pow(10.,-11.950)]
    flux_errs_up = [pow(10.,-12.180),pow(10.,-12.020),pow(10.,-12.040),pow(10.,-11.779),pow(10.,-11.700)]
    flux_errs_low = [pow(10.,-12.619),pow(10.,-12.305),pow(10.,-12.470),pow(10.,-12.267),pow(10.,-12.598)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetVeritasFluxJ2227(energy_index):
    energies = [pow(10.,0.008),pow(10.,0.200),pow(10.,0.402),pow(10.,0.721),pow(10.,1.180)]
    fluxes = [pow(10.,-11.619),pow(10.,-11.688),pow(10.,-11.709),pow(10.,-11.950),pow(10.,-11.738)]
    flux_errs = [pow(10.,-12.343),pow(10.,-12.142),pow(10.,-12.200),pow(10.,-11.959),pow(10.,-11.950)]
    flux_errs_up = [pow(10.,-11.424),pow(10.,-11.505),pow(10.,-11.590),pow(10.,-11.744),pow(10.,-11.543)]
    flux_errs_low = [pow(10.,-11.889),pow(10.,-11.901),pow(10.,-11.927),pow(10.,-12.418),pow(10.,-12.078)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = (1./1.62)*fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = (1./1.62)*0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHAWCFluxJ1825(energy_index):
    energies = [pow(10.,0.149),pow(10.,0.307),pow(10.,0.473),pow(10.,0.715),pow(10.,0.963),pow(10.,1.172),pow(10.,1.436),pow(10.,1.694)]
    fluxes = [pow(10.,-10.799),pow(10.,-10.933),pow(10.,-10.829),pow(10.,-10.843),pow(10.,-11.040),pow(10.,-11.043),pow(10.,-11.323),pow(10.,-11.601)]
    flux_errs = [pow(10.,0.149),pow(10.,0.307),pow(10.,0.473),pow(10.,0.715),pow(10.,0.963),pow(10.,1.172),pow(10.,1.436),pow(10.,1.694)]
    flux_errs_up = [pow(10.,-10.679),pow(10.,-10.873),pow(10.,-10.779),pow(10.,-10.796),pow(10.,-10.989),pow(10.,-10.993),pow(10.,-11.260),pow(10.,-11.484)]
    flux_errs_low = [pow(10.,-10.973),pow(10.,-11.003),pow(10.,-10.883),pow(10.,-10.893),pow(10.,-11.096),pow(10.,-11.093),pow(10.,-11.393),pow(10.,-11.747)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetLHAASOFluxJ1825(energy_index):
    energies = [pow(10.,1.5),pow(10.,1.70),pow(10.,1.901),pow(10.,2.102),pow(10.,2.302),pow(10.,2.500),pow(10.,2.704)]
    fluxes = [pow(10.,-10.776),pow(10.,-11.156),pow(10.,-11.256),pow(10.,-11.418),pow(10.,-12.016),pow(10.,-12.424),pow(10.,-12.715)]
    flux_errs = [pow(10.,1.5),pow(10.,1.70),pow(10.,1.901),pow(10.,2.102),pow(10.,2.302),pow(10.,2.500),pow(10.,2.704)]
    flux_errs_up = [pow(10.,-10.659),pow(10.,-11.089),pow(10.,-11.206),pow(10.,-11.357),pow(10.,-11.865),pow(10.,-12.173),pow(10.,-12.262)]
    flux_errs_low = [pow(10.,-10.949),pow(10.,-11.268),pow(10.,-11.318),pow(10.,-11.491),pow(10.,-12.206),pow(10.,-12.882),pow(10.,-13.162)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHESSFluxR0p8J1825(energy_index):
    energies = [pow(10.,-0.573),pow(10.,-0.414),pow(10.,-0.252),pow(10.,-0.084),pow(10.,0.080),pow(10.,0.248),pow(10.,0.413),pow(10.,0.578),pow(10.,0.746),pow(10.,0.914),pow(10.,1.079),pow(10.,1.247),pow(10.,1.412),pow(10.,1.58),pow(10.,1.733)]
    fluxes = [pow(10.,-10.452),pow(10.,-10.4),pow(10.,-10.5),pow(10.,-10.488),pow(10.,-10.544),pow(10.,-10.592),pow(10.,-10.616),pow(10.,-10.72),pow(10.,-10.808),pow(10.,-10.808),pow(10.,-10.9),pow(10.,-11.168),pow(10.,-11.392),pow(10.,-11.64),pow(10.,-11.732)]
    flux_errs = [pow(10.,-0.573),pow(10.,-0.414),pow(10.,-0.252),pow(10.,-0.084),pow(10.,0.080),pow(10.,0.248),pow(10.,0.413),pow(10.,0.578),pow(10.,0.746),pow(10.,0.914),pow(10.,1.079),pow(10.,1.247),pow(10.,1.412),pow(10.,1.58),pow(10.,1.733)]
    flux_errs_up = [pow(10.,-10.424),pow(10.,-10.372),pow(10.,-10.472),pow(10.,-10.46),pow(10.,-10.516),pow(10.,-10.564),pow(10.,-10.588),pow(10.,-10.692),pow(10.,-10.780),pow(10.,-10.776),pow(10.,-10.86),pow(10.,-11.096),pow(10.,-11.272),pow(10.,-11.436),pow(10.,-11.444)]
    flux_errs_low = [pow(10.,-10.48),pow(10.,-10.428),pow(10.,-10.528),pow(10.,-10.52),pow(10.,-10.572),pow(10.,-10.624),pow(10.,-10.648),pow(10.,-10.748),pow(10.,-10.844),pow(10.,-10.84),pow(10.,-10.944),pow(10.,-11.252),pow(10.,-11.556),pow(10.,-12.02),pow(10.,-12.884)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetHESSFluxR0p4J1825(energy_index):
    energies = [pow(10.,-0.615),pow(10.,-0.456),pow(10.,-0.291),pow(10.,-0.126),pow(10.,0.038),pow(10.,0.206),pow(10.,0.370),pow(10.,0.535),pow(10.,0.703),pow(10.,0.871),pow(10.,1.040),pow(10.,1.204),pow(10.,1.373),pow(10.,1.537),pow(10.,1.715)]
    fluxes = [pow(10.,-10.836),pow(10.,-10.856),pow(10.,-10.904),pow(10.,-10.92),pow(10.,-10.928),pow(10.,-10.964),pow(10.,-11.02),pow(10.,-11.096),pow(10.,-11.128),pow(10.,-11.252),pow(10.,-11.328),pow(10.,-11.436),pow(10.,-11.704),pow(10.,-11.66),pow(10.,-11.924)]
    flux_errs = [pow(10.,-0.573),pow(10.,-0.414),pow(10.,-0.252),pow(10.,-0.084),pow(10.,0.080),pow(10.,0.248),pow(10.,0.413),pow(10.,0.578),pow(10.,0.746),pow(10.,0.914),pow(10.,1.079),pow(10.,1.247),pow(10.,1.412),pow(10.,1.58),pow(10.,1.733)]
    flux_errs_up = [pow(10.,-10.812),pow(10.,-10.832),pow(10.,-10.88),pow(10.,-10.896),pow(10.,-10.9),pow(10.,-10.936),pow(10.,-10.996),pow(10.,-11.072),pow(10.,-11.1),pow(10.,-11.22),pow(10.,-11.288),pow(10.,-11.384),pow(10.,-11.608),pow(10.,-11.56),pow(10.,-11.76)]
    flux_errs_low = [pow(10.,-10.864),pow(10.,-10.884),pow(10.,-10.932),pow(10.,-10.948),pow(10.,-10.952),pow(10.,-10.988),pow(10.,-11.048),pow(10.,-11.124),pow(10.,-11.152),pow(10.,-11.284),pow(10.,-11.372),pow(10.,-11.496),pow(10.,-11.824),pow(10.,-11.775),pow(10.,-12.171)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetFermiFluxJ1825(energy_index):
    energies = [pow(10.,-2.948),pow(10.,-2.850),pow(10.,-2.748),pow(10.,-2.649),pow(10.,-2.550),pow(10.,-2.448),pow(10.,-2.349),pow(10.,-2.250),pow(10.,-2.148),pow(10.,-2.049),pow(10.,-1.950),pow(10.,-1.848),pow(10.,-1.749),pow(10.,-1.650),pow(10.,-1.551),pow(10.,-1.449),pow(10.,-1.347),pow(10.,-1.248),pow(10.,-1.149),pow(10.,-1.051),pow(10.,-0.948),pow(10.,-0.850),pow(10.,-0.748),pow(10.,-0.649)]
    fluxes = [pow(10.,-10.614),pow(10.,-10.729),pow(10.,-10.666),pow(10.,-10.576),pow(10.,-10.524),pow(10.,-10.437),pow(10.,-10.447),pow(10.,-10.371),pow(10.,-10.385),pow(10.,-10.364),pow(10.,-10.298),pow(10.,-10.215),pow(10.,-10.201),pow(10.,-10.236),pow(10.,-10.131),pow(10.,-10.218),pow(10.,-10.034),pow(10.,-10.114),pow(10.,-10.086),pow(10.,-10.270),pow(10.,-10.371),pow(10.,-10.222),pow(10.,-10.243),pow(10.,-10.038)]
    flux_errs = [pow(10.,-2.948),pow(10.,-2.850),pow(10.,-2.748),pow(10.,-2.649),pow(10.,-2.550),pow(10.,-2.448),pow(10.,-2.349),pow(10.,-2.250),pow(10.,-2.148),pow(10.,-2.049),pow(10.,-1.950),pow(10.,-1.848),pow(10.,-1.749),pow(10.,-1.650),pow(10.,-1.551),pow(10.,-1.449),pow(10.,-1.347),pow(10.,-1.248),pow(10.,-1.149),pow(10.,-1.051),pow(10.,-0.948),pow(10.,-0.850),pow(10.,-0.748),pow(10.,-0.649)]
    flux_errs_up = [pow(10.,-10.555),pow(10.,-10.656),pow(10.,-10.597),pow(10.,-10.496),pow(10.,-10.440),pow(10.,-10.392),pow(10.,-10.392),pow(10.,-10.298),pow(10.,-10.312),pow(10.,-10.312),pow(10.,-10.250),pow(10.,-10.145),pow(10.,-10.128),pow(10.,-10.177),pow(10.,-10.072),pow(10.,-10.114),pow(10.,-9.947),pow(10.,-10.041),pow(10.,-10.006),pow(10.,-10.152),pow(10.,-10.225),pow(10.,-10.086),pow(10.,-10.083),pow(10.,-9.881)]
    flux_errs_low = [pow(10.,-10.684),pow(10.,-10.819),pow(10.,-10.750),pow(10.,-10.687),pow(10.,-10.628),pow(10.,-10.500),pow(10.,-10.510),pow(10.,-10.461),pow(10.,-10.479),pow(10.,-10.430),pow(10.,-10.361),pow(10.,-10.302),pow(10.,-10.291),pow(10.,-10.312),pow(10.,-10.201),pow(10.,-10.354),pow(10.,-10.145),pow(10.,-10.211),pow(10.,-10.194),pow(10.,-10.434),pow(10.,-10.604),pow(10.,-10.430),pow(10.,-10.510),pow(10.,-10.288)]

    erg_to_TeV = 0.62
    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])*erg_to_TeV/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

    return energies, fluxes, flux_errs

def GetLHAASOFluxJ2108(energy_index):
    energies = [pow(10.,1.300),pow(10.,1.500),pow(10.,1.699),pow(10.,1.901),pow(10.,2.100),pow(10.,2.301)]
    fluxes = [pow(10.,-12.318),pow(10.,-12.357),pow(10.,-12.515),pow(10.,-12.665),pow(10.,-12.685),pow(10.,-13.521)]
    flux_errs = [pow(10.,1.300),pow(10.,1.500),pow(10.,1.699),pow(10.,1.901),pow(10.,2.100),pow(10.,2.301)]
    flux_errs_up = [pow(10.,-12.123),pow(10.,-12.252),pow(10.,-12.420),pow(10.,-12.560),pow(10.,-12.561),pow(10.,-13.132)]
    flux_errs_low = [pow(10.,-12.658),pow(10.,-12.495),pow(10.,-12.639),pow(10.,-12.803),pow(10.,-12.831),pow(10.,-14.030)]

    for entry in range(0,len(energies)):
        energies[entry] = energies[entry]*1e3
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry]/1e6)*pow(energies[entry]/1e3,energy_index)

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

def GetHAWCDiffusionFluxGeminga(energy_index):
    energies = [pow(10.,0.90),pow(10.,1.60)]
    fluxes = [pow(10.,-11.12),pow(10.,-11.36)]
    flux_errs = [0.,0.]
    flux_errs_up = [pow(10.,-11.04),pow(10.,-11.28)]
    flux_errs_low = [pow(10.,-11.21),pow(10.,-11.44)]

    for entry in range(0,len(energies)):
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        energies[entry] = energies[entry]*1e3

    return energies, fluxes, flux_errs

def GetHAWCGaussianFluxGeminga(energy_index):
    energies = [pow(10.,0.90),pow(10.,1.60)]
    fluxes = [pow(10.,-11.36),pow(10.,-11.52)]
    flux_errs = [0.,0.]
    flux_errs_up = [pow(10.,-11.28),pow(10.,-11.45)]
    flux_errs_low = [pow(10.,-11.44),pow(10.,-11.59)]

    for entry in range(0,len(energies)):
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        energies[entry] = energies[entry]*1e3

    return energies, fluxes, flux_errs

def GetHAWCDiskFluxGeminga(energy_index):
    energies = [pow(10.,0.00),pow(10.,1.70)]
    fluxes = [pow(10.,-11.42),pow(10.,-11.81)]
    flux_errs = [0.,0.]
    flux_errs_up = [pow(10.,-11.30),pow(10.,-11.68)]
    flux_errs_low = [pow(10.,-11.56),pow(10.,-11.94)]

    for entry in range(0,len(energies)):
        fluxes[entry] = fluxes[entry]/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        flux_errs[entry] = 0.5*(flux_errs_up[entry]-flux_errs_low[entry])/(energies[entry]*energies[entry])*pow(energies[entry],energy_index)
        energies[entry] = energies[entry]*1e3

    return energies, fluxes, flux_errs

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

def MaskKnownSources():

    theta_cut_min = 0.4

    target_name = []
    target_ra = []
    target_dec = []
    target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age = CommonPlotFunctions.ReadATNFTargetListFromFile('ATNF_pulsar_list.txt')
    target_name += target_psr_name
    target_ra += target_psr_ra
    target_dec += target_psr_dec
    target_snr_name, target_snr_ra, target_snr_dec = CommonPlotFunctions.ReadSNRTargetListFromCSVFile()
    target_name += target_snr_name
    target_ra += target_snr_ra
    target_dec += target_snr_dec

    if CommonPlotFunctions.doGalacticCoord:
        for src in range(0,len(target_name)):
            target_ra[src], target_dec[src] = ConvertRaDecToGalactic(target_ra[src], target_dec[src])

    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        for binx in range(0,hist_real_flux_skymap[ebin].GetNbinsX()):
            for biny in range(0,hist_real_flux_skymap[ebin].GetNbinsY()):
                bin_center_x = hist_real_flux_skymap[ebin].GetXaxis().GetBinCenter(binx+1)
                bin_center_y = hist_real_flux_skymap[ebin].GetYaxis().GetBinCenter(biny+1)
                for src in range(0,len(target_name)):
                    dist_to_src_x = bin_center_x-target_ra[src]
                    dist_to_src_y = bin_center_y-target_dec[src]
                    dist_to_src = pow(dist_to_src_x*dist_to_src_x+dist_to_src_y*dist_to_src_y,0.5)
                    if dist_to_src<theta_cut_min:
                        hist_real_flux_skymap[ebin].SetBinContent(binx+1,biny+1,0.)
                        hist_real_flux_skymap[ebin].SetBinError(binx+1,biny+1,0.)
                        hist_real_data_skymap[ebin].SetBinContent(binx+1,biny+1,0.)
                        hist_real_data_skymap[ebin].SetBinError(binx+1,biny+1,0.)
                        hist_real_bkgd_skymap[ebin].SetBinContent(binx+1,biny+1,0.)
                        hist_real_bkgd_skymap[ebin].SetBinError(binx+1,biny+1,0.)

def CleanFluxMapNoise(hist_skymap):

    mean_flux_err = 0.
    mean_flux = 0.
    rms_flux = 0.
    local_flux = 0.
    local_flux_err = 0.
    used_bins = 0.
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        mean_flux_err = 0.
        mean_flux = 0.
        rms_flux = 0.
        used_bins = 0.
        for binx in range(0,hist_skymap[ebin].GetNbinsX()):
            for biny in range(0,hist_skymap[ebin].GetNbinsY()):
                bin_content = hist_skymap[ebin].GetBinContent(binx+1,biny+1)
                bin_error = hist_skymap[ebin].GetBinError(binx+1,biny+1)
                if bin_content!=0.:
                    used_bins += 1.
                    mean_flux += bin_content
                    mean_flux_err += bin_error
        if used_bins>0.:
            mean_flux = mean_flux/used_bins
            mean_flux_err = mean_flux_err/used_bins
        used_bins = 0.
        for binx in range(0,hist_skymap[ebin].GetNbinsX()):
            for biny in range(0,hist_skymap[ebin].GetNbinsY()):
                bin_content = hist_skymap[ebin].GetBinContent(binx+1,biny+1)
                if bin_content!=0.:
                    used_bins += 1.
                    rms_flux += pow(bin_content-mean_flux,2)
        if used_bins>0.:
            rms_flux = pow(rms_flux/used_bins,0.5)

        for binx in range(0,hist_skymap[ebin].GetNbinsX()):
            for biny in range(0,hist_skymap[ebin].GetNbinsY()):
                if rms_flux==0.: continue
                local_flux = 0.
                local_flux_err = 0.
                local_flux = hist_skymap[ebin].GetBinContent(binx+1,biny+1)
                local_flux_stat_err = hist_skymap[ebin].GetBinError(binx+1,biny+1)
                local_flux_syst_err = 0.
                for imposter in range(0,n_imposters):
                    bin_content = hist_imposter_flux_skymap[imposter][ebin].GetBinContent(binx+1,biny+1) 
                    local_flux_syst_err += pow(bin_content,2)
                local_flux_syst_err = pow(local_flux_syst_err/n_imposters,0.5)
                local_flux_err = pow(local_flux_stat_err*local_flux_stat_err+local_flux_syst_err*local_flux_syst_err,0.5)
                if local_flux_err==0.: continue
                if abs(local_flux-mean_flux)/rms_flux>2.:
                    if abs(local_flux)/local_flux_err<2.:
                        hist_skymap[ebin].SetBinContent(binx+1,biny+1,0.)
                        hist_skymap[ebin].SetBinError(binx+1,biny+1,0.)
                #if local_flux_err/mean_flux_err>3.0:
                #    hist_skymap[ebin].SetBinContent(binx+1,biny+1,0.)
                #    hist_skymap[ebin].SetBinError(binx+1,biny+1,0.)


def SumFluxMap():

    hist_real_flux_skymap_sum.Reset()
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_real_flux_skymap_sum.Add(hist_real_flux_skymap[ebin])

    hist_real_flux_skymap_le.Reset()
    for ebin in range(max(0,energy_bin_cut_low),plot_energy_break_bin):
        hist_real_flux_skymap_le.Add(hist_real_flux_skymap[ebin])
    hist_real_flux_skymap_he.Reset()
    for ebin in range(max(plot_energy_break_bin,energy_bin_cut_low),6):
        hist_real_flux_skymap_he.Add(hist_real_flux_skymap[ebin])

    hist_bkgd_flux_skymap_sum.Reset()
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_bkgd_flux_skymap_sum.Add(hist_bkgd_flux_skymap[ebin])

    for imposter in range(0,n_imposters):
        hist_imposter_flux_skymap_sum[imposter].Reset()
        for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
            hist_imposter_flux_skymap_sum[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])
        hist_imposter_flux_skymap_le[imposter].Reset()
        for ebin in range(max(0,energy_bin_cut_low),plot_energy_break_bin):
            hist_imposter_flux_skymap_le[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])
        hist_imposter_flux_skymap_he[imposter].Reset()
        for ebin in range(max(plot_energy_break_bin,energy_bin_cut_low),6):
            hist_imposter_flux_skymap_he[imposter].Add(hist_imposter_flux_skymap[imposter][ebin])

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)

    hist_bkgd_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_bkgd_flux_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_bkgd_flux_skymap_reflect,hist_zscore_skymap_sum_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapBkgdFlux_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_sum)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,hist_zscore_skymap_sum_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_Sum_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    hist_zscore_skymap_le_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_le)

    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_le)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,hist_zscore_skymap_le_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_LE_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    hist_zscore_skymap_he_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_he)

    hist_real_flux_skymap_reflect = CommonPlotFunctions.reflectXaxis(hist_real_flux_skymap_he)
    CommonPlotFunctions.MatplotlibMap2D(hist_real_flux_skymap_reflect,hist_zscore_skymap_he_reflect,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapFlux_HE_%s'%(plot_tag),roi_x=region_x,roi_y=region_y,roi_r=region_r[0],rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)


    CommonPlotFunctions.SaveAsFITS(hist_real_flux_skymap_reflect,'skymap_flux_sum_%s'%(plot_tag))

def PrintSpectralDataForNaima(energy_axis,src_flux,src_flux_err,data_name):
    
    energy_mean_log = [] 
    energy_mean = [] 
    energy_edge_lo = [] 
    energy_edge_hi = [] 
    flux_mean = [] 
    flux_error = []
    ul = []
    for eb in range(0,len(energy_axis)):
        energy_mean_log += [math.log10(energy_axis[eb]/1000.)]
    for eb in range(0,len(energy_axis)):
        energy_log_delta = 0.
        if eb+1<len(energy_axis):
            energy_log_delta = energy_mean_log[eb+1]-energy_mean_log[eb]
        else:
            energy_log_delta = energy_mean_log[eb]-energy_mean_log[eb-1]
        energy_mean += [pow(10,energy_mean_log[eb])]
        energy_edge_lo += [pow(10,energy_mean_log[eb]-energy_log_delta)]
        energy_edge_hi += [pow(10,energy_mean_log[eb]+energy_log_delta)]
        flux_mean += [src_flux[eb]/((energy_axis[eb]/1000.)*(energy_axis[eb]/1000.))]
        flux_error += [src_flux_err[eb]/((energy_axis[eb]/1000.)*(energy_axis[eb]/1000.))]
    print ('=======================================================')
    print ('data_name = %s'%(data_name))
    for eb in range(0,len(energy_axis)):
        print ('%.2f %.2f %.2f %.2e %.2e %s'%(energy_mean[eb],energy_edge_lo[eb],energy_edge_hi[eb],flux_mean[eb],flux_error[eb],0))

def MakeSpectrum(roi_x,roi_y,roi_r,roi_name,excl_roi_x,excl_roi_y,excl_roi_r,excl_roi_name):

    energy_axis, energy_error, real_flux, real_flux_stat_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_real_flux_skymap,hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
    imposter_flux_list = []
    imposter_flux_err_list = []
    for imposter in range(0,n_imposters):
        energy_axis, energy_error, imposter_flux, imposter_flux_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_flux_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        imposter_flux_list += [imposter_flux]
        imposter_flux_err_list += [imposter_flux_err]

    real_flux_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_err += pow(imposter_flux_list[imposter][ebin],2)
        if correct_bias:
            syst_err = pow(syst_err/float(n_imposters-1),0.5)
        else:
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_flux_syst_err += [syst_err]

    energy_axis, energy_error, real_data, real_data_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_real_data_skymap,hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)

    imposter_data_list = []
    imposter_bkgd_list = []
    for imposter in range(0,n_imposters):
        energy_axis, energy_error, imposter_data, imposter_data_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_data_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        energy_axis, energy_error, imposter_bkgd, imposter_bkgd_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_imposter_bkgd_skymap[imposter],hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
        imposter_data_list += [imposter_data]
        imposter_bkgd_list += [imposter_bkgd]
    imposter_data_avg = []
    for ebin in range(0,len(energy_axis)):
        avg = 0.
        for imposter in range(0,n_imposters):
            avg += imposter_data_list[imposter][ebin]/n_imposters
        imposter_data_avg += [avg]
    imposter_s2b_lower_bound_list = []
    imposter_s2b_list = []
    imposter_s2b_err_list = []
    for imposter in range(0,n_imposters):
        imposter_s2b_lower_bound = []
        imposter_s2b = []
        imposter_s2b_err = []
        for ebin in range(0,len(energy_axis)):
            s2b_lower_bound = 0.
            s2b = 0.
            s2b_err = 0.
            if doImposter==1 and imposter_bkgd_list[imposter][ebin]>0. and imposter_data_avg[ebin]>0.:
                s2b_lower_bound = (imposter_data_list[imposter][ebin]-imposter_data_avg[ebin])/imposter_data_avg[ebin]
                s2b = (imposter_data_list[imposter][ebin]-imposter_bkgd_list[imposter][ebin])/imposter_bkgd_list[imposter][ebin]
                s2b_err = pow(imposter_data_list[imposter][ebin],0.5)/imposter_bkgd_list[imposter][ebin]
            imposter_s2b_lower_bound += [s2b_lower_bound]
            imposter_s2b += [s2b]
            imposter_s2b_err += [s2b_err]
        imposter_s2b_lower_bound_list += [imposter_s2b_lower_bound]
        imposter_s2b_list += [imposter_s2b]
        imposter_s2b_err_list += [imposter_s2b_err]

    real_rel_syst_lower_bound_bias = []
    real_rel_syst_lower_bound_err = []
    real_rel_syst_bias = []
    real_rel_syst_err = []
    for ebin in range(0,len(energy_axis)):
        syst_lower_bound_bias = 0.
        syst_lower_bound_err = 0.
        syst_bias = 0.
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_lower_bound_bias += imposter_s2b_lower_bound_list[imposter][ebin]
            syst_lower_bound_err += pow(imposter_s2b_lower_bound_list[imposter][ebin],2)
            syst_bias += imposter_s2b_list[imposter][ebin]
            syst_err += pow(imposter_s2b_list[imposter][ebin],2)
        if correct_bias:
            syst_lower_bound_bias = syst_lower_bound_bias/float(n_imposters-2)
            syst_lower_bound_err = pow(syst_lower_bound_err/float(n_imposters-2),0.5)
            syst_bias = syst_bias/float(n_imposters-1)
            syst_err = pow(syst_err/float(n_imposters-1),0.5)
        else:
            syst_lower_bound_bias = syst_lower_bound_bias/float(n_imposters-1)
            syst_lower_bound_err = pow(syst_lower_bound_err/float(n_imposters-1),0.5)
            syst_bias = syst_bias/float(n_imposters)
            syst_err = pow(syst_err/float(n_imposters),0.5)
        if syst_err>0:
            real_rel_syst_lower_bound_bias += [pow(n_imposters,0.5)*syst_lower_bound_bias/syst_lower_bound_err]
            real_rel_syst_lower_bound_err += [syst_lower_bound_err]
            real_rel_syst_bias += [pow(n_imposters,0.5)*syst_bias/syst_err]
            real_rel_syst_err += [syst_err]
        else:
            real_rel_syst_lower_bound_bias += [0.]
            real_rel_syst_lower_bound_err += [0.]
            real_rel_syst_bias += [0.]
            real_rel_syst_err += [0.]

    real_rel_stat_err = []
    for ebin in range(0,len(energy_axis)):
        if real_data[ebin]>0.:
            real_rel_stat_err += [1./pow(real_data[ebin],0.5)]
        else:
            real_rel_stat_err += [0.]

    vectorize_f_crab = np.vectorize(flux_crab_func)
    xdata_array = []
    for binx in range(0,len(energy_axis)):
        xdata_array += [energy_axis[binx]]
    ydata_crab_2 = pow(np.array(xdata_array)/1e3,energy_index_scale)*vectorize_f_crab(xdata_array)

    if 'Crab' in source_name:
        log_energy = np.linspace(log10(2e2),log10(1.2e4),50)
        xdata = pow(10.,log_energy)
        ydata_crab = pow(xdata/1e3,energy_index_scale)*vectorize_f_crab(xdata)
        calibration_new = []
        for binx in range(0,len(energy_axis)):
            if real_flux[binx]>0.:
                calibration_new += [ydata_crab_2[binx]/real_flux[binx]]
            else:
                calibration_new += [0.]
        print ('=======================================================================')
        print ('new flux_calibration = %s'%(calibration_new))
        print ('=======================================================================')
    if source_name=='MGRO_J1908' or source_name=='PSR_J1907_p0602':
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
        ydata_hawc = pow(xdata_ref/1e3,energy_index_scale)*vectorize_f_hawc(xdata_ref)
        # HAWC systematic uncertainty, The Astrophysical Journal 881, 134. Fig 13
        #HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1908(energy_index_scale)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHawcDiffusionFluxJ1908(energy_index_scale)
        HESS_energies, HESS_fluxes, HESS_flux_errs = GetHessFluxJ1908(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ1908(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1908(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1908(energy_index_scale)
    if source_name=='PSR_J1856_p0245':
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1856(energy_index_scale)
    if source_name=='MGRO_J2019' or source_name=='PSR_J2021_p3651':
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ2019(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ2019(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ2019(energy_index_scale)
    if source_name=='Boomerang' or source_name=='PSR_J2229_p6114':
        TAS_energies, TAS_fluxes, TAS_flux_errs = GetTibetASFluxJ2227(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ2227(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ2227(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ2226(energy_index_scale)
    if source_name=='HESS_J1825':
        HESS_0p8_energies, HESS_0p8_fluxes, HESS_0p8_flux_errs = GetHESSFluxR0p8J1825(energy_index_scale)
        HESS_0p4_energies, HESS_0p4_fluxes, HESS_0p4_flux_errs = GetHESSFluxR0p4J1825(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1825(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1825(energy_index_scale)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1825(energy_index_scale)
    if source_name=='LHAASO_J2108':
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ2108(energy_index_scale)
    if source_name=='IC443HotSpot' or source_name=='SNR_G189_p03':
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_ic443_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)
        vectorize_f_hawc = np.vectorize(flux_ic443_hawc_func)
        ydata_hawc = pow(xdata/1e3,energy_index_scale)*vectorize_f_hawc(xdata)
    if source_name=='Geminga':
        HAWC_diff_energies, HAWC_diff_fluxes, HAWC_diff_flux_errs = GetHAWCDiffusionFluxGeminga(energy_index_scale)
        HAWC_disk_energies, HAWC_disk_fluxes, HAWC_disk_flux_errs = GetHAWCDiskFluxGeminga(energy_index_scale)
        HAWC_gaus_energies, HAWC_gaus_fluxes, HAWC_gaus_flux_errs = GetHAWCGaussianFluxGeminga(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxGeminga(energy_index_scale)
        Fermi_UL_energies, Fermi_UL_fluxes, Fermi_UL_err = GetFermiUpperLimitFluxGeminga(energy_index_scale)
    if source_name=='1ES0229':
        log_energy = np.linspace(log10(2e2),log10(1.2e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_1es0229_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)
        HESS_energies, HESS_fluxes, HESS_flux_errs = GetHessFlux1ES0229(energy_index_scale)
    if source_name=='MGRO_J2019' or source_name=='PSR_J2021_p3651':
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_mgro_j2019_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)

    imposter_flux_crab_unit = []
    imposter_flux_err_crab_unit = []
    flux_crab_unit = []
    flux_stat_crab_unit = []
    flux_syst_crab_unit = []
    for imposter in range(0,n_imposters):
        imposter_flux_crab_unit_sublist = []
        imposter_flux_err_crab_unit_sublist = []
        for binx in range(0,len(energy_axis)):
            if ydata_crab_2[binx]>0.:
                imposter_flux_crab_unit_sublist += [imposter_flux_list[imposter][binx]/ydata_crab_2[binx]]
                imposter_flux_err_crab_unit_sublist += [imposter_flux_err_list[imposter][binx]/ydata_crab_2[binx]]
            else:
                imposter_flux_crab_unit_sublist += [0.]
                imposter_flux_err_crab_unit_sublist += [0.]
        imposter_flux_crab_unit += [imposter_flux_crab_unit_sublist]
        imposter_flux_err_crab_unit += [imposter_flux_err_crab_unit_sublist]
    for binx in range(0,len(energy_axis)):
        if ydata_crab_2[binx]>0.:
            flux_crab_unit += [real_flux[binx]/ydata_crab_2[binx]]
            flux_stat_crab_unit += [real_flux_stat_err[binx]/ydata_crab_2[binx]]
            flux_syst_crab_unit += [real_flux_syst_err[binx]/ydata_crab_2[binx]]
        else:
            flux_crab_unit += [0.]
            flux_stat_crab_unit += [0.]
            flux_syst_crab_unit += [0.]
    print ('=======================================================================')
    print ('Relative errors')
    for binx in range(0,len(energy_axis)):
        print ('+/- %0.3f (stat) +/- %0.3f (syst)'%(real_rel_stat_err[binx],real_rel_syst_err[binx]))
    print ('=======================================================================')
    print ('Data flux in Crab unit')
    for binx in range(0,len(energy_axis)):
        print ('%0.3f +/- %0.3f (stat) +/- %0.3f (syst)'%(flux_crab_unit[binx],flux_stat_crab_unit[binx],flux_syst_crab_unit[binx]))
    print ('=======================================================================')
    print ('Mimic flux in Crab unit')
    for imposter in range(0,n_imposters):
        for binx in range(0,len(energy_axis)):
            print ('(%s) %0.3f +/- %0.3f (stat)'%(imposter,imposter_flux_crab_unit[imposter][binx],imposter_flux_err_crab_unit[imposter][binx]))
    print ('=======================================================================')
    print ('Data flux in TeV/cm2/s')
    for binx in range(0,len(energy_axis)):
        print ('%0.3e +/- %0.3e (stat) +/- %0.3e (syst)'%(real_flux[binx],real_flux_stat_err[binx],real_flux_syst_err[binx]))
    print ('=======================================================================')
    print ('Mimic flux in TeV/cm2/s')
    for imposter in range(0,n_imposters):
        for binx in range(0,len(energy_axis)):
            print ('(%s) %0.3e +/- %0.3e (stat)'%(imposter,imposter_flux_list[imposter][binx],imposter_flux_err_list[imposter][binx]))
    print ('=======================================================================')

    energy_axis = np.array(energy_axis)
    energy_error = np.array(energy_error)
    real_flux = np.array(real_flux)
    real_flux_UL = np.array(real_flux)
    real_flux_stat_err = np.array(real_flux_stat_err)
    real_flux_syst_err = np.array(real_flux_syst_err)
    real_rel_syst_lower_bound_bias = np.array(real_rel_syst_lower_bound_bias)
    real_rel_syst_lower_bound_err = np.array(real_rel_syst_lower_bound_err)
    real_rel_syst_bias = np.array(real_rel_syst_bias)
    real_rel_syst_err = np.array(real_rel_syst_err)

    uplims = []
    zscore = []
    for eb in range(0,len(energy_axis)):
        zscore += [real_flux[eb]/pow(pow(real_flux_stat_err[eb],2)+pow(real_flux_syst_err[eb],2),0.5)]
        uplims += [0.]
        #if zscore[eb]<2.:
        #    uplims += [1.]
        #    real_flux_UL[eb] = 2.*pow(pow(real_flux_stat_err[eb],2)+pow(real_flux_syst_err[eb],2),0.5)
        #else:
        #    uplims += [0.]
        #uplims += [1.]
        #ratio_syst_2_stat_err = 20.
        #real_flux_UL[eb] = 2.*pow(pow(real_flux_stat_err[eb],2)+pow(ratio_syst_2_stat_err*real_flux_stat_err[eb],2),0.5)
    zscore = np.array(zscore)
    uplims = np.array(uplims, dtype=bool)

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_flux_list[imposter],imposter_flux_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=0.-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    plotname = 'ImposterSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_s2b_list[imposter],imposter_s2b_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_rel_syst_err, bottom=0.-real_rel_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    for ebin in range(0,len(energy_axis)):
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.2*real_rel_syst_err[len(real_rel_syst_err)-1], '%0.1f %%'%(real_rel_syst_err[ebin]*100.))
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.0*real_rel_syst_err[len(real_rel_syst_err)-1], '%0.1f $\sigma$'%(real_rel_syst_bias[ebin]))
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('relative error')
    axbig.set_xscale('log')
    plotname = 'ImposterRelError_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(energy_axis,imposter_s2b_lower_bound_list[imposter],imposter_s2b_err_list[imposter],xerr=energy_error,color=next_color,marker='_',ls='none')
    axbig.bar(energy_axis, 2.*real_rel_syst_lower_bound_err, bottom=0.-real_rel_syst_lower_bound_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
    for ebin in range(0,len(energy_axis)):
        axbig.text(energy_axis[ebin]-0.5*energy_error[ebin], 2.0*real_rel_syst_lower_bound_err[ebin], '%0.1f %%'%(real_rel_syst_lower_bound_err[ebin]*100.))
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('relative error')
    axbig.set_xscale('log')
    plotname = 'ImposterLowerBoundRelError_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    if source_name=='MGRO_J1908' or source_name=='PSR_J1907_p0602':

        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC',zorder=3)
        #axbig.errorbar(HESS_energies,HESS_fluxes,HESS_flux_errs,color='g',marker='s',ls='none',label='HESS',zorder=2)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='y',marker='s',ls='none',label='VERITAS (2014)',zorder=1)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=7)

        #axbig.plot(xdata_ref, ydata_hawc,'r-',label='1909.08609 (HAWC)',zorder=2)
        #axbig.fill_between(xdata_ref, ydata_hawc-0.15*ydata_hawc, ydata_hawc+0.15*ydata_hawc, alpha=0.2, color='r',zorder=1)

        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)

        real_flux_total_err = pow(real_flux_syst_err*real_flux_syst_err+real_flux_stat_err*real_flux_stat_err,0.5)
        PrintSpectralDataForNaima(energy_axis,real_flux_UL,real_flux_total_err,'VERITAS')
        PrintSpectralDataForNaima(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,'Fermi')
        PrintSpectralDataForNaima(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,'LHAASO')
        PrintSpectralDataForNaima(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,'HAWC')

    elif source_name=='HESS_J1825':
        axbig.errorbar(HESS_0p8_energies,HESS_0p8_fluxes,HESS_0p8_flux_errs,color='r',marker='s',ls='none',label='HESS 0.8-deg',zorder=1)
        axbig.errorbar(HESS_0p4_energies,HESS_0p4_fluxes,HESS_0p4_flux_errs,color='g',marker='s',ls='none',label='HESS 0.4-deg',zorder=2)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi',zorder=3)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='y',marker='s',ls='none',label='HAWC',zorder=5)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=6)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=7)
    elif source_name=='PSR_J1856_p0245':
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi (2209.11855)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='MGRO_J2019' or source_name=='PSR_J2021_p3651':
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='eHWC J2019+368',zorder=1)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='g',marker='s',ls='none',label='VER J2019+368',zorder=2)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi (2209.11855)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='Boomerang' or source_name=='PSR_J2229_p6114':
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.errorbar(TAS_energies,TAS_fluxes,TAS_flux_errs,color='r',marker='s',ls='none',label='Tibet AS',zorder=1)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='g',marker='s',ls='none',label='Fermi',zorder=2)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='y',marker='s',ls='none',label='VERITAS(2009)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='LHAASO_J2032':
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='LHAASO_J2108':
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='b',marker='_',ls='none',label='VERITAS (30 hrs)',uplims=uplims,zorder=6)
        #axbig.errorbar(energy_axis,real_flux_UL*(1./pow(4.,0.5)),real_flux_stat_err,xerr=energy_error,color='g',marker='_',ls='none',label='VERITAS (120 hrs)',uplims=uplims,zorder=6)
        #axbig.errorbar(energy_axis,real_flux_UL*(1./pow(12.,0.5)),real_flux_stat_err,xerr=energy_error,color='r',marker='_',ls='none',label='VERITAS (360 hrs)',uplims=uplims,zorder=6)
    elif source_name=='IC443HotSpot' or source_name=='SNR_G189_p03':
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (0905.3291)',zorder=1)
        axbig.plot(xdata, ydata_hawc,'g-',label='HAWC (2007.08582)',zorder=2)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=3)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=4)
    elif source_name=='1ES0229':
        #axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (2013ICRC)',zorder=1)
        axbig.errorbar(HESS_energies,HESS_fluxes,HESS_flux_errs,color='r',marker='s',ls='none',label='HESS',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=3)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=4)
    elif 'Crab' in source_name:
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=3)
    elif 'Geminga' in source_name:
        axbig.plot(HAWC_diff_energies, HAWC_diff_fluxes,'g-',label='HAWC diffusion')
        axbig.fill_between(HAWC_diff_energies, np.array(HAWC_diff_fluxes)-np.array(HAWC_diff_flux_errs), np.array(HAWC_diff_fluxes)+np.array(HAWC_diff_flux_errs), alpha=0.2, color='g')
        axbig.plot(HAWC_disk_energies, HAWC_disk_fluxes,'m-',label='HAWC disk')
        axbig.fill_between(HAWC_disk_energies, np.array(HAWC_disk_fluxes)-np.array(HAWC_disk_flux_errs), np.array(HAWC_disk_fluxes)+np.array(HAWC_disk_flux_errs), alpha=0.2, color='m')
        axbig.plot(HAWC_gaus_energies, HAWC_gaus_fluxes,'y-',label='HAWC gaussian')
        axbig.fill_between(HAWC_gaus_energies, np.array(HAWC_gaus_fluxes)-np.array(HAWC_gaus_flux_errs), np.array(HAWC_gaus_fluxes)+np.array(HAWC_gaus_flux_errs), alpha=0.2, color='y')
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='r',marker='_',ls='none',label='Fermi')
        fermi_uplims = np.array([1,1,1,1], dtype=bool)
        axbig.errorbar(Fermi_UL_energies,Fermi_UL_fluxes,Fermi_UL_err,color='r',marker='_',ls='none',uplims=fermi_uplims)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=3)
    else:
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS',uplims=uplims)

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'RealSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()


    energy_axis, energy_error, bkgd_flux, bkgd_flux_stat_err = CommonPlotFunctions.GetRegionSpectrum_v2(hist_bkgd_flux_skymap,hist_real_flux_syst_skymap,None,energy_bin_cut_low,energy_bin_cut_up,roi_x,roi_y,roi_r,excl_roi_x,excl_roi_y,excl_roi_r)
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.errorbar(energy_axis,bkgd_flux,bkgd_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='bkg flux')
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'BkgdSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.step(energy_axis, zscore, where='mid', color='b')
    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('Z score')
    axbig.set_xscale('log')
    axbig.legend(loc='best')
    plotname = 'ZscoreSpectrum_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

    eff_area = []
    for eb in range(0,len(hist_real_Aeff_skymap)):
        eff_area += [hist_real_Aeff_skymap[eb].GetMaximum()]
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('energy_axis = %s'%(energy_axis))
    print ('bkgd_flux = %s'%(bkgd_flux))
    print ('eff_area = %s'%(eff_area))

def MakeDiffusionSpectrum(energy_axis,energy_error,r_axis,y_axis,y_axis_stat_err,y_axis_syst_err,flux_norm,flux_norm_err):

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
            flux_syst_err += pow(brightness_syst_err*2.*3.14*radius*delta_r,2)
        #real_flux[eb] = flux
        real_flux[eb] = flux_norm[eb]
        real_flux_stat_err[eb] = pow(flux_stat_err,0.5)
        #real_flux_syst_err[eb] = pow(flux_syst_err,0.5)
        real_flux_syst_err[eb] = flux_norm_err[eb]

    vectorize_f_crab = np.vectorize(flux_crab_func)
    xdata_array = []
    for binx in range(0,len(energy_axis)):
        xdata_array += [energy_axis[binx]]
    ydata_crab_2 = pow(np.array(xdata_array)/1e3,energy_index_scale)*vectorize_f_crab(xdata_array)

    if 'Crab' in source_name:
        log_energy = np.linspace(log10(2e2),log10(1.2e4),50)
        xdata = pow(10.,log_energy)
        ydata_crab = pow(xdata/1e3,energy_index_scale)*vectorize_f_crab(xdata)
        calibration_new = []
        for binx in range(0,len(energy_axis)):
            if real_flux[binx]>0.:
                calibration_new += [ydata_crab_2[binx]/real_flux[binx]]
            else:
                calibration_new += [0.]
        print ('=======================================================================')
        print ('(diffusion spectrum) new flux_calibration = %s'%(calibration_new))
        print ('=======================================================================')
    if source_name=='MGRO_J1908' or source_name=='PSR_J1907_p0602':
        log_energy = np.linspace(log10(1e2),log10(1e5),50)
        xdata_ref = pow(10.,log_energy)
        vectorize_f_hawc = np.vectorize(flux_hawc_j1908_func)
        ydata_hawc = pow(xdata_ref/1e3,energy_index_scale)*vectorize_f_hawc(xdata_ref)
        # HAWC systematic uncertainty, The Astrophysical Journal 881, 134. Fig 13
        #HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1908(energy_index_scale)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHawcDiffusionFluxJ1908(energy_index_scale)
        HESS_energies, HESS_fluxes, HESS_flux_errs = GetHessFluxJ1908(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ1908(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1908(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1908(energy_index_scale)
    if source_name=='PSR_J1856_p0245':
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1856(energy_index_scale)
    if source_name=='MGRO_J2019' or source_name=='PSR_J2021_p3651':
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ2019(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ2019(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ2019(energy_index_scale)
    if source_name=='Boomerang' or source_name=='PSR_J2229_p6114':
        TAS_energies, TAS_fluxes, TAS_flux_errs = GetTibetASFluxJ2227(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ2227(energy_index_scale)
        OldV_energies, OldV_fluxes, OldV_flux_errs = GetVeritasFluxJ2227(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ2226(energy_index_scale)
    if source_name=='HESS_J1825':
        HESS_0p8_energies, HESS_0p8_fluxes, HESS_0p8_flux_errs = GetHESSFluxR0p8J1825(energy_index_scale)
        HESS_0p4_energies, HESS_0p4_fluxes, HESS_0p4_flux_errs = GetHESSFluxR0p4J1825(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxJ1825(energy_index_scale)
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ1825(energy_index_scale)
        HAWC_energies, HAWC_fluxes, HAWC_flux_errs = GetHAWCFluxJ1825(energy_index_scale)
    if source_name=='LHAASO_J2108':
        LHAASO_energies, LHAASO_fluxes, LHAASO_flux_errs = GetLHAASOFluxJ2108(energy_index_scale)
    if source_name=='IC443HotSpot' or source_name=='SNR_G189_p03':
        log_energy = np.linspace(log10(1e2),log10(1e4),50)
        xdata = pow(10.,log_energy)
        vectorize_f_veritas_paper = np.vectorize(flux_ic443_func)
        ydata_veritas_paper = pow(xdata/1e3,energy_index_scale)*vectorize_f_veritas_paper(xdata)
        vectorize_f_hawc = np.vectorize(flux_ic443_hawc_func)
        ydata_hawc = pow(xdata/1e3,energy_index_scale)*vectorize_f_hawc(xdata)
    if source_name=='Geminga':
        HAWC_diff_energies, HAWC_diff_fluxes, HAWC_diff_flux_errs = GetHAWCDiffusionFluxGeminga(energy_index_scale)
        HAWC_disk_energies, HAWC_disk_fluxes, HAWC_disk_flux_errs = GetHAWCDiskFluxGeminga(energy_index_scale)
        HAWC_gaus_energies, HAWC_gaus_fluxes, HAWC_gaus_flux_errs = GetHAWCGaussianFluxGeminga(energy_index_scale)
        Fermi_energies, Fermi_fluxes, Fermi_flux_errs = GetFermiFluxGeminga(energy_index_scale)
        Fermi_UL_energies, Fermi_UL_fluxes, Fermi_UL_err = GetFermiUpperLimitFluxGeminga(energy_index_scale)

    energy_axis = np.array(energy_axis)
    energy_error = np.array(energy_error)
    real_flux = np.array(real_flux)
    real_flux_stat_err = np.array(real_flux_stat_err)
    real_flux_syst_err = np.array(real_flux_syst_err)
    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    if source_name=='MGRO_J1908' or source_name=='PSR_J1907_p0602':
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='HAWC',zorder=3)
        axbig.errorbar(HESS_energies,HESS_fluxes,HESS_flux_errs,color='g',marker='s',ls='none',label='HESS',zorder=2)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='y',marker='s',ls='none',label='VERITAS (2014)',zorder=1)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=7)
        #axbig.plot(xdata_ref, ydata_hawc,'r-',label='1909.08609 (HAWC)',zorder=2)
        #axbig.fill_between(xdata_ref, ydata_hawc-0.15*ydata_hawc, ydata_hawc+0.15*ydata_hawc, alpha=0.2, color='r',zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=6)
    elif source_name=='HESS_J1825':
        axbig.errorbar(HESS_0p8_energies,HESS_0p8_fluxes,HESS_0p8_flux_errs,color='r',marker='s',ls='none',label='HESS 0.8-deg',zorder=1)
        axbig.errorbar(HESS_0p4_energies,HESS_0p4_fluxes,HESS_0p4_flux_errs,color='g',marker='s',ls='none',label='HESS 0.4-deg',zorder=2)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi',zorder=3)
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='y',marker='s',ls='none',label='HAWC',zorder=5)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=6)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=7)
    elif source_name=='PSR_J1856_p0245':
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi (2209.11855)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='MGRO_J2019' or source_name=='PSR_J2021_p3651':
        axbig.errorbar(HAWC_energies,HAWC_fluxes,HAWC_flux_errs,color='r',marker='s',ls='none',label='eHWC J2019+368',zorder=1)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='g',marker='s',ls='none',label='VER J2019+368',zorder=2)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='b',marker='s',ls='none',label='Fermi (2209.11855)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=6)
    elif source_name=='Boomerang' or source_name=='PSR_J2229_p6114':
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.errorbar(TAS_energies,TAS_fluxes,TAS_flux_errs,color='r',marker='s',ls='none',label='Tibet AS',zorder=1)
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='g',marker='s',ls='none',label='Fermi',zorder=2)
        axbig.errorbar(OldV_energies,OldV_fluxes,OldV_flux_errs,color='y',marker='s',ls='none',label='VERITAS(2009)',zorder=3)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=6)
    elif source_name=='LHAASO_J2032':
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',zorder=6)
    elif source_name=='LHAASO_J2108':
        axbig.errorbar(LHAASO_energies,LHAASO_fluxes,LHAASO_flux_errs,color='m',marker='s',ls='none',label='LHAASO',zorder=4)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2, zorder=5)
        axbig.errorbar(energy_axis,real_flux_UL,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (this work)',uplims=uplims,zorder=6)
    elif source_name=='IC443HotSpot' or source_name=='SNR_G189_p03':
        axbig.plot(xdata, ydata_veritas_paper,'r-',label='VERITAS (0905.3291)',zorder=1)
        axbig.plot(xdata, ydata_hawc,'g-',label='HAWC (2007.08582)',zorder=2)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=3)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=4)
    elif 'Crab' in source_name:
        axbig.plot(xdata, ydata_crab,'r-',label='1508.06442', zorder=1)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    elif 'Geminga' in source_name:
        axbig.plot(HAWC_diff_energies, HAWC_diff_fluxes,'g-',label='HAWC diffusion')
        axbig.fill_between(HAWC_diff_energies, np.array(HAWC_diff_fluxes)-np.array(HAWC_diff_flux_errs), np.array(HAWC_diff_fluxes)+np.array(HAWC_diff_flux_errs), alpha=0.2, color='g')
        axbig.plot(HAWC_disk_energies, HAWC_disk_fluxes,'m-',label='HAWC disk')
        axbig.fill_between(HAWC_disk_energies, np.array(HAWC_disk_fluxes)-np.array(HAWC_disk_flux_errs), np.array(HAWC_disk_fluxes)+np.array(HAWC_disk_flux_errs), alpha=0.2, color='m')
        axbig.plot(HAWC_gaus_energies, HAWC_gaus_fluxes,'y-',label='HAWC gaussian')
        axbig.fill_between(HAWC_gaus_energies, np.array(HAWC_gaus_fluxes)-np.array(HAWC_gaus_flux_errs), np.array(HAWC_gaus_fluxes)+np.array(HAWC_gaus_flux_errs), alpha=0.2, color='y')
        axbig.errorbar(Fermi_energies,Fermi_fluxes,Fermi_flux_errs,color='r',marker='_',ls='none',label='Fermi')
        uplims = np.array([1,1,1,1], dtype=bool)
        axbig.errorbar(Fermi_UL_energies,Fermi_UL_fluxes,Fermi_UL_err,color='r',marker='_',ls='none',uplims=uplims)
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2,zorder=2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS (new)',zorder=3)
    else:
        axbig.bar(energy_axis, 2.*real_flux_syst_err, bottom=real_flux-real_flux_syst_err, width=2.*energy_error, color='b', align='center', alpha=0.2)
        axbig.errorbar(energy_axis,real_flux,real_flux_stat_err,xerr=energy_error,color='k',marker='_',ls='none',label='VERITAS')

    axbig.set_xlabel('Energy [GeV]')
    axbig.set_ylabel('$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
    axbig.set_xscale('log')
    axbig.set_yscale('log')
    axbig.legend(loc='best')
    plotname = 'ProfileSpectrum'
    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    axbig.remove()

def MakeGalacticProfile(proj_type,axis_title,real_map,imposter_maps,erange_tag):

    real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindGalacticProjection_v2(real_map,proj_type=proj_type)
    imposter_profile_list = []
    imposter_profile_err_list = []
    for imposter in range(0,n_imposters):
        imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindGalacticProjection_v2(imposter_maps[imposter],proj_type=proj_type)
        imposter_profile_list += [imposter_profile]
        imposter_profile_err_list += [imposter_profile_stat_err]

    real_profile_syst_err = []
    for ubin in range(0,len(theta2)):
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_err += pow(imposter_profile_list[imposter][ubin],2)
        if correct_bias:
            syst_err = pow(syst_err/float(n_imposters-1),0.5)
        else:
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_profile_syst_err += [syst_err]

    theta2 = np.array(theta2)
    theta2_err = np.array(theta2_err)
    real_profile = np.array(real_profile)
    real_profile_stat_err = np.array(real_profile_stat_err)
    real_profile_syst_err = np.array(real_profile_syst_err)
    real_profile_total_err = []
    for ubin in range(0,len(theta2)):
        stat_err = real_profile_stat_err[ubin]
        syst_err = real_profile_syst_err[ubin]
        real_profile_total_err += [pow(stat_err*stat_err+syst_err*syst_err,0.5)]
    real_profile_total_err = np.array(real_profile_total_err)

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='%s-%s GeV'%(energy_bin[energy_bin_cut_low],energy_bin[energy_bin_cut_up]))
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{sr}^{-1}$]')
    axbig.set_xlabel(axis_title)
    #axbig.legend(loc='best')
    plotname = 'Profile%s'%(proj_type)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(theta2,imposter_profile_list[imposter],imposter_profile_err_list[imposter],color=next_color,marker='s',ls='none')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{sr}^{-1}$]')
    axbig.set_xlabel(axis_title)
    plotname = 'Profile%s_Imposter'%(proj_type)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()

def MakeExtensionProfile(roi_x,roi_y,roi_r,fit_profile,roi_name,real_map,imposter_maps,erange_tag):

    real_profile, real_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(real_map,roi_x,roi_y,3.0*roi_r[2])
    imposter_profile_list = []
    imposter_profile_err_list = []
    for imposter in range(0,n_imposters):
        imposter_profile, imposter_profile_stat_err, theta2, theta2_err = CommonPlotFunctions.FindExtension_v2(imposter_maps[imposter],roi_x,roi_y,3.0*roi_r[2])
        imposter_profile_list += [imposter_profile]
        imposter_profile_err_list += [imposter_profile_stat_err]

    real_profile_syst_err = []
    for ubin in range(0,len(theta2)):
        syst_err = 0.
        for imposter in range(0,n_imposters):
            syst_err += pow(imposter_profile_list[imposter][ubin],2)
        if correct_bias:
            syst_err = pow(syst_err/float(n_imposters-1),0.5)
        else:
            syst_err = pow(syst_err/float(n_imposters),0.5)
        real_profile_syst_err += [syst_err]

    theta2 = np.array(theta2)
    theta2_err = np.array(theta2_err)
    real_profile = np.array(real_profile)
    real_profile_stat_err = np.array(real_profile_stat_err)
    real_profile_syst_err = np.array(real_profile_syst_err)
    real_profile_total_err = []
    for ubin in range(0,len(theta2)):
        stat_err = real_profile_stat_err[ubin]
        syst_err = real_profile_syst_err[ubin]
        real_profile_total_err += [pow(stat_err*stat_err+syst_err*syst_err,0.5)]
    real_profile_total_err = np.array(real_profile_total_err)

    profile_sum = 0.
    for ubin in range(0,len(theta2)):
        profile_sum += real_profile[ubin]
    if fit_profile==1:
        start = (profile_sum, 0.5)
        popt, pcov = curve_fit(diffusion_func,theta2,real_profile,p0=start,sigma=real_profile_total_err,absolute_sigma=True,bounds=((0, 0.01), (np.inf, np.inf)))
        profile_fit = diffusion_func(theta2, *popt)
        residual = real_profile - profile_fit
        chisq = np.sum((residual/real_profile_stat_err)**2)
        dof = len(theta2)-2
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print (erange_tag)
        print ('diffusion flux = %0.2E +/- %0.2E'%(popt[0],pow(pcov[0][0],0.5)))
        print ('diffusion radius = %0.2f +/- %0.2f deg (chi2/dof = %0.2f)'%(popt[1],pow(pcov[1][1],0.5),chisq/dof))
    elif fit_profile==2:
        start = (profile_sum, 0.5)
        popt, pcov = curve_fit(gauss_func,theta2,real_profile,p0=start,sigma=real_profile_stat_err,absolute_sigma=True,bounds=((0, 0.01), (np.inf, np.inf)))
        profile_fit = gauss_func(theta2, *popt)
        residual = real_profile - profile_fit
        chisq = np.sum((residual/real_profile_stat_err)**2)
        dof = len(theta2)-2
        print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print (erange_tag)
        print ('gaussian flux = %0.2E +/- %0.2E'%(popt[0],pow(pcov[0][0],0.5)))
        print ('gaussian radius = %0.2f +/- %0.2f deg (chi2/dof = %0.2f)'%(popt[1],pow(pcov[1][1],0.5),chisq/dof))

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    axbig.errorbar(theta2,real_profile,real_profile_stat_err,color='k',marker='s',ls='none',label='ON data')
    if fit_profile!=0:
        axbig.plot(theta2,diffusion_func(theta2,*popt),color='r')
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    #axbig.set_ylim([-0.2e-11, 1.2e-11])
    plotname = 'ProfileVsTheta2_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()

    fig.clf()
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)
    axbig = fig.add_subplot()
    cycol = cycle('rgbcmy')
    axbig.bar(theta2, 2.*real_profile_syst_err, bottom=-real_profile_syst_err, width=1.*theta2_err, color='b', align='center', alpha=0.2)
    for imposter in range(0,n_imposters):
        next_color = next(cycol)
        axbig.errorbar(theta2,imposter_profile_list[imposter],imposter_profile_err_list[imposter],color=next_color,marker='s',ls='none',label='Mimic #%s'%(imposter+1))
    axbig.set_ylabel('surface brightness [$\mathrm{TeV}\ \mathrm{cm}^{-2}\mathrm{s}^{-1}\mathrm{deg}^{-2}$]')
    axbig.set_xlabel('angular distance from source [degree]')
    axbig.legend(loc='best')
    plotname = 'ProfileVsTheta2_Imposter_%s'%(roi_name)
    fig.savefig("output_plots/%s_%s_%s.png"%(plotname,erange_tag,plot_tag),bbox_inches='tight')
    axbig.remove()

    #if fit_profile:

    #    energy_axis = []
    #    energy_axis_err = []
    #    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    #        energy_axis += [0.5*(energy_bin[ebin]+energy_bin[ebin+1])]
    #        energy_axis_err += [0.5*(energy_bin[ebin+1]-energy_bin[ebin])]

    #    fig.clf()
    #    axbig = fig.add_subplot()
    #    axbig.errorbar(energy_axis,diff_radius,diff_radius_err,color='k',marker='s',ls='none')
    #    axbig.set_ylabel('diffusion radius [deg]')
    #    axbig.set_xlabel('photon energy [GeV]')
    #    axbig.set_xscale('log')
    #    plotname = 'RadiusVsEnergy_%s'%(roi_name)
    #    fig.savefig("output_plots/%s_%s.png"%(plotname,plot_tag),bbox_inches='tight')
    #    axbig.remove()

    #    MakeDiffusionSpectrum(energy_axis,energy_axis_err,r_axis_allE,y_axis_allE,y_axis_stat_err_allE,y_axis_syst_err_allE,norm_allE,norm_err_allE)



def MakeFluxMap(flux_map, data_map, bkgd_map, norm_map, aeff_map, expo_map):

    skymap_bin_size_x = data_map[0].GetXaxis().GetBinCenter(2)-data_map[0].GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = data_map[0].GetYaxis().GetBinCenter(2)-data_map[0].GetYaxis().GetBinCenter(1)
    calibration_radius = CommonPlotFunctions.calibration_radius
    for ebin in range(0,len(energy_bin)-1):
        flux_map[ebin].Reset()
        expo_content_max = expo_map[ebin].GetMaximum()
        norm_content_max = norm_map[0].GetMaximum()
        print ('energy = %0.1f GeV, expo_content_max = %0.2f hrs'%(energy_bin[ebin],expo_content_max))
        for binx in range(0,bkgd_map[ebin].GetNbinsX()):
            for biny in range(0,bkgd_map[ebin].GetNbinsY()):
                data_content = data_map[ebin].GetBinContent(binx+1,biny+1)
                data_error = data_map[ebin].GetBinError(binx+1,biny+1)
                bkgd_content = bkgd_map[ebin].GetBinContent(binx+1,biny+1)
                norm_content = norm_map[0].GetBinContent(binx+1,biny+1)
                expo_content = expo_map[ebin].GetBinContent(binx+1,biny+1)
                aeff_content = aeff_map[ebin].GetBinContent(binx+1,biny+1)

                #if data_error*data_error<10.: continue
                if norm_content==0.: continue
                if norm_content/norm_content_max<0.1: continue
                #if expo_content<5.: continue

                if UseEffectiveArea:
                    norm_content = norm_content*3600.
                    Aeff_correct = GetEffectiveAreaCorrection(ebin)
                    correction = Aeff_correct*1./(aeff_content*100.*100.)
                else:
                    correction = GetFluxCalibration(ebin)

                norm_ratio = norm_content/norm_content_max
                norm_weight = 1./(1.+np.exp(-(norm_ratio-0.3)/0.1))
                if expo_content_max<10.:
                    norm_weight = 1./(1.+np.exp(-(norm_ratio-0.7)/0.1))
                correction = correction*norm_weight

                stat_data_err = pow(max(data_content,0.),0.5)
                #stat_data_err = data_error
                flux_stat_err = max(stat_data_err,1.)/norm_content*correction*pow(energy_bin[ebin]/1e3,energy_index_scale)
                flux_content = (data_content-bkgd_content)/norm_content*correction*pow(energy_bin[ebin]/1e3,energy_index_scale)
                flux_map[ebin].SetBinContent(binx+1,biny+1,flux_content)
                flux_map[ebin].SetBinError(binx+1,biny+1,flux_stat_err)



folder_name = CommonPlotFunctions.folder_path
source_name = sys.argv[1]
source_name = source_name.split('_ON')[0]
source_name = source_name.split('_OFF')[0]
plot_tag = source_name+'_'+folder_name+'_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)
if doImposter==1:
    plot_tag += '_science'
else:
    plot_tag += '_fiction'
plot_tag += CommonPlotFunctions.additional_tag

InputFile = ROOT.TFile("output_fitting/%s_skymap_%s.root"%(sys.argv[1],folder_name))
HistName = "hist_data_skymap_0"
nbins_x = InputFile.Get(HistName).GetNbinsX()
nbins_y = InputFile.Get(HistName).GetNbinsY()
MapEdge_left = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(1)
MapEdge_right = InputFile.Get(HistName).GetXaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsX()+1)
MapEdge_lower = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(1)
MapEdge_upper = InputFile.Get(HistName).GetYaxis().GetBinLowEdge(InputFile.Get(HistName).GetNbinsY()+1)
InputFile.Close()

MapCenter_x = (MapEdge_left+MapEdge_right)/2.
MapCenter_y = (MapEdge_lower+MapEdge_upper)/2.

print ('MapEdge_left = %0.2f'%(MapEdge_left))
print ('MapEdge_right = %0.2f'%(MapEdge_right))
print ('MapEdge_lower = %0.2f'%(MapEdge_lower))
print ('MapEdge_upper = %0.2f'%(MapEdge_upper))

hist_real_elev_skymap = ROOT.TH2D("hist_real_elev_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_azim_skymap = ROOT.TH2D("hist_real_azim_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_nsb_skymap = ROOT.TH2D("hist_real_nsb_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_mjd_skymap = ROOT.TH2D("hist_real_mjd_skymap","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_zscore_skymap_sum = ROOT.TH2D("hist_real_zscore_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_zscore_skymap_le = ROOT.TH2D("hist_real_zscore_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_zscore_skymap_he = ROOT.TH2D("hist_real_zscore_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_sum = ROOT.TH2D("hist_real_flux_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_le = ROOT.TH2D("hist_real_flux_skymap_le","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_flux_skymap_he = ROOT.TH2D("hist_real_flux_skymap_he","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_bkgd_flux_skymap_sum = ROOT.TH2D("hist_bkgd_flux_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap_sum = ROOT.TH2D("hist_real_data_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_bkgd_skymap_sum = ROOT.TH2D("hist_real_bkgd_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_syst_skymap_sum = ROOT.TH2D("hist_real_syst_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_mimic_syst_skymap_sum = ROOT.TH2D("hist_mimic_syst_skymap_sum","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
hist_real_data_skymap = []
hist_real_bkgd_skymap = []
hist_real_norm_skymap = []
hist_real_Aeff_skymap = []
hist_real_expo_skymap = []
hist_real_bias_skymap = []
hist_real_flux_skymap = []
hist_real_flux_syst_skymap = []
hist_bkgd_flux_skymap = []
hist_null_skymap = []
for ebin in range(0,len(energy_bin)-1):
    hist_real_data_skymap += [ROOT.TH2D("hist_real_data_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bkgd_skymap += [ROOT.TH2D("hist_real_bkgd_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_norm_skymap += [ROOT.TH2D("hist_real_norm_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_Aeff_skymap += [ROOT.TH2D("hist_real_Aeff_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_expo_skymap += [ROOT.TH2D("hist_real_expo_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_bias_skymap += [ROOT.TH2D("hist_real_bias_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_skymap += [ROOT.TH2D("hist_real_flux_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_real_flux_syst_skymap += [ROOT.TH2D("hist_real_flux_syst_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_bkgd_flux_skymap += [ROOT.TH2D("hist_bkgd_flux_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_null_skymap += [ROOT.TH2D("hist_null_skymap_E%s"%(ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
InputFile = ROOT.TFile("output_fitting/%s_skymap_%s.root"%(sys.argv[1],folder_name))
HistName = "Hist_Data_Elev_Skymap"
hist_real_elev_skymap.Add(InputFile.Get(HistName))
HistName = "Hist_Data_Azim_Skymap"
hist_real_azim_skymap.Add(InputFile.Get(HistName))
HistName = "Hist_Data_NSB_Skymap"
hist_real_nsb_skymap.Add(InputFile.Get(HistName))
HistName = "Hist_Data_MJD_Skymap"
hist_real_mjd_skymap.Add(InputFile.Get(HistName))
for ebin in range(0,len(energy_bin)-1):
    HistName = "hist_data_skymap_%s"%(ebin)
    hist_real_data_skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "hist_bkgd_skymap_%s"%(ebin)
    if use_rfov:
        HistName = "hist_rfov_skymap_%s"%(ebin)
    hist_real_bkgd_skymap[ebin].Add(InputFile.Get(HistName))
    if UseEffectiveArea:
        HistName = "hist_expo_hour_skymap_sum"
    else:
        HistName = "hist_bkgd_skymap_smooth_%s"%(ebin)
    hist_real_norm_skymap[0].Add(InputFile.Get(HistName))
    HistName = "hist_effarea_skymap_%s"%(ebin)
    hist_real_Aeff_skymap[ebin].Add(InputFile.Get(HistName))
    #HistName = "hist_expo_hour_skymap_%s"%(ebin)
    HistName = "hist_expo_hour_skymap_sum"
    hist_real_expo_skymap[ebin].Add(InputFile.Get(HistName))
InputFile.Close()

hist_imposter_elev_skymap = []
hist_imposter_azim_skymap = []
hist_imposter_nsb_skymap = []
hist_imposter_mjd_skymap = []
hist_imposter_data_skymap = []
hist_imposter_bkgd_skymap = []
hist_imposter_norm_skymap = []
hist_imposter_Aeff_skymap = []
hist_imposter_expo_skymap = []
hist_imposter_bias_skymap = []
hist_imposter_flux_skymap = []
hist_imposter_data_skymap_sum = []
hist_imposter_bkgd_skymap_sum = []
hist_imposter_flux_skymap_sum = []
hist_imposter_flux_skymap_le = []
hist_imposter_flux_skymap_he = []
for imposter in range(0,n_imposters):
    hist_imposter_elev_skymap += [ROOT.TH2D("hist_imposter%s_elev_skymap"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_azim_skymap += [ROOT.TH2D("hist_imposter%s_azim_skymap"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_nsb_skymap += [ROOT.TH2D("hist_imposter%s_nsb_skymap"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_mjd_skymap += [ROOT.TH2D("hist_imposter%s_mjd_skymap"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sum += [ROOT.TH2D("hist_imposter%s_data_skymap_sum"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_bkgd_skymap_sum += [ROOT.TH2D("hist_imposter%s_bkgd_skymap_sum"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_sum += [ROOT.TH2D("hist_imposter%s_flux_skymap_sum"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_le += [ROOT.TH2D("hist_imposter%s_flux_skymap_le"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_flux_skymap_he += [ROOT.TH2D("hist_imposter%s_flux_skymap_he"%(imposter),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap_sublist = []
    hist_imposter_bkgd_skymap_sublist = []
    hist_imposter_norm_skymap_sublist = []
    hist_imposter_Aeff_skymap_sublist = []
    hist_imposter_expo_skymap_sublist = []
    hist_imposter_bias_skymap_sublist = []
    hist_imposter_flux_skymap_sublist = []
    for ebin in range(0,len(energy_bin)-1):
        hist_imposter_data_skymap_sublist += [ROOT.TH2D("hist_imposter%s_data_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bkgd_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bkgd_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_norm_skymap_sublist += [ROOT.TH2D("hist_imposter%s_norm_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_Aeff_skymap_sublist += [ROOT.TH2D("hist_imposter%s_Aeff_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_expo_skymap_sublist += [ROOT.TH2D("hist_imposter%s_expo_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_bias_skymap_sublist += [ROOT.TH2D("hist_imposter%s_bias_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
        hist_imposter_flux_skymap_sublist += [ROOT.TH2D("hist_imposter%s_flux_skymap_E%s"%(imposter,ebin),"",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)]
    hist_imposter_data_skymap += [hist_imposter_data_skymap_sublist]
    hist_imposter_bkgd_skymap += [hist_imposter_bkgd_skymap_sublist]
    hist_imposter_norm_skymap += [hist_imposter_norm_skymap_sublist]
    hist_imposter_Aeff_skymap += [hist_imposter_Aeff_skymap_sublist]
    hist_imposter_expo_skymap += [hist_imposter_expo_skymap_sublist]
    hist_imposter_bias_skymap += [hist_imposter_bias_skymap_sublist]
    hist_imposter_flux_skymap += [hist_imposter_flux_skymap_sublist]
if doImposter==1:
    for imposter in range(0,n_imposters):
        InputFile = ROOT.TFile("output_fitting/%s_Imposter%s_skymap_%s.root"%(source_name,imposter+1,folder_name))
        HistName = "Hist_Data_Elev_Skymap"
        hist_imposter_elev_skymap[imposter].Add(InputFile.Get(HistName))
        HistName = "Hist_Data_Azim_Skymap"
        hist_imposter_azim_skymap[imposter].Add(InputFile.Get(HistName))
        HistName = "Hist_Data_NSB_Skymap"
        hist_imposter_nsb_skymap[imposter].Add(InputFile.Get(HistName))
        HistName = "Hist_Data_MJD_Skymap"
        hist_imposter_mjd_skymap[imposter].Add(InputFile.Get(HistName))
        for ebin in range(0,len(energy_bin)-1):
            HistName = "hist_data_skymap_%s"%(ebin)
            hist_imposter_data_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            if UseEffectiveArea:
                HistName = "hist_expo_hour_skymap_sum"
            else:
                HistName = "hist_bkgd_skymap_smooth_%s"%(ebin)
            hist_imposter_norm_skymap[imposter][0].Add(InputFile.Get(HistName))
            HistName = "hist_effarea_skymap_%s"%(ebin)
            hist_imposter_Aeff_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            #HistName = "hist_expo_hour_skymap_%s"%(ebin)
            HistName = "hist_expo_hour_skymap_sum"
            hist_imposter_expo_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            HistName = "hist_bkgd_skymap_%s"%(ebin)
            if use_rfov:
                HistName = "hist_rfov_skymap_%s"%(ebin)
            hist_imposter_bkgd_skymap[imposter][ebin].Add(InputFile.Get(HistName))
            data_norm = hist_imposter_data_skymap[imposter][ebin].Integral()
            bkgd_norm = hist_imposter_bkgd_skymap[imposter][ebin].Integral()
            #hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.*data_norm/bkgd_norm)
            hist_imposter_bias_skymap[imposter][ebin].Add(InputFile.Get(HistName),-1.)
        InputFile.Close()

    exposure_scaling_factors = []
    for imposter in range(0,n_imposters):
        exposure_scaling_factors_per_imposter = []
        for ebin in range(0,len(energy_bin)-1):
            if hist_imposter_norm_skymap[imposter][ebin].Integral()==0.: continue
            expo_scale_global = hist_real_norm_skymap[ebin].Integral()/hist_imposter_norm_skymap[imposter][ebin].Integral()
            exposure_scaling_factors_per_imposter += [expo_scale_global]
            for binx in range(0,hist_real_norm_skymap[ebin].GetNbinsX()):
                for biny in range(0,hist_real_norm_skymap[ebin].GetNbinsY()):
                    real_expo = hist_real_expo_skymap[ebin].GetBinContent(binx+1,biny+1)
                    imposter_expo = hist_imposter_expo_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    if imposter_expo==0.: continue
                    expo_scale = real_expo/imposter_expo
                    old_data = hist_imposter_data_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    old_bias = hist_imposter_bias_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    old_norm = hist_imposter_norm_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    old_bkgd = hist_imposter_bkgd_skymap[imposter][ebin].GetBinContent(binx+1,biny+1)
                    hist_imposter_data_skymap[imposter][ebin].SetBinContent(binx+1,biny+1,old_data*expo_scale)
                    hist_imposter_bias_skymap[imposter][ebin].SetBinContent(binx+1,biny+1,old_bias*expo_scale)
                    hist_imposter_norm_skymap[imposter][ebin].SetBinContent(binx+1,biny+1,old_norm*expo_scale)
                    hist_imposter_bkgd_skymap[imposter][ebin].SetBinContent(binx+1,biny+1,old_bkgd*expo_scale)
        exposure_scaling_factors += [exposure_scaling_factors_per_imposter]

if correct_bias:
    for ebin in range(0,len(energy_bin)-1):
        for imposter in range(0,n_imposters):
            hist_real_bias_skymap[ebin].Add(hist_imposter_bias_skymap[imposter][ebin],1./float(n_imposters))
    for ebin in range(0,len(energy_bin)-1):
        #if ebin>=3: continue
        hist_real_bkgd_skymap[ebin].Add(hist_real_bias_skymap[ebin])
        for imposter in range(0,n_imposters):
            hist_imposter_bkgd_skymap[imposter][ebin].Add(hist_real_bias_skymap[ebin])
else:
    plot_tag += '_NC'

if use_rfov:
    plot_tag += '_rfov'

#if 'Crab' in source_name:
#    region_x = MapCenter_x
#    region_y = MapCenter_y
#    region_r = CommonPlotFunctions.calibration_radius
#    ring_r = region_r*2.
#    for ebin in range(0,len(energy_bin)-1):
#        nbins_ring = 0.
#        total_bkgd_ring = 0.
#        total_expo_ring = 0.
#        for binx in range(0,hist_real_bkgd_skymap[ebin].GetNbinsX()):
#            for biny in range(0,hist_real_bkgd_skymap[ebin].GetNbinsY()):
#                cell_x = hist_real_bkgd_skymap[ebin].GetXaxis().GetBinCenter(binx+1)
#                cell_y = hist_real_bkgd_skymap[ebin].GetYaxis().GetBinCenter(biny+1)
#                distance_to_center = pow(pow(region_x-cell_x,2)+pow(region_y-cell_y,2),0.5)
#                if distance_to_center>region_r and distance_to_center<ring_r:
#                    nbins_ring += 1.
#                    total_bkgd_ring += hist_real_bkgd_skymap[ebin].GetBinContent(binx+1,biny+1)
#                    total_expo_ring += hist_real_norm_skymap[ebin].GetBinContent(binx+1,biny+1)
#        avg_bkgd_ring = total_bkgd_ring/nbins_ring
#        avg_expo_ring = total_expo_ring/nbins_ring
#        for binx in range(0,hist_real_bkgd_skymap[ebin].GetNbinsX()):
#            for biny in range(0,hist_real_bkgd_skymap[ebin].GetNbinsY()):
#                cell_x = hist_real_bkgd_skymap[ebin].GetXaxis().GetBinCenter(binx+1)
#                cell_y = hist_real_bkgd_skymap[ebin].GetYaxis().GetBinCenter(biny+1)
#                distance_to_center = pow(pow(region_x-cell_x,2)+pow(region_y-cell_y,2),0.5)
#                if distance_to_center<region_r:
#                    hist_real_bkgd_skymap[ebin].SetBinContent(binx+1,biny+1,avg_bkgd_ring)
#                    hist_real_norm_skymap[ebin].SetBinContent(binx+1,biny+1,avg_expo_ring)


#MaskKnownSources()

hist_real_data_skymap_sum.Reset()
hist_real_bkgd_skymap_sum.Reset()
for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
    hist_real_data_skymap_sum.Add(hist_real_data_skymap[ebin])
    hist_real_bkgd_skymap_sum.Add(hist_real_bkgd_skymap[ebin])
for imposter in range(0,n_imposters):
    hist_imposter_data_skymap_sum[imposter].Reset()
    hist_imposter_bkgd_skymap_sum[imposter].Reset()
    for ebin in range(energy_bin_cut_low,energy_bin_cut_up):
        hist_imposter_data_skymap_sum[imposter].Add(hist_imposter_data_skymap[imposter][ebin])
        hist_imposter_bkgd_skymap_sum[imposter].Add(hist_imposter_bkgd_skymap[imposter][ebin])

MakeFluxMap(hist_bkgd_flux_skymap, hist_real_bkgd_skymap, hist_null_skymap, hist_real_norm_skymap, hist_real_Aeff_skymap, hist_real_expo_skymap)
MakeFluxMap(hist_real_flux_skymap, hist_real_data_skymap, hist_real_bkgd_skymap, hist_real_norm_skymap, hist_real_Aeff_skymap, hist_real_expo_skymap)

excl_region_x = MapCenter_x
excl_region_y = MapCenter_y
excl_region_r = [0.0 for element in range(len(energy_bin)-1)]
excl_region_name = 'Center'
region_x = MapCenter_x
region_y = MapCenter_y
region_r = [1.0 for element in range(len(energy_bin)-1)]
region_name = 'Center'
text_angle = 0.

prime_psr_name = None
prime_psr_ra = None
prime_psr_dec = None

MakeSignificanceMap(hist_real_data_skymap,hist_real_bkgd_skymap,hist_imposter_data_skymap,hist_imposter_bkgd_skymap)

prime_psr_name, prime_psr_ra, prime_psr_dec = FindDetectedPWN(hist_real_data_skymap,hist_real_bkgd_skymap,hist_imposter_data_skymap,hist_imposter_bkgd_skymap)

if 'Crab' in source_name:
    text_angle = 0.
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = [CommonPlotFunctions.calibration_radius for element in range(len(energy_bin)-1)]
    #region_r = [2.0 for element in range(len(energy_bin)-1)]
    region_name = 'Center'
    do_fit = 2
elif source_name=='PSR_J2021_p4026' and not CommonPlotFunctions.doGalacticCoord:
    text_angle = 0.
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = [1.0 for element in range(len(energy_bin)-1)]
    region_name = 'Center'
    do_fit = 0

    #hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)

    #Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    #Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)

    #MWL_map_file = 'MWL_maps/CGPS_1420MHz.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_title = 'H_{2} column density (1/cm^{2})'
    #Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    #Hist_mc_intensity.Scale(1.823*1e18) # optical thin
    #Hist_mc_column.Reset()
    #Hist_mc_column.Add(Hist_mc_intensity)
    #Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    #CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapCGPS_%s.png'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    ##MWL_map_file = 'MWL_maps/skv4052116888244_gammaCygni.txt' # CO intensity (K km s^{-1} deg)
    ##MWL_map_title = 'H_{2} column density (1/cm^{2})'
    ##Hist_mc_intensity.Reset()
    ##Hist_mc_column.Reset()
    ##Hist_mc_intensity = CommonPlotFunctions.GetSkyViewMap(MWL_map_file, Hist_mc_intensity, True)
    ##Hist_mc_intensity.Scale(2.*1e20) # H2 column density
    ##Hist_mc_column.Add(Hist_mc_intensity)
    ##Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    ##CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapSkv_%s.png'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    #pc_to_cm = 3.086e+18
    #CO_intensity_to_H_column_density = 2.*1e20
    #FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    #Hist_mc_intensity.Reset()
    #Hist_mc_column.Reset()
    #MWL_map_file = 'MWL_maps/DHT10_Cygnus_interp_m20_p20_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_title = 'H_{2} column density (1/cm^{2})'
    #Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    #Hist_mc_intensity.Scale(FITS_correction)
    #Hist_mc_column.Reset()
    #Hist_mc_column.Add(Hist_mc_intensity)
    #Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    #Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    #CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapDHT10_%s.png'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

elif source_name=='SS433' and not CommonPlotFunctions.doGalacticCoord:
    text_angle = 30.

    ##SS 433
    #region_x = 287.9565
    #region_y = 4.9827
    #region_r = [0.5 for element in range(len(energy_bin)-1)]
    #region_name = 'SS433'

    #SS 433 e1
    region_x = 288.404
    region_y = 4.930
    region_r = [0.25 for element in range(len(energy_bin)-1)]
    region_name = 'SS433e1'

    ##SS 433 e2
    #region_x = 288.58
    #region_y = 4.91
    #region_r = [0.5 for element in range(len(energy_bin)-1)]
    #region_name = 'SS433e2'

    ##SS 433 w1
    #287.654,5.037
    #region_x = 287.654
    #region_y = 5.037
    #region_r = [0.5 for element in range(len(energy_bin)-1)]
    #region_name = 'SS433w1'

    do_fit = 0

elif (source_name=='MGRO_J1908' or source_name=='PSR_J1907_p0602') and not CommonPlotFunctions.doGalacticCoord:
    text_angle = 30.

    ##3HWC J1908+063, 287.05, 6.39 
    #region_x = 287.05
    #region_y = 6.39
    #region_r = [1.5 for element in range(len(energy_bin)-1)]
    #region_name = '3HWC'

    #Fermi J1906+0626
    #region_x = 286.88
    #region_y = 6.29
    #region_r = 0.5
    #region_name = 'Fermi J1906+0626'

    #PSR J1907+0602 # cover more exposure hours, allow larger radius
    region_x = 286.975
    region_y = 6.03777777778
    region_r = [1.5 for element in range(len(energy_bin)-1)]
    region_name = 'PSR'

    # SNR G40.5-0.5
    #region_x = 286.79
    #region_y = 6.52
    #region_r = [0.5 for element in range(len(energy_bin)-1)]
    #region_name = 'SNR'

    ##4FGL J1906.9+0712
    excl_region_x = 286.7473
    excl_region_y = 7.2165
    excl_region_r = [0.3,0.3,0.0,0.0,0.0,0.0]
    excl_region_name = 'North'

    do_fit = 1

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)
    Hist_Tobias = ROOT.TH2D("Hist_Tobias","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    MWL_map_file = 'MWL_maps/tobias_map.fits'
    Hist_Tobias = CommonPlotFunctions.GetFITSMap(MWL_map_file, Hist_Tobias, True)
    Hist_Tobias_reflect = CommonPlotFunctions.reflectXaxis(Hist_Tobias)
    CommonPlotFunctions.MatplotlibMap2D(Hist_Tobias_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','Z score','SkymapTobias_%s'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    Hist_Fermi = ROOT.TH2D("Hist_Fermi","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    fermi_map_list = []
    fermi_map_list += ['1GeV']
    fermi_map_list += ['3GeV']
    fermi_map_list += ['5GeV']
    fermi_map_list += ['10GeV']
    fermi_map_list += ['30GeV']
    for hfile in range(0,len(fermi_map_list)):
        MWL_map_file = 'MWL_maps/TSmap_%s_pointsource_powerlaw_2.00_tsmap.fits'%(fermi_map_list[hfile])
        Hist_Fermi = CommonPlotFunctions.GetFITSMap(MWL_map_file, Hist_Fermi, True)
        Hist_Fermi = CommonPlotFunctions.Smooth2DMap(Hist_Fermi,0.1,False,False)
        Hist_Fermi_reflect = CommonPlotFunctions.reflectXaxis(Hist_Fermi)
        CommonPlotFunctions.MatplotlibMap2D(Hist_Fermi_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','TS','SkymapFermi_%s_%s'%(fermi_map_list[hfile],plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

    #Hist_Hawc = ROOT.TH2D("Hist_Hawc","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    #hawc_map_list = []
    #hawc_map_list += ['cd'] # 1-3.16 TeV
    ##hawc_map_list += ['ef'] # 3.16-10 TeV
    #hawc_map_list += ['gh'] # 10-31.6 TeV
    ##hawc_map_list += ['ij'] # 31.6-100 TeV
    #hawc_map_list += ['kl'] # 100-300 TeV
    #for hfile in range(0,len(hawc_map_list)):
    #    MWL_map_file = 'MWL_maps/%s-j1908.fits.gz'%(hawc_map_list[hfile])
    #    Hist_Hawc = CommonPlotFunctions.GetHealpixMap(MWL_map_file, Hist_Hawc, True)
    #    Hist_Hawc_reflect = CommonPlotFunctions.reflectXaxis(Hist_Hawc)
    #    CommonPlotFunctions.MatplotlibMap2D(Hist_Hawc_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','$E^{2}$ dN/dE [$\mathrm{TeV}\cdot\mathrm{cm}^{-2}\mathrm{s}^{-1}$]','SkymapHAWC_%s_%s.png'%(hawc_map_list[hfile],plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)



    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)
    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    # Dame, T. M.; Hartmann, Dap; Thaddeus, P., 2011, "Replication data for: First Quadrant, main survey (DHT08)", https://doi.org/10.7910/DVN/1PG9NV, Harvard Dataverse, V3
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PG9NV
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p10_p50_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p30_p50_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p10_p30_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p50_p70_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p50_p60_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    #MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_45_65_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    Hist_mc_intensity.Scale(FITS_correction)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapCOMap_p30p50_%s'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p10_p30_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    Hist_mc_intensity.Scale(FITS_correction)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapCOMap_p10p30_%s'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)
    MWL_map_file = 'MWL_maps/DHT08_Quad1_interp_p50_p70_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    Hist_mc_intensity.Scale(FITS_correction)
    Hist_mc_column.Reset()
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapCOMap_p50p70_%s'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)


elif source_name=='IC443HotSpot' or source_name=='SNR_G189_p03':
    text_angle = 45.
    region_x = 94.25
    region_y = 22.57
    region_r = [0.35 for element in range(len(energy_bin)-1)]
    region_name = 'G189.1+03.0'
    #region_x = 94.62
    #region_y = 22.17
    #region_r = [1.0 for element in range(len(energy_bin)-1)]
    #region_name = 'G189.6+03.3'
    do_fit = 0
elif source_name=='WComae':
    text_angle = 45.
    # 1ES 1218+304
    region_x = 185.360
    region_y = 30.191
    region_r = [0.2 for element in range(len(energy_bin)-1)]
    region_name = '1ES1218'
    do_fit = 0
elif source_name=='Geminga':
    text_angle = 45.
    region_x = MapCenter_x
    region_y = MapCenter_y
    region_r = [1.5 for element in range(len(energy_bin)-1)]
    region_name = 'Center'
    do_fit = 0
elif source_name=='HESS_J1825':
    text_angle = 45.
    region_x = 276.454
    region_y = -13.776
    region_r = [0.8 for element in range(len(energy_bin)-1)]
    region_name = '2HWC J1825-134'
    do_fit = 0
elif source_name=='Boomerang' or source_name=='PSR_J2229_p6114':
    text_angle = 45.
    region_x = 336.962
    region_y = 61.0482
    region_r = [1.0 for element in range(len(energy_bin)-1)]
    region_name = 'Boomerang'
    do_fit = 0
elif 'LHAASO_J2032' in source_name or 'PSR_J2032_p4127' in source_name:
    text_angle = 45.
    region_x = 308.042
    region_y = 41.459
    region_r = [1.0 for element in range(len(energy_bin)-1)]
    region_name = 'PSR J2032+4127'
    do_fit = 1
elif 'Tycho' in source_name:
    text_angle = 45.
    region_x = 6.340
    region_y = 64.130
    region_r = [2.0 for element in range(len(energy_bin)-1)]
    region_name = 'Tycho'
    do_fit = 0

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)

    Hist_fermi5 = ROOT.TH2D("Hist_fermi5","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_fermi5 = CommonPlotFunctions.GetSkyViewMap("MWL_maps/skv23149968874996_1.txt", Hist_fermi5, True)
    # 3-300 GeV Band 5, Atwood et al. 2009
    Hist_fermi5_reflect = CommonPlotFunctions.reflectXaxis(Hist_fermi5)
    CommonPlotFunctions.MatplotlibMap2D(Hist_fermi5_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','cnts$/s/cm^{2}/sr$','SkymapFermi5.png',rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

elif 'PSR_J2021_p3651' in source_name:
    text_angle = 45.

    region_x = 304.3
    region_y = 36.8
    region_r = [1.2 for element in range(len(energy_bin)-1)]
    region_name = 'Dragonfly'
    do_fit = 0

    #region_x = 304.85
    #region_y = 36.80
    #region_r = [1.0 for element in range(len(energy_bin)-1)]
    #region_name = 'VER J2019+368'
    #do_fit = 0

    #region_x = 305.27
    #region_y = 36.85
    #region_r = [0.14 for element in range(len(energy_bin)-1)]
    #region_name = 'PSR J2021+3651'
    #do_fit = 0

    hist_zscore_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_zscore_skymap_sum)

    Hist_mc_intensity = ROOT.TH2D("Hist_mc_intensity","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    Hist_mc_column = ROOT.TH2D("Hist_mc_column","",nbins_x,MapEdge_left,MapEdge_right,nbins_y,MapEdge_lower,MapEdge_upper)
    pc_to_cm = 3.086e+18
    CO_intensity_to_H_column_density = 2.*1e20
    FITS_correction = 1000.# the source FITS file has a mistake in velocity km/s -> m/s
    MWL_map_file = 'MWL_maps/DHT17_Quad2_bw_interp_p0_p40_0th_moment.txt' # CO intensity (K km s^{-1} deg)
    Hist_mc_intensity = CommonPlotFunctions.GetGalacticCoordMap(MWL_map_file, Hist_mc_intensity, True)
    Hist_mc_intensity.Scale(FITS_correction)
    Hist_mc_column.Add(Hist_mc_intensity)
    Hist_mc_column.Scale(CO_intensity_to_H_column_density) # H2 column density in unit of 1/cm2
    Hist_mc_column_reflect = CommonPlotFunctions.reflectXaxis(Hist_mc_column)
    CommonPlotFunctions.MatplotlibMap2D(Hist_mc_column_reflect,None,hist_zscore_skymap_sum_reflect,fig,'RA','Dec','column density [$1/cm^{2}$]','SkymapMolecularColumn_%s.png'%(plot_tag),rotation_angle=text_angle,prime_psr_name=prime_psr_name,prime_psr_ra=prime_psr_ra,prime_psr_dec=prime_psr_dec)

elif 'PSR_J1928_p1746' in source_name:
    text_angle = 45.
    region_x = 292.1115
    region_y = 18.0192
    region_r = [0.3 for element in range(len(energy_bin)-1)]
    region_name = '4FGL J1928.4+1801c'
    do_fit = 0

elif 'Segue1' in source_name:
    text_angle = 45.
    region_x = 151.83
    region_y = 16.76
    region_r = [0.3 for element in range(len(energy_bin)-1)]
    region_name = 'Star'
    do_fit = 0

else:
    text_angle = 45.
    region_x = MapCenter_x
    region_y = MapCenter_y
    #region_r = [CommonPlotFunctions.calibration_radius for element in range(len(energy_bin)-1)]
    region_r = [2.0 for element in range(len(energy_bin)-1)]
    region_name = 'Center'
    do_fit = 0

MakeDiagnisticPlots()
#CleanFluxMapNoise(hist_real_flux_skymap)
SumFluxMap()
MakeSpectrum(region_x,region_y,region_r,region_name,excl_region_x,excl_region_y,excl_region_r,excl_region_name)
MakeExtensionProfile(region_x,region_y,region_r,do_fit,region_name,hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeExtensionProfile(region_x,region_y,region_r,do_fit,region_name,hist_real_flux_skymap_le,hist_imposter_flux_skymap_le,'le')
MakeExtensionProfile(region_x,region_y,region_r,do_fit,region_name,hist_real_flux_skymap_he,hist_imposter_flux_skymap_he,'he')
MakeGalacticProfile('X','RA [deg]',hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeGalacticProfile('Y','Dec [deg]',hist_real_flux_skymap_sum,hist_imposter_flux_skymap_sum,'sum')
MakeGalacticProfile('X','RA [deg]',hist_real_flux_skymap_le,hist_imposter_flux_skymap_le,'le')
MakeGalacticProfile('Y','Dec [deg]',hist_real_flux_skymap_le,hist_imposter_flux_skymap_le,'le')
MakeGalacticProfile('X','RA [deg]',hist_real_flux_skymap_he,hist_imposter_flux_skymap_he,'he')
MakeGalacticProfile('Y','Dec [deg]',hist_real_flux_skymap_he,hist_imposter_flux_skymap_he,'he')
MakeSignificanceMap(hist_real_data_skymap,hist_real_bkgd_skymap,hist_imposter_data_skymap,hist_imposter_bkgd_skymap)

hist_real_data_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_data_skymap_sum)
hist_real_bkgd_skymap_sum_reflect = CommonPlotFunctions.reflectXaxis(hist_real_bkgd_skymap_sum)
CommonPlotFunctions.BackgroundSubtractMap(fig,hist_real_data_skymap_sum_reflect,hist_real_bkgd_skymap_sum_reflect,'RA','Dec','Count','SkymapBkgSubtraction_%s'%(plot_tag))

fit_2d_model(hist_real_data_skymap_sum, hist_real_bkgd_skymap_sum, region_x, region_y)

for esf in range(0,len(exposure_scaling_factors)):
    print ('exposure_scaling_factors = %s'%(exposure_scaling_factors[esf]))

