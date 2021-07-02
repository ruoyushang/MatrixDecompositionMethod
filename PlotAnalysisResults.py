
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
#from scipy import special
import scipy.stats as st

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

n_rebins = 2
smooth_size = 0.1

#method_tag = 'loose_mdm_default'
#method_tag = 'loose_mdm_rank3'
#method_tag = 'loose_mdm_rank5'
#method_tag = 'loose_mdm_cutoff'
#method_tag = 'loose_mdm_tikhonov'
method_tag = 'tight_mdm_default'
#method_tag = 'tight_mdm_rank3'
#method_tag = 'tight_mdm_rank5'
#method_tag = 'tight_mdm_cutoff'
#method_tag = 'tight_mdm_tikhonov'
#method_tag = 'tight_mdm_cutoff_eigen'

lowrank_tag = '_svd'
#lowrank_tag = '_eigen'
method_tag += lowrank_tag

energy_bin_cut_low = 0
energy_bin_cut_up = 6

elev_range = [45,85]
#elev_range = [25,55]
#elev_range = [45,55,65,75,85]
theta2_bins = [0,9]

ONOFF_tag = 'ON'
sample_list = []


if sys.argv[1]=='SgrA_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SgrAV6_ON']
    theta2_bins = [0,9]
    elev_range = [25,55]
    
if sys.argv[1]=='3C273_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['3C273V6_OFF']
    theta2_bins = [0,9]
if sys.argv[1]=='3C273_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['3C273V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0502_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0502V6_OFF']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES0502_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0502V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Draco_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['DracoV6_OFF']
    theta2_bins = [0,9]
if sys.argv[1]=='Draco_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['DracoV6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='OJ287_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['OJ287V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='OJ287_Model1':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model1'
    sample_list = []
    sample_list += ['OJ287V6_Model01']
    theta2_bins = [0,9]
if sys.argv[1]=='OJ287_Model2':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model2'
    sample_list = []
    sample_list += ['OJ287V6_Model02']
    theta2_bins = [0,9]
if sys.argv[1]=='OJ287_Model3':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model3'
    sample_list = []
    sample_list += ['OJ287V6_Model03']
    theta2_bins = [0,9]
    
if sys.argv[1]=='OJ287_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['OJ287V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0229_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0229V6_OFF']
    sample_list += ['1ES0229V5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0229_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0229V6_ON']
    sample_list += ['1ES0229V5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='H1426_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['H1426V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='H1426_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['H1426V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='PKS1424_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['PKS1424V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='PKS1424_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['PKS1424V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='3C264_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['3C264V6_ON']
    theta2_bins = [0,9]
if sys.argv[1]=='3C264_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['3C264V6_OFF']
    theta2_bins = [0,9]
    # The VERITAS observations of 3C 264 were taken from February through May 2017, from February through April 2018, and from January through May 2019.
    
if sys.argv[1]=='PG1553_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['PG1553V5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='PG1553_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['PG1553V5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='PG1553_Raster':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['PG1553V5_Raster']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES1011_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES1011V6_ON']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES1011_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES1011V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES1011_Model1':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model1'
    sample_list = []
    sample_list += ['1ES1011V6_Model01']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES1011_Model2':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model2'
    sample_list = []
    sample_list += ['1ES1011V6_Model02']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES1011_Model3':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model3'
    sample_list = []
    sample_list += ['1ES1011V6_Model03']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0229_Model1':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model1'
    sample_list = []
    sample_list += ['1ES0229V6_Model01']
    sample_list += ['1ES0229V5_Model01']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES0229_Model2':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model2'
    sample_list = []
    sample_list += ['1ES0229V6_Model02']
    sample_list += ['1ES0229V5_Model02']
    theta2_bins = [0,9]
if sys.argv[1]=='1ES0229_Model3':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model3'
    sample_list = []
    sample_list += ['1ES0229V6_Model03']
    sample_list += ['1ES0229V5_Model03']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Crab_Model1':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model1'
    sample_list = []
    sample_list += ['CrabV6_Model01']
    sample_list += ['CrabV5_Model01']
    sample_list += ['CrabV4_Model01']
    theta2_bins = [0,9]
if sys.argv[1]=='Crab_Model2':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model2'
    sample_list = []
    sample_list += ['CrabV6_Model02']
    sample_list += ['CrabV5_Model02']
    sample_list += ['CrabV4_Model02']
    theta2_bins = [0,9]
if sys.argv[1]=='Crab_Model3':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model3'
    sample_list = []
    sample_list += ['CrabV6_Model03']
    sample_list += ['CrabV5_Model03']
    sample_list += ['CrabV4_Model03']
    theta2_bins = [0,9]
if sys.argv[1]=='Crab_Model4':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model4'
    sample_list = []
    sample_list += ['CrabV6_Model04']
    sample_list += ['CrabV5_Model04']
    sample_list += ['CrabV4_Model04']
    theta2_bins = [0,9]
if sys.argv[1]=='Coma_Model5':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model5'
    sample_list = []
    sample_list += ['ComaV6_Model05']
    #sample_list += ['ComaV4_Model05']
    theta2_bins = [0,9]
if sys.argv[1]=='Coma_Model6':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model6'
    sample_list = []
    sample_list += ['ComaV6_Model06']
    #sample_list += ['ComaV4_Model06']
    theta2_bins = [0,9]
    
if sys.argv[1]=='RBS0413_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['RBS0413V6_OFF']
    sample_list += ['RBS0413V5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='RBS0413_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['RBS0413V6_ON']
    sample_list += ['RBS0413V5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0647_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0647V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='1ES0647_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['1ES0647V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='RGBJ0710_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['RGBJ0710V5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='RGBJ0710_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['RGBJ0710V5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Segue1_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['Segue1V6_OFF']
    sample_list += ['Segue1V5_OFF']
    theta2_bins = [0,9]
    # only V5 data published
    
if sys.argv[1]=='Segue1_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['Segue1V6_ON']
    sample_list += ['Segue1V5_ON']
    theta2_bins = [0,9]
    # only V5 data published
    
if sys.argv[1]=='BLLac_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['BLLacV6_OFF']
    sample_list += ['BLLacV5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='BLLac_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['BLLacV6_ON']
    sample_list += ['BLLacV5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='NGC1275_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['NGC1275V6_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='NGC1275_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['NGC1275V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='SNR_G150p3Plus04p5_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SNR_G150p3Plus04p5_V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='SNR_G150p3Plus04p5_Jamie_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SNR_G150p3Plus04p5_Jamie_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='CasA_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['CasAV6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='M82_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['M82V6_OFF']
    sample_list += ['M82V5_OFF']
    #sample_list += ['M82V4_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='M82_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['M82V6_ON']
    sample_list += ['M82V5_ON']
    #sample_list += ['M82V4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='M87_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['M87V6_ON']
    sample_list += ['M87V5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Crab_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['CrabV6_ON']
    sample_list += ['CrabV5_ON']
    sample_list += ['CrabV4_ON']
    theta2_bins = [0,9]
if sys.argv[1]=='Crab_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['CrabV6_OFF']
    sample_list += ['CrabV5_OFF']
    sample_list += ['CrabV4_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Mrk421_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['Mrk421V5_ON']
    theta2_bins = [0,9]
if sys.argv[1]=='Mrk421_OFF':
    ONOFF_tag = 'OFF'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['Mrk421V5_OFF']
    theta2_bins = [0,9]
    
if sys.argv[1]=='2HWC_J1930_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['2HWC_J1930V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='2HWC_J1953_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['2HWC_J1953V6_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='WComae_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['WComaeV6_ON']
    sample_list += ['WComaeV5_ON']
    sample_list += ['WComaeV4_ON']
    theta2_bins = [0,9]
    # https://arxiv.org/pdf/2002.04119.pdf VERITAS observations of 1ES 1215+303 from 2008 December to 2017 May.
    
if sys.argv[1]=='IC443HotSpot_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['IC443HotSpotV6_ON']
    sample_list += ['IC443HotSpotV5_ON']
    sample_list += ['IC443HotSpotV4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Boomerang_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['BoomerangV6_ON']
    sample_list += ['BoomerangV5_ON']
    sample_list += ['BoomerangV4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Tycho_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['TychoV6_ON']
    sample_list += ['TychoV5_ON']
    sample_list += ['TychoV4_ON']
    theta2_bins = [0,9]

if sys.argv[1]=='HESS_J1825_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['HESS_J1825_V6_ON']
    theta2_bins = [0,9]
    elev_range = [25,55]

if sys.argv[1]=='MGRO_J1908_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['MGRO_J1908_V6_ON']
    sample_list += ['MGRO_J1908_V5_ON']
    sample_list += ['MGRO_J1908_V4_ON']
    theta2_bins = [0,9]
    # this is a Tevatron
    
if sys.argv[1]=='SS433_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SS433_V6_ON']
    sample_list += ['SS433_V5_ON']
    sample_list += ['SS433_V4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='SS433Half1_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SS433Half1_V6_ON']
    sample_list += ['SS433Half1_V5_ON']
    sample_list += ['SS433Half1_V4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='SS433Half2_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['SS433Half2_V6_ON']
    sample_list += ['SS433Half2_V5_ON']
    sample_list += ['SS433Half2_V4_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='MAGIC_J1857_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['MAGIC_J1857_V6_ON']
    sample_list += ['MAGIC_J1857_V5_ON']
    #sample_list += ['MAGIC_J1857_V4_ON']
    theta2_bins = [0,9]
    # this is a Tevatron, largely extended at 400 GeV
    
if sys.argv[1]=='MGRO_J2031_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['MGRO_J2031_V6_ON']
    sample_list += ['MGRO_J2031_V5_ON']
    sample_list += ['MGRO_J2031_V4_ON']
    theta2_bins = [0,9]
    # this is a Tevatron with time-variable morphology
    
if sys.argv[1]=='Cygnus_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['CygnusV6_ON']
    sample_list += ['CygnusV5_ON']
    theta2_bins = [0,9]
    # this is a Tevatron
    
if sys.argv[1]=='GammaCygni_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['GammaCygniV4_ON']
    sample_list += ['GammaCygniV5_ON']
    sample_list += ['GammaCygniV6_ON']
    theta2_bins = [0,9]
    # this is a Tevatron
    
if sys.argv[1]=='Geminga_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['GemingaV6_ON']
    sample_list += ['GemingaV5_ON']
    theta2_bins = [0,9]
    
if sys.argv[1]=='Coma_ON':
    ONOFF_tag = 'ON'
    ONOFF_tag += '_Model0'
    sample_list = []
    sample_list += ['ComaV6_ON']
    #sample_list += ['ComaV4_ON']
    theta2_bins = [0,9]

root_file_tags = []
# all time
mjd_tag = []
mjd_tag += ['']
#mjd_tag += ['_MJD53613to55074']
#mjd_tag += ['_MJD55074to56535']
#mjd_tag += ['_MJD56535to57996']
#mjd_tag += ['_MJD57996to59457']

for elev in range(0,len(elev_range)-1):
    elev_tag = '_TelElev%sto%s'%(elev_range[elev],elev_range[elev+1])
    for u in range(0,len(theta2_bins)-1):
        theta2_tag = '_Theta2%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        #theta2_tag = '_Y%sto%s'%(theta2_bins[u],theta2_bins[u+1])
        for d in range(0,len(mjd_tag)):
            root_file_tags += [method_tag+elev_tag+theta2_tag+mjd_tag[d]+'_'+ONOFF_tag]

print 'Get %s'%(root_file_tags[0])

#for elev in range(0,len(elev_range)-1):
#    # all time
#    root_file_tags += [method_tag+'_TelElev%sto%s_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # 1ES 1215 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56372to56464_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # WComae and 1ES 1218 flare
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD54400to54700_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # Geminga
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56000to56500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58000to58700_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD58700to59000_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # MGRO J2031
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to57997_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57997to59000_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    # 3 slices of 1ES 0229
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD55000to56500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD56501to57500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]
#    #root_file_tags += [method_tag+'_TelElev%sto%s_MJD57501to59500_%s'%(elev_range[elev],elev_range[elev+1],ONOFF_tag)]

selection_tag = root_file_tags[0]

#folder_path = 'output_test'
folder_path = 'output_root'
PercentCrab = ''

selection_tag += '_%s'%(folder_path)

#N_bins_for_deconv = 16
N_bins_for_deconv = 12
gamma_hadron_dim_ratio_w = 1.
gamma_hadron_dim_ratio_l = 1.
MSCW_blind_cut = 0.5
MSCL_blind_cut = 0.5
MSCW_chi2_upper = -0.6
MSCL_chi2_upper = -0.6
MSCW_plot_lower = -0.6
MSCL_plot_lower = -0.6
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
ErecS_lower_cut = 0
ErecS_upper_cut = 0

n_good_matches = 0
exposure_hours = 0.
exposure_hours_ref = 0.
NSB_mean_data = 0.
Zenith_mean_data = 0.
NSB_RMS_data = 0.
Zenith_RMS_data = 0.
Skymap_size = 3.
Skymap_nbins = 120
source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
n_control_samples = 1
if sys.argv[1]=='SgrA':
    n_control_samples = 1
MJD_Start = 2147483647
MJD_End = 0
Data_runlist_NSB = ROOT.std.vector("double")(10)
Data_runlist_L3Rate = ROOT.std.vector("double")(10)
Data_runlist_L3Rate_all = ROOT.std.vector("double")(10)
roi_name = ROOT.std.vector("string")(10)
roi_ra = ROOT.std.vector("double")(10)
roi_dec = ROOT.std.vector("double")(10)
roi_radius = ROOT.std.vector("double")(10)
max_chi2_diff2_position = ROOT.std.vector("double")(10)
max_chi2_diff2_position_this_energy = 0.

Syst_MDM = 0.02
Syst_Init = 0.02
Syst_Redu = 0.02
Syst_Corr = 0.02
Syst_Clos = 0.02

energy_bin = []
energy_bin += [int(pow(10,2.0))]
energy_bin += [int(pow(10,2.33))]
energy_bin += [int(pow(10,2.66))]
energy_bin += [int(pow(10,3.0))]
energy_bin += [int(pow(10,3.33))]
energy_bin += [int(pow(10,3.66))]
energy_bin += [int(pow(10,4.0))]

energy_syst = []
energy_syst += [0.012]
energy_syst += [0.011]
energy_syst += [0.025]
energy_syst += [0.038]
energy_syst += [0.061]
energy_syst += [0.240]

#energy_fine_bin = []
#energy_fine_bin += [pow(10,2.0)]
#energy_fine_bin += [pow(10,2.2)]
#energy_fine_bin += [pow(10,2.4)]
#energy_fine_bin += [pow(10,2.6)]
#energy_fine_bin += [pow(10,2.8)]
#energy_fine_bin += [pow(10,3.0)]
#energy_fine_bin += [pow(10,3.2)]
#energy_fine_bin += [pow(10,3.4)]
#energy_fine_bin += [pow(10,3.6)]
#energy_fine_bin += [pow(10,3.8)]
#energy_fine_bin += [pow(10,4.0)]

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

elev_bins = [25,35,45,55,65,75,85]
MJD_bins = [53613,55074,56535,57996,59457]

other_stars = []
other_star_coord = []

bright_star_ra = []
bright_star_dec = []
bright_star_brightness = []
faint_star_ra = []
faint_star_dec = []
faint_star_brightness = []

#NRGBs = 5
#NCont = 104
#stops = [0.00,0.40,0.50,0.60,1.00]
#red =   [0.00,1.00,1.00,1.00,1.00]
#green = [0.00,1.00,1.00,1.00,0.00]
#blue =  [1.00,1.00,1.00,1.00,0.00]
#ROOT.TColor.CreateGradientColorTable(NRGBs,array('d',stops),array('d',red),array('d',green),array('d',blue),NCont)
#ROOT.gStyle.SetNumberContours(NCont)

def set_palette(name, ncontours=999):
    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    ROOT.TColor.CreateGradientColorTable(len(stops),array('d',stops),array('d',red),array('d',green),array('d',blue),ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

def GetBrightStarInfo(file_path):

    global bright_star_ra
    global bright_star_dec
    global bright_star_brightness
    global faint_star_ra
    global faint_star_dec
    global faint_star_brightness

    print 'Read file: %s'%(file_path)
    if not os.path.isfile(file_path): 
        print 'No such a file!'
        return
    InputFile = ROOT.TFile(file_path)
    StarTree = InputFile.Get("StarTree")
    print 'StarTree.GetEntries() = %s'%(StarTree.GetEntries())
    for entry in range(0,StarTree.GetEntries()):
        StarTree.GetEntry(entry)
        distance_to_center = pow(pow(source_ra-StarTree.star_ra,2)+pow(source_dec-StarTree.star_dec,2),0.5)
        bright_star_ra += [StarTree.star_ra]
        bright_star_dec += [StarTree.star_dec]
        bright_star_brightness += [StarTree.star_brightness]
    FaintStarTree = InputFile.Get("FaintStarTree")
    for entry in range(0,FaintStarTree.GetEntries()):
        FaintStarTree.GetEntry(entry)
        distance_to_center = pow(pow(source_ra-FaintStarTree.faint_star_ra,2)+pow(source_dec-FaintStarTree.faint_star_dec,2),0.5)
        faint_star_ra += [FaintStarTree.faint_star_ra]
        faint_star_dec += [FaintStarTree.faint_star_dec]
        faint_star_brightness += [FaintStarTree.faint_star_brightness]
    InputFile.Close()

def GetGammaSourceInfo():

    global other_stars
    global other_star_coord

    inputFile = open('PSR_RaDec_w_Names.txt')
    for line in inputFile:
        gamma_source_name = line.split(',')[0]
        gamma_source_ra = float(line.split(',')[1])
        gamma_source_dec = float(line.split(',')[2])
        near_a_source = False
        for entry in range(0,len(other_stars)):
            distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
            if distance<0.3*0.3:
                near_a_source = True
        if not near_a_source and not '%' in gamma_source_name:
            other_stars += [gamma_source_name]
            other_star_coord += [[gamma_source_ra,gamma_source_dec]]

    inputFile = open('TeVCat_RaDec_w_Names.txt')
    for line in inputFile:
        gamma_source_name = line.split(',')[0]
        gamma_source_ra = float(line.split(',')[1])
        gamma_source_dec = float(line.split(',')[2])
        near_a_source = False
        for entry in range(0,len(other_stars)):
            distance = pow(gamma_source_ra-other_star_coord[entry][0],2)+pow(gamma_source_dec-other_star_coord[entry][1],2)
            if distance<0.3*0.3:
                near_a_source = True
        if not near_a_source and not '%' in gamma_source_name:
            other_stars += [gamma_source_name]
            other_star_coord += [[gamma_source_ra,gamma_source_dec]]

def ResetStackedShowerHistograms():

    Hist_EffArea_Sum.Reset()

    Hist2D_OnData_Sum.Reset()
    Hist2D_OnBkgd_Sum.Reset()
    Hist2D_OnGamma_Sum.Reset()
    Hist2D_OnDark_Sum.Reset()
    Hist_OnData_MSCL_Sum.Reset()
    Hist_OnBkgd_MSCL_Sum.Reset()
    Hist_OnDark_MSCL_Sum.Reset()
    Hist_OnData_MSCW_Sum.Reset()
    Hist_OnBkgd_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_wGamma_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_woGamma_MSCW_Sum.Reset()
    Hist_OnGamma_MSCW_Sum.Reset()
    Hist_OnBkgd_Unblind_wGamma_MSCL_Sum.Reset()
    Hist_OnBkgd_Unblind_woGamma_MSCL_Sum.Reset()
    Hist_OnGamma_MSCL_Sum.Reset()
    Hist_OnDark_MSCW_Sum.Reset()

    Hist2D_Rank0_Data_Sum.Reset()
    Hist2D_Rank1_Data_Sum.Reset()
    Hist2D_Rank2_Data_Sum.Reset()
    Hist2D_Rank3_Data_Sum.Reset()
    Hist2D_Rank4_Data_Sum.Reset()
    Hist2D_Rank0_Dark_Sum.Reset()
    Hist2D_Rank1_Dark_Sum.Reset()
    Hist2D_Rank2_Dark_Sum.Reset()
    Hist2D_Rank3_Dark_Sum.Reset()
    Hist2D_Rank4_Dark_Sum.Reset()

    Hist2D_U_Proj_Sum.Reset()
    Hist2D_V_Proj_Sum.Reset()
    Hist2D_Coeff_Data_Sum.Reset()
    Hist2D_Coeff_Bkgd_Sum.Reset()

    Hist_SystErr_MSCL.Reset()
    Hist_SystErr_MSCW.Reset()
    Hist_SystErr_Energy.Reset()
    Hist_SystErr_Theta2.Reset()

    Hist_OnData_Yoff_Sum.Reset()
    Hist_OnBkgd_Yoff_Sum.Reset()
    Hist_OnBkgd_Yoff_Raw_Sum.Reset()
    Hist_OnData_Theta2_Sum.Reset()
    Hist_OnBkgd_Theta2_Sum.Reset()
    Hist_OnRFoV_Theta2_Sum.Reset()
    Hist_OnDark_Theta2_Sum.Reset()
    Hist_OnData_Energy_Sum.Reset()
    Hist_OnBkgd_Energy_Sum.Reset()
    Hist_OnRFoV_Energy_Sum.Reset()
    Hist_OnDark_Energy_Sum.Reset()
    Hist_OnData_Zenith_Sum.Reset()
    Hist_OnBkgd_Zenith_Sum.Reset()
    Hist_OnData_Energy_CamCenter_Sum.Reset()
    Hist_OnBkgd_Energy_CamCenter_Sum.Reset()

    Hist_OnData_Skymap_Sum.Reset()
    Hist_OnBkgd_Skymap_Sum.Reset()
    Hist_OnData_Skymap_Galactic_Sum.Reset()
    Hist_OnBkgd_Skymap_Galactic_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_MDM.Reset()
    Hist_OnBkgd_Skymap_Galactic_Syst_MDM.Reset()

    Hist_OnData_Skymap_ProjX_Sum.Reset()
    Hist_OnBkgd_Skymap_ProjX_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_ProjX_Sum.Reset()
    Hist_OnData_Skymap_ProjY_Sum.Reset()
    Hist_OnBkgd_Skymap_ProjY_Sum.Reset()
    Hist_OnBkgd_Skymap_Syst_ProjY_Sum.Reset()

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Reset() 
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Reset() 
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Reset()
        Hist_OnData_RoI_MJD_Sum[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Reset()

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Reset()
        Hist_OffData_MSCL_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCL_Sum[nth_sample].Reset()
        Hist_OffData_MSCW_Sum[nth_sample].Reset()
        Hist_OffBkgd_MSCW_Sum[nth_sample].Reset()
        Hist_OffData_Energy_Sum[nth_sample].Reset()
        Hist_OffBkgd_Energy_Sum[nth_sample].Reset()
        Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].Reset()
        Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].Reset()

def GetSourceInfo(file_list):

    global N_bins_for_deconv
    global MSCW_blind_cut
    global MSCL_blind_cut
    global MSCW_chi2_upper
    global MSCL_chi2_upper
    #global MSCW_plot_lower
    #global MSCL_plot_lower
    global n_good_matches
    global exposure_hours
    global exposure_hours_ref
    global NSB_mean_data
    global Zenith_mean_data
    global NSB_RMS_data
    global Zenith_RMS_data
    global Skymap_size
    global Skymap_nbins
    global source_ra
    global source_dec
    global source_l
    global source_b
    global n_control_samples
    global MJD_Start
    global MJD_End
    global Data_runlist_NSB
    global Data_runlist_L3Rate
    global Data_runlist_L3Rate_all
    global roi_name
    global roi_ra
    global roi_dec
    global roi_radius
    global max_chi2_diff2_position

    n_good_matches = 0
    exposure_hours = 0.
    exposure_hours_ref = 0.
    NSB_mean_data = 0.
    Zenith_mean_data = 0.
    NSB_RMS_data = 0.
    Zenith_RMS_data = 0.
    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.SetBranchAddress('Data_runlist_NSB',ROOT.AddressOf(Data_runlist_NSB))
        InfoTree.SetBranchAddress('Data_runlist_L3Rate',ROOT.AddressOf(Data_runlist_L3Rate))
        InfoTree.SetBranchAddress('Data_runlist_L3Rate_all',ROOT.AddressOf(Data_runlist_L3Rate_all))
        InfoTree.SetBranchAddress('roi_name',ROOT.AddressOf(roi_name))
        InfoTree.SetBranchAddress('roi_ra',ROOT.AddressOf(roi_ra))
        InfoTree.SetBranchAddress('roi_dec',ROOT.AddressOf(roi_dec))
        InfoTree.SetBranchAddress('roi_radius',ROOT.AddressOf(roi_radius))
        InfoTree.GetEntry(0)
        Hist_NSB.Reset()
        Hist_L3Rate.Reset()
        Hist_L3Rate_all.Reset()
        for entry in range(0,len(Data_runlist_NSB)):
            Hist_NSB.Fill(Data_runlist_NSB[entry])
        for entry in range(0,len(Data_runlist_L3Rate)):
            Hist_L3Rate.Fill(Data_runlist_L3Rate[entry])
        for entry in range(0,len(Data_runlist_L3Rate_all)):
            Hist_L3Rate_all.Fill(Data_runlist_L3Rate_all[entry])
        #N_bins_for_deconv = InfoTree.N_bins_for_deconv
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        MSCW_chi2_upper = InfoTree.MSCW_chi2_upper
        MSCL_chi2_upper = InfoTree.MSCL_chi2_upper
        #MSCW_plot_lower = InfoTree.MSCW_plot_lower
        #MSCL_plot_lower = InfoTree.MSCL_plot_lower
        Skymap_size = InfoTree.Skymap_size
        Skymap_nbins = InfoTree.Skymap_nbins
        Skymap_nbins = Skymap_nbins/n_rebins
        source_ra = InfoTree.mean_tele_point_ra
        source_dec = InfoTree.mean_tele_point_dec
        source_l = InfoTree.mean_tele_point_l
        source_b = InfoTree.mean_tele_point_b
        if 'Theta20' in file_list[path]:
            exposure_hours += InfoTree.exposure_hours_usable
        MJD_Start = min(InfoTree.MJD_Start,MJD_Start)
        MJD_End = max(InfoTree.MJD_End,MJD_End)

        NewInfoTree = InputFile.Get("NewInfoTree")
        NewInfoTree.SetBranchAddress('max_chi2_diff2_position',ROOT.AddressOf(max_chi2_diff2_position))
        NewInfoTree.GetEntry(0)
        NSB_mean_data = NewInfoTree.NSB_mean_data
        Zenith_mean_data = NewInfoTree.Zenith_mean_data
        NSB_RMS_data = NewInfoTree.NSB_RMS_data
        Zenith_RMS_data = NewInfoTree.Zenith_RMS_data
        #if 'Theta20' in file_list[path] and 'Up' in file_list[path]:
        #if 'Y0' in file_list[path]:

        HistName = "Hist_EffArea"
        Hist_EffArea.Reset()
        Hist_EffArea.Add(InputFile.Get(HistName))
        Hist_EffArea.Scale(InfoTree.exposure_hours_usable*3600.)
        Hist_EffArea_Sum.Add(Hist_EffArea)

        InputFile.Close()

def MergeHistogram(hist_sum,hist_input):

    bin_size_x_1 = hist_sum.GetXaxis().GetBinCenter(1)-hist_sum.GetXaxis().GetBinCenter(0)
    bin_size_y_1 = hist_sum.GetYaxis().GetBinCenter(1)-hist_sum.GetYaxis().GetBinCenter(0)
    bin_size_x_2 = hist_input.GetXaxis().GetBinCenter(1)-hist_input.GetXaxis().GetBinCenter(0)
    bin_size_y_2 = hist_input.GetYaxis().GetBinCenter(1)-hist_input.GetYaxis().GetBinCenter(0)

    bin_size_x_ratio = bin_size_x_2/bin_size_x_1
    bin_size_y_ratio = bin_size_y_2/bin_size_y_1

    hist_output = hist_sum.Clone()
    for binx in range(0,hist_sum.GetNbinsX()):
        for biny in range(0,hist_sum.GetNbinsY()):
            old_content = hist_sum.GetBinContent(binx+1,biny+1)
            old_error = pow(hist_sum.GetBinError(binx+1,biny+1),2)
            x = hist_sum.GetXaxis().GetBinCenter(binx+1)
            y = hist_sum.GetYaxis().GetBinCenter(biny+1)
            new_binx = hist_input.GetXaxis().FindBin(x)-1
            new_biny = hist_input.GetYaxis().FindBin(y)-1
            new_content = hist_input.GetBinContent(new_binx+1,new_biny+1)/(bin_size_x_ratio*bin_size_y_ratio) 
            new_error = pow(hist_input.GetBinError(new_binx+1,new_biny+1),2)/(bin_size_x_ratio*bin_size_y_ratio) 
            hist_output.SetBinContent(binx+1,biny+1,old_content+new_content)
            hist_output.SetBinError(binx+1,biny+1,pow(old_error+new_error,0.5))

    return hist_output

def GetShowerHistogramsFromFile(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    global MSCW_chi2_upper
    global MSCL_chi2_upper
    global max_chi2_diff2_position_this_energy

    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    MSCW_chi2_upper = InfoTree.MSCW_chi2_upper
    MSCL_chi2_upper = InfoTree.MSCL_chi2_upper
    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnData.Reset()
    Hist2D_OnData.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnData,InputFile.Get(HistName))
    #Hist2D_OnData.Add(HistTemp)

    HistName = "Hist_OnDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnDark.Reset()
    Hist2D_OnDark.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnDark,InputFile.Get(HistName))
    #Hist2D_OnDark.Add(HistTemp)

    HistName = "Hist_OnBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnBkgd.Reset()
    Hist2D_OnBkgd.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnBkgd,InputFile.Get(HistName))
    #Hist2D_OnBkgd.Add(HistTemp)

    HistName = "Hist_Gamma_MSCLW_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    print 'Getting histogram %s'%(HistName)
    Hist2D_OnGamma.Reset()
    Hist2D_OnGamma.Add(InputFile.Get(HistName))
    #HistTemp = MergeHistogram(Hist2D_OnGamma,InputFile.Get(HistName))
    #Hist2D_OnGamma.Add(HistTemp)

    Hist2D_Rank0_Data.Reset()
    Hist2D_Rank1_Data.Reset()
    Hist2D_Rank2_Data.Reset()
    Hist2D_Rank3_Data.Reset()
    Hist2D_Rank4_Data.Reset()
    Hist2D_Rank0_Dark.Reset()
    Hist2D_Rank1_Dark.Reset()
    Hist2D_Rank2_Dark.Reset()
    Hist2D_Rank3_Dark.Reset()
    Hist2D_Rank4_Dark.Reset()

    Hist2D_U_Proj.Reset()
    Hist2D_V_Proj.Reset()
    Hist2D_Coeff_Data.Reset()
    Hist2D_Coeff_Bkgd.Reset()

    HistName = "Hist_Rank0_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank0_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank1_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank1_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank2_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank2_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank3_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank3_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank4_MSCLW_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank4_Data.Add(InputFile.Get(HistName))

    HistName = "Hist_Rank0_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank0_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank1_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank1_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank2_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank2_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank3_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank3_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_Rank4_MSCLW_Dark_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Rank4_Dark.Add(InputFile.Get(HistName))

    HistName = "Hist_U_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_U_Proj.Add(InputFile.Get(HistName))

    HistName = "Hist_V_Proj_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_V_Proj.Add(InputFile.Get(HistName))

    HistName = "Hist_Coeff_Data_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Coeff_Data.Add(InputFile.Get(HistName))

    HistName = "Hist_Coeff_Bkgd_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist2D_Coeff_Bkgd.Add(InputFile.Get(HistName))

    #if Hist1D_Data_Rank0_LeftVector.Integral()==0:
    Hist_VVV_Eigenvalues.Reset()
    Hist_Bkgd_Chi2.Reset()
    Hist_Bkgd_Chi2_Diff.Reset()
    Hist_Bkgd_Chi2_Diff2.Reset()
    Hist_Bkgd_Optimization.Reset()
    Hist_Bkgd_Converge_Blind.Reset()
    Hist_Bkgd_Converge_Unblind.Reset()
    Hist1D_Data_Rank0_LeftVector.Reset()
    Hist1D_Data_Rank1_LeftVector.Reset()
    Hist1D_Data_Rank2_LeftVector.Reset()
    Hist1D_Data_Rank3_LeftVector.Reset()
    Hist1D_Data_Rank0_RightVector.Reset()
    Hist1D_Data_Rank1_RightVector.Reset()
    Hist1D_Data_Rank2_RightVector.Reset()
    Hist1D_Data_Rank3_RightVector.Reset()
    Hist1D_Bkgd_Rank0_LeftVector.Reset()
    Hist1D_Bkgd_Rank1_LeftVector.Reset()
    Hist1D_Bkgd_Rank2_LeftVector.Reset()
    Hist1D_Bkgd_Rank3_LeftVector.Reset()
    Hist1D_Bkgd_Rank0_RightVector.Reset()
    Hist1D_Bkgd_Rank1_RightVector.Reset()
    Hist1D_Bkgd_Rank2_RightVector.Reset()
    Hist1D_Bkgd_Rank3_RightVector.Reset()
    Hist1D_Dark_Rank0_LeftVector.Reset()
    Hist1D_Dark_Rank1_LeftVector.Reset()
    Hist1D_Dark_Rank2_LeftVector.Reset()
    Hist1D_Dark_Rank3_LeftVector.Reset()
    Hist1D_Dark_Rank0_RightVector.Reset()
    Hist1D_Dark_Rank1_RightVector.Reset()
    Hist1D_Dark_Rank2_RightVector.Reset()
    Hist1D_Dark_Rank3_RightVector.Reset()
    HistName = "Hist_VVV_Eigenvalues_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_VVV_Eigenvalues.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_Diff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2_Diff.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Chi2_Diff2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Chi2_Diff2.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Optimization_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Optimization.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Converge_Blind_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Converge_Blind.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Converge_Unblind_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_Bkgd_Converge_Unblind.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Data_Rank3_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Bkgd_Rank3_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank0_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank0_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank1_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank1_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank2_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank2_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank3_LeftVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank3_LeftVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank0_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank0_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank1_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank1_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank2_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank2_RightVector.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_Rank3_RightVector_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist1D_Dark_Rank3_RightVector.Add(InputFile.Get(HistName))

    Hist_VVV_Eigenvalues.SetMinimum(1e-2)
    #MakeOneHistPlot(Hist_VVV_Eigenvalues,'entry','eigenvalue','VVV_Eigenvalue_%s_%s'%(sample_list[0],ErecS_lower_cut_int),True)
    #MakeOneHistPlot(Hist_Bkgd_Chi2,'log10 #alpha','#chi^{2} in CR','Bkgd_Chi2_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)
    #MakeOneHistPlot(Hist_Bkgd_Chi2_Diff,'log10 #alpha','#chi^{2} 1st derivative in CR','Bkgd_Chi2_Diff_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)
    #MakeOneHistPlot(Hist_Bkgd_Chi2_Diff2,'log10 #alpha','#chi^{2} 2nd derivative in CR','Bkgd_Chi2_Diff2_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)
    #MakeOneHistPlot(Hist_Bkgd_Optimization,'log10 #alpha','abs(N_{#gamma}-N_{model})/N_{#gamma}','Bkgd_Optimization_%s_%s'%(sample_list[0],ErecS_lower_cut_int),False)

    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank0_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank0_LeftVector]
    #legends += ['inital vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank0_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{1}}','Rank0_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank1_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank1_LeftVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank1_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{2}}','Rank1_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank2_LeftVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank2_LeftVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank2_LeftVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{p_{3}}','Rank2_LeftVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank0_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank0_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank0_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{1}}','Rank0_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank1_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank1_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank1_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{2}}','Rank1_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist1D_Data_Rank2_RightVector]
    #legends += ['data vector']
    #colors += [1]
    #Hists += [Hist1D_Dark_Rank2_RightVector]
    #legends += ['initial vector']
    #colors += [2]
    #Hists += [Hist1D_Bkgd_Rank2_RightVector]
    #legends += ['predict. vector']
    #colors += [3]
    #MakeComparisonPlot(Hists,legends,colors,'entry','#vec{q_{3}}','Rank2_RightVector_%s_%s'%(sample_list[0],ErecS_lower_cut_int),0,0,False,False)

    if Hist2D_OnData.Integral()<1600. or Hist2D_OnDark.Integral()<1600.:
        Hist2D_OnData.Reset()
        Hist2D_OnDark.Reset()
        Hist2D_OnBkgd.Reset()
        Hist2D_OnGamma.Reset()
        Hist2D_Rank0_Data.Reset()
        Hist2D_Rank1_Data.Reset()
        Hist2D_Rank2_Data.Reset()
        Hist2D_Rank3_Data.Reset()
        Hist2D_Rank4_Data.Reset()
        Hist2D_Rank0_Dark.Reset()
        Hist2D_Rank1_Dark.Reset()
        Hist2D_Rank2_Dark.Reset()
        Hist2D_Rank3_Dark.Reset()
        Hist2D_Rank4_Dark.Reset()
        Hist2D_U_Proj.Reset()
        Hist2D_V_Proj.Reset()
        Hist2D_Coeff_Data.Reset()
        Hist2D_Coeff_Bkgd.Reset()

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
    Hist_OnGamma_MSCW.Reset()
    Hist_OnGamma_MSCW.Add(Hist2D_OnGamma.ProjectionY("Hist1D_OnGamma_MSCW",bin_lower_x,bin_upper_x))
    Hist_OnGamma_MSCL.Reset()
    Hist_OnGamma_MSCL.Add(Hist2D_OnGamma.ProjectionX("Hist1D_OnGamma_MSCL",bin_lower_y,bin_upper_y))

    for nth_sample in range(0,n_control_samples):

        HistName = "Hist_OffData2_MSCLW_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        print 'Getting histogram %s'%(HistName)
        Hist2D_OffData[nth_sample].Reset()
        Hist2D_OffData[nth_sample].Add(InputFile.Get(HistName))
        #HistTemp = MergeHistogram(Hist2D_OffData[nth_sample],InputFile.Get(HistName))
        #Hist2D_OffData[nth_sample].Add(HistTemp)

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist2D_OffData[nth_sample].Reset()
        if math.isnan(Hist2D_OffData[nth_sample].Integral()):
            Hist2D_OffData[nth_sample].Reset()

        Hist_OffData_MSCL[nth_sample].Reset()
        Hist_OffData_MSCL[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionX("Hist1D_OffData_MSCL_%s"%(nth_sample),bin_lower_y,bin_upper_y))

        Hist_OffData_MSCW[nth_sample].Reset()
        Hist_OffData_MSCW[nth_sample].Add(Hist2D_OffData[nth_sample].ProjectionY("Hist1D_OffData_MSCW_%s"%(nth_sample),bin_lower_x,bin_upper_x))


def StackShowerHistograms():

    Hist2D_OnData_Sum.Add(Hist2D_OnData)
    Hist2D_OnDark_Sum.Add(Hist2D_OnDark)
    Hist2D_OnBkgd_Sum.Add(Hist2D_OnBkgd)
    Hist2D_OnGamma_Sum.Add(Hist2D_OnGamma)

    if Hist2D_U_Proj_Sum.Integral()==0.:
        Hist2D_U_Proj_Sum.Add(Hist2D_U_Proj)
        Hist2D_V_Proj_Sum.Add(Hist2D_V_Proj)
        Hist2D_Coeff_Data_Sum.Add(Hist2D_Coeff_Data)
        Hist2D_Coeff_Bkgd_Sum.Add(Hist2D_Coeff_Bkgd)

    #Hist2D_Rank0_Data_Sum.Add(Hist2D_Rank0_Data)
    #Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank1_Data)
    #Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank2_Data)
    #Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank3_Data)
    Hist2D_Rank0_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    Hist2D_Rank0_Data_Sum.Add(Hist2D_OnData)
    Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    Hist2D_Rank1_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    Hist2D_Rank1_Data_Sum.Add(Hist2D_OnData)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    Hist2D_Rank2_Data_Sum.Add(Hist2D_OnData)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_Rank3_Data,-1.)
    Hist2D_Rank3_Data_Sum.Add(Hist2D_OnData)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank0_Data,-1.)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank1_Data,-1.)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank2_Data,-1.)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank3_Data,-1.)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_Rank4_Data,-1.)
    Hist2D_Rank4_Data_Sum.Add(Hist2D_OnData)

    Hist2D_Rank0_Dark_Sum.Add(Hist2D_Rank0_Data)
    Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank1_Data)
    Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank2_Data)
    Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank3_Data)
    #Hist2D_Rank0_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank0_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank1_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank2_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_Rank3_Dark,-1.)
    #Hist2D_Rank3_Dark_Sum.Add(Hist2D_OnDark)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank0_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank1_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank2_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank3_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_Rank4_Dark,-1.)
    #Hist2D_Rank4_Dark_Sum.Add(Hist2D_OnDark)

    Hist_OnData_MSCL_Sum.Add(Hist_OnData_MSCL)
    Hist_OnData_MSCW_Sum.Add(Hist_OnData_MSCW)
    Hist_OnDark_MSCL_Sum.Add(Hist_OnDark_MSCL)
    Hist_OnDark_MSCW_Sum.Add(Hist_OnDark_MSCW)
    Hist_OnBkgd_MSCL_Sum.Add(Hist_OnBkgd_MSCL)
    Hist_OnBkgd_MSCW_Sum.Add(Hist_OnBkgd_MSCW)
    Hist_OnBkgd_Unblind_wGamma_MSCW_Sum.Add(Hist_OnBkgd_Unblind_wGamma_MSCW)
    Hist_OnBkgd_Unblind_woGamma_MSCW_Sum.Add(Hist_OnBkgd_Unblind_woGamma_MSCW)
    Hist_OnGamma_MSCW_Sum.Add(Hist_OnGamma_MSCW)
    Hist_OnBkgd_Unblind_wGamma_MSCL_Sum.Add(Hist_OnBkgd_Unblind_wGamma_MSCL)
    Hist_OnBkgd_Unblind_woGamma_MSCL_Sum.Add(Hist_OnBkgd_Unblind_woGamma_MSCL)
    Hist_OnGamma_MSCL_Sum.Add(Hist_OnGamma_MSCL)

    for nth_sample in range(0,n_control_samples):

        Hist2D_OffData_Sum[nth_sample].Add(Hist2D_OffData[nth_sample])

        Hist_OffData_MSCL_Sum[nth_sample].Add(Hist_OffData_MSCL[nth_sample])
        Hist_OffData_MSCW_Sum[nth_sample].Add(Hist_OffData_MSCW[nth_sample])
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

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

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

    gap = 0.1*(max_heigh-min_heigh)
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(2)
        Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy,normalized):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.1)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

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

    gap = 0.1*(max_heigh-min_heigh)
    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        if normalized: Hists[0].DrawNormalized("E")
        else: Hists[0].Draw("E")
    else:
        if not logy:
            Hists[max_hist].SetMaximum(max_heigh+gap)
            Hists[max_hist].SetMinimum(min_heigh-gap)
        if normalized: Hists[max_hist].DrawNormalized("E")
        else: Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        #if Hists[h]!=0:
        Hists[h].SetLineColor(colors[h])
        Hists[h].SetLineWidth(2)
        if normalized: Hists[h].DrawNormalized("E same")
        else: Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

#def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
#    
#    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
#    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
#    pad3.SetBottomMargin(0.0)
#    pad3.SetTopMargin(0.03)
#    pad3.SetBorderMode(1)
#    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
#    pad1.SetBottomMargin(0.1)
#    pad1.SetTopMargin(0.0)
#    pad1.SetBorderMode(0)
#    pad1.SetGrid()
#    pad1.Draw()
#    pad3.Draw()
#
#    pad1.cd()
#    if logy: pad1.SetLogy()
#
#    min_heigh = 0
#    max_heigh = 0
#    max_hist = 0
#    mean = []
#    rms = []
#    amp = []
#    for h in range(0,len(Hists)):
#        mean += [0]
#        rms += [0]
#        amp += [0]
#        if Hists[h]!=0:
#            Hists[h].GetXaxis().SetTitleOffset(0.8)
#            Hists[h].GetXaxis().SetTitleSize(0.06)
#            Hists[h].GetXaxis().SetLabelSize(0.06)
#            Hists[h].GetYaxis().SetLabelSize(0.06)
#            Hists[h].GetYaxis().SetTitleOffset(1.2)
#            Hists[h].GetYaxis().SetTitleSize(0.06)
#            Hists[h].GetXaxis().SetTitle(title_x)
#            Hists[h].GetYaxis().SetTitle(title_y)
#            if max_heigh < Hists[h].GetMaximum(): 
#                max_heigh = Hists[h].GetMaximum()
#                max_hist = h
#            if min_heigh > Hists[h].GetMinimum(): 
#                min_heigh = Hists[h].GetMinimum()
#
#    gap = 0.1*(max_heigh-min_heigh)
#    if not y_max==0. and not y_min==0.:
#        Hists[0].SetMaximum(y_max)
#        Hists[0].SetMinimum(y_min)
#        Hists[0].DrawNormalized("E")
#    else:
#        if not logy:
#            Hists[max_hist].SetMaximum(max_heigh+gap)
#            Hists[max_hist].SetMinimum(min_heigh-gap)
#        Hists[max_hist].DrawNormalized("E")
#
#    for h in range(0,len(Hists)):
#        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
#        #if Hists[h]!=0:
#        Hists[h].SetLineColor(colors[h])
#        Hists[h].SetLineWidth(2)
#        Hists[h].DrawNormalized("E same")
#
#    pad3.cd()
#    legend = ROOT.TLegend(0.5,0.1,0.9,0.9)
#    legend.SetTextFont(42)
#    legend.SetBorderSize(0)
#    legend.SetTextSize(0.1)
#    legend.SetFillColor(0)
#    legend.SetFillStyle(0)
#    legend.SetLineColor(0)
#    legend.Clear()
#    for h in range(0,len(Hists)):
#        if Hists[h]!=0:
#            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
#    legend.Draw("SAME")
#
#    if logx: 
#        pad1.SetLogx()
#
#    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def GetCrabFlux():

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    func_crab.SetParameters(37.5,-2.467,-0.16)

    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1
    count_crab_total = 0.
    for binx in range(0,Hist_EffArea_Sum.GetNbinsX()):
        if binx+1<bin_lower: continue
        if binx+1>bin_upper: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        flux_crab_this_energy = func_crab.Eval(Hist_EffArea_Sum.GetBinCenter(binx+1))
        count_crab_total += flux_crab_this_energy*Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE
    return count_crab_total

def GetBackgroundFluxInterval(lower_cut,upper_cut):

    bin_lower = Hist_OnData_Energy.FindBin(lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(upper_cut)-1
    flux_bkgd_total = 0.
    for binx in range(0,Hist_EffArea_Sum.GetNbinsX()):
        if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
        if binx+1<bin_lower: continue
        if binx+1>bin_upper: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        bkg_this_energy, bkg_err_this_energy = IntegralAndError(Hist_OnBkgd_Energy_Sum,binx+1,binx+2)
        flux_bkgd_this_energy = bkg_this_energy
        flux_bkgd_total += flux_bkgd_this_energy
    return flux_bkgd_total

def MakeCrabUnitSpectrumPlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    n_roi = len(Hist_data)
    pads = []

    nth_pad = 0
    pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,1-0.7,1,1)]
    pads[nth_pad].SetBorderMode(1)
    pads[nth_pad].SetTopMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.0)
    pads[nth_pad].SetTopMargin(0.2)
    pads[nth_pad].Draw("same")
    nth_pad = 1
    pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,0,1,1-0.7)]
    pads[nth_pad].SetBorderMode(1)
    pads[nth_pad].SetTopMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.0)
    pads[nth_pad].SetBottomMargin(0.2)
    pads[nth_pad].Draw("same")

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    Func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    Func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_data_energy = []
    Hist_bkgd_energy = []
    Hist_ratio_energy = []
    for roi in range(0,n_roi):

        Hist_data_energy += [ROOT.TH1D("Hist_data_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            Hist_data_energy[roi].SetBinContent(binx+1,Hist_data[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_data_energy[roi].GetNbinsX()):
            Hist_data_energy[roi].SetBinError(binx+1,pow(Hist_data_energy[roi].GetBinContent(binx+1),0.5))

        Hist_bkgd_energy += [ROOT.TH1D("Hist_bkgd_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            Hist_bkgd_energy[roi].SetBinContent(binx+1,Hist_bkgd[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_energy[roi].GetNbinsX()):
            Hist_bkgd_energy[roi].SetBinError(binx+1,Hist_bkgd_energy[roi].GetBinContent(binx+1)*Syst_MDM)

        Hist_ratio_energy += [ROOT.TH1D("Hist_ratio_energy_%s"%(roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        Hist_ratio_energy[roi].Add(Hist_data_energy[roi])
        Hist_ratio_energy[roi].Add(Hist_bkgd_energy[roi],-1.)
        Hist_ratio_energy[roi].Divide(Hist_bkgd_energy[roi])
        for binx in range(0,Hist_data_energy[roi].GetNbinsX()):
            lower_cut = Hist_data[roi].GetBinLowEdge(binx+1)
            upper_cut = Hist_data[roi].GetBinLowEdge(binx+2)
            deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
            bkgd_flux = 0.
            if Hist_EffArea_Sum.GetBinContent(binx+1)!=0.:
                bkgd_flux = Hist_bkgd_energy[roi].GetBinContent(binx+1)/(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
            crab_flux = func_crab.Eval(Hist_EffArea_Sum.GetBinCenter(binx+1))
            Hist_ratio_energy[roi].SetBinContent(binx+1,Hist_ratio_energy[roi].GetBinContent(binx+1)*bkgd_flux/crab_flux)
            Hist_ratio_energy[roi].SetBinError(binx+1,Hist_ratio_energy[roi].GetBinError(binx+1)*bkgd_flux/crab_flux)

        pads[0].cd()
        Hist_ratio_energy[roi].SetLineColor(1)
        Hist_ratio_energy[roi].SetLineWidth(3)
        Hist_ratio_energy[roi].SetMinimum(-0.1)
        Hist_ratio_energy[roi].Draw("E")

        pads[1].cd()
        Hist_data_energy[roi].GetXaxis().SetTitle("MJD")
        Hist_data_energy[roi].SetLineColor(1)
        Hist_data_energy[roi].SetLineWidth(3)
        Hist_data_energy[roi].Draw("E")
        fill_color = 38
        stack = ROOT.THStack("stack", "")
        set_histStyle( Hist_bkgd_energy[roi] , fill_color)
        stack.Add( Hist_bkgd_energy[roi] )
        stack.Draw("hist same")
        Hist_data_energy[roi].Draw("E same")

        pads[0].SetLogx()
        pads[1].SetLogx()
        pads[1].SetLogy()

        c_both.SaveAs('output_plots/%s_RoI%s_%s.png'%(name,roi,selection_tag))

def MakeSpectrumBackgroundE27(hist_data,hist_bkgd,legends,title,name,syst):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    Hist_EffArea_tmp = Hist_EffArea_Sum.Clone()

    Hist_Flux = []
    Hist_Flux += [ROOT.TH1D("Hist_Flux_Data","",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_Flux += [ROOT.TH1D("Hist_Flux_Bkgd","",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_Flux[0].Add(hist_data[0])
    Hist_Flux[1].Add(hist_bkgd[0])
    for binx in range(0,Hist_Flux[0].GetNbinsX()):
        if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
        deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
        avgE = 0.5*(energy_fine_bin[binx+1]+energy_fine_bin[binx])/1000.
        scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
        Hist_Flux[0].SetBinContent(binx+1,Hist_Flux[0].GetBinContent(binx+1)*scale)
        Hist_Flux[0].SetBinError(binx+1,Hist_Flux[0].GetBinError(binx+1)*scale)
        Hist_Flux[1].SetBinContent(binx+1,Hist_Flux[1].GetBinContent(binx+1)*scale)
        Hist_Flux[1].SetBinError(binx+1,Hist_Flux[1].GetBinError(binx+1)*scale)

    Hist_Invisible = Hist_Flux[0].Clone()
    Hist_Invisible.SetLineColor(0)
    Hist_Invisible.GetXaxis().SetTitleOffset(0.8)
    Hist_Invisible.GetXaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetTitleOffset(1.2)
    Hist_Invisible.GetYaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetTitle('energy [GeV]')
    Hist_Invisible.GetYaxis().SetTitle('flux [counts/GeV/s/cm2]')
    Hist_Invisible.Draw()

    func_source = []
    func_source += [ROOT.TF1("func_source_data","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)]
    func_source[0].SetParameters(37.5,-4.)
    Hist_Flux[0].Fit("func_source_data","N")
    Hist_Flux[0].SetLineColor(1)
    Hist_Flux[0].Draw("E same")
    func_source[0].SetLineColor(1)
    func_source[0].Draw("E same")
    func_source += [ROOT.TF1("func_source_bkgd","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)]
    func_source[1].SetParameters(37.5,-4.)
    Hist_Flux[1].Fit("func_source_bkgd","N")
    Hist_Flux[1].SetLineColor(2)
    Hist_Flux[1].Draw("E same")
    func_source[1].SetLineColor(2)
    func_source[1].Draw("E same")

    #func_crab.SetLineColor(2)
    #func_crab.Draw("same")

    pad3.cd()

    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_Flux[0],"data, amp = %0.2f #pm %0.2f, index = %0.2f #pm %0.2f"%(func_source[0].GetParameter(0),func_source[0].GetParError(0),func_source[0].GetParameter(1),func_source[0].GetParError(1)),"pl")
    legend.AddEntry(Hist_Flux[1],"model, amp = %0.2f #pm %0.2f, index = %0.2f #pm %0.2f"%(func_source[1].GetParameter(0),func_source[1].GetParError(0),func_source[1].GetParameter(1),func_source[1].GetParError(1)),"pl")
    legend.Draw("SAME")

    #lumilab1 = ROOT.TLatex(0.15,0.80,'Exposure %.1f hrs'%(exposure_hours) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    pad1.SetLogy()
    pad1.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeSpectrumInNonCrabUnitE2(hist_data,hist_bkgd,legends,title,name,syst):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    Hist_EffArea_tmp = Hist_EffArea_Sum.Clone()

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    func_crab = ROOT.TF1("func_crab","pow(x/1000.,2)*[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    func_crab.SetParameters(37.5,-2.467,-0.16)

    func_source = []

    Hist_Flux = []
    for nth_roi in range(0,len(hist_data)):
        Hist_Flux += [ROOT.TH1D("Hist_Flux_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        Hist_Flux[nth_roi].Add(hist_data[nth_roi])
        Hist_Flux[nth_roi].Add(hist_bkgd[nth_roi],-1.)
        for binx in range(0,Hist_Flux[nth_roi].GetNbinsX()):
            if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
            deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
            avgE = 0.5*(energy_fine_bin[binx+1]+energy_fine_bin[binx])/1000.
            scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
            Hist_Flux[nth_roi].SetBinContent(binx+1,avgE*avgE*Hist_Flux[nth_roi].GetBinContent(binx+1)*scale)
            Hist_Flux[nth_roi].SetBinError(binx+1,avgE*avgE*Hist_Flux[nth_roi].GetBinError(binx+1)*scale)

    Hist_Invisible = Hist_Flux[0].Clone()
    ymax = 0.
    ymin = 0.
    for nth_roi in range(0,len(hist_data)):
        ymax = max(Hist_Flux[nth_roi].GetMaximum(),Hist_Invisible.GetBinContent(1))
        ymin = min(Hist_Flux[nth_roi].GetMinimum(),Hist_Invisible.GetBinContent(2))
        Hist_Invisible.SetBinContent(1,ymax)
        Hist_Invisible.SetBinContent(2,ymin)
    Hist_Invisible.SetLineColor(0)
    Hist_Invisible.GetXaxis().SetTitleOffset(0.8)
    Hist_Invisible.GetXaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetTitleOffset(1.2)
    Hist_Invisible.GetYaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetTitle('energy [GeV]')
    Hist_Invisible.GetYaxis().SetTitle('E^{2} #times flux [counts/GeV/s/cm2]')
    Hist_Invisible.Draw()

    for nth_roi in range(0,len(hist_data)):
        Hist_Flux[nth_roi].SetLineColor(nth_roi+1)
        func_source += [ROOT.TF1("func_source_%s"%(nth_roi),"pow(x/1000.,2)*[0]*pow(10,-12)*pow(x/1000.,[1])", 150, 10000)]
        func_source[nth_roi].SetParameters(37.5,-2.)
        Hist_Flux[nth_roi].Fit("func_source_%s"%(nth_roi),"N")
        Hist_Flux[nth_roi].Draw("E same")
        func_source[nth_roi].SetLineColor(nth_roi+1)
        func_source[nth_roi].Draw("E same")

    #func_crab.SetLineColor(2)
    #func_crab.Draw("same")

    # MGRO J1908
    if 'MGRO_J1908' in name:
        func_1908 = ROOT.TF1("func_1908","pow(x/1000.,2)*[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_1908.SetParameters(4.23,-2.2)
        func_1908.SetLineColor(4)
        func_1908.Draw("same")
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    if 'IC443' in name:
        print 'Compare to official IC 443 flux...'
        func_ic443 = ROOT.TF1("func_ic443","pow(x/1000.,2)*[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_ic443.SetParameters(0.838,-2.99)
        func_ic443.SetLineColor(4)
        func_ic443.Draw("same")
    # 1ES 1218 https://arxiv.org/pdf/0810.0301.pdf
    if '1ES1218' in name:
        func_1218 = ROOT.TF1("func_1218","pow(x/1000.,2)*[0]*pow(10,-12)*pow(x/500.,[1])", 200, 4000)
        func_1218.SetParameters(7.5,-3.08)
        func_1218.SetLineColor(4)
        func_1218.Draw("same")


    pad3.cd()

    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for nth_roi in range(0,len(hist_data)):
        #legend.AddEntry(Hist_Flux[nth_roi],"%s, A #times 10^{-12} #times E^{x}, A = %0.2f, x = %0.2f"%(legends[nth_roi],func_source[nth_roi].GetParameter(0),func_source[nth_roi].GetParameter(1)),"pl")
        legend.AddEntry(Hist_Flux[nth_roi],"%s, index = %0.2f #pm %0.2f"%(legends[nth_roi],func_source[nth_roi].GetParameter(1),func_source[nth_roi].GetParError(1)),"pl")
    if 'IC443' in name:
        legend.AddEntry(func_ic443,"flux from arXiv:0905.3291","pl")
    if 'MGRO_J1908' in name:
        legend.AddEntry(func_1908,"flux from arXiv:1404.7185","pl")
    legend.Draw("SAME")

    #lumilab1 = ROOT.TLatex(0.15,0.80,'Exposure %.1f hrs'%(exposure_hours) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    pad1.SetLogy()
    pad1.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeSpectrumInNonCrabUnit(hist_data,hist_bkgd,legends,title,name,syst):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    Hist_EffArea_tmp = Hist_EffArea_Sum.Clone()

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    func_crab.SetParameters(37.5,-2.467,-0.16)

    func_source = []

    Hist_Flux = []
    for nth_roi in range(0,len(hist_data)):
        Hist_Flux += [ROOT.TH1D("Hist_Flux_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
        Hist_Flux[nth_roi].Add(hist_data[nth_roi])
        Hist_Flux[nth_roi].Add(hist_bkgd[nth_roi],-1.)
        for binx in range(0,Hist_Flux[nth_roi].GetNbinsX()):
            if Hist_EffArea_Sum.GetBinContent(binx+1)==0.: continue
            deltaE = (energy_fine_bin[binx+1]-energy_fine_bin[binx])/1000.
            scale = 1./(Hist_EffArea_Sum.GetBinContent(binx+1)*10000.*deltaE)
            Hist_Flux[nth_roi].SetBinContent(binx+1,Hist_Flux[nth_roi].GetBinContent(binx+1)*scale)
            Hist_Flux[nth_roi].SetBinError(binx+1,Hist_Flux[nth_roi].GetBinError(binx+1)*scale)

    Hist_Invisible = Hist_Flux[0].Clone()
    ymax = 0.
    ymin = 0.
    for nth_roi in range(0,len(hist_data)):
        ymax = max(Hist_Flux[nth_roi].GetMaximum(),Hist_Invisible.GetBinContent(1))
        ymin = min(Hist_Flux[nth_roi].GetMinimum(),Hist_Invisible.GetBinContent(2))
        Hist_Invisible.SetBinContent(1,ymax)
        Hist_Invisible.SetBinContent(2,ymin)
    Hist_Invisible.SetLineColor(0)
    Hist_Invisible.GetXaxis().SetTitleOffset(0.8)
    Hist_Invisible.GetXaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetLabelSize(0.04)
    Hist_Invisible.GetYaxis().SetTitleOffset(1.2)
    Hist_Invisible.GetYaxis().SetTitleSize(0.06)
    Hist_Invisible.GetXaxis().SetTitle('energy [GeV]')
    Hist_Invisible.GetYaxis().SetTitle('flux [counts/GeV/s/cm2]')
    Hist_Invisible.Draw()

    for nth_roi in range(0,len(hist_data)):
        Hist_Flux[nth_roi].SetLineColor(nth_roi+1)
        func_source += [ROOT.TF1("func_source_%s"%(nth_roi),"[0]*pow(10,-12)*pow(x/1000.,[1])", 150, 10000)]
        func_source[nth_roi].SetParameters(37.5,-2.)
        Hist_Flux[nth_roi].Fit("func_source_%s"%(nth_roi),"N")
        Hist_Flux[nth_roi].Draw("E same")
        func_source[nth_roi].SetLineColor(nth_roi+1)
        func_source[nth_roi].Draw("E same")

    func_crab.SetLineColor(2)
    func_crab.Draw("same")

    # MGRO J1908
    if 'MGRO_J1908' in name:
        func_1908 = ROOT.TF1("func_1908","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_1908.SetParameters(4.23,-2.2)
        func_1908.SetLineColor(4)
        func_1908.Draw("same")
    # IC 443 https://arxiv.org/pdf/0905.3291.pdf
    if 'IC443' in name:
        print 'Compare to official IC 443 flux...'
        func_ic443 = ROOT.TF1("func_ic443","[0]*pow(10,-12)*pow(x/1000.,[1])", 200, 4000)
        func_ic443.SetParameters(0.838,-2.99)
        func_ic443.SetLineColor(4)
        func_ic443.Draw("same")
    # 1ES 1218 https://arxiv.org/pdf/0810.0301.pdf
    if '1ES1218' in name:
        func_1218 = ROOT.TF1("func_1218","[0]*pow(10,-12)*pow(x/500.,[1])", 200, 4000)
        func_1218.SetParameters(7.5,-3.08)
        func_1218.SetLineColor(4)
        func_1218.Draw("same")


    pad3.cd()

    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for nth_roi in range(0,len(hist_data)):
        #legend.AddEntry(Hist_Flux[nth_roi],"%s, A #times 10^{-12} #times E^{x}, A = %0.2f, x = %0.2f"%(legends[nth_roi],func_source[nth_roi].GetParameter(0),func_source[nth_roi].GetParameter(1)),"pl")
        legend.AddEntry(Hist_Flux[nth_roi],"%s, index = %0.2f #pm %0.2f"%(legends[nth_roi],func_source[nth_roi].GetParameter(1),func_source[nth_roi].GetParError(1)),"pl")
    if 'IC443' in name:
        legend.AddEntry(func_ic443,"flux from arXiv:0905.3291","pl")
    if 'MGRO_J1908' in name:
        legend.AddEntry(func_1908,"flux from arXiv:1404.7185","pl")
    legend.Draw("SAME")

    #lumilab1 = ROOT.TLatex(0.15,0.80,'Exposure %.1f hrs'%(exposure_hours) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.15)
    #lumilab1.Draw()

    pad1.SetLogy()
    pad1.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeSignalToBkgdRatioPlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    n_roi = len(Hist_data)

    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.3)
    pad1.SetTopMargin(0.1)
    pad1.SetBorderMode(1)
    pad1.SetGrid()
    pad1.Draw()

    pad1.cd()

    time = Time(MJD_Start, format='mjd')
    time.format = 'decimalyear'
    year_start = time.value
    time = Time(MJD_End, format='mjd')
    time.format = 'decimalyear'
    year_end = time.value

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    Func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    Func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    Hist_ratio_mjd = []
    IncValues = ROOT.TF1( "IncValues", "x", year_start , year_end )
    raLowerAxis = []
    days_per_bin = 80.
    for roi in range(0,n_roi):

        Hist_data_mjd += [ROOT.TH1D("Hist_data_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_data_mjd[roi].Fill(mjd,Hist_data[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
            Hist_data_mjd[roi].SetBinError(binx+1,pow(Hist_data_mjd[roi].GetBinContent(binx+1),0.5))

        Hist_bkgd_mjd += [ROOT.TH1D("Hist_bkgd_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_bkgd_mjd[roi].Fill(mjd,Hist_bkgd[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_mjd[roi].GetNbinsX()):
            Hist_bkgd_mjd[roi].SetBinError(binx+1,Hist_bkgd_mjd[roi].GetBinContent(binx+1)*Syst_MDM)

    #for roi in range(0,n_roi):
    #    for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
    #        content = Hist_data_mjd[2].GetBinContent(binx+1)
    #        error = Hist_data_mjd[2].GetBinError(binx+1)
    #        if content<30.:
    #            Hist_data_mjd[roi].SetBinContent(binx+1,0.)
    #            Hist_data_mjd[roi].SetBinError(binx+1,0.)
    #            Hist_bkgd_mjd[roi].SetBinContent(binx+1,0.)
    #            Hist_bkgd_mjd[roi].SetBinError(binx+1,0.)


    for roi in range(0,n_roi):
        Hist_ratio_mjd += [ROOT.TH1D("Hist_ratio_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        Hist_ratio_mjd[roi].Add(Hist_data_mjd[roi])
        Hist_ratio_mjd[roi].Add(Hist_bkgd_mjd[roi],-1.)
        Hist_ratio_mjd[roi].Divide(Hist_bkgd_mjd[roi])

    Hist_ratio_mjd[0].Draw("E same")
    for roi in range(0,n_roi):

        Hist_ratio_mjd[roi].SetLineColor(roi+1)
        Hist_ratio_mjd[roi].SetLineWidth(3)
        Hist_ratio_mjd[roi].Draw("E same")
        c_both.Update()
        x1 = ROOT.gPad.GetUxmin()
        x2 = ROOT.gPad.GetUxmax()
        y1 = ROOT.gPad.GetUymin()
        y2 = ROOT.gPad.GetUymax()
        raLowerAxis += [ROOT.TGaxis( x1, y2, x2, y2,"IncValues", 510, "-")]
        raLowerAxis[roi].SetLabelSize(Hist_ratio_mjd[roi].GetXaxis().GetLabelSize())
        raLowerAxis[roi].SetTitle("Year")
        raLowerAxis[roi].Draw()

    c_both.SaveAs('output_plots/%s_RoI%s_%s.png'%(name,roi,selection_tag))


def MakeLightCurvePlot(Hist_data,Hist_bkgd,legends,colors,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    n_roi = len(Hist_data)
    pads = []

    #step_size = 1./float(n_roi)
    #for roi in range(0,n_roi):
    #    pads += [ROOT.TPad("pad%s"%(roi),"",0,1-(roi+1)*step_size,1,1-(roi)*step_size)]
    #    pads[roi].SetBorderMode(1)
    #    pads[roi].SetTopMargin(0.0)
    #    pads[roi].SetBottomMargin(0.0)
    #    if roi==0:
    #        pads[roi].SetTopMargin(0.1/step_size)
    #    if roi==n_roi-1:
    #        pads[roi].SetBottomMargin(0.1/step_size)
    #    pads[roi].Draw("same")

    step_size = 0.5
    for nth_pad in range(0,2):
        pads += [ROOT.TPad("pad%s"%(nth_pad),"",0,1-(nth_pad+1)*step_size,1,1-(nth_pad)*step_size)]
        pads[nth_pad].SetBorderMode(1)
        pads[nth_pad].SetTopMargin(0.0)
        pads[nth_pad].SetBottomMargin(0.0)
        if nth_pad==0:
            pads[nth_pad].SetTopMargin(0.1/step_size)
        if nth_pad==1:
            pads[nth_pad].SetBottomMargin(0.1/step_size)
        pads[nth_pad].Draw("same")

    time = Time(MJD_Start, format='mjd')
    time.format = 'decimalyear'
    year_start = time.value
    time = Time(MJD_End, format='mjd')
    time.format = 'decimalyear'
    year_end = time.value

    # Crab https://arxiv.org/pdf/1508.06442.pdf
    Func_crab = ROOT.TF1("func_crab","[0]*pow(10,-12)*pow(x/1000.,[1]+[2]*log(x/1000.))", 100, 10000)
    Func_crab.SetParameters(37.5,-2.467,-0.16)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    Hist_ratio_mjd = []
    IncValues = ROOT.TF1( "IncValues", "x", year_start , year_end )
    raLowerAxis = []
    days_per_bin = 90.
    for roi in range(0,n_roi):

        Hist_data_mjd += [ROOT.TH1D("Hist_data_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_data_mjd[roi].Fill(mjd,Hist_data[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
            Hist_data_mjd[roi].SetBinError(binx+1,pow(Hist_data_mjd[roi].GetBinContent(binx+1),0.5))

        Hist_bkgd_mjd += [ROOT.TH1D("Hist_bkgd_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        for binx in range(0,Hist_data[roi].GetNbinsX()):
            mjd = Hist_data[roi].GetBinCenter(binx+1)
            Hist_bkgd_mjd[roi].Fill(mjd,Hist_bkgd[roi].GetBinContent(binx+1))
        for binx in range(0,Hist_bkgd_mjd[roi].GetNbinsX()):
            Hist_bkgd_mjd[roi].SetBinError(binx+1,Hist_bkgd_mjd[roi].GetBinContent(binx+1)*Syst_MDM)


        Hist_ratio_mjd += [ROOT.TH1D("Hist_ratio_mjd_%s"%(roi),"",int((MJD_End+1-MJD_Start)/days_per_bin),MJD_Start,MJD_End+1)]
        Hist_ratio_mjd[roi].Add(Hist_data_mjd[roi])
        Hist_ratio_mjd[roi].Add(Hist_bkgd_mjd[roi],-1.)
        Hist_ratio_mjd[roi].Divide(Hist_bkgd_mjd[roi])

        ratio_threshold = 1.0
        for binx in range(0,Hist_data_mjd[roi].GetNbinsX()):
            if Hist_ratio_mjd[roi].GetBinError(binx+1)==0.: continue
            significance = Hist_ratio_mjd[roi].GetBinContent(binx+1)/Hist_ratio_mjd[roi].GetBinError(binx+1)
            if significance<4.0: continue
            if Hist_ratio_mjd[roi].GetBinContent(binx+1)>ratio_threshold:
                print '>%s Crab, bin edge: %s, %s'%(ratio_threshold,Hist_ratio_mjd[roi].GetBinLowEdge(binx+1),Hist_ratio_mjd[roi].GetBinLowEdge(binx+2))

        pads[0].cd()
        Hist_ratio_mjd[roi].SetLineColor(1)
        Hist_ratio_mjd[roi].SetLineWidth(3)
        Hist_ratio_mjd[roi].Draw("E")
        c_both.Update()
        x1 = ROOT.gPad.GetUxmin()
        x2 = ROOT.gPad.GetUxmax()
        y1 = ROOT.gPad.GetUymin()
        y2 = ROOT.gPad.GetUymax()
        raLowerAxis += [ROOT.TGaxis( x1, y2, x2, y2,"IncValues", 510, "-")]
        raLowerAxis[roi].SetLabelSize(Hist_ratio_mjd[roi].GetXaxis().GetLabelSize())
        raLowerAxis[roi].SetTitle("Year")
        raLowerAxis[roi].Draw()

        pads[1].cd()
        Hist_data_mjd[roi].GetXaxis().SetTitle("MJD")
        Hist_data_mjd[roi].SetLineColor(1)
        Hist_data_mjd[roi].SetLineWidth(3)
        Hist_data_mjd[roi].Draw("E")
        fill_color = 38
        stack = ROOT.THStack("stack", "")
        set_histStyle( Hist_bkgd_mjd[roi] , fill_color)
        stack.Add( Hist_bkgd_mjd[roi] )
        stack.Draw("hist same")
        Hist_data_mjd[roi].Draw("E same")

        c_both.SaveAs('output_plots/%s_RoI%s_%s.png'%(name,roi,selection_tag))

def MakeStarEffectPlot(Hist_data,Hist_bkgd,legends,colors,name,range_lower,range_upper):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad = ROOT.TPad("pad","pad",0,0,1,1)
    pad.SetBottomMargin(0.1)
    pad.SetTopMargin(0.1)
    pad.SetLeftMargin(0.15)
    pad.SetBorderMode(1)
    pad.Draw()
    pad.cd()

    norm_bin_low_target = Hist_data[0].FindBin(range_lower)
    norm_bin_up_target = Hist_data[0].FindBin(range_upper)-1

    n_roi = 0
    star_brightness, star_effect, star_brightness_err, star_effect_err = array( 'd' ), array( 'd' ) , array( 'd' ) , array( 'd' )
    for roi in range(0,len(Hist_data)):
        if not 'b-mag' in legends[roi]: continue
        star_brightness.append(float(legends[roi].strip('b-mag ')))
        star_brightness_err.append(0.)
        data_count = Hist_data[roi].Integral(norm_bin_low_target,norm_bin_up_target)
        bkgd_count = Hist_bkgd[roi].Integral(norm_bin_low_target,norm_bin_up_target)
        if (bkgd_count>0.):
            star_effect.append( (data_count-bkgd_count)/bkgd_count )
            star_effect_err.append( pow(data_count,0.5)/bkgd_count )
            n_roi += 1

    gr = ROOT.TGraphErrors( n_roi, star_brightness, star_effect , star_brightness_err, star_effect_err )
    gr.SetLineColor( 1 )
    gr.SetLineWidth( 4 )
    gr.SetMarkerColor( 1 )
    gr.SetMarkerStyle( 21 )
    gr.GetXaxis().SetTitle( 'b-mag' )
    gr.GetYaxis().SetTitle( 'star syst.' )
    gr.Draw( 'AP' )

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def MakeChi2Plot(Hists,legends,colors,stack_it,title,name,doSum,range_lower,range_upper,syst):
    
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

    norm_bin_low_target = 1
    norm_bin_up_target = Hists[0].GetNbinsX()
    if 'MSCW' in name or 'MSCL' in name:
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
    #Hists[max_hist].Draw("E")
    Hists[0].Draw("E")

    Hist_Bkgd = Hists[1].Clone()
    Hist_Bkgd.Reset()
    Hist_Syst = Hists[1].Clone()
    Hist_Syst.Reset()
    found_syst_hist = False
    for h in range(0,len(Hists)):
        if stack_it[h] or 'syst' in legends[h]:
            Hist_Syst.Add(Hists[h])
        if 'syst' in legends[h]:
            found_syst_hist = True
        if colors[h]!=0 and colors[h]!=1:
            Hist_Bkgd.Add(Hists[h].Clone())
    #if not found_syst_hist:
    #    for binx in range(0,Hist_Syst.GetNbinsX()):
    #        syst_err = Syst_MDM*Hist_Bkgd.GetBinContent(binx+1)
    #        stat_err = Hist_Bkgd.GetBinError(binx+1)
    #        Hist_Syst.SetBinError(binx+1,pow(syst_err*syst_err+stat_err*stat_err,0.5))

    fill_color = [0,0,46,0,38,30]
    if doSum:
        stack = ROOT.THStack("stack", "")
        for h in range(1,len(Hists)):
            if stack_it[h]:
                set_histStyle( Hists[h] , fill_color[colors[h]])
                stack.Add( Hists[h] )
        stack.Draw("hist same")
    if found_syst_hist:
        Hist_Syst.SetFillColor(1)
        Hist_Syst.SetFillStyle(3004)
        Hist_Syst.SetMarkerSize(0)
        Hist_Syst.Draw("e2 same")

    #for h in range(0,len(Hists)):
    #    if colors[h]==0: continue
    #    Hists[h].SetLineWidth(3)
    #    Hists[h].Draw("E same")
    Hists[0].SetLineWidth(3)
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
    lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
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
    Hist_Band = Hist_Bkgd.Clone()
    if found_syst_hist:
        Hist_Band.Reset()
        Hist_Band.Add(Hist_Syst)
    Hist_Band.Divide(Hist_Bkgd)
    Hist_Band.SetFillColor(1)
    if not found_syst_hist:
        Hist_Band.SetFillColor(0)
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
    if 'MJD' in title:
        Hist_Band.SetMaximum(3.0)
        Hist_Band.SetMinimum(0.)
    Hist_Band.Draw("e2")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    Hist_Ratio = Hists[0].Clone()
    Hist_Ratio.Divide(Hist_Bkgd)
    Hist_Ratio.SetLineWidth(2)
    Hist_Ratio.SetLineColor(1)
    Hist_Ratio.Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,selection_tag))

def PlotsStackedHistograms(tag):

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjX_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Skymap_ProjX_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnBkgd_Skymap_Syst_ProjX_Sum]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Skymap_ProjX_%s'%(tag)
    title = 'RA'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Skymap_ProjY_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Skymap_ProjY_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnBkgd_Skymap_Syst_ProjY_Sum]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Skymap_ProjY_%s'%(tag)
    title = 'Dec'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0,1,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_MSCW]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_wGamma_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnGamma_MSCW_Sum]
    legends += ['#gamma fixed at 10% of bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Unblind_wGamma_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_wGamma_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_OnGamma_MSCL_Sum]
    legends += ['#gamma fixed at 10% of bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Unblind_wGamma_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_woGamma_MSCW_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Unblind_woGamma_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Unblind_woGamma_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Unblind_woGamma_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCW_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_MSCW_Sum]
    legends += ['init. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_MSCW_Dark_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCW_plot_lower,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_MSCL_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_MSCL]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_MSCL_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_MSCL_Sum]
    legends += ['init. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_MSCL_Dark_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,MSCL_plot_lower,MSCL_blind_cut,-1)

    if energy_bin[energy_bin_cut_low]>=2000.:
        Hist_OnData_Theta2_Sum.Rebin(2)
        Hist_OnBkgd_Theta2_Sum.Rebin(2)
        Hist_OnRFoV_Theta2_Sum.Rebin(2)
        Hist_OnDark_Theta2_Sum.Rebin(2)
        Hist_SystErr_Theta2.Rebin(2)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hist_OnBkgd_Yoff_Raw_Sum.Divide(Hist_OnData_Yoff_Sum)
    Hist_OnBkgd_Yoff_Sum.Divide(Hist_OnData_Yoff_Sum)
    Hist_OnData_Yoff_Sum.Divide(Hist_OnData_Yoff_Sum)
    Hists += [Hist_OnData_Yoff_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Yoff_Raw_Sum]
    legends += ['uncorrected']
    colors += [2]
    Hists += [Hist_OnBkgd_Yoff_Sum]
    legends += ['corrected']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Yoff_MDM_%s'%(tag)
    title = 'Y off'
    MakeComparisonPlot(Hists,legends,colors,title,'counts',plotname,0.5,1.5,False,False,False)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Theta2]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Theta2_MDM_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)
    plotname = 'SignificanceCurve_Theta2_MDM_%s'%(tag)
    title = 'analysis cut on #theta^{2}'
    MakeAccumulatedSignificancePlot(Hist_OnData_Theta2_Sum,Hist_OnBkgd_Theta2_Sum,Hist_SystErr_Theta2,0,title,'significance',plotname)
    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnRFoV_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Theta2]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Theta2_RFoV_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)
    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Theta2_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Theta2_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    Hists += [Hist_SystErr_Theta2]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Theta2_Dark_%s'%(tag)
    title = 'squared angle from source location #theta^{2}'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

    #if energy_bin[energy_bin_cut_low]>=300.:
    #    for nth_roi in range(0,len(roi_ra)):
    #        Hist_OnData_RoI_Theta2_Sum[nth_roi].Rebin(2)
    #        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Rebin(2)
    for nth_roi in range(0,len(roi_ra)):
        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Theta2_Sum[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_Theta2_MDM_%s'%(nth_roi,tag)
        title = 'squared angle from source location #theta^{2}'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,1.,-1)

        Hists = []
        legends = []
        colors = []
        stack_it = []
        Hists += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        legends += ['obs. data (%s)'%(roi_name[nth_roi])]
        colors += [1]
        stack_it += [False]
        Hists += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        legends += ['predict. bkg.']
        colors += [4]
        stack_it += [True]
        plotname = 'Stack_RoI%s_Energy_MDM_%s'%(nth_roi,tag)
        title = 'energy [GeV]'
        MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hist_data = []
    Hist_bkgd = []
    legends = []
    #for nth_roi in range(0,len(roi_ra)):
    for nth_roi in range(1,len(roi_ra)):
        if 'b-mag' in roi_name[nth_roi]: continue
        Hist_data += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
        Hist_bkgd += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
        legends += ['%s'%(roi_name[nth_roi])]
    title = 'energy [GeV]'
    plotname = 'Flux_MDM_%s'%(tag)
    MakeSpectrumInNonCrabUnit(Hist_data,Hist_bkgd,legends,title,plotname,-1)
    plotname = 'FluxE2_MDM_%s'%(tag)
    MakeSpectrumInNonCrabUnitE2(Hist_data,Hist_bkgd,legends,title,plotname,-1)
    Hist_data = []
    Hist_bkgd = []
    Hist_data += [Hist_OnData_RoI_Energy_Sum[1]]
    Hist_bkgd += [Hist_OnBkgd_RoI_Energy_Sum[1]]
    title = 'energy [GeV]'
    plotname = 'FluxE27_MDM_%s'%(tag)
    MakeSpectrumBackgroundE27(Hist_data,Hist_bkgd,legends,title,plotname,-1)

    Hist_data_mjd = []
    Hist_bkgd_mjd = []
    legends = []
    colors = []
    for nth_roi in range(0,len(roi_ra)):
        legends += ['%s'%(roi_name[nth_roi])]
        colors += [nth_roi+1]
        Hist_data_mjd += [Hist_OnData_RoI_MJD_Sum[nth_roi]]
        Hist_bkgd_mjd += [Hist_OnBkgd_RoI_MJD_Sum[nth_roi]]
    plotname = 'LightCurve_RoI_Year_MDM_%s'%(tag)
    title = 'Year'
    MakeLightCurvePlot(Hist_data_mjd,Hist_bkgd_mjd,legends,colors,title,plotname)
    plotname = 'LightCurve_Ratio_Year_MDM_%s'%(tag)
    MakeSignalToBkgdRatioPlot(Hist_data_mjd,Hist_bkgd_mjd,legends,colors,title,plotname)

    #Hist_data_energy = []
    #Hist_bkgd_energy = []
    #legends = []
    #colors = []
    #for nth_roi in range(0,len(roi_ra)):
    #    legends += ['%s'%(roi_name[nth_roi])]
    #    colors += [nth_roi+1]
    #    Hist_data_energy += [Hist_OnData_RoI_Energy_Sum[nth_roi]]
    #    Hist_bkgd_energy += [Hist_OnBkgd_RoI_Energy_Sum[nth_roi]]
    #plotname = 'Spectrum_RoI%s_CrabUnit_MDM_%s'%(nth_roi,tag)
    #title = 'GeV'
    #MakeCrabUnitSpectrumPlot(Hist_data_energy,Hist_bkgd_energy,legends,colors,title,plotname)

    #Hist_data_theta2 = []
    #Hist_bkgd_theta2 = []
    #legends = []
    #colors = []
    #for nth_roi in range(0,len(roi_ra)):
    #    legends += ['%s'%(roi_name[nth_roi])]
    #    colors += [nth_roi+1]
    #    Hist_data_theta2 += [Hist_OnData_RoI_Theta2_Sum[nth_roi]]
    #    Hist_bkgd_theta2 += [Hist_OnBkgd_RoI_Theta2_Sum[nth_roi]]
    #plotname = 'StarEffect_RoI%s_MDM_%s'%(nth_roi,tag)
    #MakeStarEffectPlot(Hist_data_theta2,Hist_bkgd_theta2,legends,colors,plotname,0.,0.075)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Energy]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Energy_MDM_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnRFoV_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    Hists += [Hist_SystErr_Energy]
    legends += ['syst. error']
    colors += [0]
    stack_it += [False]
    plotname = 'Stack_Energy_RFoV_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_CamCenter_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Energy_CamCenter_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Energy_CamCenter_MDM_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Energy_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnDark_Energy_Sum]
    legends += ['predict. bkg.']
    colors += [2]
    stack_it += [True]
    plotname = 'Stack_Energy_Dark_%s'%(tag)
    title = 'energy [GeV]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)

    Hists = []
    legends = []
    colors = []
    stack_it = []
    Hists += [Hist_OnData_Zenith_Sum]
    legends += ['obs. data']
    colors += [1]
    stack_it += [False]
    Hists += [Hist_OnBkgd_Zenith_Sum]
    legends += ['predict. bkg.']
    colors += [4]
    stack_it += [True]
    plotname = 'Stack_Zenith_MDM_%s'%(tag)
    title = 'Zenith [degree]'
    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,True,0.,pow(10,4.0),-1)


    #for nth_sample in range(0,n_control_samples):

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_MSCW_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_MSCW_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_MSCW_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'MSCW'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCW_plot_lower,MSCW_blind_cut,-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_MSCL_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_MSCL_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_MSCL_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'MSCL'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,MSCL_plot_lower,MSCL_blind_cut,-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    Hists += [Hist_OffData_Energy_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_Energy_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_Energy_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'energy [GeV]'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,pow(10,4.0),-1)

    #    Hists = []
    #    legends = []
    #    colors = []
    #    stack_it = []
    #    #Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
    #    #Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].GetXaxis().SetRangeUser(0,9)
    #    Hists += [Hist_OffData_CameraFoV_Theta2_Sum[nth_sample]]
    #    legends += ['obs. data']
    #    colors += [1]
    #    stack_it += [False]
    #    Hists += [Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample]]
    #    legends += ['predict. bkg.']
    #    colors += [4]
    #    stack_it += [False]
    #    plotname = 'Stack_CameraFoV_Theta2_MDM_%s_Control%s'%(tag,nth_sample)
    #    title = 'squared angle from camera center #theta^{2}'
    #    MakeChi2Plot(Hists,legends,colors,stack_it,title,plotname,False,0.,1.,-1)

def CalculateSystErrorIndividual_v3():

    Syst_MDM_single = 0.

    binx_low_target = Hist_OnData_MSCL.FindBin(MSCL_plot_lower)
    binx_up_target = Hist_OnData_MSCL.FindBin(MSCL_plot_upper)-1
    binx_blind_target = Hist_OnData_MSCL.FindBin(MSCL_blind_cut)-1
    biny_low_target = Hist_OnData_MSCW.FindBin(MSCW_plot_lower)
    biny_up_target = Hist_OnData_MSCW.FindBin(MSCW_plot_upper)-1
    biny_blind_target = Hist_OnData_MSCW.FindBin(MSCW_blind_cut)-1

    #Total_OnBkgd = Hist2D_OnBkgd.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
    Total_OnBkgd = Hist2D_OnBkgd.Integral()
    if Total_OnBkgd==0.: return Syst_MDM_single

    Hist2D_Diff = Hist2D_OnBkgd.Clone()
    #Hist2D_Diff = Hist2D_OnDark.Clone()
    Hist2D_Diff.Add(Hist2D_OnData,-1.)
    for binx in range(0,Hist2D_Diff.GetNbinsX()):
        for biny in range(0,Hist2D_Diff.GetNbinsY()):
            if binx<=binx_blind_target and biny<=biny_blind_target: continue
            #if binx>binx_blind_target and biny>biny_blind_target: continue
            Syst_MDM_single += abs(Hist2D_Diff.GetBinContent(binx+1,biny+1))
    Syst_MDM_single = 2.*(pow(Syst_MDM_single,1.0))/Total_OnBkgd
    #Syst_MDM_single = 1.*pow(max(0.,Syst_MDM_single-Hist2D_OnData.Integral()+Total_OnBkgd),0.5)/Total_OnBkgd

    return Syst_MDM_single

def CalculateSystError_v3():

    global Syst_MDM
    Syst_MDM = 0.

    binx_low_target = Hist_OnData_MSCL.FindBin(MSCL_plot_lower)
    binx_up_target = Hist_OnData_MSCL.FindBin(MSCL_plot_upper)-1
    binx_blind_target = Hist_OnData_MSCL.FindBin(MSCL_blind_cut)-1
    biny_low_target = Hist_OnData_MSCW.FindBin(MSCW_plot_lower)
    biny_up_target = Hist_OnData_MSCW.FindBin(MSCW_plot_upper)-1
    biny_blind_target = Hist_OnData_MSCW.FindBin(MSCW_blind_cut)-1

    #Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral(binx_low_target,binx_blind_target,biny_low_target,biny_blind_target)
    Total_OnBkgd = Hist2D_OnBkgd_Sum.Integral()
    if Total_OnBkgd==0.: return

    Hist2D_Diff = Hist2D_OnBkgd_Sum.Clone()
    #Hist2D_Diff = Hist2D_OnDark_Sum.Clone()
    Hist2D_Diff.Add(Hist2D_OnData_Sum,-1.)
    for binx in range(0,Hist2D_Diff.GetNbinsX()):
        for biny in range(0,Hist2D_Diff.GetNbinsY()):
            if binx<=binx_blind_target and biny<=biny_blind_target: continue
            #if binx>binx_blind_target and biny>biny_blind_target: continue
            Syst_MDM += abs(Hist2D_Diff.GetBinContent(binx+1,biny+1))
    Syst_MDM = 2.*(pow(Syst_MDM,1.0))/Total_OnBkgd
    #Syst_MDM = 1.*pow(max(0.,Syst_MDM-Hist2D_OnData_Sum.Integral()+Total_OnBkgd),0.5)/Total_OnBkgd

def Theta2HistScale(Hist,scale,scale_err):

    for b in range(1,Hist.GetNbinsX()+1):
        old_content = Hist.GetBinContent(b)
        old_error = Hist.GetBinError(b)
        new_content = old_content*scale
        new_error = 0
        if old_content>0 and scale>0:
            #new_error = new_content*(old_error/old_content)
            new_error = pow(pow(old_content*scale_err,2)+pow(old_error*scale,2),0.5)
        Hist.SetBinContent(b,new_content)

def NormalizeEnergyHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Energy.Reset()
    Hist_OnData_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Energy.Reset()
    Hist_OnBkgd_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnRFoV_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnRFoV_Energy.Reset()
    Hist_OnRFoV_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_CR_Energy_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_Energy.Reset()
    Hist_OnDark_Energy.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Energy_CamCenter_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Energy_CamCenter.Reset()
    Hist_OnData_Energy_CamCenter.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Energy_CamCenter_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Energy_CamCenter.Reset()
    Hist_OnBkgd_Energy_CamCenter.Add(InputFile.Get(HistName))

    HistName = "Hist_OnData_SR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Zenith.Reset()
    Hist_OnData_Zenith.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Zenith_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Zenith.Reset()
    Hist_OnBkgd_Zenith.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Energy[nth_roi].Reset()
        Hist_OnData_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_Energy_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
        Hist_OnBkgd_RoI_Energy[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_SR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_MJD[nth_roi].Reset()
        Hist_OnData_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_RoI_MJD_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_MJD[nth_roi].Reset()
        Hist_OnBkgd_RoI_MJD[nth_roi].Add(InputFile.Get(HistName))


    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Energy.Reset()
        Hist_OnBkgd_Energy.Reset()
        Hist_OnRFoV_Energy.Reset()
        Hist_OnData_Energy_CamCenter.Reset()
        Hist_OnBkgd_Energy_CamCenter.Reset()
        Hist_OnDark_Energy.Reset()
        Hist_OnData_Zenith.Reset()
        Hist_OnBkgd_Zenith.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Energy[nth_roi].Reset()
            Hist_OnBkgd_RoI_Energy[nth_roi].Reset()
            Hist_OnData_RoI_MJD[nth_roi].Reset()
            Hist_OnBkgd_RoI_MJD[nth_roi].Reset()

def NormalizeTheta2Histograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_Yoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Yoff.Reset()
    Hist_OnData_Yoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Yoff_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Yoff.Reset()
    Hist_OnBkgd_Yoff.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Yoff_Raw_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Yoff_Raw.Reset()
    Hist_OnBkgd_Yoff_Raw.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Theta2.Reset()
    Hist_OnData_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Theta2.Reset()
    Hist_OnBkgd_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnRFoV_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnRFoV_Theta2.Reset()
    Hist_OnRFoV_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_OnDark_CR_Skymap_Theta2_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnDark_Theta2.Reset()
    Hist_OnDark_Theta2.Add(InputFile.Get(HistName))

    for nth_roi in range(0,len(roi_ra)):
        HistName = "Hist_OnData_SR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnData_RoI_Theta2[nth_roi].Reset()
        Hist_OnData_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))
        HistName = "Hist_OnData_CR_Skymap_RoI_Theta2_V%s_ErecS%sto%s"%(nth_roi,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()
        Hist_OnBkgd_RoI_Theta2[nth_roi].Add(InputFile.Get(HistName))


    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Yoff.Reset()
        Hist_OnBkgd_Yoff.Reset()
        Hist_OnBkgd_Yoff_Raw.Reset()
        Hist_OnData_Theta2.Reset()
        Hist_OnBkgd_Theta2.Reset()
        Hist_OnRFoV_Theta2.Reset()
        Hist_OnDark_Theta2.Reset()
        for nth_roi in range(0,len(roi_ra)):
            Hist_OnData_RoI_Theta2[nth_roi].Reset()
            Hist_OnBkgd_RoI_Theta2[nth_roi].Reset()

def StackEnergyHistograms():

    Hist_OnData_Energy_Sum.Add(Hist_OnData_Energy)
    Hist_OnBkgd_Energy_Sum.Add(Hist_OnBkgd_Energy)
    Hist_OnRFoV_Energy_Sum.Add(Hist_OnRFoV_Energy)
    Hist_OnData_Energy_CamCenter_Sum.Add(Hist_OnData_Energy_CamCenter)
    Hist_OnBkgd_Energy_CamCenter_Sum.Add(Hist_OnBkgd_Energy_CamCenter)
    Hist_OnDark_Energy_Sum.Add(Hist_OnDark_Energy)
    Hist_OnData_Zenith_Sum.Add(Hist_OnData_Zenith)
    Hist_OnBkgd_Zenith_Sum.Add(Hist_OnBkgd_Zenith)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Energy_Sum[nth_roi].Add(Hist_OnData_RoI_Energy[nth_roi])
        Hist_OnBkgd_RoI_Energy_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Energy[nth_roi])
        Hist_OnData_RoI_MJD_Sum[nth_roi].Add(Hist_OnData_RoI_MJD[nth_roi])
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Add(Hist_OnBkgd_RoI_MJD[nth_roi])

def CorrectTheta2Histograms():

    Hist_OnBkgd_Theta2_Sum.Scale(1.+Syst_Corr)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Scale(1.+Syst_Corr)
        Hist_OnBkgd_RoI_MJD_Sum[nth_roi].Scale(1.+Syst_Corr)

def StackTheta2Histograms():

    Hist_OnData_Yoff_Sum.Add(Hist_OnData_Yoff)
    Hist_OnBkgd_Yoff_Sum.Add(Hist_OnBkgd_Yoff)
    Hist_OnBkgd_Yoff_Raw_Sum.Add(Hist_OnBkgd_Yoff_Raw)
    Hist_OnData_Theta2_Sum.Add(Hist_OnData_Theta2)
    Hist_OnBkgd_Theta2_Sum.Add(Hist_OnBkgd_Theta2)
    Hist_OnRFoV_Theta2_Sum.Add(Hist_OnRFoV_Theta2)
    Hist_OnDark_Theta2_Sum.Add(Hist_OnDark_Theta2)

    for nth_roi in range(0,len(roi_ra)):
        Hist_OnData_RoI_Theta2_Sum[nth_roi].Add(Hist_OnData_RoI_Theta2[nth_roi])
        Hist_OnBkgd_RoI_Theta2_Sum[nth_roi].Add(Hist_OnBkgd_RoI_Theta2[nth_roi])

def NormalizeCameraFoVHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    for nth_sample in range(0,n_control_samples):

        HistName = "Hist_OffData_SR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffData_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))
        HistName = "Hist_OffData_CR_CameraFoV_Theta2_V%s_ErecS%sto%s"%(nth_sample,ErecS_lower_cut_int,ErecS_upper_cut_int)
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
        Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Add(InputFile.Get(HistName))

        if Hist2D_OffData[nth_sample].Integral()<1600.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()

        data_total, data_err = IntegralAndError(Hist_OffData_Energy[nth_sample],bin_lower,bin_upper)
        if data_total==0.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
            return

        old_integral = Hist_OffData_CameraFoV_Theta2[nth_sample].Integral()
        data_scale = 0
        data_scale_err = 0
        if not data_total==0 and not old_integral==0:
            data_scale = data_total/old_integral
            data_scale_err = data_scale*(data_err/data_total)
        else:
            data_scale = 0
            data_scale_err = 0
        if not data_total==0:
            Theta2HistScale(Hist_OffData_CameraFoV_Theta2[nth_sample],data_scale,data_scale_err)
        else:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Scale(0)

        bkg_total, bkg_err = IntegralAndError(Hist_OffBkgd_Energy[nth_sample],bin_lower,bin_upper)
        if bkg_total==0.:
            Hist_OffData_CameraFoV_Theta2[nth_sample].Reset()
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Reset()
            return

        old_integral = Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Integral()
        bkgd_scale = 0
        bkgd_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            bkgd_scale = bkg_total/old_integral
            bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        else:
            bkgd_scale = 0
            bkgd_scale_err = 0
        if not bkg_total==0:
            Theta2HistScale(Hist_OffBkgd_CameraFoV_Theta2[nth_sample],bkgd_scale,bkgd_scale_err)
        else:
            Hist_OffBkgd_CameraFoV_Theta2[nth_sample].Scale(0)

def StackCameraFoVHistograms():

    for nth_sample in range(0,n_control_samples):
        Hist_OffData_CameraFoV_Theta2_Sum[nth_sample].Add(Hist_OffData_CameraFoV_Theta2[nth_sample])
        Hist_OffBkgd_CameraFoV_Theta2_Sum[nth_sample].Add(Hist_OffBkgd_CameraFoV_Theta2[nth_sample])

def RaDecHistScale(Hist,scale,scale_err):

    for bx in range(1,Hist.GetNbinsX()+1):
        for by in range(1,Hist.GetNbinsY()+1):
            old_content = Hist.GetBinContent(bx,by)
            old_error = Hist.GetBinError(bx,by)
            new_content = old_content*scale
            new_error = 0
            if scale>0:
                new_error = old_error*scale
            Hist.SetBinContent(bx,by,new_content)
            Hist.SetBinError(bx,by,new_error)

def NormalizeSkyMapHistograms(FilePath,ebin):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist_OnData_Energy.FindBin(ErecS_lower_cut)
    bin_upper = Hist_OnData_Energy.FindBin(ErecS_upper_cut)-1

    ErecS_lower_cut_int = int(ErecS_lower_cut)
    ErecS_upper_cut_int = int(ErecS_upper_cut)

    HistName = "Hist_OnData_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap.Reset()
    Hist_OnData_Skymap.Add(InputFile.Get(HistName))
    Hist_Data_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_SR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnData_Skymap_Galactic.Reset()
    Hist_OnData_Skymap_Galactic.Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap.Reset()
    Hist_OnBkgd_Skymap.Add(InputFile.Get(HistName))
    Hist_Bkgd_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    Hist_Syst_Energy_Skymap[ebin].Add(InputFile.Get(HistName))
    HistName = "Hist_OnData_CR_Skymap_Galactic_ErecS%sto%s"%(ErecS_lower_cut_int,ErecS_upper_cut_int)
    Hist_OnBkgd_Skymap_Galactic.Reset()
    Hist_OnBkgd_Skymap_Galactic.Add(InputFile.Get(HistName))

    if Hist2D_OnData.Integral()<1600.:
        Hist_OnData_Skymap.Reset()
        Hist_OnData_Skymap_Galactic.Reset()
        Hist_OnBkgd_Skymap.Reset()
        Hist_OnBkgd_Skymap_Galactic.Reset()

def StackSkymapHistograms(ebin):

    Hist_OnData_Skymap_Sum.Add(Hist_OnData_Skymap)
    Hist_OnData_Skymap_Galactic_Sum.Add(Hist_OnData_Skymap_Galactic)
    Hist_OnBkgd_Skymap_Sum.Add(Hist_OnBkgd_Skymap)
    Hist_OnBkgd_Skymap_Galactic_Sum.Add(Hist_OnBkgd_Skymap_Galactic)

    Hist_OnData_Skymap_ProjX_Sum.Add(Hist_OnData_Skymap.ProjectionX())
    Hist_OnBkgd_Skymap_ProjX_Sum.Add(Hist_OnBkgd_Skymap.ProjectionX())
    Hist_OnData_Skymap_ProjY_Sum.Add(Hist_OnData_Skymap.ProjectionY())
    Hist_OnBkgd_Skymap_ProjY_Sum.Add(Hist_OnBkgd_Skymap.ProjectionY())

    #Syst_MDM_single = CalculateSystErrorIndividual_v3()
    Syst_MDM_single = Syst_MDM
    print 'Syst_MDM_single = %s'%(Syst_MDM_single)

    Hist_Syst_Energy_Skymap[ebin].Scale(Syst_MDM_single)

    #Hist_OnBkgd_Skymap_Syst_MDM.Add(Hist_OnBkgd_Skymap,Syst_MDM_single)
    #Hist_OnBkgd_Skymap_Galactic_Syst_MDM.Add(Hist_OnBkgd_Skymap_Galactic,Syst_MDM_single)
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_MDM.GetNbinsX()):
        for biny in range(0,Hist_OnBkgd_Skymap_Syst_MDM.GetNbinsY()):
            old_content = Hist_OnBkgd_Skymap_Syst_MDM.GetBinContent(binx+1,biny+1)
            new_content = Hist_OnBkgd_Skymap.GetBinContent(binx+1,biny+1)*Syst_MDM_single
            Hist_OnBkgd_Skymap_Syst_MDM.SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content*new_content,0.5))
    for binx in range(0,Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetNbinsX()):
        for biny in range(0,Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetNbinsY()):
            old_content = Hist_OnBkgd_Skymap_Galactic_Syst_MDM.GetBinContent(binx+1,biny+1)
            new_content = Hist_OnBkgd_Skymap_Galactic.GetBinContent(binx+1,biny+1)*Syst_MDM_single
            Hist_OnBkgd_Skymap_Galactic_Syst_MDM.SetBinContent(binx+1,biny+1,pow(old_content*old_content+new_content*new_content,0.5))

    #for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
    #    old_err = Hist_SystErr_MSCL.GetBinError(binx+1)
    #    new_err = old_err+Hist_OnBkgd_MSCL.GetBinContent(binx+1)*Syst_MDM_single
    #    Hist_SystErr_MSCL.SetBinError(binx+1,new_err)
    #for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
    #    old_err = Hist_SystErr_MSCW.GetBinError(binx+1)
    #    new_err = old_err+Hist_OnBkgd_MSCW.GetBinContent(binx+1)*Syst_MDM_single
    #    Hist_SystErr_MSCW.SetBinError(binx+1,new_err)
    #for binx in range(0,Hist_SystErr_Energy.GetNbinsX()):
    #    old_err = Hist_SystErr_Energy.GetBinError(binx+1)
    #    new_err = old_err+Hist_OnBkgd_Energy.GetBinContent(binx+1)*Syst_MDM_single
    #    Hist_SystErr_Energy.SetBinError(binx+1,new_err)
    #for binx in range(0,Hist_SystErr_Theta2.GetNbinsX()):
    #    old_err = Hist_SystErr_Theta2.GetBinError(binx+1)
    #    new_err = old_err+Hist_OnBkgd_Theta2.GetBinContent(binx+1)*Syst_MDM_single
    #    Hist_SystErr_Theta2.SetBinError(binx+1,new_err)
    for binx in range(0,Hist_SystErr_MSCL.GetNbinsX()):
        old_err = Hist_SystErr_MSCL.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_MSCL.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_MSCL.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_MSCW.GetNbinsX()):
        old_err = Hist_SystErr_MSCW.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_MSCW.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_MSCW.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_Energy.GetNbinsX()):
        old_err = Hist_SystErr_Energy.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_Energy.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_Energy.SetBinError(binx+1,pow(new_err,0.5))
    for binx in range(0,Hist_SystErr_Theta2.GetNbinsX()):
        old_err = Hist_SystErr_Theta2.GetBinError(binx+1)
        new_err = pow(old_err,2)+pow(Hist_OnBkgd_Theta2.GetBinContent(binx+1)*Syst_MDM_single,2)
        Hist_SystErr_Theta2.SetBinError(binx+1,pow(new_err,0.5))
    Hist_OnBkgd_Skymap_Syst_ProjX.Reset()
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_ProjX_Sum.GetNbinsX()):
        old_err = Hist_OnBkgd_Skymap_Syst_ProjX_Sum.GetBinError(binx+1)
        new_err_mdm = Hist_OnBkgd_Skymap.ProjectionX().GetBinContent(binx+1)*Syst_MDM_single
        new_err_shape = Hist_OnBkgd_Skymap_Syst_ProjX.GetBinContent(binx+1)
        new_err = pow(old_err,2)+pow(new_err_mdm,2)
        Hist_OnBkgd_Skymap_Syst_ProjX_Sum.SetBinError(binx+1,pow(new_err,0.5))
    Hist_OnBkgd_Skymap_Syst_ProjY.Reset()
    for binx in range(0,Hist_OnBkgd_Skymap_Syst_ProjY_Sum.GetNbinsX()):
        old_err = Hist_OnBkgd_Skymap_Syst_ProjY_Sum.GetBinError(binx+1)
        new_err_mdm = Hist_OnBkgd_Skymap.ProjectionY().GetBinContent(binx+1)*Syst_MDM_single
        new_err_shape = Hist_OnBkgd_Skymap_Syst_ProjY.GetBinContent(binx+1)
        new_err = pow(old_err,2)+pow(new_err_mdm,2)
        Hist_OnBkgd_Skymap_Syst_ProjY_Sum.SetBinError(binx+1,pow(new_err,0.5))

def Smooth2DMap(Hist_Old,smooth_size,addLinearly):

    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(2*smooth_size/bin_size) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            bin_content = 0
            bin_error = 0
            bin_norm = 0
            locationx1 = Hist_Old.GetXaxis().GetBinCenter(bx1)
            locationy1 = Hist_Old.GetYaxis().GetBinCenter(by1)
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2>=1 and bx2<=Hist_Old.GetNbinsX():
                        if by2>=1 and by2<=Hist_Old.GetNbinsY():
                            locationx2 = Hist_Old.GetXaxis().GetBinCenter(bx2)
                            locationy2 = Hist_Old.GetYaxis().GetBinCenter(by2)
                            distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5)
                            bin_content += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinContent(bx2,by2)
                            bin_norm += ROOT.TMath.Gaus(distance,0,smooth_size)
                            if not addLinearly:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*pow(Hist_Old.GetBinError(bx2,by2),2)
                            else:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    return Hist_Smooth

def CorrectLEE(z_score,n_tests,threshold_sigma):

    if abs(z_score)<threshold_sigma:
        return z_score
    new_z_score = 0.
    if z_score>0.:
        p_value = 1.-st.norm.cdf(z_score)
        p_value = min(p_value*n_tests,1.0)
        new_z_score = st.norm.ppf(1.-p_value)
    else:
        p_value = st.norm.cdf(z_score)
        p_value = min(p_value*n_tests,1.0)
        new_z_score = st.norm.ppf(p_value)
    return min(new_z_score,1000.)

def GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,syst_method):

    Hist_Skymap = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            Shape_Err = Hist_Syst.GetBinContent(bx+1,by+1)
            #Shape_Err = 0.
            #Shape_Err = NBkg*syst_method
            Stat_Err = Hist_Bkg.GetBinError(bx+1,by+1)
            #NBkg_Err = Shape_Err
            NBkg_Err = pow(pow(Stat_Err,2)+pow(Shape_Err,2),0.5)
            if syst_method==0.: NBkg_Err = Stat_Err
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
    return Hist_Skymap

def reflectXaxis(hist):

    # taken from VPlotAnasumHistograms.cpp
	
    # temporary histogram
    hT = ROOT.TH2D( "%s_REFLECTED"%(hist.GetName()), "", hist.GetNbinsX(), -1.*hist.GetXaxis().GetXmax(), -1.*hist.GetXaxis().GetXmin(), hist.GetNbinsY(), hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax() )
    hT.SetStats( 0 )
    hT.SetXTitle( hist.GetXaxis().GetTitle() )
    hT.SetYTitle( hist.GetYaxis().GetTitle() )
	
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsX()+1):
            hT.SetBinContent( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinContent( binx, biny ) )
            hT.SetBinError( hist.GetNbinsX() + 1 - binx, biny, hist.GetBinError( binx, biny ) )
    return hT

def FillSkymapHoles(hist_map, map_resolution):

    hist_map_new = hist_map.Clone()
    bin_size = hist_map.GetXaxis().GetBinCenter(2)-hist_map.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(map_resolution/bin_size) + 1
    for bx1 in range(1,hist_map.GetNbinsX()+1):
        for by1 in range(1,hist_map.GetNbinsY()+1):
            bin_content_1 = hist_map.GetBinContent(bx1,by1)
            if bin_content_1!=0.: continue
            min_distance = 1e10
            min_distance_content = 0.
            locationx1 = hist_map.GetXaxis().GetBinCenter(bx1)
            locationy1 = hist_map.GetYaxis().GetBinCenter(by1)
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2>=1 and bx2<=hist_map.GetNbinsX():
                        if by2>=1 and by2<=hist_map.GetNbinsY():
                            bin_content_2 = hist_map.GetBinContent(bx2,by2)
                            if bin_content_2==0.: continue
                            locationx2 = hist_map.GetXaxis().GetBinCenter(bx2)
                            locationy2 = hist_map.GetYaxis().GetBinCenter(by2)
                            distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5)
                            if min_distance>distance:
                                min_distance = distance
                                min_distance_content = bin_content_2
            hist_map_new.SetBinContent(bx1,by1,min_distance_content)
    return hist_map_new

def GetSkyViewMap(map_file, hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open(map_file)
    for line in inputFile:
        sig = float(line.split(' ')[2])
        l = float(line.split(' ')[0])
        b = float(line.split(' ')[1])
        if not isRaDec: 
            l, b = ConvertRaDecToGalactic(l,b)
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,0.))

    map_resolution = 0.1
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def GetHawcSkymap(hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open('MWL_maps/hawc_map.txt')
    for line in inputFile:
        sig = float(line.split(' ')[0])
        l = float(line.split(' ')[1])
        b = float(line.split(' ')[2])
        if not isRaDec: 
            l, b = ConvertRaDecToGalactic(l,b)
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,0.))

    map_resolution = 0.1
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def GetCOSkymap(hist_map, isRaDec):

    hist_map.Reset()
    inputFile = open('MWL_maps/CO_skymap_v2.txt')
    if isRaDec: 
        inputFile = open('MWL_maps/CO_skymap_RaDec_v2.txt')
    for line in inputFile:
        sig = float(line.split(' ')[0])
        l = float(line.split(' ')[1])
        b = float(line.split(' ')[2])
        binx = hist_map.GetXaxis().FindBin(l)
        biny = hist_map.GetYaxis().FindBin(b)
        old_sig = hist_map.GetBinContent(binx+1,biny+1)
        hist_map.SetBinContent(binx+1,biny+1,max(sig,0.))

    map_resolution = 0.2
    #map_resolution = 0.4
    hist_map_new = FillSkymapHoles(hist_map, map_resolution)

    return hist_map_new

def OutputToTxtFile(hist_input,name):

    file1 = open(name,"w")
    for binx in range(0,hist_input.GetNbinsX()):
        for biny in range(0,hist_input.GetNbinsY()):
            file1.write('x: %s y: %s sigma: %s \n'%(hist_input.GetXaxis().GetBinCenter(binx+1),hist_input.GetYaxis().GetBinCenter(biny+1),hist_input.GetBinContent(binx+1,biny+1)))
    file1.close()

def Count5SigmaBins(Hist_sigma,Hist_excess,threshold_sigma):

    skymap_bin_size_x = Hist_sigma.GetXaxis().GetBinCenter(2)-Hist_sigma.GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = Hist_sigma.GetYaxis().GetBinCenter(2)-Hist_sigma.GetYaxis().GetBinCenter(1)
    unit_area = skymap_bin_size_x*skymap_bin_size_y
    n_tests = Hist_sigma.GetNbinsX()*Hist_sigma.GetNbinsY()
    global_p_value = 1.
    nbins = 0
    for binx in range(Hist_sigma.GetNbinsX()):
        for biny in range(Hist_sigma.GetNbinsY()):
            Sig = Hist_sigma.GetBinContent(binx+1,biny+1)
            if Sig>threshold_sigma:
                p_value = 1.-st.norm.cdf(Sig)
                p_value_correct = min(p_value*n_tests,1.0)
                global_p_value = global_p_value*p_value_correct
                nbins += 1
    if nbins==0:
        return Hist_sigma.GetMaximum()
    return st.norm.ppf(1.-global_p_value)

def MakeSpectrumIndexSkymap(hist_exposure,hist_data,hist_bkgd,hist_syst,hist_contour,title_x,title_y,name,nbins,zoomin_scale):

    max_nbins = 60
    best_nbins = nbins

    MapEdge_left = hist_contour[0].GetXaxis().GetBinLowEdge(1)
    MapEdge_right = hist_contour[0].GetXaxis().GetBinLowEdge(hist_contour[0].GetNbinsX()+1)
    MapEdge_lower = hist_contour[0].GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = hist_contour[0].GetYaxis().GetBinLowEdge(hist_contour[0].GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.

    hist_data_skymap = []
    hist_diff_skymap = []
    hist_syst_skymap = []
    hist_expo_skymap = []
    for ebin in range(0,len(energy_bin)-1):
        hist_diff_skymap += [ROOT.TH2D("hist_diff_skymap_%s"%(ebin),"",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_data_skymap += [ROOT.TH2D("hist_data_skymap_%s"%(ebin),"",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_syst_skymap += [ROOT.TH2D("hist_syst_skymap_%s"%(ebin),"",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
        hist_expo_skymap += [ROOT.TH2D("hist_expo_skymap_%s"%(ebin),"",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)]
    for ebin in range(0,len(energy_bin)-1):
        for bx in range(0,hist_contour[0].GetNbinsX()):
            for by in range(0,hist_contour[0].GetNbinsY()):
                bx_center = hist_contour[0].GetXaxis().GetBinCenter(bx+1)
                by_center = hist_contour[0].GetYaxis().GetBinCenter(by+1)
                bx2 = hist_syst_skymap[ebin].GetXaxis().FindBin(bx_center)
                by2 = hist_syst_skymap[ebin].GetYaxis().FindBin(by_center)
                new_content = hist_syst[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_syst_skymap[ebin].GetBinContent(bx2,by2)
                hist_syst_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_data[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_data_skymap[ebin].GetBinContent(bx2,by2)
                hist_data_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_data[ebin].GetBinContent(bx+1,by+1)-hist_bkgd[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_diff_skymap[ebin].GetBinContent(bx2,by2)
                hist_diff_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
                new_content = hist_exposure[ebin].GetBinContent(bx+1,by+1)
                old_content = hist_expo_skymap[ebin].GetBinContent(bx2,by2)
                hist_expo_skymap[ebin].SetBinContent(bx2,by2,old_content+new_content)
    for ebin in range(0,len(energy_bin)-1):
        hist_syst_skymap[ebin].Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
        hist_data_skymap[ebin].Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
        hist_diff_skymap[ebin].Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
        hist_expo_skymap[ebin].Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    for ebin in range(0,len(energy_bin)-1):
        max_exposure = hist_expo_skymap[ebin].GetMaximum()
        if max_exposure>0.:
            hist_expo_skymap[ebin].Scale(1./max_exposure)
        else:
            hist_expo_skymap[ebin].Scale(0.)


    hist_index_skymap = ROOT.TH2D("hist_index_skymap","",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_index_skymap.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    hist_index_skymap_full = ROOT.TH2D("hist_index_skymap_full","",Skymap_nbins/zoomin_scale,MapCenter_x-MapSize_x/zoomin_scale,MapCenter_x+MapSize_x/zoomin_scale,Skymap_nbins/zoomin_scale,MapCenter_y-MapSize_y/zoomin_scale,MapCenter_y+MapSize_y/zoomin_scale)
    hist_index_skymap_full.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    pad1.SetLogx()
    pad1.SetLogy()

    Hist_Flux = ROOT.TH1D("Hist_Flux","",len(energy_bin)-1,array('d',energy_bin))
    func_source = []
    nth_roi = 0
    min_index = 0.
    for binx in range(0,hist_index_skymap.GetNbinsX()):
        for biny in range(0,hist_index_skymap.GetNbinsY()):
            bx_center = hist_index_skymap.GetXaxis().GetBinCenter(binx+1)
            by_center = hist_index_skymap.GetYaxis().GetBinCenter(biny+1)
            binx2 = hist_contour[0].GetXaxis().FindBin(bx_center)
            biny2 = hist_contour[0].GetYaxis().FindBin(by_center)
            nth_roi += 1
            func_source += [ROOT.TF1("func_source_%s"%(nth_roi-1),"[0]*pow(10,-12)*pow(x/1000.,[1])", 150, 10000)]
            Hist_Flux.Reset()
            max_z_score = 0.
            min_z_score = 99.
            z_score = 0.
            n_2sigma = 0
            n_3sigma = 0
            n_4sigma = 0
            n_2sigma_connect = 0
            n_3sigma_connect = 0
            n_4sigma_connect = 0
            n_2sigma_connect_max = 0
            n_3sigma_connect_max = 0
            n_4sigma_connect_max = 0
            for ebin in range(0,len(energy_bin)-1):
                acceptance = hist_expo_skymap[ebin].GetBinContent(binx+1,biny+1)
                data_count = hist_data_skymap[ebin].GetBinContent(binx+1,biny+1)
                diff_count = hist_diff_skymap[ebin].GetBinContent(binx+1,biny+1)
                syst_error = hist_syst_skymap[ebin].GetBinContent(binx+1,biny+1)
                z_score = CalculateSignificance(diff_count,data_count-diff_count,syst_error)
                if max_z_score<z_score: max_z_score = z_score
                if min_z_score>z_score: min_z_score = z_score
                if z_score>2.: 
                    n_2sigma += 1
                    n_2sigma_connect += 1
                else:
                    n_2sigma_connect = 0
                if z_score>2.5: 
                    n_3sigma += 1
                    n_3sigma_connect += 1
                else:
                    n_3sigma_connect = 0
                if z_score>3.0: 
                    n_4sigma += 1
                    n_4sigma_connect += 1
                else:
                    n_4sigma_connect = 0
                n_2sigma_connect_max = max(n_2sigma_connect_max,n_2sigma_connect)
                n_3sigma_connect_max = max(n_3sigma_connect_max,n_3sigma_connect)
                n_4sigma_connect_max = max(n_4sigma_connect_max,n_4sigma_connect)
                energy = Hist_Flux.GetBinCenter(ebin+1)
                ebin_fine = Hist_EffArea_Sum.FindBin(energy)
                if Hist_EffArea_Sum.GetBinContent(ebin_fine)==0.: continue
                deltaE = (energy_bin[ebin+1]-energy_bin[ebin])/1000.
                if acceptance>0.05:
                    scale = 1./(Hist_EffArea_Sum.GetBinContent(ebin_fine)*10000.*deltaE*acceptance)
                    Hist_Flux.SetBinContent(ebin+1,diff_count*scale)
                    Hist_Flux.SetBinError(ebin+1,pow(data_count+syst_error*syst_error,0.5)*scale)
                else:
                    Hist_Flux.SetBinContent(ebin+1,0.)
                    Hist_Flux.SetBinError(ebin+1,0.)
            func_source[nth_roi-1].SetParameters(37.5/(nbins*nbins),-2.)
            Hist_Flux.Fit("func_source_%s"%(nth_roi-1),"N","",ErecS_lower_cut,ErecS_upper_cut)
            good_fit = True
            good_data = True
            if (max_nbins/best_nbins)>=4:
                if n_3sigma_connect_max<4 and n_4sigma_connect_max<3: good_data = False
            else:
                if n_2sigma_connect_max<3 and n_3sigma_connect_max<2: good_data = False
            if abs(func_source[nth_roi-1].GetParameter(1))<5.*func_source[nth_roi-1].GetParError(1): good_fit = False
            if good_fit and good_data:
                hist_index_skymap.SetBinContent(binx+1,biny+1,func_source[nth_roi-1].GetParameter(1))
                hist_index_skymap.SetBinError(binx+1,biny+1,func_source[nth_roi-1].GetParError(1))
                if func_source[nth_roi-1].GetParameter(1)<min_index:
                    min_index = func_source[nth_roi-1].GetParameter(1)
                #Hist_Flux.Draw("E")
                #func_source[nth_roi-1].Draw("same")
                #c_both.SaveAs('output_plots/Fit_x%s_y%s.png'%(binx,biny))
            else:
                hist_index_skymap.SetBinContent(binx+1,biny+1,-99.)
                hist_index_skymap.SetBinError(binx+1,biny+1,-99.)
                #Hist_Flux.Draw("E")
                #func_source[nth_roi-1].Draw("same")
                #c_both.SaveAs('output_plots/Fit_x%s_y%s.png'%(binx,biny))
            if good_fit:
                hist_index_skymap_full.SetBinContent(binx+1,biny+1,func_source[nth_roi-1].GetParameter(1))
                hist_index_skymap_full.SetBinError(binx+1,biny+1,func_source[nth_roi-1].GetParError(1))
            else:
                hist_index_skymap_full.SetBinContent(binx+1,biny+1,-99.)
                hist_index_skymap_full.SetBinError(binx+1,biny+1,-99.)

    hist_index_skymap_reflect = reflectXaxis(hist_index_skymap)
    hist_index_skymap_reflect.SetMinimum(min_index)
    hist_index_skymap_full_reflect = reflectXaxis(hist_index_skymap_full)
    hist_index_skymap_full_reflect.SetMaximum(hist_index_skymap_reflect.GetMaximum())
    hist_index_skymap_full_reflect.SetMinimum(min_index)

    hist_contour_reflect = []
    for ebin in range(0,len(hist_contour)):
        hist_contour_reflect += [reflectXaxis(hist_contour[ebin])]
        hist_contour_reflect[ebin].SetLineColor(1)
        hist_contour_reflect[ebin].SetContour(3)
        hist_contour_reflect[ebin].SetContourLevel(0,3)
        hist_contour_reflect[ebin].SetContourLevel(1,4)
        hist_contour_reflect[ebin].SetContourLevel(2,5)

    Hist_CO = hist_data[0].Clone()
    ContourLevel_1 = 0.
    ContourLevel_2 = 0.
    ContourLevel_3 = 0.
    if sys.argv[1]=='GammaCygni_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv19487154506678_1420MHz_GammaCygni.txt", Hist_CO, True)
        ContourLevel_1 = 11000.
        ContourLevel_2 = 13000.
        ContourLevel_3 = 15000.
    if sys.argv[1]=='MGRO_J1908_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv22338613766121_SwiftBat14keV_J1908.txt", Hist_CO, True)
        ContourLevel_1 = 1.0
        ContourLevel_2 = 2.0
        ContourLevel_3 = 4.0
    if sys.argv[1]=='MAGIC_J1857_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv19487414377965_1420MHz_J1857.txt", Hist_CO, True)
        ContourLevel_1 = 7000.
        ContourLevel_2 = 9000.
        ContourLevel_3 = 11000.
    Hist_CO = reflectXaxis(Hist_CO)
    Hist_CO.SetContour(3)
    Hist_CO.SetContourLevel(0,ContourLevel_1)
    Hist_CO.SetContourLevel(1,ContourLevel_2)
    Hist_CO.SetContourLevel(2,ContourLevel_3)
    Hist_CO.SetLineColor(2)

    other_star_labels = []
    other_star_markers = []
    other_star_names = []
    bright_star_labels = []
    bright_star_markers = []
    faint_star_labels = []
    faint_star_markers = []
    star_range = 2.5
    if title_x=="RA":
        for star in range(0,len(other_stars)):
            if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
            other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
            other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
            other_star_names += [other_stars[star]]
        for star in range(0,len(bright_star_ra)):
            bright_star_markers += [ROOT.TMarker(-bright_star_ra[star],bright_star_dec[star],30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-bright_star_ra[star]-0.15,bright_star_dec[star]+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
        for star in range(0,len(faint_star_ra)):
            faint_star_markers += [ROOT.TMarker(-faint_star_ra[star],faint_star_dec[star],30)]
            faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
            faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
            faint_star_labels += [ROOT.TLatex(-faint_star_ra[star]-0.15,faint_star_dec[star]+0.15,'b-mag %s'%(faint_star_brightness[star]))]
            faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
            faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)
    else:
        for star in range(0,len(other_stars)):
            gal_l, gal_b = ConvertRaDecToGalactic(other_star_coord[star][0],other_star_coord[star][1])
            if pow(source_l-gal_l,2)+pow(source_b-gal_b,2)>star_range*star_range: continue
            other_star_markers += [ROOT.TMarker(-gal_l,gal_b,2)]
            other_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
            other_star_names += [other_stars[star]]
        for star in range(0,len(bright_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(bright_star_ra[star],bright_star_dec[star])
            bright_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
        for star in range(0,len(faint_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(faint_star_ra[star],faint_star_dec[star])
            faint_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
            faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
            faint_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(faint_star_brightness[star]))]
            faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
            faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)


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

    pad1.cd()
    hist_index_skymap_reflect.Draw("COL4Z")
    #hist_index_skymap_full_reflect.Draw("COL4Z")
    #Hist_CO.Draw("CONT3 same")
    #hist_index_skymap_reflect.Draw("TEXT45 same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    hist_index_skymap_reflect.GetXaxis().SetLabelOffset(999)
    hist_index_skymap_reflect.GetXaxis().SetTickLength(0)
    x1 = hist_index_skymap_reflect.GetXaxis().GetXmin()
    x2 = hist_index_skymap_reflect.GetXaxis().GetXmax()
    y1 = hist_index_skymap_reflect.GetYaxis().GetXmin()
    y2 = hist_index_skymap_reflect.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(hist_index_skymap_reflect.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    pad3.cd()
    lumilab3 = ROOT.TLatex(0.15,0.50,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.30,'MJD %s-%s'%(MJD_Start,MJD_End) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    canvas.SaveAs('output_plots/SkymapIndex_%s_%s.png'%(name,selection_tag))

def VariableSkymapBins(syst_method,Hist_SR_input,Hist_Bkg_input,Hist_Syst_input,Hist_RBM_input,Hist_Exposure_input,xtitle,ytitle,name, input_nbins):

    threshold_sigma = 2.5
    max_nbins = 60
    #list_nbins = [60,40,30,20,15,12,10,8,6,5,4,3,2,1]
    list_nbins = [30,20,15,12,10,6,5,4,3,2,1]

    hist_maxsig = ROOT.TH1D("hist_maxsig","",len(list_nbins),0,len(list_nbins))
    list_maxsig = []
    list_binsize = []

    max_nbins_5sigma = 0
    best_nbins = 10
    for nbins in list_nbins:
        Hist_SR = reflectXaxis(Hist_SR_input)
        Hist_Bkg = reflectXaxis(Hist_Bkg_input)
        Hist_Syst = reflectXaxis(Hist_Syst_input)
        Hist_RBM = reflectXaxis(Hist_RBM_input)
        Hist_Exposure = reflectXaxis(Hist_Exposure_input)
        Hist_SR.Rebin2D(max_nbins/nbins,max_nbins/nbins)
        Hist_Bkg.Rebin2D(max_nbins/nbins,max_nbins/nbins)
        Hist_Syst.Rebin2D(max_nbins/nbins,max_nbins/nbins)
        Hist_RBM.Rebin2D(max_nbins/nbins,max_nbins/nbins)
        Hist_Exposure.Rebin2D(max_nbins/nbins,max_nbins/nbins)
        Hist_Skymap = GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,syst_method)
        nbins_5sigma = Count5SigmaBins(Hist_Skymap,Hist_SR,threshold_sigma)
        max_sigma = Hist_Skymap.GetMaximum()
        max_sigma = CorrectLEE(max_sigma,Hist_Skymap.GetNbinsX()*Hist_Skymap.GetNbinsY(),threshold_sigma)
        list_maxsig += [max_sigma]
        skymap_bin_size_x = Hist_SR.GetXaxis().GetBinCenter(2)-Hist_SR.GetXaxis().GetBinCenter(1)
        skymap_bin_size_y = Hist_SR.GetYaxis().GetBinCenter(2)-Hist_SR.GetYaxis().GetBinCenter(1)
        list_binsize += [skymap_bin_size_x]
        print 'syst_method = %s, nbins = %s, nbins_5sigma = %s, max_sigma = %0.2f'%(syst_method,nbins,nbins_5sigma,max_sigma)
        if nbins==1:
            total_data = Hist_SR.GetBinContent(1,1)
            total_bkgd = Hist_Bkg.GetBinContent(1,1)
            total_syst = Hist_Syst.GetBinContent(1,1)
            print 'total_data = %0.2f, total_bkgd = %0.2f, total_syst = %0.2f'%(total_data,total_bkgd,total_syst)
        if max_nbins_5sigma<nbins_5sigma:
            max_nbins_5sigma = nbins_5sigma
            best_nbins = nbins

    for binx in range(0,len(list_nbins)):
        hist_maxsig.SetBinContent(binx+1,list_maxsig[binx])
        hist_maxsig.GetXaxis().SetBinLabel(binx+1,'%0.1f'%(list_binsize[binx]))
    hist_maxsig.SetMinimum(0.)
    hist_maxsig.SetMaximum(10.)
    MakeOneHistPlot(hist_maxsig,'bin size','max z score','MaxSigma_%s_%s'%(name,selection_tag),False)

    best_nbins = input_nbins

    Hist_SR = reflectXaxis(Hist_SR_input)
    Hist_Bkg = reflectXaxis(Hist_Bkg_input)
    Hist_Syst = reflectXaxis(Hist_Syst_input)
    Hist_RBM = reflectXaxis(Hist_RBM_input)
    Hist_Exposure = reflectXaxis(Hist_Exposure_input)
    Hist_SR.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Bkg.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Syst.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_RBM.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Exposure.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Skymap = GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,syst_method)
    Hist_Skymap.GetYaxis().SetTitle(ytitle)
    Hist_Skymap.GetXaxis().SetTitle(xtitle)
    Hist_Skymap.GetZaxis().SetTitle('Significance')
    Hist_Skymap.GetZaxis().SetTitleOffset(1.1)
    Hist_Skymap.SetMaximum(5)
    Hist_Skymap.SetMinimum(-5)
    #Hist_Skymap.SetMaximum(3)
    #Hist_Skymap.SetMinimum(-3)

    skymap_bin_size_x = Hist_Skymap.GetXaxis().GetBinCenter(2)-Hist_Skymap.GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = Hist_Skymap.GetYaxis().GetBinCenter(2)-Hist_Skymap.GetYaxis().GetBinCenter(1)
    unit_area = skymap_bin_size_x*skymap_bin_size_y

    Hist_Contour = Hist_Skymap.Clone()
    Hist_Contour.Reset()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            Hist_Contour.SetBinContent(bx+1,by+1,(Hist_Skymap.GetBinContent(bx+1,by+1)))
    Hist_Contour.SetContour(1)
    Hist_Contour.SetContourLevel(0,3.)

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

    pad1.cd()
    Hist_Skymap.Draw("COL4Z")
    Hist_Contour.Draw("CONT3 same")
    Hist_Skymap.GetXaxis().SetLabelOffset(999)
    Hist_Skymap.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap.GetXaxis().GetXmin()
    x2 = Hist_Skymap.GetXaxis().GetXmax()
    y1 = Hist_Skymap.GetYaxis().GetXmin()
    y2 = Hist_Skymap.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    pad3.cd()
    lumilab3 = ROOT.TLatex(0.15,0.50,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.30,'MJD %s-%s'%(MJD_Start,MJD_End) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    canvas.SaveAs('output_plots/SkymapSig_%s_%s.png'%(name,selection_tag))

def Make2DSignificancePlot(syst_method,Hist_SR_input,Hist_Bkg_input,Hist_Syst_input,Hist_RBM_input,Hist_Exposure_input,xtitle,ytitle,name):

    Hist_SR = reflectXaxis(Hist_SR_input)
    Hist_Bkg = reflectXaxis(Hist_Bkg_input)
    Hist_Syst = reflectXaxis(Hist_Syst_input)
    Hist_RBM = reflectXaxis(Hist_RBM_input)
    Hist_Exposure = reflectXaxis(Hist_Exposure_input)

    Hist_Skymap = GetSignificanceMap(Hist_SR,Hist_Bkg,Hist_Syst,syst_method)
    max_sig = Hist_Skymap.GetMaximum()

    isRaDec = False
    if 'RaDec' in name: isRaDec = True

    Hist_HAWC = Hist_SR_input.Clone()
    Hist_HAWC = GetHawcSkymap(Hist_HAWC, isRaDec)
    Hist_HAWC = reflectXaxis(Hist_HAWC)
    Hist_HAWC.SetLineColor(0)
    Hist_HAWC.SetContour(3)
    Hist_HAWC.SetContourLevel(0,5)
    Hist_HAWC.SetContourLevel(1,10)
    Hist_HAWC.SetContourLevel(2,15)

    Hist_CO = Hist_SR_input.Clone()
    ContourLevel_1 = 0.
    ContourLevel_2 = 0.
    ContourLevel_3 = 0.
    #Hist_CO = GetCOSkymap(Hist_CO, isRaDec)
    #Hist_CO = reflectXaxis(Hist_CO)
    #Hist_CO.SetLineColor(0)
    #Hist_CO.SetContour(3)
    #Hist_CO.SetContourLevel(0,5)
    #Hist_CO.SetContourLevel(1,20)
    #Hist_CO.SetContourLevel(2,80)
    if sys.argv[1]=='GammaCygni_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv19487154506678_1420MHz_GammaCygni.txt", Hist_CO, isRaDec)
        ContourLevel_1 = 11000.
        ContourLevel_2 = 13000.
        ContourLevel_3 = 15000.
    if sys.argv[1]=='MGRO_J1908_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv22338613766121_SwiftBat14keV_J1908.txt", Hist_CO, isRaDec)
        ContourLevel_1 = 1.0
        ContourLevel_2 = 2.0
        ContourLevel_3 = 4.0
    if sys.argv[1]=='MAGIC_J1857_ON':
        Hist_CO = GetSkyViewMap("MWL_maps/skv19487414377965_1420MHz_J1857.txt", Hist_CO, isRaDec)
        ContourLevel_1 = 7000.
        ContourLevel_2 = 9000.
        ContourLevel_3 = 11000.
    Hist_CO = reflectXaxis(Hist_CO)
    Hist_CO.SetContour(3)
    Hist_CO.SetContourLevel(0,ContourLevel_1)
    Hist_CO.SetContourLevel(1,ContourLevel_2)
    Hist_CO.SetContourLevel(2,ContourLevel_3)
    Hist_CO.SetLineColor(0)

    other_star_labels = []
    other_star_markers = []
    other_star_names = []
    other_star_significance = []
    bright_star_labels = []
    bright_star_markers = []
    faint_star_labels = []
    faint_star_markers = []
    star_range = 2.5
    if xtitle=="RA":
        for star in range(0,len(other_stars)):
            if pow(source_ra-other_star_coord[star][0],2)+pow(source_dec-other_star_coord[star][1],2)>star_range*star_range: continue
            star_significance = FindLocalMaximum(Hist_Skymap, -other_star_coord[star][0], other_star_coord[star][1])
            #if star_significance<2.0: continue
            other_star_markers += [ROOT.TMarker(-other_star_coord[star][0],other_star_coord[star][1],2)]
            other_star_labels += [ROOT.TLatex(-other_star_coord[star][0]-0.15,other_star_coord[star][1]+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
            other_star_names += [other_stars[star]]
            other_star_significance += [star_significance]
        for star in range(0,len(bright_star_ra)):
            bright_star_markers += [ROOT.TMarker(-bright_star_ra[star],bright_star_dec[star],30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-bright_star_ra[star]-0.15,bright_star_dec[star]+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
        for star in range(0,len(faint_star_ra)):
            faint_star_markers += [ROOT.TMarker(-faint_star_ra[star],faint_star_dec[star],30)]
            faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
            faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
            faint_star_labels += [ROOT.TLatex(-faint_star_ra[star]-0.15,faint_star_dec[star]+0.15,'b-mag %s'%(faint_star_brightness[star]))]
            faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
            faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)
    else:
        for star in range(0,len(other_stars)):
            gal_l, gal_b = ConvertRaDecToGalactic(other_star_coord[star][0],other_star_coord[star][1])
            if pow(source_l-gal_l,2)+pow(source_b-gal_b,2)>star_range*star_range: continue
            star_significance = FindLocalMaximum(Hist_Skymap, -gal_l, gal_b)
            #if star_significance<2.0: continue
            other_star_markers += [ROOT.TMarker(-gal_l,gal_b,2)]
            other_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,other_stars[star])]
            other_star_markers[len(other_star_markers)-1].SetMarkerSize(1.5)
            other_star_labels[len(other_star_labels)-1].SetTextSize(0.03)
            other_star_names += [other_stars[star]]
            other_star_significance += [star_significance]
        for star in range(0,len(bright_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(bright_star_ra[star],bright_star_dec[star])
            bright_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            bright_star_markers[len(bright_star_markers)-1].SetMarkerColor(2)
            bright_star_markers[len(bright_star_markers)-1].SetMarkerSize(1.5)
            bright_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(bright_star_brightness[star]))]
            bright_star_labels[len(bright_star_labels)-1].SetTextColor(2)
            bright_star_labels[len(bright_star_labels)-1].SetTextSize(0.02)
        for star in range(0,len(faint_star_ra)):
            gal_l, gal_b = ConvertRaDecToGalactic(faint_star_ra[star],faint_star_dec[star])
            faint_star_markers += [ROOT.TMarker(-gal_l,gal_b,30)]
            faint_star_markers[len(faint_star_markers)-1].SetMarkerColor(3)
            faint_star_markers[len(faint_star_markers)-1].SetMarkerSize(1.5)
            faint_star_labels += [ROOT.TLatex(-gal_l-0.15,gal_b+0.15,'b-mag %s'%(faint_star_brightness[star]))]
            faint_star_labels[len(faint_star_labels)-1].SetTextColor(3)
            faint_star_labels[len(faint_star_labels)-1].SetTextSize(0.02)

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    #pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    #pad1.SetTopMargin(0.15)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for star in range(0,len(other_star_names)):
        legend.AddEntry(other_stars[star],'%s (%0.1f#sigma)'%(other_star_names[star],other_star_significance[star]))

    Hist_Contour = Hist_Skymap.Clone()
    Hist_Contour.Reset()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            Hist_Contour.SetBinContent(bx+1,by+1,(Hist_Skymap.GetBinContent(bx+1,by+1)))
    Hist_Contour.SetContour(3)
    Hist_Contour.SetContourLevel(0,3)
    Hist_Contour.SetContourLevel(1,4)
    Hist_Contour.SetContourLevel(2,5)

    Hist_Skymap.GetYaxis().SetTitle(ytitle)
    Hist_Skymap.GetXaxis().SetTitle(xtitle)
    Hist_Skymap.GetZaxis().SetTitle('Significance')
    Hist_Skymap.GetZaxis().SetTitleOffset(1.1)
    Hist_Skymap.SetMaximum(5)
    Hist_Skymap.SetMinimum(-5)
    #Hist_Skymap.SetMaximum(3)
    #Hist_Skymap.SetMinimum(-3)

    pad1.cd()
    Hist_Skymap.Draw("COL4Z")
    Hist_Contour.Draw("CONT3 same")
    #Hist_HAWC.Draw("CONT3 same")
    #Hist_CO.Draw("CONT3 same")
    Hist_Skymap.GetXaxis().SetLabelOffset(999)
    Hist_Skymap.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap.GetXaxis().GetXmin()
    x2 = Hist_Skymap.GetXaxis().GetXmax()
    y1 = Hist_Skymap.GetYaxis().GetXmin()
    y2 = Hist_Skymap.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        if nth_roi==0: continue
        if nth_roi==1: continue
        mycircles[nth_roi].Draw("same")
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'max. %0.1f#sigma (syst = %0.1f%%)'%(max_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.50,'E >%0.1f GeV (%.1f hrs)'%(energy_bin[energy_bin_cut_low],exposure_hours) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.30,'MJD %s-%s'%(MJD_Start,MJD_End) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/SkymapSig_%s_%s.png'%(name,selection_tag))
    OutputToTxtFile(Hist_Skymap,'output_plots/TxtSkymapSig_%s_%s.txt'%(name,selection_tag))


    Hist_Skymap_Excess = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            Hist_Skymap_Excess.SetBinContent(bx+1,by+1,NSR-NBkg)

    pad1.cd()
    Hist_Skymap_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_Excess.Draw("COL4Z")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        if nth_roi==0: continue
        if nth_roi==1: continue
        mycircles[nth_roi].Draw("same")
    Hist_Skymap_Excess.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Excess.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_Excess.GetXaxis().GetXmin()
    x2 = Hist_Skymap_Excess.GetXaxis().GetXmax()
    y1 = Hist_Skymap_Excess.GetYaxis().GetXmin()
    y2 = Hist_Skymap_Excess.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_Excess.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapExcess_%s_%s.png'%(name,selection_tag))


    pad1.cd()
    Hist_Exposure.GetYaxis().SetTitle(ytitle)
    Hist_Exposure.GetXaxis().SetTitle(xtitle)
    Hist_Exposure.Draw("COL4Z")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        if nth_roi==0: continue
        if nth_roi==1: continue
        mycircles[nth_roi].Draw("same")
    Hist_Exposure.GetXaxis().SetLabelOffset(999)
    Hist_Exposure.GetXaxis().SetTickLength(0)
    x1 = Hist_Exposure.GetXaxis().GetXmin()
    x2 = Hist_Exposure.GetXaxis().GetXmax()
    y1 = Hist_Exposure.GetYaxis().GetXmin()
    y2 = Hist_Exposure.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Exposure.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/Skymap_Exposure_%s_%s.png'%(name,selection_tag))

    Hist_Skymap_Ratio = Hist_SR.Clone()
    Hist_Skymap_Ratio.Add(Hist_Bkg,-1.)
    Hist_Skymap_Ratio.Divide(Hist_Exposure)
    for bx in range(0,Hist_Skymap_Ratio.GetNbinsX()):
        for by in range(0,Hist_Skymap_Ratio.GetNbinsY()):
            #if Hist_Bkg.GetBinContent(bx+1,by+1)<1.5 and Hist_Skymap.GetBinContent(bx+1,by+1)<4.0:
            if Hist_Exposure.GetBinContent(bx+1,by+1)<1.0:
                Hist_Skymap_Ratio.SetBinContent(bx+1,by+1,0.)

    max_nbins = 60
    best_nbins = 5
    Hist_SR_reflect = reflectXaxis(Hist_SR_input)
    Hist_Bkg_reflect = reflectXaxis(Hist_Bkg_input)
    Hist_Syst_reflect = reflectXaxis(Hist_Syst_input)
    Hist_RBM_reflect = reflectXaxis(Hist_RBM_input)
    Hist_SR_reflect.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Bkg_reflect.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Syst_reflect.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_RBM_reflect.Rebin2D(max_nbins/best_nbins,max_nbins/best_nbins)
    Hist_Skymap_reflect = GetSignificanceMap(Hist_SR_reflect,Hist_Bkg_reflect,Hist_Syst_reflect,syst_method)
    Hist_Skymap_ImportantBins = Hist_Skymap_reflect.Clone()
    for bx in range(0,Hist_Skymap_ImportantBins.GetNbinsX()):
        for by in range(0,Hist_Skymap_ImportantBins.GetNbinsY()):
            if Hist_Skymap_reflect.GetBinContent(bx+1,by+1)>5.0:
                Hist_Skymap_ImportantBins.SetBinContent(bx+1,by+1,1.)
            else:
                Hist_Skymap_ImportantBins.SetBinContent(bx+1,by+1,0.)

    pad1.cd()
    Hist_Skymap_Ratio.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_Ratio.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_Ratio.GetZaxis().SetTitle('Rate [evt/hour]')
    Hist_Skymap_Ratio.GetZaxis().SetTitleOffset(1.1)
    Hist_Skymap_Ratio_ImportantBins = Hist_Skymap_Ratio.Clone()
    Hist_Skymap_Ratio.SetMinimum(0.)
    Hist_Skymap_Ratio_ImportantBins.SetMinimum(0.)
    for bx in range(0,Hist_Skymap_Ratio.GetNbinsX()):
        for by in range(0,Hist_Skymap_Ratio.GetNbinsY()):
            bin_location_x = Hist_Skymap_Ratio.GetXaxis().GetBinCenter(bx+1)
            bin_location_y = Hist_Skymap_Ratio.GetYaxis().GetBinCenter(by+1)
            big_binx = Hist_Skymap_ImportantBins.GetXaxis().FindBin(bin_location_x)
            big_biny = Hist_Skymap_ImportantBins.GetYaxis().FindBin(bin_location_y)
            if (Hist_Skymap_ImportantBins.GetBinContent(big_binx,big_biny)==0.):
                Hist_Skymap_Ratio_ImportantBins.SetBinContent(bx+1,by+1,0.)
    Hist_Skymap_Ratio.SetMaximum(Hist_Skymap_Ratio_ImportantBins.GetMaximum())
    Hist_Skymap_Ratio.Draw("CONT1Z")
    Hist_Skymap_Ratio_ImportantBins.Draw("COL4Z same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        if nth_roi==0: continue
        if nth_roi==1: continue
        mycircles[nth_roi].Draw("same")
    Hist_Skymap_Ratio.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Ratio.GetXaxis().SetTickLength(0)
    Hist_Skymap_Ratio_ImportantBins.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_Ratio_ImportantBins.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_Ratio.GetXaxis().GetXmin()
    x2 = Hist_Skymap_Ratio.GetXaxis().GetXmax()
    y1 = Hist_Skymap_Ratio.GetYaxis().GetXmin()
    y2 = Hist_Skymap_Ratio.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_Ratio.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapRatio_%s_%s.png'%(name,selection_tag))

    MapEdge_left = Hist_Skymap_Excess.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Skymap_Excess.GetXaxis().GetBinLowEdge(Hist_Skymap_Excess.GetNbinsX()+1)
    MapEdge_lower = Hist_Skymap_Excess.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Skymap_Excess.GetYaxis().GetBinLowEdge(Hist_Skymap_Excess.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapSize_x = (MapEdge_right-MapEdge_left)/2.
    MapSize_y = (MapEdge_upper-MapEdge_lower)/2.
    Hist_Skymap_zoomin = ROOT.TH2D("Hist_Skymap_zoomin","",Skymap_nbins/3,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,Skymap_nbins/3,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_Contour_zoomin = ROOT.TH2D("Hist_Contour_zoomin","",Skymap_nbins/3,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,Skymap_nbins/3,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_HAWC_zoomin = ROOT.TH2D("Hist_HAWC_zoomin","",Skymap_nbins/3,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,Skymap_nbins/3,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    Hist_CO_zoomin = ROOT.TH2D("Hist_CO_zoomin","",Skymap_nbins/3,MapCenter_x-MapSize_x/3.,MapCenter_x+MapSize_x/3.,Skymap_nbins/3,MapCenter_y-MapSize_y/3.,MapCenter_y+MapSize_y/3)
    for bx in range(0,Hist_Skymap_Excess.GetNbinsX()):
        for by in range(0,Hist_Skymap_Excess.GetNbinsY()):
            bx_center = Hist_Skymap_Excess.GetXaxis().GetBinCenter(bx+1)
            by_center = Hist_Skymap_Excess.GetYaxis().GetBinCenter(by+1)
            bx2 = Hist_Skymap_zoomin.GetXaxis().FindBin(bx_center)
            by2 = Hist_Skymap_zoomin.GetYaxis().FindBin(by_center)
            #Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Skymap_Excess.GetBinContent(bx+1,by+1))
            Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Skymap_Ratio.GetBinContent(bx+1,by+1))
            #Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Skymap.GetBinContent(bx+1,by+1))
            #Hist_Skymap_zoomin.SetBinContent(bx2,by2,Hist_Bkg.GetBinContent(bx+1,by+1))
            Hist_Contour_zoomin.SetBinContent(bx2,by2,Hist_Contour.GetBinContent(bx+1,by+1))
            Hist_HAWC_zoomin.SetBinContent(bx2,by2,Hist_HAWC.GetBinContent(bx+1,by+1))
            Hist_CO_zoomin.SetBinContent(bx2,by2,Hist_CO.GetBinContent(bx+1,by+1))
    Hist_Contour_zoomin.SetContour(3)
    Hist_Contour_zoomin.SetContourLevel(0,3)
    Hist_Contour_zoomin.SetContourLevel(1,4)
    Hist_Contour_zoomin.SetContourLevel(2,5)
    #Hist_Skymap_zoomin.SetMaximum(0.6)
    #Hist_Skymap_zoomin.SetMinimum(-0.2)
    #Hist_Skymap_zoomin.SetMaximum(5)
    #Hist_Skymap_zoomin.SetMinimum(-5)

    pad1.cd()
    Hist_Skymap_zoomin.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_zoomin.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_zoomin.GetZaxis().SetTitle('Rate [evt/hour]')
    #Hist_Skymap_zoomin.GetZaxis().SetTitle('S/B ratio')
    Hist_Skymap_zoomin.GetZaxis().SetTitleOffset(1.1)
    Hist_Skymap_zoomin.Draw("COL4Z")
    Hist_Contour_zoomin.Draw("CONT3 same")
    Hist_HAWC_zoomin.SetLineColor(0)
    Hist_HAWC_zoomin.SetContour(3)
    Hist_HAWC_zoomin.SetContourLevel(0,5)
    Hist_HAWC_zoomin.SetContourLevel(1,10)
    Hist_HAWC_zoomin.SetContourLevel(2,15)
    #Hist_HAWC_zoomin.Draw("CONT3 same")
    Hist_CO_zoomin.SetLineColor(2)
    Hist_CO_zoomin.SetContour(3)
    Hist_CO_zoomin.SetContourLevel(0,ContourLevel_1)
    Hist_CO_zoomin.SetContourLevel(1,ContourLevel_2)
    Hist_CO_zoomin.SetContourLevel(2,ContourLevel_3)
    Hist_CO_zoomin.Draw("CONT3 same")
    for star in range(0,len(other_star_markers)):
        other_star_markers[star].Draw("same")
        other_star_labels[star].Draw("same")
    for star in range(0,len(bright_star_markers)):
        bright_star_markers[star].Draw("same")
        bright_star_labels[star].Draw("same")
    for star in range(0,len(faint_star_markers)):
        faint_star_markers[star].Draw("same")
        faint_star_labels[star].Draw("same")
    mycircles = []
    for nth_roi in range(0,len(roi_ra)):
        mycircles += [ROOT.TEllipse(-1.*roi_ra[nth_roi],roi_dec[nth_roi],roi_radius[nth_roi])]
        mycircles[nth_roi].SetFillStyle(0)
        mycircles[nth_roi].SetLineColor(2)
        if nth_roi==0: continue
        if nth_roi==1: continue
        mycircles[nth_roi].Draw("same")
    Hist_Skymap_zoomin.GetXaxis().SetLabelOffset(999)
    Hist_Skymap_zoomin.GetXaxis().SetTickLength(0)
    x1 = Hist_Skymap_zoomin.GetXaxis().GetXmin()
    x2 = Hist_Skymap_zoomin.GetXaxis().GetXmax()
    y1 = Hist_Skymap_zoomin.GetYaxis().GetXmin()
    y2 = Hist_Skymap_zoomin.GetYaxis().GetXmax()
    IncValues = ROOT.TF1( "IncValues", "-x", -x2, -x1 )
    raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
    raLowerAxis.SetLabelSize(Hist_Skymap_zoomin.GetXaxis().GetLabelSize())
    raLowerAxis.Draw()
    canvas.SaveAs('output_plots/SkymapZoomin_%s_%s.png'%(name,selection_tag))

def FindLocalMaximum(Hist_sig, init_x, init_y):

    max_sig = 0.
    for bx in range(0,Hist_sig.GetNbinsX()):
        for by in range(0,Hist_sig.GetNbinsY()):
            bin_x = Hist_sig.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_sig.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_x-init_x,2)+pow(bin_y-init_y,2),0.5)
            if distance < 0.1:
                if max_sig < Hist_sig.GetBinContent(bx+1,by+1):
                    max_sig = Hist_sig.GetBinContent(bx+1,by+1)
    return max_sig

def GetExtention(Hist_data, Hist_bkgd, Hist_sig, Hist_exposure, highlight_threshold, init_x, init_y):

    Hist_Excess = Hist_data.Clone()
    Hist_Excess.Add(Hist_bkgd,-1.)
    Hist_Excess.Divide(Hist_exposure)

    max_exposure = Hist_exposure.GetMaximum()
    if max_exposure==0.:
        return 0., 0., 0.
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_x = Hist_Excess.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_Excess.GetYaxis().GetBinCenter(by+1)
            distance = pow(pow(bin_x-init_x,2)+pow(bin_y-init_y,2),0.5)
            exposure = Hist_exposure.GetBinContent(bx+1,by+1)
            if distance > 1.0:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
            if exposure/max_exposure<0.05:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
            if Hist_Excess.GetBinContent(bx+1,by+1)<0.: 
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
            if not Hist_sig.GetBinContent(bx+1,by+1)>=highlight_threshold: 
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
    xx, yy, zz = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
    maxbin = Hist_Excess.GetMaximumBin()
    Hist_Excess.GetBinXYZ(maxbin, xx, yy, zz)
    print 'max excess at %s, %s'%(Hist_Excess.GetXaxis().GetBinCenter(xx),Hist_Excess.GetYaxis().GetBinCenter(yy))
    excess_center_x_init = Hist_Excess.GetXaxis().GetBinCenter(xx)
    excess_center_y_init = Hist_Excess.GetYaxis().GetBinCenter(yy)

    total_mass = 0.
    total_mass_distance = 0.
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_x = Hist_Excess.GetXaxis().GetBinCenter(bx+1)
            bin_y = Hist_Excess.GetYaxis().GetBinCenter(by+1)
            mass = Hist_Excess.GetBinContent(bx+1,by+1)
            mass_distance = mass*mass*(pow(bin_x-excess_center_x_init,2)+pow(bin_y-excess_center_y_init,2))
            total_mass += mass*mass
            total_mass_distance += mass_distance
    if total_mass>0.: total_mass_distance = pow(total_mass_distance/total_mass,0.5)
    else: total_mass_distance = 0.
    excess_radius_init = total_mass_distance
    print 'excess_radius_init = %s'%(excess_radius_init)

    Hist_Excess.Reset()
    Hist_Excess.Add(Hist_data)
    Hist_Excess.Add(Hist_bkgd,-1.)
    Hist_Excess.Divide(Hist_exposure)
    for bx in range(0,Hist_Excess.GetNbinsX()):
        for by in range(0,Hist_Excess.GetNbinsY()):
            bin_delta_ra = Hist_Excess.GetXaxis().GetBinCenter(bx+1)-excess_center_x_init
            bin_delta_dec = Hist_Excess.GetYaxis().GetBinCenter(by+1)-excess_center_y_init
            if pow(bin_delta_ra*bin_delta_ra+bin_delta_dec*bin_delta_dec,0.5)>(3.*excess_radius_init):
                Hist_Excess.SetBinContent(bx+1,by+1,0.)
            exposure = Hist_exposure.GetBinContent(bx+1,by+1)
            if exposure/max_exposure<0.05:
                Hist_Excess.SetBinContent(bx+1,by+1,0.)

    MapEdge_left = Hist_Excess.GetXaxis().GetBinLowEdge(1)
    MapEdge_right = Hist_Excess.GetXaxis().GetBinLowEdge(Hist_Excess.GetNbinsX()+1)
    MapEdge_lower = Hist_Excess.GetYaxis().GetBinLowEdge(1)
    MapEdge_upper = Hist_Excess.GetYaxis().GetBinLowEdge(Hist_Excess.GetNbinsY()+1)
    MapCenter_x = (MapEdge_right+MapEdge_left)/2.
    MapCenter_y = (MapEdge_upper+MapEdge_lower)/2.
    MapBinSize = Hist_Excess.GetXaxis().GetBinLowEdge(2)-Hist_Excess.GetXaxis().GetBinLowEdge(1)
    Func_Gauss2D = ROOT.TF2("Func_Gauss2D","[0]*TMath::Gaus(x,[1],[3])*TMath::Gaus(y,[2],[3])",MapEdge_left,MapEdge_right,MapEdge_lower,MapEdge_upper)
    excess_radius = excess_radius_init
    excess_center_x = excess_center_x_init
    excess_center_y = excess_center_y_init
    amplitude = Hist_Excess.GetMaximum()
    Func_Gauss2D.SetParameters(amplitude,excess_center_x,excess_center_y,excess_radius)
    Hist_Excess.Fit("Func_Gauss2D","N")
    amplitude = Func_Gauss2D.GetParameter(0)
    amplitude_err = Func_Gauss2D.GetParError(0)
    excess_center_x = Func_Gauss2D.GetParameter(1)
    excess_center_x_err = Func_Gauss2D.GetParError(1)
    excess_center_y = Func_Gauss2D.GetParameter(2)
    excess_center_y_err = Func_Gauss2D.GetParError(2)
    excess_radius = Func_Gauss2D.GetParameter(3)
    excess_radius_err = Func_Gauss2D.GetParError(3)

    return excess_center_x, excess_center_y, excess_radius

def Make2DProjectionPlot(Hist_Data,xtitle,ytitle,name,doProj):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    #pad1.SetGrid()
    pad1.Draw()
    pad1.cd()
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.Draw("COL4Z")
    #Hist_Data.Draw("CONT3 same")
    bins = []
    for b in range(0,Hist_Data.GetNbinsX()+1):
        bins += [Hist_Data.GetXaxis().GetBinLowEdge(b+1)]
    Hist_1D = ROOT.TH1D("Hist_1D","",len(bins)-1,array('d',bins))
    for b in range(0,Hist_Data.GetNbinsX()):
        hist_temp = Hist_Data.ProjectionY("hist_temp",b+1,b+1)
        Hist_1D.SetBinContent(b+1,hist_temp.GetMean())
        Hist_1D.SetBinError(b+1,hist_temp.GetRMS())
    Hist_1D.SetLineColor(2)
    if doProj: Hist_1D.Draw("E same")

    canvas.SaveAs('output_plots/%s.png'%(name))

def MatrixDecompositionDemo(name):

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

    line1 = ROOT.TLine(MSCL_plot_lower,MSCW_blind_cut,MSCL_blind_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line2 = ROOT.TLine(MSCL_blind_cut,MSCW_plot_lower,MSCL_blind_cut,MSCW_blind_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line3 = ROOT.TLine(MSCL_plot_lower,MSCW_chi2_upper,MSCL_chi2_upper,MSCW_chi2_upper)
    line3.SetLineStyle(10)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line4 = ROOT.TLine(MSCL_chi2_upper,MSCW_plot_lower,MSCL_chi2_upper,MSCW_chi2_upper)
    line4.SetLineStyle(10)
    line4.SetLineColor(2)
    line4.SetLineWidth(2)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnData_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_OnData_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_OnData_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/OnData_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{OFF}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_OnDark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_OnDark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_OnDark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/OnDark_%s_%s.png'%(name,selection_tag))

    entry = 1
    for entry in range(0,3):
        Hists = []
        legends = []
        colors = []
        Hists += [Hist2D_Coeff_Data_Sum.ProjectionX("hist_temp_x",Hist2D_Coeff_Data_Sum.GetNbinsY()-entry,Hist2D_Coeff_Data_Sum.GetNbinsY()-entry)]
        legends += ['C_{i,1}']
        colors += [1]
        Hists += [Hist2D_Coeff_Data_Sum.ProjectionY("hist_temp_y",Hist2D_Coeff_Data_Sum.GetNbinsX()-entry,Hist2D_Coeff_Data_Sum.GetNbinsX()-entry)]
        legends += ['C_{1,j}']
        colors += [2]
        MakeMultiplePlot(Hists,legends,colors,'entry','C_{ij}','Coeff_C_entry%s_%s'%(entry,name),0,0,False,False)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'t_{ij} coefficents' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.20,'(Zenith = %0.1f #pm %0.1f, NSB = %0.1f #pm %0.1f)'%(Zenith_mean_data,Zenith_RMS_data,NSB_mean_data,NSB_RMS_data) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.2)
    lumilab2.Draw()
    pad1.cd()
    #pad1.SetLogz()
    Hist2D_Coeff_Bkgd_Sum.GetYaxis().SetTitle('n')
    Hist2D_Coeff_Bkgd_Sum.GetXaxis().SetTitle('k')
    Hist2D_Coeff_Bkgd_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_C_bkgd_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.70,'t_{ij} coefficents')
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.20,'(Zenith = %0.1f #pm %0.1f, NSB = %0.1f #pm %0.1f)'%(Zenith_mean_data,Zenith_RMS_data,NSB_mean_data,NSB_RMS_data) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.2)
    lumilab2.Draw()
    pad1.cd()
    #pad1.SetLogz()
    #Hist2D_Coeff_Data_Sum.SetMaximum(Hist2D_Coeff_Bkgd_Sum.GetMaximum());
    #Hist2D_Coeff_Data_Sum.SetMinimum(Hist2D_Coeff_Bkgd_Sum.GetMinimum());
    Hist2D_Coeff_Data_Sum.GetYaxis().SetTitle('n')
    Hist2D_Coeff_Data_Sum.GetXaxis().SetTitle('k')
    Hist2D_Coeff_Data_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_C_data_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#delta H_{ij} coefficents' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    pad1.SetLogz()
    #Hist2D_U_Proj_Sum.SetMaximum(1e-1);
    #Hist2D_U_Proj_Sum.SetMinimum(1e-4);
    Hist2D_U_Proj_Sum.GetYaxis().SetTitle('y')
    Hist2D_U_Proj_Sum.GetXaxis().SetTitle('x')
    Hist2D_U_Proj_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_U_proj_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#delta H_{ij} coefficents' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    pad1.SetLogz()
    #Hist2D_V_Proj_Sum.SetMaximum(1e-1);
    #Hist2D_V_Proj_Sum.SetMinimum(1e-4);
    Hist2D_V_Proj_Sum.GetYaxis().SetTitle('y')
    Hist2D_V_Proj_Sum.GetXaxis().SetTitle('x')
    Hist2D_V_Proj_Sum.Draw("COL4Z")
    canvas.SaveAs('output_plots/Coeff_V_proj_%s_%s.png'%(name,selection_tag))
    pad1.SetLogz(0)

    pad3.Clear()
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON} - #sum_{k=1}^{1} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank0_Data_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank0_Data_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank0_Data_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank0_Data_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sigma_{1} u_{1} v_{1}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank0_Dark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank0_Dark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank0_Dark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank0_Dark_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON} - #sum_{k=1}^{2} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank1_Data_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank1_Data_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank1_Data_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank1_Data_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sigma_{2} u_{2} v_{2}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank1_Dark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank1_Dark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank1_Dark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank1_Dark_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON} - #sum_{k=1}^{3} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank2_Data_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank2_Data_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank2_Data_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank2_Data_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sigma_{3} u_{3} v_{3}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank2_Dark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank2_Dark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank2_Dark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank2_Dark_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON} - #sum_{k=1}^{4} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank3_Data_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank3_Data_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank3_Data_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank3_Data_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sigma_{4} u_{4} v_{4}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank3_Dark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank3_Dark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank3_Dark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank3_Dark_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON} - #sum_{k=1}^{5} #sigma_{k} u_{k} v_{k}^{T}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank4_Data_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank4_Data_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank4_Data_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank4_Data_%s_%s.png'%(name,selection_tag))

    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'#sigma_{5} u_{5} v_{5}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_Rank4_Dark_Sum.GetYaxis().SetTitle('MSCW')
    Hist2D_Rank4_Dark_Sum.GetXaxis().SetTitle('MSCL')
    Hist2D_Rank4_Dark_Sum.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    canvas.SaveAs('output_plots/Rank4_Dark_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrDark = ROOT.TH2D("Hist2D_ErrDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrDark.Add(Hist2D_OnData_Sum)
    Hist2D_ErrDark.Add(Hist2D_OnDark_Sum,-1)
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrDark.GetNbinsX()):
        for biny in range (0,Hist2D_ErrDark.GetNbinsY()):
            if Hist2D_ErrDark.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrDark.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnDark_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    error_avg = sum_error/sum_bins
    sum_var = 0.
    for binx in range (0,Hist2D_ErrDark.GetNbinsX()):
        for biny in range (0,Hist2D_ErrDark.GetNbinsY()):
            if Hist2D_ErrDark.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrDark.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnDark_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd-error_avg,2)
    error_var = sum_var/sum_bins
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON}-M^{OFF}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrDark.GetYaxis().SetTitle('MSCW')
    Hist2D_ErrDark.GetXaxis().SetTitle('MSCL')
    Hist2D_ErrDark.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    #line3.Draw("same")
    #line4.Draw("same")
    canvas.SaveAs('output_plots/ErrorDark_%s_%s.png'%(name,selection_tag))

    Hist2D_ErrBkgd = ROOT.TH2D("Hist2D_ErrBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_ErrBkgd.Add(Hist2D_OnData_Sum)
    Hist2D_ErrBkgd.Add(Hist2D_OnBkgd_Sum,-1)
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    error_avg = sum_error/sum_bins
    sum_var = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)<MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)<MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd-error_avg,2)
    error_var = sum_var/sum_bins
    sum_bins = 0.
    sum_error = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)>MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)>MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_bins += 1.
                    sum_error += (n_data-n_bkgd)
    sum_var = 0.
    for binx in range (0,Hist2D_ErrBkgd.GetNbinsX()):
        for biny in range (0,Hist2D_ErrBkgd.GetNbinsY()):
            if Hist2D_ErrBkgd.GetXaxis().GetBinCenter(binx+1)>MSCL_blind_cut:
                if Hist2D_ErrBkgd.GetYaxis().GetBinCenter(biny+1)>MSCW_blind_cut:
                    n_data = Hist2D_OnData_Sum.GetBinContent(binx+1,biny+1)
                    n_bkgd = Hist2D_OnBkgd_Sum.GetBinContent(binx+1,biny+1)
                    sum_var += pow(n_data-n_bkgd,2)
    pad3.cd()
    lumilab1 = ROOT.TLatex(0.15,0.50,'M^{ON}-#tilde{M}^{ON}' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.3)
    lumilab1.Draw()
    pad1.cd()
    Hist2D_ErrBkgd.SetMaximum(Hist2D_ErrDark.GetMaximum());
    Hist2D_ErrBkgd.SetMinimum(Hist2D_ErrDark.GetMinimum());
    Hist2D_ErrBkgd.GetYaxis().SetTitle('MSCW')
    Hist2D_ErrBkgd.GetXaxis().SetTitle('MSCL')
    Hist2D_ErrBkgd.Draw("COL4Z")
    line1.Draw("same")
    line2.Draw("same")
    #line3.Draw("same")
    #line4.Draw("same")
    canvas.SaveAs('output_plots/ErrorBkgd_%s_%s.png'%(name,selection_tag))

def MakeOneHistPlot(Hist,title_x,title_y,name,logy):
    
    global max_chi2_diff2_position_this_energy

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
    #if logy: pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad3.cd()
    #lumilab2 = ROOT.TLatex(0.15,0.50,'#epsilon = #sum_{#gamma region} (M_{data} - #sum_{k=1}^{k=n} #lambda_{k} q_{k} p_{k}^{T})  / #sum_{#gamma region} M_{data}' )
    #lumilab2 = ROOT.TLatex(0.15,0.50,'#epsilon = #sum_{#gamma region} (#lambda_{n} q_{n} p_{n}^{T})  / #sum_{#gamma region} M_{data}' )
    #lumilab2.SetNDC()
    #lumilab2.SetTextSize(0.2)
    #lumilab2.Draw()

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

    #print 'plot max_chi2_diff2_position_this_energy = %0.3f'%(max_chi2_diff2_position_this_energy)
    #myline = ROOT.TLine(max_chi2_diff2_position_this_energy,0.,max_chi2_diff2_position_this_energy,100.)
    #myline.SetLineStyle(1)
    #myline.SetLineColor(2)
    #myline.SetLineWidth(2)
    #myline.Draw("same")

    if Hist.GetXaxis().GetLabelOffset()==999:
        x1 = Hist.GetXaxis().GetXmin()
        x2 = Hist.GetXaxis().GetXmax()
        y1 = Hist.GetYaxis().GetXmin()
        y2 = Hist.GetYaxis().GetXmax()
        IncValues = ROOT.TF1( "IncValues", "x", 0, 256 )
        raLowerAxis = ROOT.TGaxis( x1, y1, x2, y1,"IncValues", 510, "+")
        raLowerAxis.SetLabelSize(Hist.GetXaxis().GetLabelSize())
        raLowerAxis.Draw()

    pad3.cd()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeAccumulatedSignificancePlot(hist_data,hist_bkgd,hist_syst,direction,title_x,title_y,name):

    hist_sig = hist_data.Clone()
    hist_syst_content = hist_data.Clone()
    for binx in range(0,hist_data.GetNbinsX()):
        hist_syst_content.SetBinContent(binx+1,hist_syst.GetBinError(binx+1))
    for binx in range(0,hist_data.GetNbinsX()):
        data = 0.
        bkgd = 0.
        syst = 0.
        if direction==0:
            data = hist_data.Integral(1,binx+1)
            bkgd = hist_bkgd.Integral(1,binx+1)
            syst = hist_syst_content.Integral(1,binx+1)
        else:
            data = hist_data.Integral(binx+1,hist_data.GetNbinsX())
            bkgd = hist_bkgd.Integral(binx+1,hist_data.GetNbinsX())
            syst = hist_syst_content.Integral(binx+1,hist_data.GetNbinsX())
        sig = CalculateSignificance(data-bkgd,bkgd,syst)
        print 'data = %0.2f, bkgd = %0.2f, syst = %0.2f, sig = %0.2f'%(data,bkgd,syst,sig)
        hist_sig.SetBinContent(binx+1,sig)
        hist_sig.SetBinError(binx+1,0.)
    MakeOneHistPlot(hist_sig,title_x,title_y,'%s_%s'%(name,selection_tag),False)

def MakeRankResidualPlots(name):

    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

    data_integral = Hist2D_OnData_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)
    data_rank_integral = []
    data_rank_integral += [Hist2D_Rank0_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank1_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank2_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank3_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]
    data_rank_integral += [Hist2D_Rank4_Data_Sum.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)]

    Hist_Residual_Gamma = ROOT.TH1D("Hist_Residual_Gamma","",5,0,5)
    for binx in range(0,len(data_rank_integral)):
        bkgd_integral = 0.
        #for rank in range (0,binx+1):
        #    bkgd_integral += data_rank_integral[rank]
        bkgd_integral = data_rank_integral[binx]
        #Hist_Residual_Gamma.SetBinContent(binx+1,abs(1.-bkgd_integral/data_integral))
        Hist_Residual_Gamma.SetBinContent(binx+1,abs(bkgd_integral/data_integral))
        ratio_err = bkgd_integral/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Gamma.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Gamma,'n ranks','Residual in gamma region','Residual_Gamma_%s_%s'%(name,selection_tag),True)

    data_integral = Hist2D_OnData_Sum.Integral()
    data_rank_integral = []
    data_rank_integral += [Hist2D_Rank0_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank1_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank2_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank3_Data_Sum.Integral()]
    data_rank_integral += [Hist2D_Rank4_Data_Sum.Integral()]

    Hist_Residual_Full = ROOT.TH1D("Hist_Residual_Full","",5,0,5)
    for binx in range(0,len(data_rank_integral)):
        bkgd_integral = 0.
        #for rank in range (0,binx+1):
        #    bkgd_integral += data_rank_integral[rank]
        bkgd_integral = data_rank_integral[binx]
        #Hist_Residual_Full.SetBinContent(binx+1,abs(1.-bkgd_integral/data_integral))
        Hist_Residual_Full.SetBinContent(binx+1,abs(bkgd_integral/data_integral))
        ratio_err = bkgd_integral/data_integral*1./pow(data_integral,0.5)
        Hist_Residual_Full.SetBinError(binx+1,ratio_err)

    MakeOneHistPlot(Hist_Residual_Full,'n ranks','Residual in gamma region','Residual_Full_%s_%s'%(name,selection_tag),True)

def MakeSystematicPlots(hists, colors, fillstyles, legends, exposure, name):

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
    lumilab1 = ROOT.TLatex(0.15,0.70,'Total exposure = %0.1f h'%(exposure) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(hists)):
        legend.AddEntry(hists[h],legends[h])
    legend.Draw("same")

    pad1.cd()

    hists[0].GetXaxis().SetTitleOffset(0.8)
    hists[0].GetXaxis().SetTitleSize(0.06)
    hists[0].GetXaxis().SetLabelSize(0.06)
    hists[0].GetYaxis().SetLabelSize(0.06)
    hists[0].GetYaxis().SetTitleOffset(1.2)
    hists[0].GetYaxis().SetTitleSize(0.06)
    hists[0].GetXaxis().SetTitle('relative error')
    hists[0].GetYaxis().SetTitle('number of measurements')
    for h in range(0,len(hists)):
        hists[h].SetFillStyle(fillstyles[h])
        hists[h].SetFillColor(colors[h])
        hists[h].SetLineColor(colors[h])
        hists[h].SetLineWidth(2)
        hists[h].Draw("hist same")
    pad3.cd()
    c_both.SaveAs('output_plots/%s.png'%(name))

def SystematicAnalysis():

    FilePath_List = []
    for elev in range(0,len(root_file_tags)):
        FilePath = "%s/Netflix_"%(folder_path)+"Syst_%s"%(root_file_tags[elev])+".root"
        FilePath_List += [FilePath]
        print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
        if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
        InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
        bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
        bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
        bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1

        Hist2D_OnData_AllGroups = []
        Hist2D_OnDark_AllGroups = []
        Hist2D_OnBkgd_AllGroups = []
        for nth_group in range(0,InfoTree.number_groups):
            Hist2D_OnData_AllGroups += [ROOT.TH2D("Hist2D_OnData_%s"%(nth_group),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnDark_AllGroups += [ROOT.TH2D("Hist2D_OnDark_%s"%(nth_group),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnBkgd_AllIterations = []
            for nth_iteration in range(0,20):
                Hist2D_OnBkgd_AllIterations += [ROOT.TH2D("Hist2D_OnBkgd_%s_%s"%(nth_group,nth_iteration),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
            Hist2D_OnBkgd_AllGroups += [Hist2D_OnBkgd_AllIterations]

        for e in range(0,len(energy_bin)-1):
            ErecS_lower_cut = energy_bin[e]
            ErecS_upper_cut = energy_bin[e+1]
            if ErecS_upper_cut<=energy_bin[energy_bin_cut_low]: continue
            if ErecS_lower_cut>=energy_bin[energy_bin_cut_up]: continue
            ErecS_lower_cut_int = int(ErecS_lower_cut)
            ErecS_upper_cut_int = int(ErecS_upper_cut)

            for nth_group in range(0,InfoTree.number_groups):
                HistName = "Hist_NthGroup_Data_MSCLW_V%s_ErecS%sto%s"%(nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                Hist2D_OnData_AllGroups[nth_group].Add(InputFile.Get(HistName))
                for nth_iteration in range(0,20):
                    HistName = "Hist_NthGroup_Bkgd_MSCLW_I%s_V%s_ErecS%sto%s"%(nth_iteration,nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                    Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration].Add(InputFile.Get(HistName))
                HistName = "Hist_NthGroup_Dark_MSCLW_V%s_ErecS%sto%s"%(nth_group,ErecS_lower_cut_int,ErecS_upper_cut_int)
                Hist2D_OnDark_AllGroups[nth_group].Add(InputFile.Get(HistName))

        Hist2D_Converge = ROOT.TH2D("Hist_Converge","",20,0,20,20,0.,0.01)
        for nth_iteration in range(0,20):
            bkgd_chi2 = 0.
            bkgd_syst_err = 0.
            dark_syst_err = 0.
            Hist_Syst_Bkgd = ROOT.TH1D("Hist_Syst_Bkgd","",19,-0.2,0.2)
            Hist_Syst_Dark = ROOT.TH1D("Hist_Syst_Dark","",19,-0.2,0.2)
            for nth_group in range(0,InfoTree.number_groups):
                Total_Data = Hist2D_OnData_AllGroups[nth_group].Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)
                if Total_Data == 0.: continue

                Hist2D_Diff = Hist2D_OnData_AllGroups[nth_group].Clone()
                Hist2D_Diff.Add(Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration],-1.)
                Hist2D_Chi2Diff = Hist2D_OnBkgd_AllGroups[nth_group][19].Clone()
                Hist2D_Chi2Diff.Add(Hist2D_OnBkgd_AllGroups[nth_group][nth_iteration],-1.)
                syst_this = Hist2D_Diff.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)/Total_Data
                Hist_Syst_Bkgd.Fill(syst_this)
                bkgd_syst_err += pow(syst_this,2)/InfoTree.number_groups
                chi2_this = 0.
                for binx in range(0,Hist2D_Chi2Diff.GetNbinsX()):
                    for biny in range(0,Hist2D_Chi2Diff.GetNbinsY()):
                        if binx>bin_upper_x : continue
                        if biny>bin_upper_y : continue
                        chi2_this += pow(Hist2D_Chi2Diff.GetBinContent(binx+1,biny+1),2)
                chi2_this = pow(chi2_this,0.5)/Total_Data
                bkgd_chi2 += pow(chi2_this,2)/InfoTree.number_groups
                Hist2D_Converge.Fill(nth_iteration,chi2_this)

                Hist2D_Diff = Hist2D_OnData_AllGroups[nth_group].Clone()
                Hist2D_Diff.Add(Hist2D_OnDark_AllGroups[nth_group],-1.)
                syst_this = Hist2D_Diff.Integral(bin_lower_x,bin_upper_x,bin_lower_y,bin_upper_y)/Total_Data
                Hist_Syst_Dark.Fill(syst_this)
                dark_syst_err += pow(syst_this,2)/InfoTree.number_groups

            bkgd_chi2 = pow(bkgd_chi2,0.5)
            bkgd_syst_err = pow(bkgd_syst_err,0.5)
            dark_syst_err = pow(dark_syst_err,0.5)

            hists = []
            colors = []
            fillstyles = []
            legends = []
            hists += [Hist_Syst_Bkgd]
            colors += [38]
            fillstyles += [3325]
            legends += ['MDM syst. = %0.3f'%(bkgd_syst_err)]
            hists += [Hist_Syst_Dark]
            colors += [46]
            fillstyles += [3352]
            legends += ['Initial syst. = %0.3f'%(dark_syst_err)]
            name = "SystAnalysis_I%s"%(nth_iteration)
            MakeSystematicPlots(hists, colors, fillstyles, legends, InfoTree.exposure_hours_usable, name)

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
        pad3.cd()
        lumilab1 = ROOT.TLatex(0.15,0.50,'#chi^{2} = #sum_{ij} (N^{nth}_{ij}-N^{20th}_{ij})^{2} / (#sum_{ij} N^{data}_{ij})^{2}' )
        lumilab1.SetNDC()
        lumilab1.SetTextSize(0.3)
        lumilab1.Draw()
        pad1.cd()
        Hist2D_Converge.GetYaxis().SetTitle('#sqrt{#chi^{2}}')
        Hist2D_Converge.GetXaxis().SetTitle('iterations')
        Hist2D_Converge.Draw("COL4Z")
        canvas.SaveAs('output_plots/Converge_%s.png'%(selection_tag))

def MakeGaussComparisonPlot(Hists,legends,colors,title,name):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.15)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

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
            Hists[h].GetXaxis().SetTitle(title)
            #Hists[h].GetXaxis().SetRangeUser(0,8)
            Hists[h].GetXaxis().SetRangeUser(-5,8)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #Hists[max_hist].SetMinimum(0)
    #Hists[max_hist].SetMaximum(1)
    Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            if h==0:
                Hists[h].Draw("hist same")
            else:
                Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],'%s, mean = %.2f, RMS = %.2f'%(legends[h],Hists[h].GetMean(),Hists[h].GetRMS()),"pl")
            #legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    pad1.SetLogy()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeSignificanceDistribution(Hist2D_Sig,Hist2D_SR,name):

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig_%s"%(name),"",65,-5,8)
    for bx in range(0,Hist2D_Sig.GetNbinsX()):
        for by in range(0,Hist2D_Sig.GetNbinsY()):
            if Hist2D_SR.GetBinContent(bx+1,by+1)==0: continue
            content = Hist2D_Sig.GetBinContent(bx+1,by+1)
            Hist_Sig.Fill(content)
    Hist_Model = ROOT.TH1D("Hist_Model_%s"%(name),"",65,-5,8)
    Hist_Model.FillRandom("func",10000*int(Hist_Sig.GetEntries()))
    Hist_Model.Scale(1./10000.)
    Hist_Model.SetMinimum(0.5)
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Ref. Gaussian']
    color_list += [2]
    Hist_list += [Hist_Sig]
    legend_list += ['Data']
    color_list += [4]
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))

    OutputFilePath = "output_plots/SigDist.root"
    OutputFile = ROOT.TFile(OutputFilePath,"update")
    Hist_Model.Write()
    Hist_Sig.Write()
    OutputFile.Close()

    file1 = open("output_plots/source_exposure.txt","a")
    file1.write('%s, (%0.1f hours) \n'%(name,exposure_hours))
    file1.close()

def GetCRcounts(name):

    bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
    bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
    bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
    for bx in range(1,Hist2D_OnData_Sum.GetNbinsX()+1):
        for by in range(1,Hist2D_OnData_Sum.GetNbinsY()+1):
            if bx>bin_upper_x and by<=bin_upper_y:
                counts = Hist2D_OnData_Sum.GetBinContent(bx,by)
                counts = max(counts,1.)
                Hist_CR_Counts_MSCL.Fill(min(5,math.log10(counts)))
            if bx<=bin_upper_x and by>bin_upper_y:
                counts = Hist2D_OnData_Sum.GetBinContent(bx,by)
                counts = max(counts,1.)
                Hist_CR_Counts_MSCW.Fill(min(5,math.log10(counts)))

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_CR_Counts_MSCW]
    legends += ['CR MSCW']
    colors += [1]
    Hists += [Hist_CR_Counts_MSCL]
    legends += ['CR MSCL']
    colors += [2]
    MakeMultiplePlot(Hists,legends,colors,'log10 counts per bin','N bins','CR_Counts_%s'%(name),0,0,False,False)

def SingleSourceAnalysis(source_list,doMap,doSmooth,e_low,e_up):

    global ErecS_lower_cut
    global ErecS_upper_cut
    global Syst_MDM
    global energy_bin_cut_low
    global energy_bin_cut_up
    global selection_tag
    global max_chi2_diff2_position_this_energy

    energy_bin_cut_low = e_low
    energy_bin_cut_up = e_up
    selection_tag = root_file_tags[0]
    selection_tag += '_E%sto%s'%(energy_bin_cut_low,energy_bin_cut_up)

    FilePath_List = []
    ResetStackedShowerHistograms()
    FilePath = "%s/Netflix_"%(folder_path)+source_list[0]+"_%s"%(root_file_tags[0])+".root"
    GetBrightStarInfo(FilePath)
    print 'len(bright_star_ra) = %s'%(len(bright_star_ra))
    for source in range(0,len(source_list)):
        source_name = source_list[source]
        for elev in range(0,len(root_file_tags)):
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+"_%s"%(root_file_tags[elev])+".root"
            FilePath_List += [FilePath]
            if not os.path.isfile(FilePath_List[len(FilePath_List)-1]):continue
            print 'Reading file %s'%(FilePath_List[len(FilePath_List)-1])
            InputFile = ROOT.TFile(FilePath_List[len(FilePath_List)-1])
            InfoTree = InputFile.Get("InfoTree")
            InfoTree.GetEntry(0)

            MSCW_blind_cut = InfoTree.MSCW_cut_blind
            MSCL_blind_cut = InfoTree.MSCL_cut_blind
            bin_lower_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_plot_lower)
            bin_upper_x = Hist2D_OnData.GetXaxis().FindBin(MSCL_blind_cut)-1
            bin_lower_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_plot_lower)
            bin_upper_y = Hist2D_OnData.GetYaxis().FindBin(MSCW_blind_cut)-1
            for e in range(0,len(energy_bin)-1):
                max_chi2_diff2_position_this_energy = max_chi2_diff2_position[e]
                ErecS_lower_cut = energy_bin[e]
                ErecS_upper_cut = energy_bin[e+1]
                print 'max_chi2_diff2_position_this_energy = %0.3f'%(max_chi2_diff2_position_this_energy)
                if ErecS_upper_cut<=energy_bin[energy_bin_cut_low]: continue
                if ErecS_lower_cut>=energy_bin[energy_bin_cut_up]: continue
                Syst_MDM = energy_syst[e]
                GetShowerHistogramsFromFile(FilePath_List[len(FilePath_List)-1])
                StackShowerHistograms()
                NormalizeEnergyHistograms(FilePath_List[len(FilePath_List)-1])
                StackEnergyHistograms()
                NormalizeTheta2Histograms(FilePath_List[len(FilePath_List)-1])
                StackTheta2Histograms()
                NormalizeSkyMapHistograms(FilePath_List[len(FilePath_List)-1],e)
                StackSkymapHistograms(e)

    source_name = sys.argv[1]
    MatrixDecompositionDemo(source_name)
    GetCRcounts(source_name)

    GetSourceInfo(FilePath_List)

    Syst_MDM = energy_syst[energy_bin_cut_low]
    #CalculateSystError_v3()
    #CorrectTheta2Histograms()

    PlotsStackedHistograms('%s%s'%(source_list[0],PercentCrab))

    #MakeRankResidualPlots('%s%s'%(source_list[0],PercentCrab))

    if not doMap: return

    Hist_OnData_Skymap_smooth = Hist_OnData_Skymap_Sum
    Hist_OnBkgd_Skymap_smooth = Hist_OnBkgd_Skymap_Sum
    Hist_OnBkgd_Skymap_Syst_MDM_smooth = Hist_OnBkgd_Skymap_Syst_MDM
    Hist_OnBkgd_Skymap_Syst_RBM_smooth = Hist_OnBkgd_Skymap_Syst_RBM
    Hist_OnData_Skymap_Galactic_smooth = Hist_OnData_Skymap_Galactic_Sum
    Hist_OnBkgd_Skymap_Galactic_smooth = Hist_OnBkgd_Skymap_Galactic_Sum
    Hist_OnBkgd_Skymap_Galactic_Syst_MDM_smooth = Hist_OnBkgd_Skymap_Galactic_Syst_MDM
    Hist_OnBkgd_Skymap_Galactic_Syst_RBM_smooth = Hist_OnBkgd_Skymap_Galactic_Syst_RBM
    if doSmooth:
        Hist_OnData_Skymap_smooth = Smooth2DMap(Hist_OnData_Skymap_Sum,smooth_size,False)
        Hist_OnBkgd_Skymap_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Sum,smooth_size,False)
        Hist_OnBkgd_Skymap_Syst_MDM_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Syst_MDM,smooth_size,False)
        Hist_OnBkgd_Skymap_Syst_RBM_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Syst_RBM,smooth_size,False)
        Hist_OnData_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnData_Skymap_Galactic_Sum,smooth_size,False)
        Hist_OnBkgd_Skymap_Galactic_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Sum,smooth_size,False)
        Hist_OnBkgd_Skymap_Galactic_Syst_MDM_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Syst_MDM,smooth_size,False)
        Hist_OnBkgd_Skymap_Galactic_Syst_RBM_smooth = Smooth2DMap(Hist_OnBkgd_Skymap_Galactic_Syst_RBM,smooth_size,False)

    skymap_bin_size_x = Hist_OnData_Skymap_Sum.GetXaxis().GetBinCenter(2)-Hist_OnData_Skymap_Sum.GetXaxis().GetBinCenter(1)
    skymap_bin_size_y = Hist_OnData_Skymap_Sum.GetYaxis().GetBinCenter(2)-Hist_OnData_Skymap_Sum.GetYaxis().GetBinCenter(1)
    event_rate = Hist_OnBkgd_Energy_CamCenter_Sum.Integral()/exposure_hours*(skymap_bin_size_x*skymap_bin_size_y)/(3.14*0.2*0.2)
    Hist_Exposure_Skymap = Hist_OnData_Skymap_Sum.Clone()
    for bx in range(1,Hist_Exposure_Skymap.GetNbinsX()+1):
        for by in range(1,Hist_Exposure_Skymap.GetNbinsY()+1):
            Hist_Exposure_Skymap.SetBinContent(bx,by,event_rate)
    Hist_Exposure_Skymap_smooth = Smooth2DMap(Hist_Exposure_Skymap,smooth_size,False)
    center_binx = int(Hist_Exposure_Skymap.GetNbinsX()/2.)
    center_biny = int(Hist_Exposure_Skymap.GetNbinsY()/2.)
    event_rate_smooth = Hist_Exposure_Skymap_smooth.GetBinContent(center_binx,center_biny)
    Hist_Exposure_Skymap_smooth = Hist_OnBkgd_Skymap_smooth.Clone()
    Hist_Exposure_Skymap_Galactic_smooth = Hist_OnBkgd_Skymap_Galactic_smooth.Clone()
    Hist_Exposure_Skymap_smooth.Scale(1./event_rate_smooth)
    Hist_Exposure_Skymap_Galactic_smooth.Scale(1./event_rate_smooth)

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Hist_OnBkgd_Skymap_Syst_RBM,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_RaDec_MDM_%s%s'%(source_name,PercentCrab))
    VariableSkymapBins(Syst_MDM,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Hist_OnBkgd_Skymap_Syst_RBM,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_RaDec_OpimizeNbins15_%s%s'%(source_name,PercentCrab),15)
    VariableSkymapBins(0.,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Hist_OnBkgd_Skymap_Syst_RBM,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_RaDec_OpimizeNbins15_NoSyst_%s%s'%(source_name,PercentCrab),15)
    VariableSkymapBins(Syst_MDM,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Hist_OnBkgd_Skymap_Syst_RBM,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_RaDec_OpimizeNbins5_%s%s'%(source_name,PercentCrab),5)
    VariableSkymapBins(0.,Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Hist_OnBkgd_Skymap_Syst_RBM,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_RaDec_OpimizeNbins5_NoSyst_%s%s'%(source_name,PercentCrab),5)

    Hist_Significance_Skymap = GetSignificanceMap(Hist_OnData_Skymap_Sum,Hist_OnBkgd_Skymap_Sum,Hist_OnBkgd_Skymap_Syst_MDM,Syst_MDM)
    MakeSignificanceDistribution(Hist_Significance_Skymap,Hist_OnData_Skymap_Sum,'SigDist_MDM_%s%s'%(source_name,PercentCrab))

    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_smooth,Hist_OnBkgd_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_MDM_smooth,Hist_OnBkgd_Skymap_Syst_RBM_smooth,Hist_Exposure_Skymap_smooth,'RA','Dec','Skymap_Smooth_RaDec_MDM_%s%s'%(source_name,PercentCrab))
    Make2DSignificancePlot(Syst_MDM,Hist_OnData_Skymap_Galactic_smooth,Hist_OnBkgd_Skymap_Galactic_smooth,Hist_OnBkgd_Skymap_Galactic_Syst_MDM_smooth,Hist_OnBkgd_Skymap_Galactic_Syst_RBM_smooth,Hist_Exposure_Skymap_Galactic_smooth,'gal. l.','gal. b.','Skymap_Smooth_Galactic_MDM_%s%s'%(source_name,PercentCrab))

    Hist_Significance_Skymap_smooth = GetSignificanceMap(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth,Hist_OnBkgd_Skymap_Syst_MDM_smooth,Syst_MDM)
    ErecS_lower_cut = energy_bin[energy_bin_cut_low]
    ErecS_upper_cut = energy_bin[energy_bin_cut_up]
    if doSmooth:
        Hist_Data_Energy_Skymap_smooth = []
        Hist_Bkgd_Energy_Skymap_smooth = []
        Hist_Zscore_Energy_Skymap_smooth = []
        Hist_Syst_Energy_Skymap_smooth = []
        for ebin in range(0,len(energy_bin)-1):
            Hist_Zscore_Energy_Skymap_smooth += [ROOT.TH2D("Hist_Zscore_Energy_Skymap_smooth_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
            Hist_Syst_Energy_Skymap_smooth += [ROOT.TH2D("Hist_Syst_Energy_Skymap_smooth_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
            Hist_Data_Energy_Skymap_smooth += [ROOT.TH2D("Hist_Data_Energy_Skymap_smooth_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
            Hist_Bkgd_Energy_Skymap_smooth += [ROOT.TH2D("Hist_Bkgd_Energy_Skymap_smooth_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
        for ebin in range(0,len(energy_bin)-1):
            Hist_Data_Energy_Skymap_smooth[ebin] = Smooth2DMap(Hist_Data_Energy_Skymap[ebin],smooth_size,False)
            Hist_Bkgd_Energy_Skymap_smooth[ebin] = Smooth2DMap(Hist_Bkgd_Energy_Skymap[ebin],smooth_size,False)
            Hist_Syst_Energy_Skymap_smooth[ebin] = Smooth2DMap(Hist_Syst_Energy_Skymap[ebin],smooth_size,False)
            Hist_Zscore_Energy_Skymap_smooth[ebin] = GetSignificanceMap(Hist_Data_Energy_Skymap_smooth[ebin],Hist_Bkgd_Energy_Skymap_smooth[ebin],Hist_Syst_Energy_Skymap_smooth[ebin],Syst_MDM)
            Hist_Expo_Energy_Skymap[ebin] = Hist_Bkgd_Energy_Skymap_smooth[ebin].Clone()
            Hist_Expo_Energy_Skymap[ebin].Scale(1./event_rate_smooth)
        MakeSpectrumIndexSkymap(Hist_Expo_Energy_Skymap,Hist_Data_Energy_Skymap_smooth,Hist_Bkgd_Energy_Skymap_smooth,Hist_Syst_Energy_Skymap_smooth,Hist_Zscore_Energy_Skymap_smooth,'RA','Dec','%s%s'%(source_name,PercentCrab),30,1)
        MakeSpectrumIndexSkymap(Hist_Expo_Energy_Skymap,Hist_Data_Energy_Skymap_smooth,Hist_Bkgd_Energy_Skymap_smooth,Hist_Syst_Energy_Skymap_smooth,Hist_Zscore_Energy_Skymap_smooth,'RA','Dec','%s%s_zoomin'%(source_name,PercentCrab),60,2)
    else:
        MakeSpectrumIndexSkymap(Hist_Expo_Energy_Skymap,Hist_Data_Energy_Skymap,Hist_Bkgd_Energy_Skymap,Hist_Syst_Energy_Skymap,Hist_Zscore_Energy_Skymap_smooth,'RA','Dec','%s%s'%(source_name,PercentCrab),12,1)
        MakeSpectrumIndexSkymap(Hist_Expo_Energy_Skymap,Hist_Data_Energy_Skymap,Hist_Bkgd_Energy_Skymap,Hist_Syst_Energy_Skymap,Hist_Zscore_Energy_Skymap_smooth,'RA','Dec','%s%s_zoomin'%(source_name,PercentCrab),30,3)

    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    init_x = source_ra
    init_y = source_dec
    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_OnData_Skymap_smooth, Hist_OnBkgd_Skymap_smooth, Hist_Significance_Skymap_smooth,Hist_Exposure_Skymap_smooth,2,init_x,init_y)
    print 'Excess (2 sigma) center RA = %0.3f'%(excess_center_x)
    print 'Excess (2 sigma) center Dec = %0.3f'%(excess_center_y)
    print 'Excess (2 sigma) radius = %0.3f'%(excess_radius)
    #for ebin in range(0,len(energy_bin)-1):
    #    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    #    init_x = source_ra
    #    init_y = source_dec
    #    excess_center_x, excess_center_y, excess_radius = GetExtention(Hist_Data_Energy_Skymap_smooth[ebin],Hist_Bkgd_Energy_Skymap_smooth[ebin],Hist_Significance_Skymap_smooth,Hist_Expo_Energy_Skymap[ebin],2,init_x,init_y)
    #    print 'Energy %0.1f GeV'%(energy_bin[ebin])
    #    print 'Excess (2 sigma) center RA = %0.3f'%(excess_center_x)
    #    print 'Excess (2 sigma) center Dec = %0.3f'%(excess_center_y)
    #    print 'Excess (2 sigma) radius = %0.3f'%(excess_radius)

def FindSourceIndex(source_name):
    for source in range(0,len(sample_list)):
        if source_name==sample_list[source]:
            return source
    return 0

Hist_NSB = ROOT.TH1D("Hist_NSB","",20,0,8)
Hist_L3Rate = ROOT.TH1D("Hist_L3Rate","",30,0,600)
Hist_L3Rate_all = ROOT.TH1D("Hist_L3Rate_all","",30,0,600)
Hist_EffArea = ROOT.TH1D("Hist_EffArea","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_EffArea_Sum = ROOT.TH1D("Hist_EffArea_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))

source_ra = 0.
source_dec = 0.
source_l = 0.
source_b = 0.
for source in range(0,len(sample_list)):
    source_idx = FindSourceIndex(sample_list[source])
    FilePath_Folder = []
    for elev in range(0,len(root_file_tags)):
        SourceFilePath = "%s/Netflix_"%(folder_path)+sample_list[source_idx]+"_%s"%(root_file_tags[elev])+".root"
        FilePath_Folder += [SourceFilePath]
        print 'Get %s...'%(FilePath_Folder[elev])
        if not os.path.isfile(FilePath_Folder[elev]): 
            print 'Found no file!!'
            continue
        else:
            print 'Found a file.'
            GetSourceInfo(FilePath_Folder)
    MakeOneHistPlot(Hist_NSB,'NSB','number of runs','NSB_%s'%(sample_list[source]),False)
    MakeOneHistPlot(Hist_L3Rate,'L3 rate','number of runs','L3Rate_%s'%(sample_list[source]),False)
    MakeOneHistPlot(Hist_L3Rate_all,'L3 rate','number of runs','L3Rate_all_%s'%(sample_list[source]),False)

print 'analysis cut: MSCL = %s, MSCW = %s'%(MSCL_blind_cut,MSCW_blind_cut)
MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_blind_cut-MSCW_plot_lower)+MSCW_blind_cut
MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_blind_cut-MSCL_plot_lower)+MSCL_blind_cut
print 'plot range: MSCL = %s, MSCW = %s'%(MSCL_plot_upper,MSCW_plot_upper)

Hist2D_OnData_Sum = ROOT.TH2D("Hist2D_OnData_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd_Sum = ROOT.TH2D("Hist2D_OnBkgd_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnGamma_Sum = ROOT.TH2D("Hist2D_OnGamma_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark_Sum = ROOT.TH2D("Hist2D_OnDark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnData = ROOT.TH2D("Hist2D_OnData","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnBkgd = ROOT.TH2D("Hist2D_OnBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnGamma = ROOT.TH2D("Hist2D_OnGamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_OnDark = ROOT.TH2D("Hist2D_OnDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCL_Sum = ROOT.TH1D("Hist_OnData_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL_Sum = ROOT.TH1D("Hist_OnDark_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCL = ROOT.TH1D("Hist_OnData_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_MSCL = ROOT.TH1D("Hist_OnBkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCL = ROOT.TH1D("Hist_OnDark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnData_MSCW_Sum = ROOT.TH1D("Hist_OnData_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCW_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnGamma_MSCW_Sum = ROOT.TH1D("Hist_OnGamma_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCL_Sum = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnGamma_MSCL_Sum = ROOT.TH1D("Hist_OnGamma_MSCL_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCW_Sum = ROOT.TH1D("Hist_OnDark_MSCW_Sum","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnData_MSCW = ROOT.TH1D("Hist_OnData_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_MSCW = ROOT.TH1D("Hist_OnBkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCW = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCW = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnGamma_MSCW = ROOT.TH1D("Hist_OnGamma_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_OnBkgd_Unblind_wGamma_MSCL = ROOT.TH1D("Hist_OnBkgd_Unblind_wGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnBkgd_Unblind_woGamma_MSCL = ROOT.TH1D("Hist_OnBkgd_Unblind_woGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnGamma_MSCL = ROOT.TH1D("Hist_OnGamma_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_OnDark_MSCW = ROOT.TH1D("Hist_OnDark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

Hist2D_U_Proj = ROOT.TH2D("Hist2D_U_Proj","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_V_Proj = ROOT.TH2D("Hist2D_V_Proj","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_U_Proj_Sum = ROOT.TH2D("Hist2D_U_Proj_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_V_Proj_Sum = ROOT.TH2D("Hist2D_V_Proj_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Data = ROOT.TH2D("Hist2D_Coeff_Data","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Bkgd = ROOT.TH2D("Hist2D_Coeff_Bkgd","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Data_Sum = ROOT.TH2D("Hist2D_Coeff_Data_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)
Hist2D_Coeff_Bkgd_Sum = ROOT.TH2D("Hist2D_Coeff_Bkgd_Sum","",N_bins_for_deconv,0,N_bins_for_deconv,N_bins_for_deconv,0,N_bins_for_deconv)

Hist2D_Rank0_Data = ROOT.TH2D("Hist2D_Rank0_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Data = ROOT.TH2D("Hist2D_Rank1_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Data = ROOT.TH2D("Hist2D_Rank2_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Data = ROOT.TH2D("Hist2D_Rank3_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Data = ROOT.TH2D("Hist2D_Rank4_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Data_Sum = ROOT.TH2D("Hist2D_Rank0_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Data_Sum = ROOT.TH2D("Hist2D_Rank1_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Data_Sum = ROOT.TH2D("Hist2D_Rank2_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Data_Sum = ROOT.TH2D("Hist2D_Rank3_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Data_Sum = ROOT.TH2D("Hist2D_Rank4_Data_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Dark = ROOT.TH2D("Hist2D_Rank0_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Dark = ROOT.TH2D("Hist2D_Rank1_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Dark = ROOT.TH2D("Hist2D_Rank2_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Dark = ROOT.TH2D("Hist2D_Rank3_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Dark = ROOT.TH2D("Hist2D_Rank4_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank0_Dark_Sum = ROOT.TH2D("Hist2D_Rank0_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank1_Dark_Sum = ROOT.TH2D("Hist2D_Rank1_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank2_Dark_Sum = ROOT.TH2D("Hist2D_Rank2_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank3_Dark_Sum = ROOT.TH2D("Hist2D_Rank3_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Rank4_Dark_Sum = ROOT.TH2D("Hist2D_Rank4_Dark_Sum","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

n_iterations = 100
optimiz_lower = -6.
optimiz_upper = -1.
Hist_Bkgd_Chi2 = ROOT.TH1D("Hist_Bkgd_Chi2","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Chi2_Diff = ROOT.TH1D("Hist_Bkgd_Chi2_Diff","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Chi2_Diff2 = ROOT.TH1D("Hist_Bkgd_Chi2_Diff2","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Optimization = ROOT.TH1D("Hist_Bkgd_Optimization","",50,optimiz_lower,optimiz_upper)
Hist_Bkgd_Converge_Blind = ROOT.TH1D("Hist_Bkgd_Converge_Blind","",n_iterations,0,n_iterations)
Hist_Bkgd_Converge_Unblind = ROOT.TH1D("Hist_Bkgd_Converge_Unblind","",n_iterations,0,n_iterations)
Hist_VVV_Eigenvalues = ROOT.TH1D("Hist_VVV_Eigenvalues","",N_bins_for_deconv*N_bins_for_deconv,0,N_bins_for_deconv*N_bins_for_deconv)

Hist1D_Data_Rank0_LeftVector = ROOT.TH1D("Hist1D_Data_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank1_LeftVector = ROOT.TH1D("Hist1D_Data_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank2_LeftVector = ROOT.TH1D("Hist1D_Data_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank3_LeftVector = ROOT.TH1D("Hist1D_Data_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank0_RightVector = ROOT.TH1D("Hist1D_Data_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank1_RightVector = ROOT.TH1D("Hist1D_Data_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank2_RightVector = ROOT.TH1D("Hist1D_Data_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Data_Rank3_RightVector = ROOT.TH1D("Hist1D_Data_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank0_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank1_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank2_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank3_LeftVector = ROOT.TH1D("Hist1D_Bkgd_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank0_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank1_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank2_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Bkgd_Rank3_RightVector = ROOT.TH1D("Hist1D_Bkgd_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank0_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank0_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank1_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank1_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank2_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank2_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank3_LeftVector = ROOT.TH1D("Hist1D_Dark_Rank3_LeftVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank0_RightVector = ROOT.TH1D("Hist1D_Dark_Rank0_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank1_RightVector = ROOT.TH1D("Hist1D_Dark_Rank1_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank2_RightVector = ROOT.TH1D("Hist1D_Dark_Rank2_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist1D_Dark_Rank3_RightVector = ROOT.TH1D("Hist1D_Dark_Rank3_RightVector","",N_bins_for_deconv,0,N_bins_for_deconv)

Hist_OnData_Yoff_Sum = ROOT.TH1D("Hist_OnData_Yoff_Sum","",30,-3,3)
Hist_OnBkgd_Yoff_Sum = ROOT.TH1D("Hist_OnBkgd_Yoff_Sum","",30,-3,3)
Hist_OnBkgd_Yoff_Raw_Sum = ROOT.TH1D("Hist_OnBkgd_Yoff_Raw_Sum","",30,-3,3)
Hist_OnData_Yoff = ROOT.TH1D("Hist_OnData_Yoff","",30,-3,3)
Hist_OnBkgd_Yoff = ROOT.TH1D("Hist_OnBkgd_Yoff","",30,-3,3)
Hist_OnBkgd_Yoff_Raw = ROOT.TH1D("Hist_OnBkgd_Yoff_Raw","",30,-3,3)
Hist_OnData_Theta2_Sum = ROOT.TH1D("Hist_OnData_Theta2_Sum","",50,0,10)
Hist_OnBkgd_Theta2_Sum = ROOT.TH1D("Hist_OnBkgd_Theta2_Sum","",50,0,10)
Hist_OnRFoV_Theta2_Sum = ROOT.TH1D("Hist_OnRFoV_Theta2_Sum","",50,0,10)
Hist_OnDark_Theta2_Sum = ROOT.TH1D("Hist_OnDark_Theta2_Sum","",50,0,10)
Hist_OnData_Theta2 = ROOT.TH1D("Hist_OnData_Theta2","",50,0,10)
Hist_OnBkgd_Theta2 = ROOT.TH1D("Hist_OnBkgd_Theta2","",50,0,10)
Hist_OnRFoV_Theta2 = ROOT.TH1D("Hist_OnRFoV_Theta2","",50,0,10)
Hist_OnDark_Theta2 = ROOT.TH1D("Hist_OnDark_Theta2","",50,0,10)
Hist_OnData_Energy_Sum = ROOT.TH1D("Hist_OnData_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnRFoV_Energy_Sum = ROOT.TH1D("Hist_OnRFoV_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy_CamCenter_Sum = ROOT.TH1D("Hist_OnData_Energy_CamCenter_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_CamCenter_Sum = ROOT.TH1D("Hist_OnBkgd_Energy_CamCenter_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnDark_Energy_Sum = ROOT.TH1D("Hist_OnDark_Energy_Sum","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy = ROOT.TH1D("Hist_OnData_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy = ROOT.TH1D("Hist_OnBkgd_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnRFoV_Energy = ROOT.TH1D("Hist_OnRFoV_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Energy_CamCenter = ROOT.TH1D("Hist_OnData_Energy_CamCenter","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnBkgd_Energy_CamCenter = ROOT.TH1D("Hist_OnBkgd_Energy_CamCenter","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnDark_Energy = ROOT.TH1D("Hist_OnDark_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_OnData_Zenith_Sum = ROOT.TH1D("Hist_OnData_Zenith_Sum","",45,0,90)
Hist_OnBkgd_Zenith_Sum = ROOT.TH1D("Hist_OnBkgd_Zenith_Sum","",45,0,90)
Hist_OnData_Zenith = ROOT.TH1D("Hist_OnData_Zenith","",45,0,90)
Hist_OnBkgd_Zenith = ROOT.TH1D("Hist_OnBkgd_Zenith","",45,0,90)

Hist_Syst_Energy_Skymap = []
Hist_Data_Energy_Skymap = []
Hist_Bkgd_Energy_Skymap = []
Hist_Expo_Energy_Skymap = []
for ebin in range(0,len(energy_bin)-1):
    Hist_Syst_Energy_Skymap += [ROOT.TH2D("Hist_Syst_Energy_Skymap_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
    Hist_Data_Energy_Skymap += [ROOT.TH2D("Hist_Data_Energy_Skymap_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
    Hist_Bkgd_Energy_Skymap += [ROOT.TH2D("Hist_Bkgd_Energy_Skymap_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]
    Hist_Expo_Energy_Skymap += [ROOT.TH2D("Hist_Expo_Energy_Skymap_%s"%(ebin),"",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)]

Hist_OnData_Skymap = ROOT.TH2D("Hist_OnData_Skymap","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap = ROOT.TH2D("Hist_OnBkgd_Skymap","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_Skymap_Sum = ROOT.TH2D("Hist_OnData_Skymap_Sum","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Sum","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst_MDM = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_MDM","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst_RBM = ROOT.TH2D("Hist_OnBkgd_Skymap_Syst_RBM","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size,Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnData_Skymap_Galactic = ROOT.TH2D("Hist_OnData_Skymap_Galactic","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnData_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnData_Skymap_Galactic_Sum","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Sum = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Sum","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Syst_MDM = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst_MDM","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)
Hist_OnBkgd_Skymap_Galactic_Syst_RBM = ROOT.TH2D("Hist_OnBkgd_Skymap_Galactic_Syst_RBM","",Skymap_nbins,source_l-Skymap_size,source_l+Skymap_size,Skymap_nbins,source_b-Skymap_size,source_b+Skymap_size)

Hist_OnData_Skymap_ProjX_Sum = ROOT.TH1D("Hist_OnData_Skymap_ProjX_Sum","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size)
Hist_OnBkgd_Skymap_ProjX_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_ProjX_Sum","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size)
Hist_OnBkgd_Skymap_Syst_ProjX = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjX","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size)
Hist_OnBkgd_Skymap_Syst_ProjX_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjX_Sum","",Skymap_nbins,source_ra-Skymap_size,source_ra+Skymap_size)
Hist_OnData_Skymap_ProjY_Sum = ROOT.TH1D("Hist_OnData_Skymap_ProjY_Sum","",Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_ProjY_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_ProjY_Sum","",Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst_ProjY = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjY","",Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)
Hist_OnBkgd_Skymap_Syst_ProjY_Sum = ROOT.TH1D("Hist_OnBkgd_Skymap_Syst_ProjY_Sum","",Skymap_nbins,source_dec-Skymap_size,source_dec+Skymap_size)

Hist_SystErr_MSCL = ROOT.TH1D("Hist_SystErr_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_SystErr_MSCW = ROOT.TH1D("Hist_SystErr_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_SystErr_Energy = ROOT.TH1D("Hist_SystErr_Energy","",len(energy_fine_bin)-1,array('d',energy_fine_bin))
Hist_SystErr_Theta2 = ROOT.TH1D("Hist_SystErr_Theta2","",50,0,10)

Hist_CR_Counts_MSCW = ROOT.TH1D("Hist_CR_Counts_MSCW","",6,0,6)
Hist_CR_Counts_MSCL = ROOT.TH1D("Hist_CR_Counts_MSCL","",6,0,6)

print 'MJD_Start = %s'%(MJD_Start)
time = Time(MJD_Start, format='mjd')
time.format = 'decimalyear'
year_start = time.value
print 'MJD_End = %s'%(MJD_End)
time = Time(MJD_End, format='mjd')
time.format = 'decimalyear'
year_end = time.value
print 'year_start = %s'%(year_start)
print 'year_end = %s'%(year_end)
print 'roi_ra = %s'%(roi_ra[0])
print 'roi_dec = %s'%(roi_dec[0])
Hist_OnData_RoI_Energy_Sum = [] 
Hist_OnBkgd_RoI_Energy_Sum = [] 
Hist_OnData_RoI_Energy = [] 
Hist_OnBkgd_RoI_Energy = [] 
Hist_OnData_RoI_Theta2_Sum = []
Hist_OnBkgd_RoI_Theta2_Sum = []
Hist_OnData_RoI_Theta2 = []
Hist_OnBkgd_RoI_Theta2 = []
Hist_OnData_RoI_MJD = [] 
Hist_OnBkgd_RoI_MJD = [] 
Hist_OnData_RoI_MJD_Sum = [] 
Hist_OnBkgd_RoI_MJD_Sum = [] 
for nth_roi in range(0,len(roi_ra)):
    Hist_OnData_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnData_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_Sum_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Energy += [ROOT.TH1D("Hist_OnData_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnBkgd_RoI_Energy += [ROOT.TH1D("Hist_OnBkgd_RoI_Energy_%s"%(nth_roi),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OnData_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnData_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_Sum_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_Theta2 += [ROOT.TH1D("Hist_OnData_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnBkgd_RoI_Theta2 += [ROOT.TH1D("Hist_OnBkgd_RoI_Theta2_%s"%(nth_roi),"",20,0,0.5)]
    Hist_OnData_RoI_MJD += [ROOT.TH1D("Hist_OnData_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnData_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnData_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]
    Hist_OnBkgd_RoI_MJD_Sum += [ROOT.TH1D("Hist_OnBkgd_RoI_MJD_Sum_%s"%(nth_roi),"",800,56200-4000,56200+4000)]


Hist2D_OffData = []
Hist2D_OffData_Sum = []
Hist_OffData_MSCL = []
Hist_OffData_MSCL_Sum = []
Hist_OffBkgd_MSCL = []
Hist_OffBkgd_MSCL_Sum = []
Hist_OffData_MSCW = []
Hist_OffData_MSCW_Sum = []
Hist_OffBkgd_MSCW = []
Hist_OffBkgd_MSCW_Sum = []
Hist_OffData_Energy = []
Hist_OffData_Energy_Sum = []
Hist_OffBkgd_Energy = []
Hist_OffBkgd_Energy_Sum = []
Hist_OffData_CameraFoV_Theta2 = []
Hist_OffData_CameraFoV_Theta2_Sum = []
Hist_OffBkgd_CameraFoV_Theta2 = []
Hist_OffBkgd_CameraFoV_Theta2_Sum = []
for nth_sample in range(0,n_control_samples):
    Hist2D_OffData += [ROOT.TH2D("Hist2D_OffData_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist2D_OffData_Sum += [ROOT.TH2D("Hist2D_OffData_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCL += [ROOT.TH1D("Hist_OnData_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCL_Sum += [ROOT.TH1D("Hist_OnData_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffData_MSCW += [ROOT.TH1D("Hist_OnData_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_MSCW_Sum += [ROOT.TH1D("Hist_OnData_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCL += [ROOT.TH1D("Hist_OnBkgd_MSCL_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCL_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCL_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)]
    Hist_OffBkgd_MSCW += [ROOT.TH1D("Hist_OnBkgd_MSCW_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffBkgd_MSCW_Sum += [ROOT.TH1D("Hist_OnBkgd_MSCW_Sum_%s"%(nth_sample),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
    Hist_OffData_Energy += [ROOT.TH1D("Hist_OffData_Energy_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffData_Energy_Sum += [ROOT.TH1D("Hist_OffData_Energy_Sum_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffBkgd_Energy += [ROOT.TH1D("Hist_OffBkgd_Energy_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffBkgd_Energy_Sum += [ROOT.TH1D("Hist_OffBkgd_Energy_Sum_%s"%(nth_sample),"",len(energy_fine_bin)-1,array('d',energy_fine_bin))]
    Hist_OffData_CameraFoV_Theta2 += [ROOT.TH1D("Hist_OffData_CameraFoV_Theta2_%s"%(nth_sample),"",50,0,10)]
    Hist_OffData_CameraFoV_Theta2_Sum += [ROOT.TH1D("Hist_OffData_CameraFoV_Theta2_Sum_%s"%(nth_sample),"",50,0,10)]
    Hist_OffBkgd_CameraFoV_Theta2 += [ROOT.TH1D("Hist_OffBkgd_CameraFoV_Theta2_%s"%(nth_sample),"",50,0,10)]
    Hist_OffBkgd_CameraFoV_Theta2_Sum += [ROOT.TH1D("Hist_OffBkgd_CameraFoV_Theta2_Sum_%s"%(nth_sample),"",50,0,10)]

GetGammaSourceInfo()

#SystematicAnalysis()

#drawMap = False
drawMap = True
#Smoothing = False
Smoothing = True

#set_palette('default')
#set_palette('gray')

#SingleSourceAnalysis(sample_list,drawMap,Smoothing,0,6)
SingleSourceAnalysis(sample_list,drawMap,Smoothing,1,6)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,2,6)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,3,6)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,0,1)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,1,2)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,2,3)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,3,4)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,4,6)
#SingleSourceAnalysis(sample_list,drawMap,Smoothing,5,6)

print 'n_good_matches = %s'%(n_good_matches)

print "Syst_MDM = %s"%(Syst_MDM) 
print "Syst_Init = %s"%(Syst_Init) 
print "Syst_Redu = %s"%(Syst_Redu) 
print "Syst_Clos = %s"%(Syst_Clos) 
