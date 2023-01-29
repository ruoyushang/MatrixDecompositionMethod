


from prettytable import PrettyTable
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as my_unit
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle

from bs4 import BeautifulSoup

def powerlaw_func(x,A,r):
    return A+r*x
def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split(':')]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split(':')]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)           
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

def FindSourceVisibility(psr_name,psr_ra,psr_dec,obs='vts'):

    current_time = '2022-04-01 23:00:00'
    #psr_coord = SkyCoord.from_name(psr_name)
    psr_coord = SkyCoord(psr_ra, psr_dec, unit="deg")

    obs_site = EarthLocation(lat=31.6751*u.deg, lon=-110.952*u.deg, height=1268*u.m) # veritas
    if obs=='hwc':
        obs_site = EarthLocation(lat=18.994722*u.deg, lon=-97.3085*u.deg, height=4100*u.m) # HAWC
    utcoffset = -7*u.hour  # AZ time
    time = Time(current_time) - utcoffset
    midnight = Time(current_time) - utcoffset

    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_whole_night = midnight + delta_midnight
    frame_whole_night = AltAz(obstime=times_whole_night, location=obs_site)

    psr_altazs_whole_night = psr_coord.transform_to(frame_whole_night)
    max_alt = np.max(psr_altazs_whole_night.alt.deg)
    #print ('max alt = %0.1f deg'%(max_alt))
    return max_alt

def FindFermiSource(psr_name,psr_ra,psr_dec,fermi_name,fermi_ra,fermi_dec,search_radis):

    found_fermi_name = ''
    found_fermi_idx = 0
    for fermi in range(0,len(fermi_name)):
        delta_ra = psr_ra - fermi_ra[fermi]
        delta_dec = psr_dec - fermi_dec[fermi]
        distance = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5)
        #if distance<0.5:
        if distance<search_radis:
            found_fermi_name = fermi_name[fermi]
            found_fermi_idx = fermi
    return found_fermi_name, found_fermi_idx

def ReadFermiCatelog():
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open('gll_psc_v26.xml')
    target_name = ''
    target_type = ''
    target_info = ''
    target_flux = ''
    target_ra = ''
    target_dec = ''
    for line in inputFile:
        if line.split(' ')[0]=='<source':
            for block in range(0,len(line.split(' '))):
                if 'Unc_' in line.split(' ')[block]: continue
                if 'name=' in line.split(' ')[block]:
                    target_name = line.split('name="')[1].split('"')[0]
                if 'type=' in line.split(' ')[block]:
                    target_type = line.split('type="')[1].split('"')[0]
                if 'Flux1000=' in line.split(' ')[block]:
                    target_flux = line.split('Flux1000="')[1].split('"')[0]
                if 'Energy_Flux100=' in line.split(' ')[block]:
                    target_info = line.split(' ')[block]
                    target_info = target_info.strip('>\n')
                    target_info = target_info.strip('"')
                    target_info = target_info.lstrip('Energy_Flux100="')
        if '<parameter' in line and 'name="RA"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_ra = line.split(' ')[block].split('"')[1]
        if '<parameter' in line and 'name="DEC"' in line:
            for block in range(0,len(line.split(' '))):
                if 'value=' in line.split(' ')[block]:
                    target_dec = line.split(' ')[block].split('"')[1]
        if 'source>' in line:
            if target_ra=='': 
                target_name = ''
                target_type = ''
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #if target_type=='PointSource':
            #    target_name = ''
            #    target_type = ''
            #    target_info = ''
            #    target_ra = ''
            #    target_dec = ''
            #    continue
            #print ('target_name = %s'%(target_name))
            #print ('target_ra = %s'%(target_ra))
            #print ('target_dec = %s'%(target_dec))
            #source_name += [target_name]
            #source_name += ['%0.2e'%(float(target_info))]
            source_name += ['%0.2e'%(float(target_flux))]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_type = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

def FindVeritasExposure(psr_ra,psr_dec):

    Total_Livetime = 0.

    for entry in range(0,len(List_RunNumber)):
    
        RunNumber = List_RunNumber[entry]
        Name = List_Name[entry]
        offset = List_Offset[entry]
        MJD = List_MJD[entry]
        Elev = List_Elev[entry]
        Azim = List_Azim[entry]
        T1_RA = List_T1_RA[entry]
        T1_Dec = List_T1_Dec[entry]
        PedVar_DC = List_PedVar_DC[entry]
        PedVar_PE = List_PedVar_PE[entry]
        FIR_Mean = List_FIR_Mean[entry]
        FIR_RMS = List_FIR_RMS[entry]
        L3_rate = List_L3_rate[entry]
        Livetime = List_Livetime[entry]
    
        epoch = 'V4'
        if int(RunNumber)<46642:
            continue
        if int(RunNumber)>=46642 and int(RunNumber)<63373:
            epoch = 'V5'
        if int(RunNumber)>=63373:
            epoch = 'V6'
    
        if Elev<45.: continue
        #if Elev<Search_Range_Elev[0]: continue
        #if Elev>Search_Range_Elev[1]: continue
        #if Azim<Search_Range_Azim[0]: continue
        #if Azim>Search_Range_Azim[1]: continue
        #if PedVar_DC<Search_Range_PedVar_DC[0]: continue
        #if PedVar_DC>Search_Range_PedVar_DC[1]: continue

        #if MJD<59480: continue

        if epoch=='V6':
            if L3_rate<250.: continue
        else:
            if L3_rate<150.: continue
        if Livetime<5.: continue
    
        if abs(float(T1_RA)-psr_ra)>1.0: continue
        if abs(float(T1_Dec)-psr_dec)>1.0: continue

        Total_Livetime += Livetime/60.

    return Total_Livetime

def ConvertGalacticToRaDec(l, b):
    my_sky = SkyCoord(l*my_unit.deg, b*my_unit.deg, frame='galactic')
    return my_sky.icrs.ra.deg, my_sky.icrs.dec.deg

def ConvertRaDecToGalactic(ra, dec):
    my_sky = SkyCoord(ra*my_unit.deg, dec*my_unit.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg

def ReadSNRTargetListFromCSVFile():
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_size = []
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '0'
            target_size = row[15]
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            source_dist += [float(target_min_dist)]
            source_size += [float(target_size)]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            #print(row)
    return source_name, source_ra, source_dec, source_dist, source_size

def ReadHessTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        if '%' in line: continue
        target_name = line.split()[0]
        target_ra = line.split()[1]
        target_dec = line.split()[2]
        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        #print('source_name = %s'%(source_name[len(source_name)-1]))
        #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
        #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
    return source_name, source_ra, source_dec

def ReadVeritasTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        target_name = line.split(';')[0]
        target_ra = line.split(';')[2]
        target_dec = line.split(';')[3]
        source_name += [target_name]
        source_ra += [float(target_ra)]
        source_dec += [float(target_dec)]
    return source_name, source_ra, source_dec

def ReadHAWCTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_flux = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        if '- name:' in line:
            target_name = line.lstrip('   - name: ')
            target_name = target_name.strip('\n')
        if 'RA:' in line:
            target_ra = line.lstrip('     RA: ')
        if 'Dec:' in line:
            target_dec = line.lstrip('     Dec: ')
        if 'flux:' in line:
            target_flux = line.lstrip('          flux: ')
        if 'index systematic uncertainty down:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            source_flux += [float(target_flux)/2.34204e-13]
            target_name = ''
            target_ra = ''
            target_dec = ''
            target_flux = ''
    return source_name, source_ra, source_dec, source_flux

def ReadTeVCatTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        target_name = line.split(',')[0].strip(" ")
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        source_name += [target_name]
        source_ra += [float(target_ra)]
        source_dec += [float(target_dec)]
    return source_name, source_ra, source_dec

def FindGlobularClusterSumEdot(cluster_name):

    cluster_name = cluster_name.strip()
    found_cluster = False
    cluster_vars = []
    pulsar_vars = []
    inputFile = open('Pulsars_in_Globular_Clusters.html')
    Sum_Edot = 0.
    for line in inputFile:
        if '<a name=' in line:
            name_in_line = line.split('"')[1].strip()
            if cluster_name == name_in_line:
                found_cluster = True
            else:
                found_cluster = False
        if found_cluster:
            if '<td><b>' in line:
                pulsar_vars = []
                pulsar_vars += [line.lstrip('<td>').rstrip('</td>')]
            elif '<td>' in line:
                pulsar_vars += [line.lstrip('<td>').rstrip('</td>')]
            if '</tr>' in line:
                cluster_vars += [pulsar_vars]
    for psr in range(1,len(cluster_vars)):
        psr_name = cluster_vars[psr][0]
        psr_name = BeautifulSoup(psr_name).get_text().lstrip('psr_name=').strip()
        P0 = cluster_vars[psr][2]
        P0 = BeautifulSoup(P0).get_text().strip().split('(')[0]
        if P0=='*': continue
        if P0=='': continue
        if '×' in P0: continue
        P0 = float(P0)*1e-3
        Pdot = cluster_vars[psr][3]
        Pdot = BeautifulSoup(Pdot).get_text().strip().split('(')[0]
        if Pdot=='*': continue
        if Pdot=='': continue
        if '×' in Pdot: continue
        Pdot = float(Pdot)*1e-20
        Inertia = 1e45
        Edot = 4.*pow(math.pi,2)*Inertia*Pdot/pow(P0,3)
        Sum_Edot += Edot

    #print ('Sum_Edot = %0.2e erg/s'%(Sum_Edot))
    return Sum_Edot

def ReadGlobularClusterListFromFile():
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_npsr = []
    source_flux = []
    inputFile = open('Baumgardt_globular_Cat.txt')
    for line in inputFile:
        if line[0]=="%": continue
        target_name = line.split('|')[0].rstrip(' ')
        target_ra = float(line.split('|')[1].strip(' '))
        target_dec = float(line.split('|')[2].strip(' '))
        target_dist = float(line.split('|')[3].strip(' '))
        target_log10_rho_c = float(line.split('|')[13].strip(' '))
        target_rho_c = pow(10,target_log10_rho_c)
        target_MLratio = float(line.split('|')[8].strip(' '))
        target_r_c = float(line.split('|')[9].strip(' '))
        target_sig_0 = float(line.split('|')[20].strip(' '))
        target_gamma_c = pow(target_rho_c/target_MLratio,2)*pow(target_r_c,3)/target_sig_0
        target_ln_npsr = -1.1+1.5*math.log10(target_gamma_c)
        alpha = 0.01 # arbitrary factor of 0.01 to match the HESS Ter 5 flux
        target_npsr = math.exp(target_ln_npsr)
        target_flux = alpha*1.6*1e-10*target_npsr/1000./(target_dist*target_dist) #/cm2/s at 1 TeV
        Crab_flux = 3.8*1e-11 #/cm2/s at 1 TeV
        source_name += [target_name]
        source_ra += [target_ra]
        source_dec += [target_dec]
        source_dist += [target_dist]
        source_npsr += [target_npsr]
        source_flux += [target_flux/Crab_flux]

    return source_name, source_ra, source_dec, source_dist, source_npsr, source_flux

def ReadGlobularClusterListFromFile_old(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_npsr = []
    source_nmsp = []
    source_edot = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="%": continue
        target_name = line.split('|')[1].lstrip(' ').strip("\t")
        #print (target_name)
        target_l = float(line.split('|')[4].strip(" ").strip("\t"))
        target_b = float(line.split('|')[5].strip(" ").strip("\t"))
        target_dist = float(line.split('|')[6].strip(" ").strip("\t"))
        target_npsr = float(line.split('|')[7].strip(">").strip(" ").strip("\t"))
        target_nmsp = float(line.split('|')[8].strip(">").strip(" ").strip("\t"))
        target_ra, target_dec = ConvertGalacticToRaDec(target_l,target_b)
        target_edot = FindGlobularClusterSumEdot(target_name)
        source_name += [target_name]
        source_ra += [target_ra]
        source_dec += [target_dec]
        source_dist += [target_dist]
        source_npsr += [target_npsr]
        source_nmsp += [target_nmsp]
        source_edot += [target_edot]
    return source_name, source_ra, source_dec, source_dist, source_npsr, source_nmsp, source_edot

def ReadATNFTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    source_flux = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        if target_name=="\n": continue
        #print ('target_name = %s'%(target_name))
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        #print ('target_ra = %s'%(target_ra))
        #print ('target_dec = %s'%(target_dec))
        target_dist = line.split(',')[3].strip(" ")
        target_age = line.split(',')[4].strip(" ")
        target_edot = line.split(',')[5].strip(" ")
        if target_dist=='*': continue
        if target_age=='*': continue
        if target_edot=='*': continue
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        #if float(target_brightness)<1e33 and float(target_edot)<1e34: continue
        #if float(target_dist)>6.: continue
        #if float(target_age)<1e4: continue
        #if float(target_age)>1e5: continue

        #ra_deg = float(HMS2deg(target_ra,target_dec)[0])
        #dec_deg = float(HMS2deg(target_ra,target_dec)[1])
        #gal_l, gal_b = ConvertRaDecToGalactic(ra_deg,dec_deg)
        #if abs(gal_b)<5.: continue

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [float(target_dist)]
        source_age += [float(target_age)]
        source_edot += [float(target_edot)]
        source_flux += [float(target_edot)/pow(float(target_dist),2)]
    return source_name, source_ra, source_dec, source_dist, source_age, source_edot, source_flux

List_RunNumber = []
List_Name = []
List_Offset = []
List_MJD = []
List_Elev = []
List_Azim = []
List_PedVar_DC = []
List_PedVar_PE = []
List_T1_RA = []
List_T1_Dec = []
List_T2_RA = []
List_T2_Dec = []
List_T3_RA = []
List_T3_Dec = []
List_T4_RA = []
List_T4_Dec = []
List_FIR_Mean = []
List_FIR_RMS = []
List_L3_rate = []
List_Livetime = []

Source_RunNumber = []
Source_PedVar_DC = []
Source_Elev = []
Source_Azim = []
Source_Livetime = []

inputFile = open('diagnostics.txt')
for line in inputFile:
    if line.split(' ')[0]=="#": 
        #print 'this is a comment line'
        continue
    if len(line.split(' '))<40: continue
    V0 = False
    V1 = False
    V2 = False
    if line.split(' ')[1]=="0": 
        #print 'this is a V0 line'
        V0 = True
    if line.split(' ')[1]=="1": 
        #print 'this is a V1 line'
        V1 = True
    if line.split(' ')[1]=="2": 
        #print 'this is a V2 line'
        V2 = True
    RunNumber = line.split(' ')[1-1]
    Name = line.split(' ')[6-1]
    if RunNumber=='': continue
    List_RunNumber += [int(RunNumber)]
    List_Name += [Name]
    offset = line.split(' ')[7-1]
    List_Offset += [offset]
    MJD = line.split(' ')[5-1]
    Elev = line.split(' ')[8-1]
    Azim = line.split(' ')[9-1]
    List_MJD += [int(float(MJD))]
    List_Elev += [float(Elev)]
    List_Azim += [float(Azim)]
    L3_rate = line.split(' ')[12-1]
    List_L3_rate += [float(L3_rate)]
    Livetime = line.split(' ')[11-1]
    List_Livetime += [float(Livetime)]
    if V1 or V2:
        T1_RA = line.split(' ')[46-1]
        T1_Dec = line.split(' ')[47-1]
        T2_RA = line.split(' ')[46-1+14]
        T2_Dec = line.split(' ')[47-1+14]
        T3_RA = line.split(' ')[46-1+28]
        T3_Dec = line.split(' ')[47-1+28]
        T4_RA = line.split(' ')[46-1+42]
        T4_Dec = line.split(' ')[47-1+42]
        List_T1_RA += [float(T1_RA)]
        List_T1_Dec += [float(T1_Dec)]
        List_T2_RA += [float(T2_RA)]
        List_T2_Dec += [float(T2_Dec)]
        List_T3_RA += [float(T3_RA)]
        List_T3_Dec += [float(T3_Dec)]
        List_T4_RA += [float(T4_RA)]
        List_T4_Dec += [float(T4_Dec)]
    else:
        List_T1_RA += [0.]
        List_T1_Dec += [0.]
        List_T2_RA += [0.]
        List_T2_Dec += [0.]
        List_T3_RA += [0.]
        List_T3_Dec += [0.]
        List_T4_RA += [0.]
        List_T4_Dec += [0.]
    if V2:
        PedVar_DC = line.split(' ')[104-1]
        PedVar_PE = line.split(' ')[105-1]
        List_PedVar_DC += [float(PedVar_DC)]
        List_PedVar_PE += [float(PedVar_PE)]
        FIR_Mean = line.split(' ')[106-1]
        FIR_RMS = line.split(' ')[107-1]
        List_FIR_Mean += [float(FIR_Mean)]
        List_FIR_RMS += [float(FIR_RMS)]
    else:
        List_PedVar_DC += [0.]
        List_PedVar_PE += [0.]
        List_FIR_Mean += [0.]
        List_FIR_RMS += [0.]
    for entry in range(0,len(Source_RunNumber)):
        if Source_RunNumber[entry]==int(RunNumber):
            Source_PedVar_DC[entry] = float(PedVar_DC)
            Source_Elev[entry] = float(Elev)
            Source_Azim[entry] = float(Azim)
            Source_Livetime[entry] = float(Livetime)


fermi_name, fermi_ra, fermi_dec = ReadFermiCatelog()
target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age, target_psr_edot, target_psr_flux = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
target_glb_name, target_glb_ra, target_glb_dec, target_glb_dist, target_glb_npsr, target_glb_flux = ReadGlobularClusterListFromFile()
target_snr_name, target_snr_ra, target_snr_dec, target_snr_dist, target_snr_size = ReadSNRTargetListFromCSVFile()
target_tev_name, target_tev_ra, target_tev_dec = ReadTeVCatTargetListFromFile('TeVCat_RaDec_w_Names.txt')
target_hwc_name, target_hwc_ra, target_hwc_dec, target_hwc_flux = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
target_vts_name, target_vts_ra, target_vts_dec = ReadVeritasTargetListFromFile('VTS_Cat.txt')
target_hss_name, target_hss_ra, target_hss_dec = ReadHessTargetListFromFile('HESS_Cat.txt')

include_ra_bands = []
include_ra_bands += ['00']
include_ra_bands += ['02']
include_ra_bands += ['04']
include_ra_bands += ['06']
include_ra_bands += ['08']
include_ra_bands += ['10']
include_ra_bands += ['12']
include_ra_bands += ['14']
include_ra_bands += ['16']
include_ra_bands += ['18']
include_ra_bands += ['20']
include_ra_bands += ['22']

select_psr = True
select_msp = False
select_gcs = True
select_snr = True
#select_new_target = False
#alt_cut = 30.
select_new_target = True
alt_cut = 60.

R = "\033[0;31m" #RED
G = "\033[0;32m" # GREEN
Y = "\033[0;33m" # Yellow
B = "\033[0;34m" # Blue
N = "\033[0m" # Reset

nongamma_psr_age = []
nongamma_psr_dist = []
nongamma_psr_edot = []
nongamma_psr_flux = []

vtsobs_psr_age = []
vtsobs_psr_dist = []
vtsobs_psr_edot = []
vtsobs_psr_flux = []
nonvtsobs_psr_age = []
nonvtsobs_psr_dist = []
nonvtsobs_psr_edot = []
nonvtsobs_psr_flux = []

vts_psr_age = []
vts_psr_dist = []
vts_psr_edot = []
vts_psr_flux = []
nonvts_psr_age = []
nonvts_psr_dist = []
nonvts_psr_edot = []
nonvts_psr_flux = []

hss_psr_age = []
hss_psr_dist = []
hss_psr_edot = []
hss_psr_flux = []

hwc_psr_age = []
hwc_psr_dist = []
hwc_psr_edot = []
hwc_psr_flux = []
nonhwc_psr_age = []
nonhwc_psr_dist = []
nonhwc_psr_edot = []
nonhwc_psr_flux = []

fermi_psr_age = []
fermi_psr_dist = []
fermi_psr_edot = []
fermi_psr_flux = []
fermi_psr_gflu = []
nonfermi_psr_age = []
nonfermi_psr_dist = []
nonfermi_psr_edot = []
nonfermi_psr_flux = []

fermihwc_psr_age = []
fermihwc_psr_dist = []
fermihwc_psr_edot = []
fermihwc_psr_flux = []
fermihwc_psr_gflu = []
fermivts_psr_age = []
fermivts_psr_dist = []
fermivts_psr_edot = []
fermivts_psr_flux = []
fermivts_psr_gflu = []
fermihss_psr_age = []
fermihss_psr_dist = []
fermihss_psr_edot = []
fermihss_psr_flux = []
fermihss_psr_gflu = []

lines_to_write = []

for ra in range(0,len(include_ra_bands)):

    if not select_psr and not select_msp: continue

    my_psr_table = PrettyTable(['%s band'%(include_ra_bands[ra]), 'RA/Dec (deg)', 'l/b (deg)', 'TeV source', 'Age (kyr)', 'Dist (kpc)', 'Edot (erg/s)', 'Est. CU', 'HAWC CU', 'VTS expo (hr)', 'VTS alt (deg)'])
    for psr in range(0,len(target_psr_name)):
    
        myangle_ra = Angle(target_psr_ra[psr], my_unit.deg)
        include_this_psr = False
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+2.:
            include_this_psr = True
        if not include_this_psr: continue
    
        vts_alt = FindSourceVisibility('PSR %s'%(target_psr_name[psr]),target_psr_ra[psr],target_psr_dec[psr],obs='vts')
        hwc_alt = FindSourceVisibility('PSR %s'%(target_psr_name[psr]),target_psr_ra[psr],target_psr_dec[psr],obs='hwc')
    
        found_fermi_name, frm_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],fermi_name,fermi_ra,fermi_dec,0.2)
        found_tev_name  , tev_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_tev_name,target_tev_ra,target_tev_dec,0.2)
        found_hwc_name  , hwc_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_hwc_name,target_hwc_ra,target_hwc_dec,0.3)
        found_vts_name  , vts_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_vts_name,target_vts_ra,target_vts_dec,0.3)
        found_hss_name  , hss_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_hss_name,target_hss_ra,target_hss_dec,0.3)
        gal_l, gal_b = ConvertRaDecToGalactic(target_psr_ra[psr],target_psr_dec[psr])
    
        crab_edot = 4.5e+38 
        crab_dist = 2.000 # kpc
        crab_pwn_bfield = 125.0 # muG
        halo_bfield = 3. # muG
        if select_psr:
            halo_bfield = pow(10.,powerlaw_func(math.log10(target_psr_age[psr]/1000.), 2.21e+00, -8.51e-01)) # muG
            tev_flux_in_cu = target_psr_edot[psr]/crab_edot*pow(crab_pwn_bfield/halo_bfield,2)*pow(crab_dist/target_psr_dist[psr],2)
            target_psr_flux[psr] = tev_flux_in_cu
        if select_msp:
            alpha = 0.01 # arbitrary factor of 0.01 to match the HESS Ter 5 flux
            target_npsr = 1.
            rel_Edot = target_psr_edot[psr]/(3.2*1e34)
            target_flux = alpha*1.6*1e-10*rel_Edot*target_npsr/1000./(target_psr_dist[psr]*target_psr_dist[psr]) #/cm2/s at 1 TeV
            Crab_flux = 3.8*1e-11 #/cm2/s at 1 TeV
            tev_flux_in_cu = target_flux/Crab_flux
            target_psr_flux[psr] = tev_flux_in_cu

        if select_psr:
            if target_psr_age[psr]/1000.>pow(10,5): continue
            if tev_flux_in_cu<1e-3: continue
        else:
            if target_psr_age[psr]/1000.<pow(10,5): continue
            if target_psr_dist[psr]>2.0: continue
            if target_psr_edot[psr]<1e34: continue

        #if target_psr_age[psr]/1000.<10.: continue
    

        if hwc_alt>45.:
            if not found_hwc_name=='':
                hwc_psr_age += [target_psr_age[psr]/1000.]
                hwc_psr_dist += [target_psr_dist[psr]]
                hwc_psr_edot += [target_psr_edot[psr]]
                hwc_psr_flux += [target_psr_flux[psr]]
                if not found_fermi_name=='':
                    fermihwc_psr_age += [target_psr_age[psr]/1000.]
                    fermihwc_psr_dist += [target_psr_dist[psr]]
                    fermihwc_psr_edot += [target_psr_edot[psr]]
                    fermihwc_psr_flux += [target_psr_flux[psr]]
                    fermihwc_psr_gflu += [float(found_fermi_name)]
            else:
                nonhwc_psr_age += [target_psr_age[psr]/1000.]
                nonhwc_psr_dist += [target_psr_dist[psr]]
                nonhwc_psr_edot += [target_psr_edot[psr]]
                nonhwc_psr_flux += [target_psr_flux[psr]]
        if not found_fermi_name=='':
            fermi_psr_age += [target_psr_age[psr]/1000.]
            fermi_psr_dist += [target_psr_dist[psr]]
            fermi_psr_edot += [target_psr_edot[psr]]
            fermi_psr_flux += [target_psr_flux[psr]]
            fermi_psr_gflu += [float(found_fermi_name)]
        else:
            nonfermi_psr_age += [target_psr_age[psr]/1000.]
            nonfermi_psr_dist += [target_psr_dist[psr]]
            nonfermi_psr_edot += [target_psr_edot[psr]]
            nonfermi_psr_flux += [target_psr_flux[psr]]

        if found_fermi_name=='' and found_vts_name=='' and found_hss_name=='' and found_hwc_name=='':
            nongamma_psr_age += [target_psr_age[psr]/1000.]
            nongamma_psr_dist += [target_psr_dist[psr]]
            nongamma_psr_edot += [target_psr_edot[psr]]
            nongamma_psr_flux += [target_psr_flux[psr]]
    
        exposure_time = FindVeritasExposure(target_psr_ra[psr],target_psr_dec[psr])

        if vts_alt>45.:
            if exposure_time>10.:
                vtsobs_psr_age += [target_psr_age[psr]/1000.]
                vtsobs_psr_dist += [target_psr_dist[psr]]
                vtsobs_psr_edot += [target_psr_edot[psr]]
                vtsobs_psr_flux += [target_psr_flux[psr]]
            else:
                nonvtsobs_psr_age += [target_psr_age[psr]/1000.]
                nonvtsobs_psr_dist += [target_psr_dist[psr]]
                nonvtsobs_psr_edot += [target_psr_edot[psr]]
                nonvtsobs_psr_flux += [target_psr_flux[psr]]

        if vts_alt>45. and exposure_time>10.:
            if not found_vts_name=='':
                vts_psr_age += [target_psr_age[psr]/1000.]
                vts_psr_dist += [target_psr_dist[psr]]
                vts_psr_edot += [target_psr_edot[psr]]
                vts_psr_flux += [target_psr_flux[psr]]
                if not found_fermi_name=='':
                    fermivts_psr_age += [target_psr_age[psr]/1000.]
                    fermivts_psr_dist += [target_psr_dist[psr]]
                    fermivts_psr_edot += [target_psr_edot[psr]]
                    fermivts_psr_flux += [target_psr_flux[psr]]
                    fermivts_psr_gflu += [float(found_fermi_name)]
            else:
                nonvts_psr_age += [target_psr_age[psr]/1000.]
                nonvts_psr_dist += [target_psr_dist[psr]]
                nonvts_psr_edot += [target_psr_edot[psr]]
                nonvts_psr_flux += [target_psr_flux[psr]]

        if not found_hss_name=='':
            hss_psr_age += [target_psr_age[psr]/1000.]
            hss_psr_dist += [target_psr_dist[psr]]
            hss_psr_edot += [target_psr_edot[psr]]
            hss_psr_flux += [target_psr_flux[psr]]
            if not found_fermi_name=='':
                fermihss_psr_age += [target_psr_age[psr]/1000.]
                fermihss_psr_dist += [target_psr_dist[psr]]
                fermihss_psr_edot += [target_psr_edot[psr]]
                fermihss_psr_flux += [target_psr_flux[psr]]
                fermihss_psr_gflu += [float(found_fermi_name)]

        if select_new_target:
            if exposure_time>20.: continue
        else:
            if exposure_time<20.: continue

        if vts_alt<alt_cut: continue

        txt_ra = '%0.2f'%(target_psr_ra[psr])
        txt_dec = '%0.2f'%(target_psr_dec[psr])
        txt_l = '%0.2f'%(gal_l)
        txt_b = '%0.2f'%(gal_b)
        txt_dist = '%0.1f'%(target_psr_dist[psr])
        txt_age = '%0.1e'%(target_psr_age[psr]/1000.)
        txt_edot = '%0.1e'%(target_psr_edot[psr])
        txt_cu = '%0.1e'%(tev_flux_in_cu)
        txt_flux = '%0.1e'%(target_psr_flux[psr])
        txt_expo = '%0.1f'%(exposure_time)
        txt_alt = '%0.1f'%(vts_alt)
        txt_hwc_flux = ''
    
        found_counterpart = ''
        if not found_hwc_name=='':
            found_counterpart += '%s'%(found_hwc_name)
            txt_hwc_flux = '%0.1e'%(target_hwc_flux[hwc_idx])
        if not found_vts_name=='':
            found_counterpart += '\n'+'%s'%(found_vts_name)
        if not found_hss_name=='':
            found_counterpart += '\n'+'%s'%(found_hss_name)
        if not found_tev_name=='':
            found_counterpart += '\n'+'%s'%(found_tev_name)

        is_good_target = False
        if found_counterpart!='':
            is_good_target = True
    
        if is_good_target:
            my_psr_table.add_row([R+target_psr_name[psr]+N,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_counterpart,txt_age,txt_dist,txt_edot,txt_cu,txt_hwc_flux,txt_expo,txt_alt])
        else:
            my_psr_table.add_row([target_psr_name[psr],txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_counterpart,txt_age,txt_dist,txt_edot,txt_cu,txt_hwc_flux,txt_expo,txt_alt])
    
    print (my_psr_table)
    #lines_to_write += [my_psr_table]
    my_psr_table.clear_rows()

#with open('output_plots/result_pwn_search.txt', 'w') as f:
#    for line in lines_to_write:
#        f.write(line)
#        f.write('\n')

for ra in range(0,len(include_ra_bands)):

    if not select_snr: continue

    my_snr_table = PrettyTable(['%s band'%(include_ra_bands[ra]), 'RA/Dec (deg)', 'l/b (deg)', 'Fermi', 'HAWC', 'VERITAS', 'Dist (kpc)', 'Size (arcmin)', 'VTS expo (hr)', 'VTS alt (deg)'])
    for snr in range(0,len(target_snr_name)):
    
        myangle_ra = Angle(target_snr_ra[snr], my_unit.deg)
        include_this_snr = False
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+2.:
            include_this_snr = True
        if not include_this_snr: continue
    
        vts_alt = FindSourceVisibility('SNR %s'%(target_snr_name[snr]),target_snr_ra[snr],target_snr_dec[snr],obs='vts')
    
        found_fermi_name, frm_idx = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],fermi_name,fermi_ra,fermi_dec,0.2)
        found_tev_name  , tev_idx = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_tev_name,target_tev_ra,target_tev_dec,0.2)
        found_hwc_name  , hwc_idx = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_hwc_name,target_hwc_ra,target_hwc_dec,0.5)
        found_vts_name  , vts_idx = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_vts_name,target_vts_ra,target_vts_dec,0.3)
        found_hss_name  , hss_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_hss_name,target_hss_ra,target_hss_dec,0.3)
        gal_l, gal_b = ConvertRaDecToGalactic(target_snr_ra[snr],target_snr_dec[snr])
    
        exposure_time = FindVeritasExposure(target_snr_ra[snr],target_snr_dec[snr])
    
        if vts_alt<alt_cut: continue
        if select_new_target:
            if exposure_time>20.: continue
        else:
            if exposure_time<20.: continue

        #if myangle_ra.hour>17.:
        #    if vts_alt<70.: continue
        #    if found_hwc_name=='': continue
        #else:
        #    if vts_alt<50.: continue
        #    if found_hwc_name=='': continue
    
        txt_ra = '%0.2f'%(target_snr_ra[snr])
        txt_dec = '%0.2f'%(target_snr_dec[snr])
        txt_l = '%0.2f'%(gal_l)
        txt_b = '%0.2f'%(gal_b)
        txt_dist = '%0.1f'%(target_snr_dist[snr])
        txt_size = '%0.1f'%(target_snr_size[snr])
        txt_expo = '%0.1f'%(exposure_time)
        txt_alt = '%0.1f'%(vts_alt)
    
        is_good_target = False
        if found_hwc_name!='' and found_fermi_name!='':
            is_good_target = True
    
        if is_good_target:
            my_snr_table.add_row([R+target_snr_name[snr]+N,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_vts_name,txt_dist,txt_size,txt_expo,txt_alt])
        else:
            my_snr_table.add_row([target_snr_name[snr],txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_vts_name,txt_dist,txt_size,txt_expo,txt_alt])
    
    print (my_snr_table)
    my_snr_table.clear_rows()

print ('Search for Globular clusters. Find info about the table here: http://www.naic.edu/~pfreire/GCpsr.html')
for ra in range(0,len(include_ra_bands)):

    if not select_gcs: continue

    my_glb_table = PrettyTable(['%s band'%(include_ra_bands[ra]), 'RA/Dec (deg)', 'l/b (deg)', 'Dist (kpc)', 'N PSRs', 'Est. CU', 'VTS expo (hr)', 'VTS alt (deg)'])
    for glb in range(0,len(target_glb_name)):
    
        myangle_ra = Angle(target_glb_ra[glb], my_unit.deg)
        include_this_glb = False
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+2.:
            include_this_glb = True
        if not include_this_glb: continue
    
        vts_alt = FindSourceVisibility('%s'%(target_glb_name[glb]),target_glb_ra[glb],target_glb_dec[glb],obs='vts')
        hwc_alt = FindSourceVisibility('%s'%(target_glb_name[glb]),target_glb_ra[glb],target_glb_dec[glb],obs='hwc')
    
        found_fermi_name, frm_idx = FindFermiSource(target_glb_name[glb],target_glb_ra[glb],target_glb_dec[glb],fermi_name,fermi_ra,fermi_dec,0.2)
        found_tev_name  , tev_idx = FindFermiSource(target_glb_name[glb],target_glb_ra[glb],target_glb_dec[glb],target_tev_name,target_tev_ra,target_tev_dec,0.2)
        found_hwc_name  , hwc_idx = FindFermiSource(target_glb_name[glb],target_glb_ra[glb],target_glb_dec[glb],target_hwc_name,target_hwc_ra,target_hwc_dec,0.5)
        found_vts_name  , vts_idx = FindFermiSource(target_glb_name[glb],target_glb_ra[glb],target_glb_dec[glb],target_vts_name,target_vts_ra,target_vts_dec,0.3)
        found_hss_name  , hss_idx = FindFermiSource(target_glb_name[glb],target_glb_ra[glb],target_glb_dec[glb],target_hss_name,target_hss_ra,target_hss_dec,0.3)
        gal_l, gal_b = ConvertRaDecToGalactic(target_glb_ra[glb],target_glb_dec[glb])
    
        exposure_time = FindVeritasExposure(target_glb_ra[glb],target_glb_dec[glb])

        #if select_new_target:
        #    if exposure_time>20.: continue
        #else:
        #    if exposure_time<20.: continue

        if vts_alt<alt_cut: continue
    
        crab_edot = 4.5e+38 
        crab_dist = 2.000 # kpc
        crab_pwn_bfield = 125.0 # muG
        halo_bfield = 3. # muG
        tev_flux_in_cu = target_glb_flux[glb]

        txt_name = target_glb_name[glb].strip(' ')
        txt_ra = '%0.2f'%(target_glb_ra[glb])
        txt_dec = '%0.2f'%(target_glb_dec[glb])
        txt_l = '%0.2f'%(gal_l)
        txt_b = '%0.2f'%(gal_b)
        txt_dist = '%0.1f'%(target_glb_dist[glb])
        txt_npsr = '%0.1f'%(target_glb_npsr[glb])
        txt_cu = '%0.1e'%(tev_flux_in_cu)
        txt_expo = '%0.1f'%(exposure_time)
        txt_alt = '%0.1f'%(vts_alt)
    
        is_good_target = False
        if found_fermi_name!='' and found_hwc_name!='':
            is_good_target = True
    
        if is_good_target:
            my_glb_table.add_row([R+txt_name+N,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,txt_dist,txt_npsr,txt_cu,txt_expo,txt_alt])
        else:
            my_glb_table.add_row([txt_name,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,txt_dist,txt_npsr,txt_cu,txt_expo,txt_alt])
    
    print (my_glb_table)
    #lines_to_write += [my_glb_table]
    my_glb_table.clear_rows()


pulsar_type = 'PSR'

fig, ax = plt.subplots()
figsize_x = 6.4
figsize_y = 4.8
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

a_nongamma_psr_age = np.array(nongamma_psr_age)
a_nongamma_psr_dist = np.array(nongamma_psr_dist)
a_nongamma_psr_edot = np.array(nongamma_psr_edot)
a_nongamma_psr_flux = np.array(nongamma_psr_flux)
a_vtsobs_psr_age = np.array(vtsobs_psr_age)
a_vtsobs_psr_dist = np.array(vtsobs_psr_dist)
a_vtsobs_psr_edot = np.array(vtsobs_psr_edot)
a_vtsobs_psr_flux = np.array(vtsobs_psr_flux)
a_nonvtsobs_psr_age = np.array(nonvtsobs_psr_age)
a_nonvtsobs_psr_dist = np.array(nonvtsobs_psr_dist)
a_nonvtsobs_psr_edot = np.array(nonvtsobs_psr_edot)
a_nonvtsobs_psr_flux = np.array(nonvtsobs_psr_flux)
a_hss_psr_age = np.array(hss_psr_age)
a_hss_psr_dist = np.array(hss_psr_dist)
a_hss_psr_edot = np.array(hss_psr_edot)
a_hss_psr_flux = np.array(hss_psr_flux)
a_vts_psr_age = np.array(vts_psr_age)
a_vts_psr_dist = np.array(vts_psr_dist)
a_vts_psr_edot = np.array(vts_psr_edot)
a_vts_psr_flux = np.array(vts_psr_flux)
a_nonvts_psr_age = np.array(nonvts_psr_age)
a_nonvts_psr_dist = np.array(nonvts_psr_dist)
a_nonvts_psr_edot = np.array(nonvts_psr_edot)
a_nonvts_psr_flux = np.array(nonvts_psr_flux)
a_hwc_psr_age = np.array(hwc_psr_age)
a_hwc_psr_dist = np.array(hwc_psr_dist)
a_hwc_psr_edot = np.array(hwc_psr_edot)
a_hwc_psr_flux = np.array(hwc_psr_flux)
a_nonhwc_psr_age = np.array(nonhwc_psr_age)
a_nonhwc_psr_dist = np.array(nonhwc_psr_dist)
a_nonhwc_psr_edot = np.array(nonhwc_psr_edot)
a_nonhwc_psr_flux = np.array(nonhwc_psr_flux)
a_fermi_psr_age = np.array(fermi_psr_age)
a_fermi_psr_dist = np.array(fermi_psr_dist)
a_fermi_psr_edot = np.array(fermi_psr_edot)
a_fermi_psr_flux = np.array(fermi_psr_flux)
a_fermi_psr_gflu = np.array(fermi_psr_gflu)
a_nonfermi_psr_age = np.array(nonfermi_psr_age)
a_nonfermi_psr_dist = np.array(nonfermi_psr_dist)
a_nonfermi_psr_edot = np.array(nonfermi_psr_edot)
a_nonfermi_psr_flux = np.array(nonfermi_psr_flux)
a_fermihwc_psr_age = np.array(fermihwc_psr_age)
a_fermihwc_psr_dist = np.array(fermihwc_psr_dist)
a_fermihwc_psr_edot = np.array(fermihwc_psr_edot)
a_fermihwc_psr_flux = np.array(fermihwc_psr_flux)
a_fermihwc_psr_gflu = np.array(fermihwc_psr_gflu)
a_fermivts_psr_age = np.array(fermivts_psr_age)
a_fermivts_psr_dist = np.array(fermivts_psr_dist)
a_fermivts_psr_edot = np.array(fermivts_psr_edot)
a_fermivts_psr_flux = np.array(fermivts_psr_flux)
a_fermivts_psr_gflu = np.array(fermivts_psr_gflu)
a_fermihss_psr_age = np.array(fermihss_psr_age)
a_fermihss_psr_dist = np.array(fermihss_psr_dist)
a_fermihss_psr_edot = np.array(fermihss_psr_edot)
a_fermihss_psr_flux = np.array(fermihss_psr_flux)
a_fermihss_psr_gflu = np.array(fermihss_psr_gflu)
a_all_snr_size = np.array(target_snr_size)

edot_bins = [1e31,1e32,1e33,1e34,1e35,1e36,1e37,1e38,1e39]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonhwc_psr_edot, bins = edot_bins, label='Undetected',alpha=0.5)
plt.hist(a_hwc_psr_edot, bins = edot_bins, label='HAWC detection',alpha=0.5)
axbig.set_xscale('log')
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('$\dot{E}$ [erg/s]')
plotname = 'Hist_%s_Edot'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

#flux_bins = [1e31,1e32,1e33,1e34,1e35,1e36,1e37,1e38,1e39]
flux_bins = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.,1e1,1e2,1e3,1e4]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonhwc_psr_flux, bins = flux_bins, label='Undetected',alpha=0.5)
plt.hist(a_hwc_psr_flux, bins = flux_bins, label='HAWC detection',alpha=0.5)
plt.hist(a_vts_psr_flux, bins = flux_bins, label='VERITAS detection',alpha=0.5)
plt.hist(a_hss_psr_flux, bins = flux_bins, label='HESS detection',alpha=0.5)
axbig.set_xscale('log')
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
plotname = 'Hist_%s_Flux'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

age_bins = [pow(10,-0.5),pow(10,0.),pow(10,0.5),pow(10,1.),pow(10,1.5),pow(10,2.),pow(10,2.5),pow(10,3.),pow(10,3.5),pow(10,4.)]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonhwc_psr_age, bins = age_bins, label='Undetected',alpha=0.5)
plt.hist(a_hwc_psr_age, bins = age_bins, label='HAWC detection',alpha=0.5)
axbig.set_xscale('log')
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('Age [kyr]')
plotname = 'Hist_%s_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

age_bins = [pow(10,-0.5),pow(10,0.),pow(10,0.5),pow(10,1.),pow(10,1.5),pow(10,2.),pow(10,2.5),pow(10,3.),pow(10,3.5),pow(10,4.)]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonvtsobs_psr_age, bins = age_bins, label='Not observed',alpha=0.5)
plt.hist(a_vtsobs_psr_age, bins = age_bins, label='Observed',alpha=0.5)
plt.hist(a_vts_psr_age, bins = age_bins, label='Detected',alpha=0.5)
axbig.set_xscale('log')
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('Age [kyr]')
plotname = 'Hist_VTS_%s_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

dist_bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonhwc_psr_dist, bins = dist_bins, label='Undetected',alpha=0.5)
plt.hist(a_hwc_psr_dist, bins = dist_bins, label='HAWC detection',alpha=0.5)
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('Distance [kpc]')
plotname = 'Hist_%s_Dist'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

dist_bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_nonvtsobs_psr_dist, bins = dist_bins, label='Not observed',alpha=0.5)
plt.hist(a_vtsobs_psr_dist, bins = dist_bins, label='Observed',alpha=0.5)
plt.hist(a_vts_psr_dist, bins = dist_bins, label='Detected',alpha=0.5)
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('Distance [kpc]')
plotname = 'Hist_VTS_%s_Dist'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

size_bins = [0,5,10,15,20,25,30,35,40,45,50]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_all_snr_size, bins = size_bins, label='All SNRs',alpha=0.5)
axbig.set_yscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of SNRs')
axbig.set_xlabel('SNR size [arcmin]')
plotname = 'Hist_SNR_Size'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nonhwc_psr_flux, a_nonhwc_psr_age, marker='o', color='blue',alpha=0.2, label='Not a HAWC source')
plt.scatter(a_hwc_psr_flux, a_hwc_psr_age, marker='o', color='yellow',alpha=0.3, label='HAWC source')
plt.scatter(a_vts_psr_flux, a_vts_psr_age, marker='o', color='orange',alpha=0.3, label='VERITAS source')
plt.scatter(a_hss_psr_flux, a_hss_psr_age, marker='o', color='red',alpha=0.3, label='HESS source')
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_HAWC_%s_Flux_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nonvtsobs_psr_flux, a_nonvtsobs_psr_age, marker='o', color='blue',alpha=0.2, label='Not observed source')
plt.scatter(a_vtsobs_psr_flux, a_vtsobs_psr_age, marker='o', color='red',alpha=0.3, label='Observed source')
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_VTSObs_%s_Flux_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nonvts_psr_flux, a_nonvts_psr_age, marker='o', color='blue',alpha=0.2, label='Not a VERITAS source')
plt.scatter(a_vts_psr_flux, a_vts_psr_age, marker='o', color='red',alpha=0.3, label='VERITAS source')
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_VTS_%s_Flux_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nonfermi_psr_flux, a_nonfermi_psr_age, marker='o', color='blue',alpha=0.2, label='Not a Fermi source')
plt.scatter(a_fermi_psr_flux, a_fermi_psr_age, marker='o', color='red',alpha=0.3, label='Fermi source')
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_Fermi_%s_Flux_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nongamma_psr_flux, a_nongamma_psr_age, marker='o', color='blue',alpha=0.2, label='Not a GeV/TeV source')
plt.scatter(a_fermi_psr_flux, a_fermi_psr_age, marker='o', color='green',alpha=0.3, label='Fermi source')
plt.scatter(a_vts_psr_flux, a_vts_psr_age, marker='o', color='red',alpha=0.3, label='VERITAS/HESS/HAWC')
plt.scatter(a_hss_psr_flux, a_hss_psr_age, marker='o', color='red',alpha=0.3)
plt.scatter(a_hwc_psr_flux, a_hwc_psr_age, marker='o', color='red',alpha=0.3)
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_TeV_%s_Flux_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_nongamma_psr_edot, a_nongamma_psr_age, marker='o', color='blue',alpha=0.2, label='Not a GeV/TeV source')
plt.scatter(a_fermi_psr_edot, a_fermi_psr_age, marker='o', color='green',alpha=0.3, label='Fermi source')
plt.scatter(a_vts_psr_edot, a_vts_psr_age, marker='o', color='red',alpha=0.3, label='VERITAS/HESS/HAWC')
plt.scatter(a_hss_psr_edot, a_hss_psr_age, marker='o', color='red',alpha=0.3)
plt.scatter(a_hwc_psr_edot, a_hwc_psr_age, marker='o', color='red',alpha=0.3)
axbig.set_ylabel('Age [kyr]')
axbig.set_xlabel('$\dot{E}$ [erg/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_TeV_%s_Edot_Age'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_fermi_psr_flux, a_fermi_psr_gflu, marker='o', color='blue',alpha=0.3, label='Fermi-only source')
plt.scatter(a_fermihwc_psr_flux, a_fermihwc_psr_gflu, marker='o', color='red',alpha=0.3, label='HAWC source')
plt.scatter(a_fermivts_psr_flux, a_fermivts_psr_gflu, marker='o', color='green',alpha=0.3, label='VERTIAS source')
axbig.set_ylabel('Fermi flux [erg/$\mathrm{cm}^{2}$/s]')
axbig.set_xlabel('$\dot{E}/d^{2}$ [erg/$\mathrm{kpc}^{2}$/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
axbig.legend(loc='best')
plotname = 'Scatter_FermiHAWC_%s_Flux_gFlux'%(pulsar_type)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

