


from prettytable import PrettyTable
import csv
import numpy as np
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

def FindSourceVisibility(psr_name,psr_ra,psr_dec):

    current_time = '2022-04-01 23:00:00'
    #psr_coord = SkyCoord.from_name(psr_name)
    psr_coord = SkyCoord(psr_ra, psr_dec, unit="deg")

    veritas_site = EarthLocation(lat=31.6751*u.deg, lon=-110.952*u.deg, height=1268*u.m)
    utcoffset = -7*u.hour  # AZ time
    time = Time(current_time) - utcoffset
    midnight = Time(current_time) - utcoffset

    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_whole_night = midnight + delta_midnight
    frame_whole_night = AltAz(obstime=times_whole_night, location=veritas_site)

    psr_altazs_whole_night = psr_coord.transform_to(frame_whole_night)
    max_alt = np.max(psr_altazs_whole_night.alt.deg)
    #print ('max alt = %0.1f deg'%(max_alt))
    return max_alt

    #psr_altaz = psr_coord.transform_to(AltAz(obstime=time,location=veritas_site))
    #print(f"Pulsar's Altitude = {psr_altaz.alt:.2}")

def FindFermiSource(psr_name,psr_ra,psr_dec,fermi_name,fermi_ra,fermi_dec):

    found_fermi_name = ''
    for fermi in range(0,len(fermi_name)):
        delta_ra = psr_ra - fermi_ra[fermi]
        delta_dec = psr_dec - fermi_dec[fermi]
        distance = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5)
        if distance<0.2:
            found_fermi_name = fermi_name[fermi]
    return found_fermi_name

def ReadFermiCatelog():
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open('gll_psc_v26.xml')
    target_name = ''
    target_info = ''
    target_ra = ''
    target_dec = ''
    for line in inputFile:
        if line.split(' ')[0]=='<source':
            for block in range(0,len(line.split(' '))):
                if 'Unc_' in line.split(' ')[block]: continue
                if 'name=' in line.split(' ')[block]:
                    target_name = line.split('name="')[1].split('"')[0]
                if 'Flux1000=' in line.split(' ')[block]:
                    target_info = line.split(' ')[block]
                    target_info = target_info.strip('>\n')
                    target_info = target_info.strip('"')
                    target_info = target_info.lstrip('Flux1000="')
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
                target_info = ''
                target_ra = ''
                target_dec = ''
                continue
            #print ('target_name = %s'%(target_name))
            #print ('target_ra = %s'%(target_ra))
            #print ('target_dec = %s'%(target_dec))
            #source_name += [target_name]
            source_name += ['%0.2e'%(float(target_info))]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
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
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '0'
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            source_dist += [float(target_min_dist)]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            #print(row)
    return source_name, source_ra, source_dec, source_dist

def ReadATNFTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
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

        if float(target_brightness)<1e33 and float(target_edot)<1e34: continue
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
    return source_name, source_ra, source_dec, source_dist, source_age, source_edot

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
target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age, target_psr_edot = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
target_snr_name, target_snr_ra, target_snr_dec, target_snr_dist = ReadSNRTargetListFromCSVFile()

include_ra_bands = []
#include_ra_bands += ['00','01']
#include_ra_bands += ['02','03']
#include_ra_bands += ['04','05']
#include_ra_bands += ['06','07']
#include_ra_bands += ['08','09']
#include_ra_bands += ['10','11']
#include_ra_bands += ['12','13']
#include_ra_bands += ['14','15']
#include_ra_bands += ['16','17']
#include_ra_bands += ['18','19']
#include_ra_bands += ['20','21']
include_ra_bands += ['22','23']

my_psr_table = PrettyTable(['%s-%s band'%(include_ra_bands[0],include_ra_bands[1]), 'RA (deg)', 'Dec (deg)', 'l (deg)', 'b (deg)', 'Fermi source', 'Dist (kpc)', 'Age (kyr)', 'Edot (erg/s)', 'VTS expo (hr)', 'VTS max alt (deg)'])
for psr in range(0,len(target_psr_name)):

    myangle_ra = Angle(target_psr_ra[psr], my_unit.deg)
    include_this_psr = False
    for ra in range(0,len(include_ra_bands)):
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+1.:
            include_this_psr = True
    if not include_this_psr: continue

    #if target_psr_age[psr]>1e6: continue

    max_alt = FindSourceVisibility('PSR %s'%(target_psr_name[psr]),target_psr_ra[psr],target_psr_dec[psr])
    if max_alt<45.: continue
    #if max_alt<75.: continue

    found_fermi_name = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],fermi_name,fermi_ra,fermi_dec)
    gal_l, gal_b = ConvertRaDecToGalactic(target_psr_ra[psr],target_psr_dec[psr])

    exposure_time = FindVeritasExposure(target_psr_ra[psr],target_psr_dec[psr])
    if exposure_time<5.: continue
    #if exposure_time>30.: continue

    txt_ra = '%0.2f'%(target_psr_ra[psr])
    txt_dec = '%0.2f'%(target_psr_dec[psr])
    txt_l = '%0.2f'%(gal_l)
    txt_b = '%0.2f'%(gal_b)
    txt_dist = '%0.1f'%(target_psr_dist[psr])
    txt_age = '%0.1e'%(target_psr_age[psr]/1000.)
    txt_edot = '%0.1e'%(target_psr_edot[psr])
    txt_expo = '%0.1f'%(exposure_time)
    txt_alt = '%0.1f'%(max_alt)
    my_psr_table.add_row([target_psr_name[psr],txt_ra,txt_dec,txt_l,txt_b,found_fermi_name,txt_dist,txt_age,txt_edot,txt_expo,txt_alt])
print (my_psr_table)


my_snr_table = PrettyTable(['%s-%s band'%(include_ra_bands[0],include_ra_bands[1]), 'RA (deg)', 'Dec (deg)', 'l (deg)', 'b (deg)', 'Fermi source', 'Dist (kpc)', 'VTS expo (hr)', 'VTS max alt (deg)'])
for snr in range(0,len(target_snr_name)):

    myangle_ra = Angle(target_snr_ra[snr], my_unit.deg)
    include_this_snr = False
    for ra in range(0,len(include_ra_bands)):
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+1.:
            include_this_snr = True
    if not include_this_snr: continue

    max_alt = FindSourceVisibility('SNR %s'%(target_snr_name[snr]),target_snr_ra[snr],target_snr_dec[snr])
    if max_alt<45.: continue
    #if max_alt<75.: continue

    found_fermi_name = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],fermi_name,fermi_ra,fermi_dec)
    gal_l, gal_b = ConvertRaDecToGalactic(target_snr_ra[snr],target_snr_dec[snr])

    exposure_time = FindVeritasExposure(target_snr_ra[snr],target_snr_dec[snr])
    if exposure_time<5.: continue
    #if exposure_time>30.: continue

    txt_ra = '%0.2f'%(target_snr_ra[snr])
    txt_dec = '%0.2f'%(target_snr_dec[snr])
    txt_l = '%0.2f'%(gal_l)
    txt_b = '%0.2f'%(gal_b)
    txt_dist = '%0.1f'%(target_snr_dist[snr])
    txt_expo = '%0.1f'%(exposure_time)
    txt_alt = '%0.1f'%(max_alt)
    my_snr_table.add_row([target_snr_name[snr],txt_ra,txt_dec,txt_l,txt_b,found_fermi_name,txt_dist,txt_expo,txt_alt])
print (my_snr_table)

