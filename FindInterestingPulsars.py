


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
        if distance<0.5:
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
        if 'flux measurements:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            target_name = ''
            target_ra = ''
            target_dec = ''
    return source_name, source_ra, source_dec

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
target_tev_name, target_tev_ra, target_tev_dec = ReadTeVCatTargetListFromFile('TeVCat_RaDec_w_Names.txt')
target_hwc_name, target_hwc_ra, target_hwc_dec = ReadHAWCTargetListFromFile('Cat_3HWC.txt')
target_vts_name, target_vts_ra, target_vts_dec = ReadVeritasTargetListFromFile('VTS_Cat.txt')

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

edot_cut = 5e34
#edot_cut = 5e35

R = "\033[0;31m" #RED
G = "\033[0;32m" # GREEN
Y = "\033[0;33m" # Yellow
B = "\033[0;34m" # Blue
N = "\033[0m" # Reset

detected_psr_edot = []
undetected_psr_edot = []
unknown_psr_edot = []
detected_psr_flux = []
undetected_psr_flux = []
unknown_psr_flux = []
detected_snr_flux = []
undetected_snr_flux = []


for ra in range(0,len(include_ra_bands)):
    my_psr_table = PrettyTable(['%s band'%(include_ra_bands[ra]), 'RA/Dec (deg)', 'l/b (deg)', 'Fermi', 'HAWC', 'TeVCat', 'Age (kyr)', 'Dist (kpc)', 'Edot (erg/s)', 'VTS expo (hr)', 'VTS alt (deg)'])
    for psr in range(0,len(target_psr_name)):
    
        myangle_ra = Angle(target_psr_ra[psr], my_unit.deg)
        include_this_psr = False
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+2.:
            include_this_psr = True
        if not include_this_psr: continue
    
        #if target_psr_age[psr]>1e6: continue
    
        max_alt = FindSourceVisibility('PSR %s'%(target_psr_name[psr]),target_psr_ra[psr],target_psr_dec[psr])
    
        found_fermi_name = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],fermi_name,fermi_ra,fermi_dec)
        found_tev_name = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_tev_name,target_tev_ra,target_tev_dec)
        found_hwc_name = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_hwc_name,target_hwc_ra,target_hwc_dec)
        found_vts_name = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_vts_name,target_vts_ra,target_vts_dec)
        gal_l, gal_b = ConvertRaDecToGalactic(target_psr_ra[psr],target_psr_dec[psr])
    
        exposure_time = FindVeritasExposure(target_psr_ra[psr],target_psr_dec[psr])
    
        #if exposure_time>20.:
        #    if not found_vts_name=='':
        #        if found_fermi_name!='':
        #            detected_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
        #            detected_psr_flux += [float(found_fermi_name)]
        #    else:
        #        if found_fermi_name!='':
        #            undetected_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
        #            undetected_psr_flux += [float(found_fermi_name)]
        #else:
        #    if found_fermi_name!='' and max_alt>50.:
        #        unknown_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
        #        unknown_psr_flux += [float(found_fermi_name)]
        if not found_hwc_name=='':
            if found_fermi_name!='':
                detected_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
                detected_psr_flux += [float(found_fermi_name)]
        else:
            if found_fermi_name!='':
                undetected_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
                undetected_psr_flux += [float(found_fermi_name)]
        if exposure_time<20.:
            if found_fermi_name!='' and max_alt>50.:
                unknown_psr_edot += [target_psr_edot[psr]/pow(target_psr_dist[psr],0)]
                unknown_psr_flux += [float(found_fermi_name)]
    
        if exposure_time>20.: continue
        if myangle_ra.hour>17.:
            if max_alt<70.: continue
            #if found_hwc_name=='': continue
            if found_fermi_name=='' and found_hwc_name=='': continue
            if found_fermi_name!='' and found_hwc_name=='':
                if float(found_fermi_name)<1e-5 and target_psr_edot[psr]<5e35: continue
        else:
            if max_alt<50.: continue
            #if found_hwc_name=='': continue
            if found_fermi_name=='' and found_hwc_name=='': continue
            if found_fermi_name!='' and found_hwc_name=='':
                if float(found_fermi_name)<1e-5 and target_psr_edot[psr]<5e34: continue
    
        txt_ra = '%0.2f'%(target_psr_ra[psr])
        txt_dec = '%0.2f'%(target_psr_dec[psr])
        txt_l = '%0.2f'%(gal_l)
        txt_b = '%0.2f'%(gal_b)
        txt_dist = '%0.1f'%(target_psr_dist[psr])
        txt_age = '%0.1e'%(target_psr_age[psr]/1000.)
        txt_edot = '%0.1e'%(target_psr_edot[psr])
        txt_expo = '%0.1f'%(exposure_time)
        txt_alt = '%0.1f'%(max_alt)
    
        is_good_target = False
        if found_fermi_name!='' and found_hwc_name!='':
            is_good_target = True
    
        if is_good_target:
            my_psr_table.add_row([R+target_psr_name[psr]+N,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_tev_name,txt_age,txt_dist,txt_edot,txt_expo,txt_alt])
        else:
            my_psr_table.add_row([target_psr_name[psr],txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_tev_name,txt_age,txt_dist,txt_edot,txt_expo,txt_alt])
    
    print (my_psr_table)
    my_psr_table.clear_rows()


for ra in range(0,len(include_ra_bands)):
    my_snr_table = PrettyTable(['%s band'%(include_ra_bands[ra]), 'RA/Dec (deg)', 'l/b (deg)', 'Fermi', 'HAWC', 'TeVCat', 'Dist (kpc)', 'VTS expo (hr)', 'VTS alt (deg)'])
    for snr in range(0,len(target_snr_name)):
    
        myangle_ra = Angle(target_snr_ra[snr], my_unit.deg)
        include_this_snr = False
        if myangle_ra.hour>float(include_ra_bands[ra]) and myangle_ra.hour<float(include_ra_bands[ra])+2.:
            include_this_snr = True
        if not include_this_snr: continue
    
        max_alt = FindSourceVisibility('SNR %s'%(target_snr_name[snr]),target_snr_ra[snr],target_snr_dec[snr])
    
        found_fermi_name = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],fermi_name,fermi_ra,fermi_dec)
        found_tev_name = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_tev_name,target_tev_ra,target_tev_dec)
        found_hwc_name = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_hwc_name,target_hwc_ra,target_hwc_dec)
        found_vts_name = FindFermiSource(target_snr_name[snr],target_snr_ra[snr],target_snr_dec[snr],target_vts_name,target_vts_ra,target_vts_dec)
        gal_l, gal_b = ConvertRaDecToGalactic(target_snr_ra[snr],target_snr_dec[snr])
    
        exposure_time = FindVeritasExposure(target_snr_ra[snr],target_snr_dec[snr])
    
        if exposure_time>20. and found_fermi_name!='':
            if not found_vts_name=='':
                detected_snr_flux += [float(found_fermi_name)]
            else:
                undetected_snr_flux += [float(found_fermi_name)]
    
        if exposure_time>20.: continue
        if myangle_ra.hour>17.:
            if max_alt<70.: continue
            if found_hwc_name=='': continue
        else:
            if max_alt<50.: continue
            if found_hwc_name=='': continue
    
        txt_ra = '%0.2f'%(target_snr_ra[snr])
        txt_dec = '%0.2f'%(target_snr_dec[snr])
        txt_l = '%0.2f'%(gal_l)
        txt_b = '%0.2f'%(gal_b)
        txt_dist = '%0.1f'%(target_snr_dist[snr])
        txt_expo = '%0.1f'%(exposure_time)
        txt_alt = '%0.1f'%(max_alt)
    
        is_good_target = False
        if found_hwc_name!='' and found_fermi_name!='':
            is_good_target = True
    
        if is_good_target:
            my_snr_table.add_row([R+target_snr_name[snr]+N,txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_tev_name,txt_dist,txt_expo,txt_alt])
        else:
            my_snr_table.add_row([target_snr_name[snr],txt_ra+'\n'+txt_dec,txt_l+'\n'+txt_b,found_fermi_name,found_hwc_name,found_tev_name,txt_dist,txt_expo,txt_alt])
    
    print (my_snr_table)
    my_snr_table.clear_rows()

fig, ax = plt.subplots()
figsize_x = 6.4
figsize_y = 4.8
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

a_detected_psr_edot = np.array(detected_psr_edot)
a_undetected_psr_edot = np.array(undetected_psr_edot)
a_unknown_psr_edot = np.array(unknown_psr_edot)
a_detected_psr_flux = np.array(detected_psr_flux)
a_undetected_psr_flux = np.array(undetected_psr_flux)
a_unknown_psr_flux = np.array(unknown_psr_flux)

edot_bins = [1e31,1e32,1e33,1e34,1e35,1e36,1e37,1e38,1e39]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_undetected_psr_edot, bins = edot_bins, label='Undetected',alpha=0.5)
plt.hist(a_detected_psr_edot, bins = edot_bins, label='VTS detection',alpha=0.5)
axbig.set_xscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('$\dot{E}$ [erg/s]')
plotname = 'Hist_PSR_Edot'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

flux_bins = [pow(10,-11),pow(10,-10.5),pow(10,-10),pow(10,-9.5),pow(10,-9),pow(10,-8.5),pow(10,-8),pow(10,-7.5),pow(10,-7),pow(10,-6.5),pow(10,-6)]
fig.clf()
axbig = fig.add_subplot()
plt.hist(a_undetected_psr_flux, bins = flux_bins, label='Undetected',alpha=0.5)
plt.hist(a_detected_psr_flux, bins = flux_bins, label='VTS detection',alpha=0.5)
axbig.set_xscale('log')
axbig.legend(loc='best')
axbig.set_ylabel('number of pulsars')
axbig.set_xlabel('Fermi flux at 1 GeV [/cm2/s]')
plotname = 'Hist_PSR_Flux'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

fig.clf()
axbig = fig.add_subplot()
plt.scatter(a_undetected_psr_edot, a_undetected_psr_flux, marker='o', color='blue',alpha=0.5)
plt.scatter(a_detected_psr_edot, a_detected_psr_flux, marker='o', color='red',alpha=0.5)
#plt.scatter(a_unknown_psr_edot, a_unknown_psr_flux, marker='o', color='green',alpha=0.5)
axbig.set_ylabel('Fermi flux at 1 GeV [/cm2/s]')
axbig.set_xlabel('$\dot{E}$ [erg/s]')
axbig.set_yscale('log')
axbig.set_xscale('log')
plotname = 'Scatter_PSR_Edot_Flux'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
