


from prettytable import PrettyTable
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


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


target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age, target_psr_edot = ReadATNFTargetListFromFile('ATNF_pulsar_list.txt')

#include_ra_bands = ['00','01']
#include_ra_bands = ['02','03']
#include_ra_bands = ['04','05']
#include_ra_bands = ['06','07']
#include_ra_bands = ['08','09']
#include_ra_bands = ['10','11']
#include_ra_bands = ['12','13']
#include_ra_bands = ['14','15']
#include_ra_bands = ['16','17']
#include_ra_bands = ['18','19']
#include_ra_bands = ['20','21']
include_ra_bands = ['22','23']

mytable = PrettyTable(['%s-%s band'%(include_ra_bands[0],include_ra_bands[1]), 'Dist (kpc)', 'Age (kyr)', 'Edot (erg/s)', 'VTS expo (hr)', 'VTS max alt (deg)'])
mytable.float_format["Dist (kpc)"] = "0.1f"
mytable.float_format["Age (kyr)"] = "0.1e"
mytable.float_format["Edot (erg/s)"] = "0.1e"
mytable.float_format["VTS expo (hr)"] = "0.1f"
mytable.float_format["VTS max alt (deg)"] = "0.1f"
for psr in range(0,len(target_psr_name)):

    include_this_psr = False
    for ra in range(0,len(include_ra_bands)):
        if 'B%s'%(include_ra_bands[ra]) in target_psr_name[psr]:
            include_this_psr = True
        if 'J%s'%(include_ra_bands[ra]) in target_psr_name[psr]:
            include_this_psr = True
    if not include_this_psr: continue

    #if target_psr_age[psr]>1e6: continue

    max_alt = FindSourceVisibility('PSR %s'%(target_psr_name[psr]),target_psr_ra[psr],target_psr_dec[psr])
    if max_alt<45.: continue
    #if max_alt<75.: continue

    exposure_time = FindVeritasExposure(target_psr_ra[psr],target_psr_dec[psr])
    mytable.add_row([target_psr_name[psr],target_psr_dist[psr],target_psr_age[psr]/1000.,target_psr_edot[psr],exposure_time,max_alt])
print (mytable)


