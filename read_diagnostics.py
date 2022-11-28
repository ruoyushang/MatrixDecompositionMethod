
import os
import sys,ROOT
import array
import math
import csv
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from scipy import special

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

def gaussian_2d(x, y, x0, y0, xsig, ysig):
    return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))

def GetDistance(x1,y1,x2,y2):
    distance = pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5)
    return distance

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b;
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

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

def ReadSNRTargetListFromCSVFile(snr_name=''):
    source_name = []
    source_ra = []
    source_dec = []
    with open('SNRcat20221001-SNR.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row)==0: continue
            if '#' in row[0]: continue
            target_name = row[0]
            if not snr_name=='':
                if not target_name==snr_name:
                    continue
                else:
                    print ('target_name = %s'%(target_name))
            target_min_dist = row[13]
            if target_min_dist=='':
                target_min_dist = '0'
            #if float(target_min_dist)>6.: continue
            target_ra = row[19]
            target_dec = row[20]
            source_name += [target_name]
            source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
            source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
            #print('target_min_dist = %s'%(target_min_dist))
            #print('source_name = %s'%(source_name[len(source_name)-1]))
            #print('source_ra = %0.2f'%(source_ra[len(source_ra)-1]))
            #print('source_dec = %0.2f'%(source_dec[len(source_dec)-1]))
            print(row)
    return source_name, source_ra, source_dec

def ReadSNRTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_l = line.split(',')[0].strip(" ")
        target_b = line.split(',')[1].strip(" ")
        target_name = "G"+target_l+target_b
        target_ra = line.split(',')[2].strip(" ")
        target_dec = line.split(',')[3].strip(" ")
        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
    return source_name, source_ra, source_dec

def ReadATNFTargetListFromFileForPlot(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        #print ('target_ra = %s'%(target_ra))
        #print ('target_dec = %s'%(target_dec))
        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
    return source_name, source_ra, source_dec

def ReadATNFTargetListFromFile(psr_name=''):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    file_path = 'ATNF_pulsar_full_list.txt'
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        if not 'J' in line and not 'B' in line: continue
        target_name = line.split(',')[0].strip(" ")
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        if not psr_name=='':
            if not target_name==psr_name:
                continue
            else:
                print ('target_name = %s'%(target_name))
        #print ('target_ra = %s'%(target_ra))
        #print ('target_dec = %s'%(target_dec))
        target_dist = float(line.split(',')[3].strip(" "))
        target_age = float(line.split(',')[4].strip(" "))
        target_edot = float(line.split(',')[5].strip(" "))
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        #if target_brightness<1e33 and target_edot<1e34: continue
        #if target_dist>6.0: continue
        #if target_dist>2.0: continue
        #if target_age<1e4: continue

        #ra_deg = float(HMS2deg(target_ra,target_dec)[0])
        #dec_deg = float(HMS2deg(target_ra,target_dec)[1])
        #gal_l, gal_b = ConvertRaDecToGalactic(ra_deg,dec_deg)
        #if abs(gal_b)<5.: continue

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [target_dist]
        source_age += [target_age]
        source_edot += [target_edot]
    return source_name, source_ra, source_dec


target_name = []
target_ra = []
target_dec = []

#target_name += ["Crab"]
#target_ra += [83.6332083333]
#target_dec += [22.0144722222]

#target_name += ["Geminga"]
#target_ra += [98.117]
#target_dec += [17.367]
#target_name += ["VER J2019+407"]
#target_ra += [305.020]
#target_dec += [40.757]
#target_name += ["MGRO J1908+06"]
#target_ra += [286.975]
#target_dec += [6.269]
#target_name += ["Boomerang"]
#target_ra += [337.183]
#target_dec += [61.167]
#target_name += ["Cas A"]
#target_ra += [350.8075000]
#target_dec += [58.8072222]
#target_name += ["Ursa Minor"]
#target_ra += [227.2854167]
#target_dec += [67.2225000]
#target_name += ["RGB J0710+591"]
#target_ra += [107.6100000]
#target_dec += [59.1500000]
#target_name += ["LHAASO J2032+4102"]
#target_ra += [308.0500000]
#target_dec += [41.0500000]
#target_name += ["LHAASO J0621+3755"]
#target_ra += [95.47]
#target_dec += [37.92]

target_psr_name = []
target_psr_ra = []
target_psr_dec = []
#target_psr_name, target_psr_ra, target_psr_dec = ReadATNFTargetListFromFile('J2238+5903')

target_snr_name = []
target_snr_ra = []
target_snr_dec = []
target_snr_name, target_snr_ra, target_snr_dec = ReadSNRTargetListFromCSVFile('G111.7-02.1')

target_name = []
target_ra = []
target_dec = []
target_name += target_psr_name
target_ra += target_psr_ra
target_dec += target_psr_dec
target_name += target_snr_name
target_ra += target_snr_ra
target_dec += target_snr_dec

range_ra = 1.5
range_dec = 1.5

#epoch = 'V4'
#epoch = 'V5'
epoch = 'V6'


search_for_on_data = True
#search_for_on_data = False

#search_for_rhv = True
search_for_rhv = False

plot_tag = 'OFF'
if search_for_on_data:
    plot_tag = 'ON'

gal_b_cut = 5.
Galactic_observation = True
#Galactic_observation = False

Search_Range_RA = [0.,360.]
Search_Range_Dec = [-90.,90.]

Search_Range_Gal_l = [0.,360.]
Search_Range_Gal_b = [-90.,90.]
#Search_Range_Gal_l = [150.,200.]
#Search_Range_Gal_b = [-90.,90.]

search_north = 0 # north+south
#search_north = 1 # north
#search_north = 2 # south
Search_Range_Elev = [0.,90.]
Search_Range_Azim = [0.,360.]
Search_Range_PedVar_DC = [3.,10.]
#Search_Range_PedVar_DC = [0.,5.]
#Search_Range_PedVar_DC = [5.,6.]
#Search_Range_PedVar_DC = [6.,10.]
if search_for_rhv:
    Search_Range_PedVar_DC = [0.,4.5] #RHV

#Search_Range_Elev = [0.,35.]
#Search_Range_RA = [0.,50.]
#Search_Range_Dec = [0,90.]


#Search_Range_RA = [Search_Center_RA-4.,Search_Center_RA+4.]
#Search_Range_Dec = [Search_Center_Dec-4.,Search_Center_Dec+4.]

RunNumber = 0
Elev = 0
Azim = 0
PedVar_DC = 0
PedVar_PE = 0
T1_RA = 0
T1_Dec = 0
T2_RA = 0
T2_Dec = 0
T3_RA = 0
T3_Dec = 0
T4_RA = 0
T4_Dec = 0
FIR_Mean = 0.
FIR_RMS = 0.
L3_rate = 0.
Livetime = 0.

Total_Livetime = 0.

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
#if not search_for_on_data:
#    #Search_Range_RA = [60.,75.]
#    #Search_Range_Dec = [72.,76.]
#    #Search_Range_RA = [75.,80.]
#    #Search_Range_Dec = [66.,70.]
#    sourceFile = open('runlist_backup/RunList_MGRO_J1908_V5.txt')
#    for line in sourceFile:
#        Source_RunNumber += [int(linei.split(',')[1])]
#        Source_PedVar_DC += [0.]
#        Source_Elev += [0.]
#        Source_Azim += [0.]
#        Source_Livetime += [0.]

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


List_Used = []
List_Used_Exposure = []
List_Used_RA = []
List_Used_Dec = []
List_Used_Gal_l = []
List_Used_Gal_b = []
List_Used_Elev = []
List_Used_Azim = []
List_Used_NSB = []

if search_for_on_data:

    for target in range(0,len(target_ra)):

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
        
            if int(RunNumber)<46642:
                if not epoch=='V4': continue
            if int(RunNumber)>=46642 and int(RunNumber)<63373:
                if not epoch=='V5': continue
            if int(RunNumber)>=63373:
                if not epoch=='V6': continue
        
            if Elev<Search_Range_Elev[0]: continue
            if Elev>Search_Range_Elev[1]: continue
            if Azim<Search_Range_Azim[0]: continue
            if Azim>Search_Range_Azim[1]: continue
            if PedVar_DC<Search_Range_PedVar_DC[0]: continue
            if PedVar_DC>Search_Range_PedVar_DC[1]: continue

            #if MJD<59480: continue

            if epoch=='V6':
                if search_for_rhv: 
                    if L3_rate>280.: continue # RHV
                    #if L3_rate<50.: continue # RHV
                else:
                    if L3_rate<250.: continue
            else:
                if L3_rate<150.: continue
            #if L3_rate>450.: continue
            if Livetime<5.: continue
            #if Livetime<20.: continue
        
            if abs(float(T1_RA)-target_ra[target])>range_ra: continue
            if abs(float(T1_Dec)-target_dec[target])>range_dec: continue
            distance_to_source = pow(pow(float(T1_RA)-target_ra[target],2)+pow(float(T1_Dec)-target_dec[target],2),0.5)
            #if distance_to_source<0.4: continue
            #if distance_to_source<0.8: continue
            #if distance_to_source<1.1: continue
            #if distance_to_source>2.0: continue
            #if 'wobble1' in offset: continue
            #if not 'wobble1' in offset: continue
            #if not 'geminga-' in Name: continue
            #if int(RunNumber)<=99358: continue

            DoNotUse = False
            for entry2 in range(0,len(Source_RunNumber)):
                if int(Source_RunNumber[entry2])==int(RunNumber): DoNotUse = True
            if DoNotUse: continue
        
            print ('RunNumber %s, L3_rate %s, Livetime %s, Elev %s, RA %s, Dec %s, offset %s'%(RunNumber,L3_rate,Livetime,Elev,T1_RA,T1_Dec, offset))
        
            List_Used += [RunNumber]
            List_Used_Exposure += [Livetime/60.]
            List_Used_RA += [T1_RA]
            List_Used_Dec += [T1_Dec]
            gal_l, gal_b = ConvertRaDecToGalactic(T1_RA,T1_Dec)
            List_Used_Gal_l += [gal_l]
            List_Used_Gal_b += [gal_b]
            List_Used_Elev += [Elev]
            List_Used_Azim += [Azim]
            List_Used_NSB += [PedVar_DC]

            Total_Livetime += Livetime/60.


if not search_for_on_data:

    List_Produced = []

    for entry in range(0,len(List_RunNumber)):

        RunNumber = List_RunNumber[entry]
        Elev = List_Elev[entry]
        Azim = List_Azim[entry]
        T1_RA = List_T1_RA[entry]
        T1_Dec = List_T1_Dec[entry]
        T2_RA = List_T2_RA[entry]
        T2_Dec = List_T2_Dec[entry]
        T3_RA = List_T3_RA[entry]
        T3_Dec = List_T3_Dec[entry]
        T4_RA = List_T4_RA[entry]
        T4_Dec = List_T4_Dec[entry]
        PedVar_DC = List_PedVar_DC[entry]
        PedVar_PE = List_PedVar_PE[entry]
        L3_rate = List_L3_rate[entry]
        Livetime = List_Livetime[entry]

        near_a_source = False
        for target in range(0,len(target_ra)):
            if abs(float(T1_RA)-target_ra[target])<range_ra and abs(float(T1_Dec)-target_dec[target])<range_dec:
                near_a_source = True
        if near_a_source: continue

        if epoch=='V6':
            if search_for_rhv: 
                if L3_rate>280.: continue # RHV
                #if L3_rate<50.: continue # RHV
            else:
                if L3_rate<250.: continue
        else:
            if L3_rate<150.: continue
        if Livetime<15.: continue
        if int(RunNumber)<46642:
            if not epoch=='V4': continue
        if int(RunNumber)>=46642 and int(RunNumber)<63373:
            if not epoch=='V5': continue
        if int(RunNumber)>=63373:
            if not epoch=='V6': continue

        if not (T1_RA==0. and T1_Dec==0.):
            gal_l, gal_b = ConvertRaDecToGalactic(T1_RA,T1_Dec)
            if Galactic_observation:
                if abs(gal_b)>gal_b_cut: continue
            else:
                if abs(gal_b)<gal_b_cut: continue
            if abs(gal_l)<Search_Range_Gal_l[0]: continue
            if abs(gal_l)>Search_Range_Gal_l[1]: continue
            if abs(gal_b)<Search_Range_Gal_b[0]: continue
            if abs(gal_b)>Search_Range_Gal_b[1]: continue
            if GetDistance(T1_RA,T1_Dec,83.633,22.014)<3.: continue  # Crab
            if GetDistance(T1_RA,T1_Dec,166.079,38.195)<3.: continue  # Mrk 421
            if GetDistance(T1_RA,T1_Dec,98.117,17.367)<3.: continue  # Geminga
        elif not (T2_RA==0. and T2_Dec==0.):
            gal_l, gal_b = ConvertRaDecToGalactic(T2_RA,T2_Dec)
            if Galactic_observation:
                if abs(gal_b)>gal_b_cut: continue
            else:
                if abs(gal_b)<gal_b_cut: continue
            if abs(gal_l)<Search_Range_Gal_l[0]: continue
            if abs(gal_l)>Search_Range_Gal_l[1]: continue
            if abs(gal_b)<Search_Range_Gal_b[0]: continue
            if abs(gal_b)>Search_Range_Gal_b[1]: continue
            if GetDistance(T2_RA,T2_Dec,83.633,22.014)<3.: continue  # Crab
            if GetDistance(T2_RA,T2_Dec,166.079,38.195)<3.: continue  # Mrk 421
            if GetDistance(T2_RA,T2_Dec,98.117,17.367)<3.: continue  # Geminga
        elif not (T3_RA==0. and T3_Dec==0.):
            gal_l, gal_b = ConvertRaDecToGalactic(T3_RA,T3_Dec)
            if Galactic_observation:
                if abs(gal_b)>gal_b_cut: continue
            else:
                if abs(gal_b)<gal_b_cut: continue
            if abs(gal_l)<Search_Range_Gal_l[0]: continue
            if abs(gal_l)>Search_Range_Gal_l[1]: continue
            if abs(gal_b)<Search_Range_Gal_b[0]: continue
            if abs(gal_b)>Search_Range_Gal_b[1]: continue
            if GetDistance(T3_RA,T3_Dec,83.633,22.014)<3.: continue  # Crab
            if GetDistance(T3_RA,T3_Dec,166.079,38.195)<3.: continue  # Mrk 421
            if GetDistance(T3_RA,T3_Dec,98.117,17.367)<3.: continue  # Geminga
        elif not (T4_RA==0. and T4_Dec==0.):
            gal_l, gal_b = ConvertRaDecToGalactic(T4_RA,T4_Dec)
            if Galactic_observation:
                if abs(gal_b)>gal_b_cut: continue
            else:
                if abs(gal_b)<gal_b_cut: continue
            if abs(gal_l)<Search_Range_Gal_l[0]: continue
            if abs(gal_l)>Search_Range_Gal_l[1]: continue
            if abs(gal_b)<Search_Range_Gal_b[0]: continue
            if abs(gal_b)>Search_Range_Gal_b[1]: continue
            if GetDistance(T4_RA,T4_Dec,83.633,22.014)<3.: continue  # Crab
            if GetDistance(T4_RA,T4_Dec,166.079,38.195)<3.: continue  # Mrk 421
            if GetDistance(T4_RA,T4_Dec,98.117,17.367)<3.: continue  # Geminga
        else: continue # there is no pointing info

        if (T1_RA<Search_Range_RA[0]): continue
        if (T1_RA>Search_Range_RA[1]): continue
        if (T1_Dec<Search_Range_Dec[0]): continue
        if (T1_Dec>Search_Range_Dec[1]): continue
        if (T2_RA<Search_Range_RA[0]): continue
        if (T2_RA>Search_Range_RA[1]): continue
        if (T2_Dec<Search_Range_Dec[0]): continue
        if (T2_Dec>Search_Range_Dec[1]): continue
        if (T3_RA<Search_Range_RA[0]): continue
        if (T3_RA>Search_Range_RA[1]): continue
        if (T3_Dec<Search_Range_Dec[0]): continue
        if (T3_Dec>Search_Range_Dec[1]): continue
        if (T4_RA<Search_Range_RA[0]): continue
        if (T4_RA>Search_Range_RA[1]): continue
        if (T4_Dec<Search_Range_Dec[0]): continue
        if (T4_Dec>Search_Range_Dec[1]): continue

        if Elev<Search_Range_Elev[0]: continue
        if Elev>Search_Range_Elev[1]: continue
        if Azim<Search_Range_Azim[0]: continue
        if Azim>Search_Range_Azim[1]: continue
        if search_north==1:
            if Azim>180.-90. and Azim<180.+90.: continue
        if search_north==2:
            if Azim<180.-90. or Azim>180.+90.: continue
        if PedVar_DC<Search_Range_PedVar_DC[0]: continue
        if PedVar_DC>Search_Range_PedVar_DC[1]: continue
        already_used = False
        for entry2 in range(0,len(List_Used)):
            if int(List_Used[entry2])==int(RunNumber): 
                #print 'RunNumber %s is in the List_Produced'%(RunNumber)
                already_used = True
        if already_used: continue
        List_Used += [RunNumber]
        List_Used_Exposure += [Livetime/60.]
        List_Used_RA += [T1_RA]
        List_Used_Dec += [T1_Dec]
        gal_l, gal_b = ConvertRaDecToGalactic(T1_RA,T1_Dec)
        List_Used_Gal_l += [gal_l]
        List_Used_Gal_b += [gal_b]
        List_Used_Elev += [Elev]
        List_Used_Azim += [Azim]
        List_Used_NSB += [PedVar_DC]

        Total_Livetime += Livetime/60.

        
for entry in range(0,len(List_Used)):
    print (List_Used[entry])
print ('total %s runs.'%(len(List_Used)))
print ('total %s hours.'%(Total_Livetime))

plt.clf()
fig, ax = plt.subplots()
w = 4
n = math.ceil((max(List_Used_RA) - min(List_Used_RA))/w)
n = max(n,1)
plt.hist(List_Used_RA, bins=n)
ax.axis('on')
ax.set_xlabel('RA')
ax.set_ylabel('counts')
plt.savefig("output_plots/RunRA_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
w = 4
n = math.ceil((max(List_Used_Dec) - min(List_Used_Dec))/w)
n = max(n,1)
plt.hist(List_Used_Dec, bins=n)
ax.axis('on')
ax.set_xlabel('Dec')
ax.set_ylabel('counts')
plt.savefig("output_plots/RunDec_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
w = 4
n = math.ceil((max(List_Used_Elev) - min(List_Used_Elev))/w)
n = max(n,1)
plt.hist(List_Used_Elev, bins=n)
ax.axis('on')
ax.set_xlabel('Elev')
ax.set_ylabel('counts')
plt.savefig("output_plots/RunElev_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
w = 4
n = math.ceil((max(List_Used_Azim) - min(List_Used_Azim))/w)
n = max(n,1)
plt.hist(List_Used_Azim, bins=n)
ax.axis('on')
ax.set_xlabel('Azim')
ax.set_ylabel('counts')
plt.savefig("output_plots/RunAzim_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
w = 0.2
n = math.ceil((max(List_Used_NSB) - min(List_Used_NSB))/w)
n = max(n,1)
plt.hist(List_Used_NSB, bins=n)
ax.axis('on')
ax.set_xlabel('NSB')
ax.set_ylabel('counts')
plt.savefig("output_plots/RunNSB_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
hist, xbins, ybins, im = plt.hist2d(List_Used_RA, List_Used_Dec, range = [[0,360], [-90,90]], bins=(90, 45), cmap=plt.cm.Greys)
plt.colorbar()
max_cnt = 0.
max_ra = 0.
max_dec = 0.
for i in range(len(xbins)-1):
    for j in range(len(ybins)-1):
        if hist[i,j]>max_cnt:
            max_cnt = hist[i,j]
            max_ra = xbins[i]
            max_dec = ybins[j]
min_dist = 1e10;
plot_ra = 0.
plot_dec = 0.
plot_name = ''
inputFile = open('TeVCat_RaDec_w_Names.txt')
for line in inputFile:
    gamma_source_name = line.split(',')[0]
    gamma_source_ra = float(line.split(',')[1])
    gamma_source_dec = float(line.split(',')[2])
    distance = pow(gamma_source_ra-max_ra,2)+pow(gamma_source_dec-max_dec,2)
    if distance<min_dist:
        min_dist = distance
        plot_name = gamma_source_name
        plot_ra = gamma_source_ra
        plot_dec = gamma_source_dec
plt.text(plot_ra,plot_dec,plot_name,color='red')
print ('source with deepest exposure: %s, RA = %0.2f, Dec = %0.2f'%(plot_name,plot_ra,plot_dec))
for i in range(len(xbins)-1):
    for j in range(len(ybins)-1):
        if hist[i,j]>30:
            print ("counts %s, RA %s, Dec %s"%(hist[i,j],xbins[i],ybins[j]))
            #ax.text(xbins[i]+3.,ybins[j]+3.,hist[i,j],color="black", ha="center", va="center", fontweight="bold")
ax.axis('on')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
plt.savefig("output_plots/RunRADec_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
hist, xbins, ybins, im = plt.hist2d(List_Used_Gal_l, List_Used_Gal_b, range = [[0,360], [-10,10]], bins=(90, 20), cmap=plt.cm.Greys)
plt.colorbar()
ax.axis('on')
ax.set_xlabel('gal. l')
ax.set_ylabel('gal. b')
plt.savefig("output_plots/RunGalCoord_%s.png"%(plot_tag))
plt.close(fig)

elev_bins = 7
azim_bins = 9
plt.clf()
fig, ax = plt.subplots()
hist, xbins, ybins, im = plt.hist2d(List_Used_Azim, List_Used_Elev, weights=List_Used_Exposure, range = [[0,360], [20,90]], bins=(azim_bins,elev_bins), cmap=plt.cm.Greys)
plt.colorbar()
for i in range(len(xbins)-1):
    for j in range(len(ybins)-1):
        if hist[i,j]>30:
            print ("counts %s, Elev %s, Azim %s"%(hist[i,j],xbins[i],ybins[j]))
ax.axis('on')
ax.set_xlabel('Azim')
ax.set_ylabel('Elev')
plt.savefig("output_plots/RunElevAzim_%s.png"%(plot_tag))
plt.close(fig)

plt.clf()
fig, ax = plt.subplots()
hist, xbins, ybins, im = plt.hist2d(List_Used_Elev, List_Used_NSB, range = [[0,90], [3,8]], bins=(18,10), cmap=plt.cm.Greys)
plt.colorbar()
for i in range(len(xbins)-1):
    for j in range(len(ybins)-1):
        if hist[i,j]>30:
            print ("counts %s, Elev %s, NSB %s"%(hist[i,j],xbins[i],ybins[j]))
ax.axis('on')
ax.set_xlabel('Elev')
ax.set_ylabel('NSB')
plt.savefig("output_plots/RunElevNSB_%s.png"%(plot_tag))
plt.close(fig)

if search_for_on_data:

    for target in range(0,len(target_ra)):

        other_stars = []
        other_star_ra = []
        other_star_dec = []
        other_stars += [target_name[target]]
        other_star_ra += [target_ra[target]]
        other_star_dec += [target_dec[target]]

        inputFile = open('TeVCat_RaDec_w_Names.txt')
        for line in inputFile:
            gamma_source_name = line.split(',')[0]
            gamma_source_ra = float(line.split(',')[1])
            gamma_source_dec = float(line.split(',')[2])
            distance = pow(gamma_source_ra-target_ra[target],2)+pow(gamma_source_dec-target_dec[target],2)
            if distance>2.*2.: continue
            near_a_source = False
            for entry in range(0,len(other_stars)):
                distance = pow(gamma_source_ra-other_star_ra[entry],2)+pow(gamma_source_dec-other_star_dec[entry],2)
                if distance<0.3*0.3:
                    near_a_source = True
            if not near_a_source and not '%' in gamma_source_name:
                other_stars += [gamma_source_name]
                other_star_ra += [gamma_source_ra]
                other_star_dec += [gamma_source_dec]

        other_star_ra_flipped = []
        for entry in range(0,len(other_star_ra)):
            other_star_ra_flipped += [-1.*other_star_ra[entry]]

        List_Proposed_RA = []
        List_Proposed_Dec = []
        List_Proposed_Exposure = []
        # winter
        #List_Proposed_RA += [286.6]
        #List_Proposed_Dec += [6.9]
        #List_Proposed_Exposure += [10.]
        #List_Proposed_RA += [286.6]
        #List_Proposed_Dec += [6.2]
        #List_Proposed_Exposure += [5.]
        ## summer
        #List_Proposed_RA += [286.6]
        #List_Proposed_Dec += [6.9]
        #List_Proposed_Exposure += [5.]
        #List_Proposed_RA += [288.]
        #List_Proposed_Dec += [6.2]
        #List_Proposed_Exposure += [5.]
        #List_Proposed_RA += [287.3]
        #List_Proposed_Dec += [5.5]
        #List_Proposed_Exposure += [0.]
        #List_Proposed_RA += [287.3]
        #List_Proposed_Dec += [6.9]
        #List_Proposed_Exposure += [5.]
        #List_Proposed_RA += [286.6]
        #List_Proposed_Dec += [6.2]
        #List_Proposed_Exposure += [5.]

        List_Used_RA_Flipped = []
        for entry in range(0,len(List_Used_RA)):
            List_Used_RA_Flipped += [-1.*List_Used_RA[entry]]
        List_Proposed_RA_Flipped = []
        for entry in range(0,len(List_Proposed_RA)):
            List_Proposed_RA_Flipped += [-1.*List_Proposed_RA[entry]]

        delta = 0.02
        x = np.arange(-target_ra[target]-3.0, -target_ra[target]+3.0, delta)
        y = np.arange(target_dec[target]-3.0, target_dec[target]+3.0, delta)
        X, Y = np.meshgrid(x, y)
        Z = 0.*gaussian_2d(X, Y, 0., 0., 1., 1.)
        for run in range(0,len(List_Used_RA_Flipped)):
            Z += List_Used_Exposure[run]*gaussian_2d(X, Y, List_Used_RA_Flipped[run], List_Used_Dec[run], 1., 1.)
        for run in range(0,len(List_Proposed_RA_Flipped)):
            Z += List_Proposed_Exposure[run]*gaussian_2d(X, Y, List_Proposed_RA_Flipped[run], List_Proposed_Dec[run], 1., 1.)
        
        # Create a contour plot with labels using default colors.  The
        # inline argument to clabel will control whether the labels are draw
        # over the line segments of the contour, removing the lines beneath
        # the label
        plt.clf()
        fig, ax = plt.subplots()
        CS = plt.contour(X, Y, Z)
        plt.scatter(other_star_ra_flipped, other_star_dec, s=80, c='black', marker="+")
        for star in range(0,len(other_stars)):
            plt.text(other_star_ra_flipped[star],other_star_dec[star]+0.2,other_stars[star])
        mytext = plt.text(other_star_ra_flipped[0],other_star_dec[0]+0.2,other_stars[0],color='r')
        mytext.set_bbox(dict(facecolor='white', alpha=0.7))
        plt.clabel(CS, inline=1, fontsize=10)
        #plt.title('Simplest default with labels')
        plt.savefig("output_plots/ExposureMap_%s_%s.png"%(plot_tag,target+1))
        plt.close(fig)

        plt.clf()
        fig, ax = plt.subplots()
        hist, xbins, ybins, im = plt.hist2d(List_Used_RA_Flipped, List_Used_Dec, weights=List_Used_Exposure, range=[[-target_ra[target]-3.0,-target_ra[target]+3.0],[target_dec[target]-3.0, target_dec[target]+3.0]], bins=(60,60), cmap=plt.cm.Reds)
        plt.colorbar()
        ax.axis('on')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        plt.scatter(other_star_ra_flipped, other_star_dec, s=80, c='black', marker="+")
        plt.scatter(List_Proposed_RA_Flipped, List_Proposed_Dec, s=80, c='red', marker="o")
        circle1 = plt.Circle((-286.65, 7.05), 0.3, color='r', fill=False)
        ax.add_patch(circle1)
        for star in range(0,len(other_stars)):
            plt.text(-1.*other_star_ra[star],other_star_dec[star]+0.2,other_stars[star])
        plt.clabel(CS, inline=1, fontsize=10)
        plt.savefig("output_plots/PointingMap_%s_%s.png"%(plot_tag,target+1))
        plt.close(fig)
        for i in range(len(xbins)-1):
            for j in range(len(ybins)-1):
                if hist[i,j]>5.:
                    print ("counts %s, RA %s, Dec %s"%(hist[i,j],xbins[i],ybins[j]))
