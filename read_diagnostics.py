
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


## galactic center
#target_ra = 266.415
#target_dec = -29.006
#range_ra = 2.0
#range_dec = 2.0

# W Comae
#target_ra = 185.382083333
#target_dec = 28.2330555556
#range_ra = 1.0
#range_dec = 1.0

## MGRO J1908
#target_ra = 286.975
#target_dec = 6.269
#range_ra = 3.0
#range_dec = 3.0

# IC 443
#target_ra = 94.511
#target_dec = 22.660
#range_ra = 1.0
#range_dec = 1.0

# Segue 1
#target_ra = 151.917
#target_dec = 16.767
#range_ra = 1.0
#range_dec = 1.0

# OJ 287
#target_ra = 133.704583333
#target_dec = 20.0996916667
#range_ra = 1.0
#range_dec = 1.0

## 2HWC J1928
#target_ra = 292.15
#target_dec = 17.78
#range_ra = 2.0
#range_dec = 2.0

## 2HWC J1953
#target_ra = 298.26
#target_dec = 29.48
#range_ra = 2.0
#range_dec = 2.0

## Cygnus
#target_ra = 304.645958333
#target_dec = 36.8333333333
#range_ra = 1.0
#range_dec = 1.0

# Crab
#target_ra = 83.6332083333
#target_dec = 22.0144722222
#range_ra = 1.0
#range_dec = 1.0

# Geminga
#target_ra = 98.1166666667
#target_dec = 17.3666666667
#range_ra = 2.0
#range_dec = 2.0

# MGRO J2031+41
target_ra = 307.18
target_dec = 41.31
range_ra = 4.0
range_dec = 4.0

# 1ES 0229
#target_ra = 38.2216666667
#target_dec = 20.2725
#range_ra = 1.0
#range_dec = 1.0

# H 1426+428
#target_ra = 217.135870833
#target_dec = 42.6725138889
#range_ra = 1.0
#range_dec = 1.0

# PKS 1424
#target_ra = 216.75
#target_dec = 23.7833333333
#range_ra = 1.0
#range_dec = 1.0

# 3C 264
#target_ra = 176.270870833
#target_dec = 19.60631666673
#range_ra = 1.0
#range_dec = 1.0

# PG 1553
#target_ra = 238.93625
#target_dec = 11.1947222222
#range_ra = 1.0
#range_dec = 1.0

# 1ES 1011
#target_ra = 153.767245833
#target_dec = 49.4335305556
#range_ra = 1.0
#range_dec = 1.0

# RBS 0413
#target_ra = 49.9458333333
#target_dec = 18.7616666667
#range_ra = 1.0
#range_dec = 1.0

# 1ES 0647
#target_ra = 102.693708333
#target_dec = 25.0498944444
#range_ra = 1.0
#range_dec = 1.0

# 1ES 0806+524
#target_ra = 122.495833333
#target_dec = 52.3166666667
#range_ra = 1.0
#range_dec = 1.0

# S5 0716+714
#target_ra = 110.4725
#target_dec = 71.3433333333
#range_ra = 1.0
#range_dec = 1.0

# Markarian 180
#target_ra = 174.11
#target_dec = 70.1575
#range_ra = 1.0
#range_dec = 1.0

# 1ES 1727+502
#target_ra = 262.0775
#target_dec = 50.2194444444
#range_ra = 1.0
#range_dec = 1.0

# Mrk 421
#target_ra = 166.079166667
#target_dec = 38.1947222222
#range_ra = 1.0
#range_dec = 1.0

# 3C 66A
#target_ra = 35.6733333333
#target_dec = 43.0431944444
#range_ra = 1.0
#range_dec = 1.0

# RGB J0710  # deep in V5
#target_ra = 107.61
#target_dec = 59.15
#range_ra = 1.0
#range_dec = 1.0


# VER J0521+211
#target_ra = 80.4375
#target_dec = 21.2142777778
#range_ra = 1.0
#range_dec = 1.0

# 1ES 0502+675
#target_ra = 76.9841666667
#target_dec = 67.6233333333
#range_ra = 1.0
#range_dec = 1.0

# RX J0648.7+1516
#target_ra = 102.19
#target_dec = 15.27
#range_ra = 1.0
#range_dec = 1.0

# 1ES 1727+502
#target_ra = 262.0775
#target_dec = 50.2194444444
#range_ra = 1.0
#range_dec = 1.0

# Coma
#target_ra = 194.952916667
#target_dec = 27.9805555556
#range_ra = 4.0
#range_dec = 4.0

# M 82
#target_ra = 148.969583333
#target_dec = 69.6794444444
#range_ra = 1.0
#range_dec = 1.0

# Boomerang
#target_ra = 337.183333333
#target_dec = 61.1666666667
#range_ra = 2.0
#range_dec = 2.0

# NGC 1275
#target_ra = 49.9504166667
#target_dec = 41.5116666667
#range_ra = 2.0
#range_dec = 2.0

# CTA1
#target_ra = 1.60833333333
#target_dec = 72.9836111111
#range_ra = 2.0
#range_dec = 2.0

# BL Lac
#target_ra = 330.680416667
#target_dec = 42.2777777778
#range_ra = 2.0
#range_dec = 2.0

# VER J2019+407
#target_ra = 305.02
#target_dec = 40.7572222222
#range_ra = 4.0
#range_dec = 4.0

# Perseus cluster
#target_ra = 49.9466666667
#target_dec = 41.5130555556
#range_ra = 2.0
#range_dec = 2.0

# SNR G150.3+4.5
#target_ra = 66.5
#target_dec = 55.0
#range_ra = 4.0
#range_dec = 4.0

search_for_on_data = True

V4 = True
V5 = False
V6 = False


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

List_RunNumber = []
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
#sourceFile = open('../data/output_list/MGRO_J2031_V5_runlist.txt')
#for line in sourceFile:
#    Source_RunNumber += [int(line)]
#    Source_PedVar_DC += [0.]
#    Source_Elev += [0.]
#    Source_Azim += [0.]
#    Source_Livetime += [0.]

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
    List_RunNumber += [int(RunNumber)]
    Elev = line.split(' ')[8-1]
    Azim = line.split(' ')[9-1]
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

if search_for_on_data:

    for entry in range(0,len(List_RunNumber)):
    
        RunNumber = List_RunNumber[entry]
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
    
        if V4:
            if int(RunNumber)>=46642: continue
        if V5:
            if int(RunNumber)<46642: continue
            if int(RunNumber)>=63373: continue
        if V6:
            if int(RunNumber)<63373: continue
    
        if L3_rate<150.: continue
        #if L3_rate>450.: continue
        if Livetime<5.: continue
        #if Livetime<20.: continue
    
        if abs(float(T1_RA)-target_ra)>range_ra: continue
        if abs(float(T1_Dec)-target_dec)>range_dec: continue

        DoNotUse = False
        for entry2 in range(0,len(Source_RunNumber)):
            if int(Source_RunNumber[entry2])==int(RunNumber): DoNotUse = True
        if DoNotUse: continue
    
        print 'RunNumber %s, L3_rate %s, Livetime %s, Elev %s, RA %s, Dec %s'%(RunNumber,L3_rate,Livetime,Elev,T1_RA,T1_Dec)
    
        List_Used += [RunNumber]

else: 

    List_Produced = []
    inputFile = open('temp_output.txt')
    for line in inputFile:
        if not 'root' in line: continue
        line_segments = line.split('.')
        RunNumber = line_segments[len(line_segments)-2]
        List_Produced += [RunNumber]
    
    n_matches = 5
    for nth_sample in range(0,n_matches):
        for entry2 in range(0,len(Source_RunNumber)):
            found_matches = 0
            entry_in_list = 0
            if int(Source_RunNumber[entry2])<=36398: continue
            for entry in range(0,len(List_RunNumber)):
                RunNumber = List_RunNumber[entry]
                if int(Source_RunNumber[entry2])==int(RunNumber):
                    entry_in_list = entry
            for entry in range(0,len(List_RunNumber)):
                if found_matches>=1: continue
                for direction in range(0,2):
                    entry_plus = 0
                    if direction==0: entry_plus = entry_in_list+entry+1
                    else: entry_plus = entry_in_list-entry-1
                    if entry_plus<0: continue
                    if entry_plus>=len(List_RunNumber): continue
                    RunNumber = List_RunNumber[entry_plus]
                    Elev = List_Elev[entry_plus]
                    Azim = List_Azim[entry_plus]
                    T1_RA = List_T1_RA[entry_plus]
                    T1_Dec = List_T1_Dec[entry_plus]
                    T2_RA = List_T2_RA[entry_plus]
                    T2_Dec = List_T2_Dec[entry_plus]
                    T3_RA = List_T3_RA[entry_plus]
                    T3_Dec = List_T3_Dec[entry_plus]
                    T4_RA = List_T4_RA[entry_plus]
                    T4_Dec = List_T4_Dec[entry_plus]
                    PedVar_DC = List_PedVar_DC[entry_plus]
                    PedVar_PE = List_PedVar_PE[entry_plus]
                    L3_rate = List_L3_rate[entry_plus]
                    Livetime = List_Livetime[entry_plus]
                    if L3_rate<150.: continue
                    #if L3_rate>450.: continue
                    if Livetime<15.: continue
                    if V4:
                        if int(RunNumber)<=36398: continue
                        if int(RunNumber)>=46642: continue
                    if V5:
                        if int(RunNumber)<46642: continue
                        if int(RunNumber)>=63373: continue
                    if V6:
                        if int(RunNumber)<63373: continue
                    if abs(float(T1_RA)-target_ra)<10. and abs(float(T1_Dec)-target_dec)<10.: continue
                    if abs(float(T2_RA)-target_ra)<10. and abs(float(T2_Dec)-target_dec)<10.: continue
                    if abs(float(T3_RA)-target_ra)<10. and abs(float(T3_Dec)-target_dec)<10.: continue
                    if abs(float(T4_RA)-target_ra)<10. and abs(float(T4_Dec)-target_dec)<10.: continue
                    if not (T1_RA==0. and T1_Dec==0.):
                        gal_l, gal_b = ConvertRaDecToGalactic(T1_RA,T1_Dec)
                        if abs(gal_b)<10.: continue
                        #if abs(gal_b)<10. and abs(gal_l-180.)>20.: continue
                        if GetDistance(T1_RA,T1_Dec,83.633,22.014)<3.: continue  # Crab
                        if GetDistance(T1_RA,T1_Dec,166.079,38.195)<3.: continue  # Mrk 421
                        if GetDistance(T1_RA,T1_Dec,98.117,17.367)<3.: continue  # Geminga
                    elif not (T2_RA==0. and T2_Dec==0.):
                        gal_l, gal_b = ConvertRaDecToGalactic(T2_RA,T2_Dec)
                        if abs(gal_b)<10.: continue
                        #if abs(gal_b)<10. and abs(gal_l-180.)>20.: continue
                        if GetDistance(T2_RA,T2_Dec,83.633,22.014)<3.: continue  # Crab
                        if GetDistance(T2_RA,T2_Dec,166.079,38.195)<3.: continue  # Mrk 421
                        if GetDistance(T2_RA,T2_Dec,98.117,17.367)<3.: continue  # Geminga
                    elif not (T3_RA==0. and T3_Dec==0.):
                        gal_l, gal_b = ConvertRaDecToGalactic(T3_RA,T3_Dec)
                        if abs(gal_b)<10.: continue
                        #if abs(gal_b)<10. and abs(gal_l-180.)>20.: continue
                        if GetDistance(T3_RA,T3_Dec,83.633,22.014)<3.: continue  # Crab
                        if GetDistance(T3_RA,T3_Dec,166.079,38.195)<3.: continue  # Mrk 421
                        if GetDistance(T3_RA,T3_Dec,98.117,17.367)<3.: continue  # Geminga
                    elif not (T4_RA==0. and T4_Dec==0.):
                        gal_l, gal_b = ConvertRaDecToGalactic(T4_RA,T4_Dec)
                        if abs(gal_b)<10.: continue
                        #if abs(gal_b)<10. and abs(gal_l-180.)>20.: continue
                        if GetDistance(T4_RA,T4_Dec,83.633,22.014)<3.: continue  # Crab
                        if GetDistance(T4_RA,T4_Dec,166.079,38.195)<3.: continue  # Mrk 421
                        if GetDistance(T4_RA,T4_Dec,98.117,17.367)<3.: continue  # Geminga
                    else: continue # there is no pointing info
                    #if abs(Livetime-Source_Livetime[entry2])/Source_Livetime[entry2]>0.5: continue
                    if abs(Elev-Source_Elev[entry2])>5.: continue
                    if abs(PedVar_DC-Source_PedVar_DC[entry2])>1.0: continue
                    already_used = False
                    for entry3 in range(0,len(List_Produced)):
                        if int(List_Produced[entry3])==int(RunNumber): 
                            #print 'RunNumber %s is in the List_Produced'%(RunNumber)
                            already_used = True
                    for entry3 in range(0,len(List_Used)):
                        if int(List_Used[entry3])==int(RunNumber): 
                            #print 'RunNumber %s is in the List_Produced'%(RunNumber)
                            already_used = True
                    if already_used: continue
                    List_Used += [RunNumber]
                    found_matches += 1
                    print 'Size of the list: %s / %s'%(len(List_Used),n_matches*len(Source_RunNumber))

        
for entry in range(0,len(List_Used)):
    print List_Used[entry]
print 'total %s runs.'%(len(List_Used))
