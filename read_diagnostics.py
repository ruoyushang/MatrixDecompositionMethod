
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


def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ASin(sin_l_NCP_m_l))*180./ROOT.TMath.Pi()
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
target_ra = 286.975
target_dec = 6.269
range_ra = 1.0
range_dec = 1.0

V4 = False
V5 = True
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
#sourceFile = open('../data/output_list/WComaeV5_runlist.txt')
sourceFile = open('../data/output_list/MGRO_J1908_V5_runlist.txt')
for line in sourceFile:
    Source_RunNumber += [int(line)]
    Source_PedVar_DC += [0.]
    Source_Elev += [0.]
    Source_Azim += [0.]
    Source_Livetime += [0.]

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

#for entry in range(0,len(List_RunNumber)):
#
#    RunNumber = List_RunNumber[entry]
#    Elev = List_Elev[entry]
#    Azim = List_Azim[entry]
#    T1_RA = List_T1_RA[entry]
#    T1_Dec = List_T1_Dec[entry]
#    PedVar_DC = List_PedVar_DC[entry]
#    PedVar_PE = List_PedVar_PE[entry]
#    FIR_Mean = List_FIR_Mean[entry]
#    FIR_RMS = List_FIR_RMS[entry]
#    L3_rate = List_L3_rate[entry]
#    Livetime = List_Livetime[entry]
#
#    if V4:
#        if int(RunNumber)>=46642: continue
#    if V5:
#        if int(RunNumber)<46642: continue
#        if int(RunNumber)>=63373: continue
#    if V6:
#        if int(RunNumber)<63373: continue
#
#    if L3_rate<150.: continue
#    if L3_rate>450.: continue
#
#    if abs(float(T1_RA)-target_ra)>range_ra: continue
#    if abs(float(T1_Dec)-target_dec)>range_dec: continue
#
#    print 'RunNumber %s, L3_rate %s, Livetime %s'%(RunNumber,L3_rate,Livetime)
#
#    List_Used += [RunNumber]

List_Produced = []
inputFile = open('temp_output.txt')
for line in inputFile:
    if not 'root' in line: continue
    line_segments = line.split('.')
    RunNumber = line_segments[len(line_segments)-2]
    List_Produced += [RunNumber]

n_matches = 5
for entry2 in range(0,len(Source_RunNumber)):
    found_matches = 0
    entry_in_list = 0
    for entry in range(0,len(List_RunNumber)):
        RunNumber = List_RunNumber[entry]
        if int(Source_RunNumber[entry2])==int(RunNumber):
            entry_in_list = entry
    for entry in range(0,len(List_RunNumber)):
        if found_matches>=n_matches: continue
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
            if L3_rate>450.: continue
            if Livetime<15.: continue
            if V4:
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
            elif not (T2_RA==0. and T2_Dec==0.):
                gal_l, gal_b = ConvertRaDecToGalactic(T2_RA,T2_Dec)
                if abs(gal_b)<10.: continue
            elif not (T3_RA==0. and T3_Dec==0.):
                gal_l, gal_b = ConvertRaDecToGalactic(T3_RA,T3_Dec)
                if abs(gal_b)<10.: continue
            elif not (T4_RA==0. and T4_Dec==0.):
                gal_l, gal_b = ConvertRaDecToGalactic(T4_RA,T4_Dec)
                if abs(gal_b)<10.: continue
            else: continue # there is no pointing info
            #if abs(Livetime-Source_Livetime[entry2])/Source_Livetime[entry2]>0.5: continue
            if abs(Elev-Source_Elev[entry2])>5.: continue
            if abs(PedVar_DC-Source_PedVar_DC[entry2])>0.5: continue
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
