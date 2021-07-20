
import os.path
import sys,ROOT
from ROOT import *

inputname = 'runlist_backup/RunList_1ES0229V6.txt'
#inputname = 'runlist_backup/RunList_OJ287V6.txt'
#inputname = 'runlist_backup/RunList_Segue1V6.txt'
#inputname = 'runlist_backup/RunList_1ES1011V6.txt'
#inputname = 'runlist_backup/RunList_H1426V6.txt'
#inputname = 'runlist_backup/RunList_M82V6.txt'
#inputname = 'runlist_backup/RunList_SNR_G150p3Plus04p5_V6.txt'
#inputname = 'runlist_backup/RunList_MGRO_J1908_V6.txt'
#inputname = 'runlist_backup/RunList_BLLacV6.txt'
#inputname = 'runlist_backup/RunList_CasAV6.txt'
#inputname = 'runlist_backup/RunList_3C273V6.txt'
#inputname = 'runlist_backup/RunList_1ES0502V6.txt'
#inputname = 'runlist_backup/RunList_DracoV6.txt'
#inputname = 'runlist_backup/RunList_MGRO_J1908_V5.txt'
#inputname = 'runlist_backup/RunList_Segue1V5.txt'
#inputname = 'runlist_backup/RunList_M82V5.txt'
#inputname = 'runlist_backup/RunList_PG1553V5.txt'
#inputname = 'runlist_backup/RunList_3C273V5.txt'
#inputname = 'runlist_backup/RunList_1ES0502V5.txt'
#inputname = 'runlist_backup/RunList_DracoV5.txt'
#inputname = 'runlist_backup/RunList_CrabV5.txt'
#inputname = 'runlist_backup/RunList_CrabV6.txt'
#inputname = 'runlist_backup/RunList_LHAASO_J2108_V6.txt'
#inputname = 'runlist_backup/RunList_BLLacV5.txt'
#inputname = 'runlist_backup/RunList_1ES0229V5.txt'
#inputname = 'runlist_backup/RunList_1ES0414V5.txt'
#inputname = 'runlist_backup/RunList_RGBJ0710V5.txt'
#inputname = 'runlist_backup/RunList_PKS1424V5.txt'
#inputname = 'runlist_backup/RunList_PG1553V6.txt'
#inputname = 'runlist_backup/RunList_RGBJ0710V5.txt'
#inputname = 'runlist_backup/RunList_PG1553V6.txt'

inputList = open(inputname)
broken_files = []
for line_list in inputList:
    runnumber = line_list.split(',')[1]
    runnumber = runnumber.strip('\n')
    #inputLog = open('/scratch/rshang/analysis/Results/v483/%s.anasum.log'%(runnumber))
    #for line_log in inputLog:
    #    if '*** Break *** segmentation violation' in line_log:
    #        broken_files += [runnumber]
    #        break
    #    if 'terminate called after throwing an instance' in line_log:
    #        broken_files += [runnumber]
    #        break
    #    if 'fatal error: file not found' in line_log:
    #        broken_files += [runnumber]
    #        break
    #    if 'Error opening file with effective areas' in line_log:
    #        broken_files += [runnumber]
    #        break
    file_path = '/veritas/userspace/rshang/analysis/Results/v483/%s.anasum.root'%(runnumber)
    if os.path.exists(file_path):
        InputFile = ROOT.TFile(file_path)
        if InputFile.IsZombie():
            print '%s, something very wrong, cannot use this file'%(runnumber)
        elif InputFile.TestBit(ROOT.TFile.kRecovered):
            print '%s, the Recover procedure has been run when opening the file'%(runnumber)
        else:
            print '%s, looks okay.'%(runnumber)
        InfoTree = InputFile.Get("run_"+runnumber+"/stereo/data_on")
        print 'run %s, InfoTree.GetEntries() = %s'%(runnumber,InfoTree.GetEntries())
        for entry in range(0,InfoTree.GetEntries()):
            InfoTree.GetEntry(entry)
            try:
                energy = InfoTree.ErecS
                #print 'run %s,energy = %s'%(runnumber,InfoTree.ErecS)
            except AttributeError:
                print 'corrupted run %s'%(runnumber)
        InputFile.Close()
    else:
        print '%s, does not exist.'%(runnumber)

for entry in range(0,len(broken_files)):
    print broken_files[entry]
