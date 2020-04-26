
## galactic center
#target_ra = 266.415
#target_dec = -29.006
#range_ra = 2.0
#range_dec = 2.0

# W Comae
target_ra = 185.382083333
target_dec = 28.2330555556
range_ra = 1.0
range_dec = 1.0

RunNumber = 0
Elev = 0
Azim = 0
PedVar_DC = 0
PedVar_PE = 0
T1_RA = 0
T1_Dec = 0

V4 = False
V5 = True
V6 = False

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
    Elev = line.split(' ')[8-1]
    Azim = line.split(' ')[9-1]
    if V1 or V2:
        T1_RA = line.split(' ')[46-1]
        T1_Dec = line.split(' ')[47-1]
    if V2:
        PedVar_DC = line.split(' ')[104-1]
        PedVar_PE = line.split(' ')[105-1]

    if V4:
        if int(RunNumber)>=46642: continue
    if V5:
        if int(RunNumber)<46642: continue
        if int(RunNumber)>=63373: continue
    if V6:
        if int(RunNumber)<63373: continue

    if abs(float(T1_RA)-target_ra)>range_ra: continue
    if abs(float(T1_Dec)-target_dec)>range_dec: continue

    #if abs(float(T1_RA)-target_ra)<2.*range_ra and abs(float(T1_Dec)-target_dec)<2.*range_dec: continue
    #if float(Elev)>30.: continue
    #if float(Elev)<20.: continue

    print 'RunNumber %s, Elev %s, Azim %s, T1 RA %s, T1 Dec %s'%(RunNumber,Elev,Azim,T1_RA,T1_Dec)
    #print RunNumber
