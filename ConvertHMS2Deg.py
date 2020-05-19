def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)           
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC

source = []
sky_coord = []

source += ['Crab']
sky_coord += [['05 34 31.97','+22 00 52.1']]
source += ['Mrk421']
sky_coord += [['11 04 19','+38 11 41']]
source += ['H1426']
sky_coord += [['14 28 32.609','+42 40 21.05']]
source += ['1ES0229']
sky_coord += [['02 32 53.2','+20 16 21']]
source += ['PKS1424']
sky_coord += [['14 27 00','+23 47 00']]
source += ['3C264']
sky_coord += [['11 45 5.009','+19 36 22.74']]
source += ['OJ287V6']
sky_coord += [['08 54 49.1','+20 05 58.89']]
source += ['RBS0413V6']
sky_coord += [['03 19 47','+18 45 42']]
source += ['PG1553V6']
sky_coord += [['15 55 44.7','+11 11 41']]
source += ['Segue1V6']
sky_coord += [['10 07 04','+16 04 55']]
source += ['ComaV6']
sky_coord += [['12 59 48.7','+27 58 50']]
source += ['1ES1011V6']
sky_coord += [['10 15 4.139','+49 26 0.71']]
source += ['NGC1275V6']
sky_coord += [['03 19 48.1','+41 30 42']]
source += ['1ES0647V6']
sky_coord += [['06 50 46.490','+25 02 59.62']]
source += ['1ES1440V6']
sky_coord += [['14 42 48.277','+12 00 40.37']]
source += ['1ES1741V6']
sky_coord += [['17 44 01.2','+19 32 47']]
source += ['IC443HotSpot']
sky_coord += [['06 16 51','+22 30 11']]
source += ['RGBJ0710']
sky_coord += [['07 10 26.4','+59 09 00']]
source += ['CasA']
sky_coord += [['23 23 13.8','+58 48 26']]
source += ['M82']
sky_coord += [['09 55 52.7','+69 40 46']]
source += ['G079']
sky_coord += [['20 32 28.56','+40 19 41.52']]
source += ['WComaeV6']
sky_coord += [['12 21 31.7','+28 13 59']]
source += ['1ES1218V6']
sky_coord += [['12 21 26.3','+30 11 29']]
source += ['MGRO_J1908_V6']
sky_coord += [['19 07 54','+06 16 07']]
source += ['MGRO_J1908_V5']
sky_coord += [['19 07 54','+06 16 07']]
source += ['Segue1V5']
sky_coord += [['10 07 04','+16 04 55']]
source += ['S3_1227_V6']
sky_coord += [['12 30 14.1','+25 18 07']]
source += ['MS1221V6']
sky_coord += [['12 24 24.2','+24 36 24']]
source += ['PKS1441V6']
sky_coord += [['14 43 56.9','+25 01 44']]
source += ['GemingaV6']
sky_coord += [['06 32 28','+17 22 00']]
source += ['GemingaV5']
sky_coord += [['06 32 28','+17 22 00']]
source += ['SgrAV6']
sky_coord += [['17 45 39.6','-29 00 22']]
source += ['1ES1215']
sky_coord += [['12 17 48.5','+30 06 06']]
source += ['Tycho']
sky_coord += [['00 25 21.6','+64 07 48']]
source += ['CTA1']
sky_coord += [['00 06 26','+72 59 01.0']]
source += ['ARGO J1910+0720']
sky_coord += [['19 10 36','+07 21 00']]
source += ['2HWC_J1953V6']
sky_coord += [['19 53 02.4','+29 28 48']]
source += ['2HWC_J1930V6']
sky_coord += [['19 30 32','+18 52 12']]
source += ['2HWC_J1928V6']
sky_coord += [['19 28 36','+17 46 48']]
source += ['CygnusV6']
sky_coord += [['20 18 35.03','+36 50 00.0']]

for s in range(0,len(source)):
    print 'if (source_name=="%s")'%(source[s])
    print '{'
    print '    Source_RA = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]))
    print '    Source_Dec = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1]))
    print '}'
for s in range(0,len(source)):
    print 'other_stars += ["%s"]'%(source[s])
    print 'other_star_coord += [[%s,%s]]'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]),float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1]))

