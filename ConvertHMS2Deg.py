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
source += ['MGRO J2031+41']
sky_coord += [['20 28 43.2','+41 18 36']]
source += ['PSR J2032+4127']
sky_coord += [['20 32 10','+41 27 34']]
source += ['1ES 0806+524']
sky_coord += [['08 09 59','+52 19 00']]
source += ['S5 0716+714']
sky_coord += [['07 21 53.4','+71 20 36']]
source += ['Markarian 180']
sky_coord += [['11 36 26.4','+70 09 27']]
source += ['1ES 1727+502']
sky_coord += [['17 28 18.6','+50 13 10']]
source += ['3C66A']
sky_coord += [['02 22 41.6','+43 02 35.5']]
source += ['VER J0521+211']
sky_coord += [['05 21 45','+21 12 51.4']]
source += ['1ES 0502+675']
sky_coord += [['05 07 56.2','+67 37 24']]
source += ['RX J0648.7+1516']
sky_coord += [['06 48 45.6','+15 16 12']]
source += ['1ES 1727+502']
sky_coord += [['17 28 18.6','+50 13 10']]
source += ['Boomerang']
sky_coord += [['22 28 44','+61 10 00']]
source += ['NGC1275']
sky_coord += [['03 19 48.1','+41 30 42']]
source += ['TeV J2227+608']
sky_coord += [['22 27 59','+60 52 37']]
source += ['CTA1']
sky_coord += [['00 06 26','+72 59 01.0']]
source += ['BL Lac']
sky_coord += [['22 02 43.3','+42 16 40']]
source += ['VER J2019+407']
sky_coord += [['20 20 04.8','+40 45 26']]
source += ['Perseus cluster']
sky_coord += [['03 19 47.2','+41 30 47']]
source += ['SS 433']
sky_coord += [['19 11 49.56','+04 58 57.8']]
source += ['2FGL J0423.3+5612']
sky_coord += [['04 23 27.0','+56 12 24']]
source += ['1FHL J0432.2+5555']
sky_coord += [['04 32 14.3','+55 55 41.16']]
source += ['2FHL J0431.2+5553e']
sky_coord += [['04 31 16.79','+55 53 24.0']]
source += ['XTE J0421+560']
sky_coord += [['04 19 42.1353986687','+55 59 57.708146070']]
source += ['2FHL J0431.2+5342']
sky_coord += [['04 31 17.1','+53 42 15.8']]
source += ['4FGL J0427.2+5533e']
sky_coord += [['04 27 15.6','+55 33 14']]
source += ['VER J2019+368']
sky_coord += [['20 19 25','+36 48 14']]
source += ['HESS J1857+026']
sky_coord += [['18 57 11','+02 40 00']]
source += ['HESS J1858+020']
sky_coord += [['18 58 20','+02 05 24']]
source += ['Tycho']
sky_coord += [['00 25 21.6','+64 07 48']]
source += ['G40.5-0.5']
sky_coord += [['19 07 11.9','+06 15 35']]
source += ['HESS J1844-030']
sky_coord += [['18 44 41.22','-03 05 34.6']]
source += ['HESS J1912+101']
sky_coord += [['19 12 49','+10 09 06']]
source += ['HESS J1852-000']
sky_coord += [['18 52 14','00 05 56']]
source += ['PSR J1907+0602']
sky_coord += [['19 07 54','+06 02 16']]
source += ['1ES1215']
sky_coord += [['12 17 48.5','+30 06 06']]
source += ['PSR J0030+0451']
sky_coord += [['00 30 27.4277','+04 51 39.710']]
source += ['PSR B0114+58']
sky_coord += [['01 17 38.661','+59 14 38.39']]
source += ['M 31']
sky_coord += [['00 42 44.330','+41 16 07.50']]
source += ['Virgo Cluster (M 87)']
sky_coord += [['12 30 49.42338230','+12 23 28.0438581']]
source += ['PSR J1930+1852']
sky_coord += [['19 30 30.13','+18 52 14.1']]
source += ['PSR J2021+3651']
sky_coord += [['20 21 05.40','+36 51 04.5']]
source += ['PSR J2229+6114']
sky_coord += [['22 29 05.262','+61 14 08.48']]
source += ['PSR J1826-1334']
sky_coord += [['18 26 13.06','-13 34 48.1']]
source += ['PSR J1856+0245']
sky_coord += [['18 56 50.80','+02 45 50.2']]
source += ['PSR J2021+4026']
sky_coord += [['20 21 30.496','+40 26 53.5']]
source += ['SNR G150.3+4.5']
sky_coord += [['04 31 16.8','+55 53 24']]
source += ['Draco dSph']
sky_coord += [['17 20 14.335','+57 55 16.39']]
source += ['3C 273']
sky_coord += [['12 29 06.6996828061','+02 03 08.598846466']]
source += ['1ES 0502+675']
sky_coord += [['05 07 56.1461168342','+67 37 24.302255568']]

for s in range(0,len(source)):
    print 'if (source_name=="%s")'%(source[s])
    print '{'
    print '    Source_RA = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]))
    print '    Source_Dec = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1]))
    print '}'
for s in range(0,len(source)):
    print 'other_stars += ["%s"]'%(source[s])
    print 'other_star_coord += [[%s,%s]]'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]),float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1]))

