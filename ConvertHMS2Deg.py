
import sys,ROOT

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

def ConvertRaDecToGalactic(ra, dec):
    delta = dec*ROOT.TMath.Pi()/180.
    delta_G = 27.12825*ROOT.TMath.Pi()/180.
    alpha = ra*ROOT.TMath.Pi()/180.
    alpha_G = 192.85948*ROOT.TMath.Pi()/180.
    l_NCP = 122.93192*ROOT.TMath.Pi()/180.
    sin_b = ROOT.TMath.Sin(delta)*ROOT.TMath.Sin(delta_G)+ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(delta_G)*ROOT.TMath.Cos(alpha-alpha_G)
    cos_b = ROOT.TMath.Cos(ROOT.TMath.ASin(sin_b))
    sin_l_NCP_m_l = ROOT.TMath.Cos(delta)*ROOT.TMath.Sin(alpha-alpha_G)/cos_b
    cos_l_NCP_m_l = (ROOT.TMath.Cos(delta_G)*ROOT.TMath.Sin(delta)-ROOT.TMath.Sin(delta_G)*ROOT.TMath.Cos(delta)*ROOT.TMath.Cos(alpha-alpha_G))/cos_b
    b = (ROOT.TMath.ASin(sin_b))*180./ROOT.TMath.Pi()
    l = (l_NCP-ROOT.TMath.ATan2(sin_l_NCP_m_l,cos_l_NCP_m_l))*180./ROOT.TMath.Pi()
    return l, b

source = []
sky_coord = []

source += ['PSR J0007+7303']
sky_coord += [['00 07 01.7','+73 03 07.4']]
source += ['PSR J0023+0923']
sky_coord += [['00 23 16.877498','+09 23 23.8604']]

for s in range(0,len(source)):
    print ('if (source_name=="%s")'%(source[s]))
    print ('{')
    print ('    Source_RA = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0])))
    print ('    Source_Dec = %0.3f;'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1])))
    print ('}')
    gal_l, gal_b = ConvertRaDecToGalactic(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]),float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1]))
    print ('gal_l = %0.1f, gal_b = %0.1f'%(gal_l,gal_b))
for s in range(0,len(source)):
    print ('other_stars += ["%s"]'%(source[s]))
    print ('other_star_coord += [[%s,%s]]'%(float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[0]),float(HMS2deg(sky_coord[s][0],sky_coord[s][1])[1])))

