
import csv
import numpy as np
import math
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
from scipy.optimize import curve_fit

from bs4 import BeautifulSoup

def powerlaw_func(x,A,r):
    return A+r*x

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

def ReadHAWCTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_flux = []
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
        if 'flux:' in line:
            target_flux = line.lstrip('          flux: ')
        if 'index systematic uncertainty down:' in line:
            source_name += [target_name]
            source_ra += [float(target_ra)]
            source_dec += [float(target_dec)]
            source_flux += [float(target_flux)/2.34204e-13]
            target_name = ''
            target_ra = ''
            target_dec = ''
            target_flux = ''
    return source_name, source_ra, source_dec, source_flux

def ReadATNFTargetListFromFile(file_path):
    source_name = []
    source_ra = []
    source_dec = []
    source_dist = []
    source_age = []
    source_edot = []
    source_flux = []
    inputFile = open(file_path)
    for line in inputFile:
        if line[0]=="#": continue
        target_name = line.split(',')[0].strip(" ")
        if target_name=="\n": continue
        target_ra = line.split(',')[1].strip(" ")
        target_dec = line.split(',')[2].strip(" ")
        target_dist = line.split(',')[3].strip(" ")
        target_age = line.split(',')[4].strip(" ")
        target_edot = line.split(',')[5].strip(" ")
        if target_dist=='*': continue
        if target_age=='*': continue
        if target_edot=='*': continue
        target_brightness = float(target_edot)/pow(float(target_dist),2)

        source_name += [target_name]
        source_ra += [float(HMS2deg(target_ra,target_dec)[0])]
        source_dec += [float(HMS2deg(target_ra,target_dec)[1])]
        source_dist += [float(target_dist)]
        source_age += [float(target_age)]
        source_edot += [float(target_edot)]
        source_flux += [float(target_edot)/pow(float(target_dist),2)]
    return source_name, source_ra, source_dec, source_dist, source_age, source_edot, source_flux

def FindFermiSource(psr_name,psr_ra,psr_dec,fermi_name,fermi_ra,fermi_dec,search_radis):

    found_fermi_name = ''
    found_fermi_idx = 0
    for fermi in range(0,len(fermi_name)):
        delta_ra = psr_ra - fermi_ra[fermi]
        delta_dec = psr_dec - fermi_dec[fermi]
        distance = pow(delta_ra*delta_ra+delta_dec*delta_dec,0.5)
        #if distance<0.5:
        if distance<search_radis:
            found_fermi_name = fermi_name[fermi]
            found_fermi_idx = fermi
    return found_fermi_name, found_fermi_idx


target_psr_name, target_psr_ra, target_psr_dec, target_psr_dist, target_psr_age, target_psr_edot, target_psr_flux = ReadATNFTargetListFromFile('ATNF_pulsar_full_list.txt')
target_hwc_name, target_hwc_ra, target_hwc_dec, target_hwc_flux = ReadHAWCTargetListFromFile('Cat_3HWC.txt')

source_age = []
source_bfield = []
source_flux = []
fit_source_age = []
fit_source_bfield = []
for psr in range(0,len(target_psr_name)):

    found_hwc_name  , hwc_idx = FindFermiSource(target_psr_name[psr],target_psr_ra[psr],target_psr_dec[psr],target_hwc_name,target_hwc_ra,target_hwc_dec,0.5)

    crab_edot = 4.5e+38 
    crab_dist = 2.000 # kpc
    crab_pwn_bfield = 125.0 # muG

    if not found_hwc_name=='':
        est_bfield = pow(target_psr_edot[psr]/crab_edot*pow(crab_pwn_bfield,2)*pow(crab_dist/target_psr_dist[psr],2)/target_hwc_flux[hwc_idx],0.5)
        source_age += [math.log10(target_psr_age[psr]/1000.)]
        source_bfield += [math.log10(est_bfield)]
        source_flux += [math.log10(target_hwc_flux[hwc_idx])]
        if target_psr_age[psr]/1000.<1e4:
            fit_source_age += [math.log10(target_psr_age[psr]/1000.)]
            fit_source_bfield += [math.log10(est_bfield)]

source_age = np.array(source_age)
source_bfield = np.array(source_bfield)
source_flux = np.array(source_flux)
fit_source_age = np.array(fit_source_age)
fit_source_bfield = np.array(fit_source_bfield)
start = (1., -0.5)
popt, pcov = curve_fit(powerlaw_func,fit_source_age,fit_source_bfield,p0=start,bounds=((-5., -10.), (5., 10.)))

log_age = np.linspace(0.,5.,50)
func_fit = powerlaw_func(log_age, *popt)
print ('Amp = %0.2e'%(popt[0]))
print ('Idx = %0.2e'%(popt[1]))

psr_age = 100. # kyr
halo_bfield = pow(10.,powerlaw_func(math.log10(psr_age), 2.21e+00, -8.51e-01)) # muG
print ('psr_age = %0.1e kry, halo_bfield = %0.1e muG'%(psr_age,halo_bfield))

fig, ax = plt.subplots()
figsize_x = 10
figsize_y = 6
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)
fig.clf()
axbig = fig.add_subplot()
axbig.scatter(source_age,source_bfield)
axbig.plot(log_age,func_fit)
axbig.legend(loc='best')
axbig.set_xlabel('Log10 age [kyr]')
axbig.set_ylabel('Log10 magnetic field [$\mu$G]')
plotname = 'MagneticFieldVsAge'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()
