import sys,ROOT
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import astropy.utils as utils
from astropy.nddata import Cutout2D
import numpy as np

# https://docs.astropy.org/en/stable/visualization/wcsaxes/slicing_datacubes.html

def read_pwn_model_file(file_path,center_x,center_y,size_x,size_y):

    pwn_gal_l = []
    pwn_gal_b = []
    pwn_sigma = []
    pwn_file = open(file_path, 'r')
    n_pwn = 0
    pwn_x = 0.
    pwn_y = 0.
    pwn_r = 0.
    for line in pwn_file:
        if 'RadialGaussian' in line:
            n_pwn += 1
            pwn_x = 0.
            pwn_y = 0.
            pwn_r = 0.
        if 'GLON' in line:
            line_blocks = line.split(' ')
            for block in line_blocks:
                if 'value' in block:
                    coord = block.split('"')[1]
                    pwn_x = float(coord)
        if 'GLAT' in line:
            line_blocks = line.split(' ')
            for block in line_blocks:
                if 'value' in block:
                    coord = block.split('"')[1]
                    pwn_y = float(coord)
        if 'Sigma' in line:
            line_blocks = line.split(' ')
            for block in line_blocks:
                if 'value' in block:
                    coord = block.split('"')[1]
                    pwn_r = float(coord)
            delta_gal_l = abs(pwn_x-center_x)
            if delta_gal_l>180.: 
                delta_gal_l = abs(360.-delta_gal_l)
            delta_gal_b = abs(pwn_y-center_y)
            if delta_gal_l<size_x and delta_gal_b<size_y:
                pwn_gal_l += [pwn_x]
                pwn_gal_b += [pwn_y]
                pwn_sigma += [pwn_r]
    return pwn_gal_l, pwn_gal_b, pwn_sigma


fig, ax = plt.subplots()
figsize_x = 10
figsize_y = 5
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

#layout = 'm4'
#offset = '5deg'
#folder_name = 'CTAsim_%s_80hr_%s_30pc'%(layout,offset)
#bkg_folder_name = 'CTAsim_%s_80hr_%s_0pc'%(layout,offset)

#layout = 'f4'
layout = 'm4'
exposure_hour = '0'
#exposure_hour = '1'
#exposure_hour = '2'
#exposure_hour = '4'
folder_name = 'CTAsim_gps_%s_%shr_pwn'%(layout,exposure_hour)
bkg_folder_name = 'CTAsim_gps_%s_%shr_bkg'%(layout,exposure_hour)

plotname = 'Plot2D_%s'%(folder_name)
filename = 'output_tehanu/%s/cntcube_emin=0.20_emax=10.00_side=45.00.fits'%(folder_name)
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
image_data = hdu.data

bkg_plotname = 'Plot2D_%s'%(bkg_folder_name)
bkg_filename = 'output_tehanu/%s/cntcube_emin=0.20_emax=10.00_side=45.00.fits'%(bkg_folder_name)
bkg_hdu = fits.open(bkg_filename)[0]
bkg_wcs = WCS(bkg_hdu.header)
bkg_image_data = bkg_hdu.data

# This is a three-dimensional dataset which you can check by looking at the header information by:
print (hdu.header)
print (wcs)
#Number of WCS axes: 3
#CTYPE : 'RA---CAR'  'DEC--CAR'  ''
#CRVAL : 83.633331299  22.014444351  0.0
#CRPIX : 100.5  100.5  0.0
#NAXIS : 200  200  25


emin = 0
emax = 50
#center_sky_x = 266.0367
#center_sky_y = -26.9782
center_sky_x = 0.
center_sky_y = 0.
delta_sky_x = 45./2.
delta_sky_y = 6./2.

print (wcs.all_world2pix(center_sky_x,center_sky_y,0.0,1))
#
pixs_start = wcs.all_world2pix(center_sky_x+delta_sky_x,center_sky_y-delta_sky_y,emin,1)
pixs_end = wcs.all_world2pix(center_sky_x-delta_sky_x,center_sky_y+delta_sky_y,emax,1)
print (pixs_start)
print (pixs_end)
start_pix_x = int(pixs_start[0])
start_pix_y = int(pixs_start[1])
start_pix_e = emin
end_pix_x = int(pixs_end[0])
end_pix_y = int(pixs_end[1])
end_pix_e = emax
center_pix_x = int((pixs_start[0]+pixs_end[0])/2.)
length_pix_x = int((pixs_start[0]-pixs_end[0])/2.)
center_pix_y = int((pixs_start[1]+pixs_end[1])/2.)
length_pix_y = int((pixs_start[1]-pixs_end[1])/2.)
#
print ((image_data[start_pix_e, :, :].shape))
print ('start_pix_e = %s'%(start_pix_e))
print ('end_pix_e = %s'%(end_pix_e))
image_data_reduced = np.full((image_data[start_pix_e, :, :].shape),0.)
for pix in range(start_pix_e,end_pix_e):
    image_data_reduced += image_data[pix, :, :]
bkg_image_data_reduced = np.full((bkg_image_data[start_pix_e, :, :].shape),0.)
for pix in range(start_pix_e,end_pix_e):
    bkg_image_data_reduced += bkg_image_data[pix, :, :]
position = (center_pix_x,center_pix_y)
size = (length_pix_x,length_pix_y)
print (wcs.dropaxis(2))
#image_data_cutout = Cutout2D(image_data_reduced, position=position, size=size, wcs=wcs.dropaxis(2))
#ax = plt.subplot(projection=image_data_cutout.wcs)
#ax.imshow(image_data_cutout.data[:,:])
ax = plt.subplot(projection=wcs.dropaxis(2))
ax.imshow(image_data_reduced[:,:])
plt.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')

#
## If we wanted to plot the spectral axes for one pixel we can do this by slicing down to one dimension.
#ax = plt.subplot(projection=wcs, slices=('x', 277, 57))
#print (image_data[gal_y_idx_start, gal_x_idx_start, :].shape)
#image_data_reduced = np.full((image_data[gal_y_idx_start, gal_x_idx_start, :].shape),0.)
#print ('gal_x_idx_start = %s'%(gal_x_idx_start))
#print ('gal_x_idx_end = %s'%(gal_x_idx_end))
#for x_idx in range(gal_x_idx_start,gal_x_idx_end):
#    for y_idx in range(gal_y_idx_start,gal_y_idx_end):
#        image_data_reduced += image_data[y_idx, x_idx, :]
#ax.plot(image_data_reduced[:])
#ax.set_xlim([vel_idx_start, vel_idx_end])
#ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
## As this is still a WCSAxes plot, we can set the display units for the x-axis
##vel, gal_lon, gal_lat = ax.coords
##vel.set_format_unit(u.m/u.s)
#plotname = 'Plot1D'
#plt.savefig("output/%s.png"%(plotname),bbox_inches='tight')

nbins_x = end_pix_x - start_pix_x
nbins_y = end_pix_y - start_pix_y
hist_skymap_2d = ROOT.TH2D("hist_skymap_2d","",nbins_x,center_sky_x-delta_sky_x,center_sky_x+delta_sky_x,nbins_y,center_sky_y-delta_sky_y,center_sky_y+delta_sky_y)
hist_skymap_2d_bkg = ROOT.TH2D("hist_skymap_2d_bkg","",nbins_x,center_sky_x-delta_sky_x,center_sky_x+delta_sky_x,nbins_y,center_sky_y-delta_sky_y,center_sky_y+delta_sky_y)
hist_skymap_2d_zscore = ROOT.TH2D("hist_skymap_2d_zscore","",nbins_x,center_sky_x-delta_sky_x,center_sky_x+delta_sky_x,nbins_y,center_sky_y-delta_sky_y,center_sky_y+delta_sky_y)
profile_skymap_1d = []
profile_skymap_1d_bkg = []
axis_skymap_1d = []
for binx in range(0,nbins_x):
    cnt_total = 0.
    cnt_bkg = 0.
    for biny in range(0,nbins_y):
        cnt_total += image_data_reduced[biny,binx]
        cnt_bkg += bkg_image_data_reduced[biny,binx]
        hist_skymap_2d.SetBinContent(binx+1,biny+1,image_data_reduced[biny,binx])
        hist_skymap_2d_bkg.SetBinContent(binx+1,biny+1,bkg_image_data_reduced[biny,binx])
    profile_skymap_1d += [cnt_total]
    profile_skymap_1d_bkg += [cnt_bkg]
    axis_skymap_1d += [hist_skymap_2d.GetXaxis().GetBinCenter(binx+1)]

image_data_zscore = np.full((image_data[start_pix_e, :, :].shape),0.)
for binx in range(0,nbins_x):
    for biny in range(0,nbins_y):
        cnt_total = image_data_reduced[biny,binx]
        cnt_bkg = bkg_image_data_reduced[biny,binx]
        zscore = 0.
        if cnt_bkg>0.:
            zscore = (cnt_total-cnt_bkg)/pow(cnt_bkg,0.5)
        image_data_zscore[biny,binx] = pow(zscore,2)
        hist_skymap_2d_zscore.SetBinContent(binx+1,biny+1,zscore)
ax = plt.subplot(projection=wcs.dropaxis(2))
max_z = 25.
min_z = 0.
im = ax.imshow(image_data_zscore[:,:],vmin=min_z,vmax=max_z)
cbar = fig.colorbar(im,orientation="horizontal")
cbar.set_label('Test Statistics')
plt.savefig("output_plots/%s_zscore.png"%(plotname),bbox_inches='tight')

file_path = 'output_tehanu/model_galactic_pwn.xml'
src_gal_l, src_gal_b, src_sigma = read_pwn_model_file(file_path,center_sky_x,center_sky_y,delta_sky_x,delta_sky_y)
n_detect_src = 0.
for src in range(0,len(src_gal_l)):
    src_data_cnt = 0.
    src_bkgd_cnt = 0.
    for binx in range(0,nbins_x):
        for biny in range(0,nbins_y):
            cell_x = hist_skymap_2d_zscore.GetXaxis().GetBinCenter(binx+1)
            cell_y = hist_skymap_2d_zscore.GetYaxis().GetBinCenter(biny+1)
            delta_x = abs(cell_x-src_gal_l[src])
            if delta_x>180.:
                delta_x = 360.-delta_x
            delta_y = abs(cell_y-src_gal_b[src])
            dist_to_pwn = pow(pow(delta_x,2)+pow(delta_y,2),0.5)
            if dist_to_pwn>2.*max(src_sigma[src],0.1): continue
            src_data_cnt += hist_skymap_2d.GetBinContent(binx+1,biny+1)
            src_bkgd_cnt += hist_skymap_2d_bkg.GetBinContent(binx+1,biny+1)
    src_zscore = 0.
    if src_bkgd_cnt>0.:
        src_zscore = (src_data_cnt-src_bkgd_cnt)/pow(src_bkgd_cnt,0.5)
    if src_zscore>5.:
        n_detect_src += 1.
print ('n candidate = %s'%(len(src_gal_l)))
print ('n_detect_src = %s'%(n_detect_src))

figsize_x = 6.4
figsize_y = 4.8
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

fig.clf()
axbig = fig.add_subplot()
axbig.plot(axis_skymap_1d,profile_skymap_1d,label='signal+bkg')
axbig.plot(axis_skymap_1d,profile_skymap_1d_bkg,label='bkg')
axbig.set_xlabel('Gal. l [deg]')
axbig.set_ylabel('count')
axbig.legend(loc='best')
plotname = 'EmissionProfile_%s'%(folder_name)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

#log_file = 'ctlike_emin=0.20_emax=10.00_side=8.00_Crab_gauss.log'
#input_txt = open('output_tehanu/'+folder_name+'/'+log_file, 'r')
#for line in input_txt:
#    if 'Test Statistic' in line:
#        TS = line.split(' ')[5]
#        print ('TS = %s'%(TS))
#    if 'Sigma' in line:
#        Sigma = line.split(' ')[5]
#        Sigma_err = line.split(' ')[7]
#        print ('Sigma = %s +/- %s'%(Sigma,Sigma_err))
