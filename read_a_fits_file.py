#import sys,ROOT
import sys
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import astropy.utils as utils
from astropy.nddata import Cutout2D
import numpy as np
from operator import itemgetter, attrgetter
from scipy.optimize import curve_fit

# https://docs.astropy.org/en/stable/visualization/wcsaxes/slicing_datacubes.html

def gauss_func(x,A,d):
    return A*np.exp(-x*x/(2.*d*d))
    #return A*1.22/(pow(3.14,1.5)*d*(x+0.06*d))*np.exp(-x*x/(d*d))

# Our function to fit is going to be a sum of two-dimensional Gaussians
def gaussian(x, y, x0, y0, sigma, A):
    #return A * np.exp( -((x-x0)/(2.*sigma))**2 -((y-y0)/(2.*sigma))**2)
    return A * np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma*sigma))/(2*np.pi*sigma*sigma)
# This is the callable that is passed to curve_fit. M is a (2,N) array
# where N is the total number of data points in Z, which will be ravelled
# to one dimension.
def _gaussian(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(len(args)//4):
       arr += gaussian(x, y, *args[i*4:i*4+4])
    return arr

def read_pwn_model_file(file_path,center_x,center_y,size_x,size_y):

    pwn_gal_l = []
    pwn_gal_b = []
    pwn_sigma = []
    pwn_norm = []
    pwn_file = open(file_path, 'r')
    pwn_x = 0.
    pwn_y = 0.
    pwn_r = 0.
    pwn_a = 0.
    for line in pwn_file:
        if '<source' in line and 'ExtendedSource' in line:
            pwn_x = 0.
            pwn_y = 0.
            pwn_r = 0.
            pwn_a = 0.
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
        if 'Prefactor' in line:
            line_blocks = line.split(' ')
            for block in line_blocks:
                if 'value' in block:
                    coord = block.split('"')[1]
                    pwn_a = float(coord)
        if '</source>' in line:
            delta_gal_l = abs(pwn_x-center_x)
            if delta_gal_l>180.: 
                delta_gal_l = abs(360.-delta_gal_l)
            delta_gal_b = abs(pwn_y-center_y)
            if delta_gal_l<size_x and delta_gal_b<size_y:
                pwn_gal_l += [pwn_x]
                pwn_gal_b += [pwn_y]
                pwn_sigma += [pwn_r]
                pwn_norm += [pwn_a]
    return pwn_gal_l, pwn_gal_b, pwn_sigma, pwn_norm

lines_to_write = []

center_sky_x = 0.
center_sky_y = 0.
delta_sky_x = 12./2.
delta_sky_y = 3./2.

fig, ax = plt.subplots()
figsize_x = delta_sky_x*2
figsize_y = delta_sky_y*2
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

#layout = 'm4'
#offset = '5deg'
#folder_name = 'CTAsim_%s_80hr_%s_30pc'%(layout,offset)
#bkg_folder_name = 'CTAsim_%s_80hr_%s_0pc'%(layout,offset)

#layout = 'mult3_zd20_m4interp'
#layout = 'mult3_zd20_f4'
layout = 'mult3_zd20_m2'
#layout = 'mult3_zd20_c0'
#layout = 'mult3_zd20_c1'
#exposure_hour = '0'
#exposure_hour = '1'
exposure_hour = '2'
#exposure_hour = '50'
#energy_threshold = 'le'
energy_threshold = 'he'
folder_name = 'CTAsim_gps_%s_%shr_pwn'%(layout,exposure_hour)
bkg_folder_name = 'CTAsim_gps_%s_%shr_bkg'%(layout,exposure_hour)

plotname = 'Plot2D_%s_%s'%(folder_name,energy_threshold)
filename = 'output_tehanu/%s/cntcube_emin=0.20_emax=10.00_side=12.00.fits'%(folder_name)
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
image_data = hdu.data

bkg_plotname = 'Plot2D_%s'%(bkg_folder_name)
bkg_filename = 'output_tehanu/%s/cntcube_emin=0.20_emax=10.00_side=12.00.fits'%(bkg_folder_name)
bkg_hdu = fits.open(bkg_filename)[0]
bkg_wcs = WCS(bkg_hdu.header)
bkg_image_data = bkg_hdu.data

# This is a three-dimensional dataset which you can check by looking at the header information by:
print ("hdu.header")
print (hdu.header)
print ("wcs")
print (wcs)

emin = 0
if energy_threshold=='he':
    emin = 17
emax = 50

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
print ("wcs.dropaxis(2)")
print (wcs.dropaxis(2))

axbig = plt.subplot(projection=wcs.dropaxis(2))
axbig.imshow(image_data_reduced[:,:])
plt.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()


nbins_x = end_pix_x - start_pix_x
nbins_y = end_pix_y - start_pix_y
profile_skymap_1d = []
profile_skymap_1d_bkg = []
axis_skymap_1d = []
for binx in range(0,nbins_x):
    cnt_total = 0.
    cnt_bkg = 0.
    for biny in range(0,nbins_y):
        cnt_total += image_data_reduced[biny,binx]
        cnt_bkg += bkg_image_data_reduced[biny,binx]
    lon, lat, energy = wcs.all_pix2world(binx,0 , 0, 1)
    if lon>180.:
        lon = lon-360.
    profile_skymap_1d += [cnt_total]
    profile_skymap_1d_bkg += [cnt_bkg]
    axis_skymap_1d += [-lon]

figsize_x = delta_sky_x*2
figsize_y = delta_sky_y*2
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

fig.clf()
axbig = fig.add_subplot()
axbig.plot(axis_skymap_1d,profile_skymap_1d,label='signal+bkg')
axbig.plot(axis_skymap_1d,profile_skymap_1d_bkg,label='bkg')
axbig.set_xlabel('Gal. l [deg]')
axbig.set_ylabel('count')
axbig.legend(loc='best')
plotname = 'EmissionProfile_%s_%s'%(folder_name,energy_threshold)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()


image_data_zscore = np.full((image_data[start_pix_e, :, :].shape),0.)
for binx in range(0,nbins_x):
    for biny in range(0,nbins_y):
        cnt_total = image_data_reduced[biny,binx]
        cnt_bkg = bkg_image_data_reduced[biny,binx]
        zscore = 0.
        if cnt_bkg>0.:
            zscore = (cnt_total-cnt_bkg)/pow(cnt_bkg,0.5)
        image_data_zscore[biny,binx] = zscore

file_path = 'model_galactic_pwn_unispec.xml'
#file_path = 'model_galactic_pwn.xml'
src_gal_l, src_gal_b, src_sigma, src_norm = read_pwn_model_file(file_path,center_sky_x,center_sky_y,delta_sky_x,delta_sky_y)
src_zscore = []
n_detect_src = 0.
for src in range(0,len(src_gal_l)):
    src_data_cnt = 0.
    src_bkgd_cnt = 0.
    for binx in range(0,nbins_x):
        for biny in range(0,nbins_y):
            lon, lat, energy = wcs.all_pix2world(binx, biny, 0, 1)
            delta_x = abs(lon-src_gal_l[src])
            if delta_x>180.:
                delta_x = 360.-delta_x
            delta_y = abs(lat-src_gal_b[src])
            dist_to_pwn = pow(pow(delta_x,2)+pow(delta_y,2),0.5)
            #if dist_to_pwn>2.*max(src_sigma[src],0.1): continue
            #if dist_to_pwn>1.*max(src_sigma[src],0.1): continue
            if dist_to_pwn>0.06: continue
            src_data_cnt += image_data_reduced[biny,binx]
            src_bkgd_cnt += bkg_image_data_reduced[biny,binx]
    tmp_src_zscore = 0.
    if src_bkgd_cnt>0.:
        print ('+++++++++++++++++++++++++++++++++++++++++')
        print ('src_id = %s'%(src))
        print ('src_gal_l = %0.2f deg'%(src_gal_l[src]))
        print ('radius cut = %0.2f deg'%(2.*max(src_sigma[src],0.06)))
        print ('src_data_cnt-src_bkgd_cnt = %0.1f'%(src_data_cnt-src_bkgd_cnt))
        print ('src_bkgd_cnt = %0.1f'%(src_bkgd_cnt))
        tmp_src_zscore = (src_data_cnt-src_bkgd_cnt)/pow(src_bkgd_cnt,0.5)
    src_zscore += [tmp_src_zscore]
    if tmp_src_zscore>5.:
        n_detect_src += 1.
print ('n candidate = %s'%(len(src_gal_l)))
print ('n_detect_src = %s'%(n_detect_src))
lines_to_write += ['n candidate = %s'%(len(src_gal_l))]
lines_to_write += ['n_detect_src = %s'%(n_detect_src)]

sorted_src = []
for src in range(0,len(src_gal_l)):
    sorted_src += [(src_gal_l[src],src_gal_b[src],src_sigma[src],src_zscore[src])]
#sorted_src = sorted(sorted_src, key=itemgetter(3),reverse=True)

pwn_index = []
pwn_marker = []
for src in range(0,len(sorted_src)):
    pwn_marker += [(sorted_src[src][0],sorted_src[src][1])]
    pwn_index += ['%s'%(src)]

figsize_x = delta_sky_x*2
figsize_y = delta_sky_y*2
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)
fig.clf()
axbig = fig.add_subplot()
axbig = plt.subplot(projection=wcs.dropaxis(2))
max_z = 5.
min_z = 0.
im = axbig.imshow(image_data_zscore[:,:],vmin=min_z,vmax=max_z)
for star in range(0,len(pwn_marker)):
    print (pwn_marker[star])
    pwn_pixs = wcs.all_world2pix(pwn_marker[star][0],pwn_marker[star][1],0.0,1)
    axbig.scatter(pwn_pixs[0], pwn_pixs[1], c='r', marker='+', label=pwn_index[star])
    text_offset_y = 3
    text_offset_x = 3
    plt.annotate(pwn_index[star], (pwn_pixs[0]+text_offset_x, pwn_pixs[1]+text_offset_y), fontsize=10, color='r')
#cbar = fig.colorbar(im,orientation="horizontal")
#cbar.set_label('Test Statistics')
plotname = 'Plot2D_%s_%s'%(folder_name,energy_threshold)
plt.savefig("output_plots/%s_zscore.png"%(plotname),bbox_inches='tight')
axbig.remove()

#map_nbins = hist_skymap_2d_zscore.GetNbinsX()
#MapEdge_left = hist_skymap_2d_zscore.GetXaxis().GetBinLowEdge(1)
#MapEdge_right = hist_skymap_2d_zscore.GetXaxis().GetBinLowEdge(hist_skymap_2d_zscore.GetNbinsX()+1)
#MapEdge_lower = hist_skymap_2d_zscore.GetYaxis().GetBinLowEdge(1)
#MapEdge_upper = hist_skymap_2d_zscore.GetYaxis().GetBinLowEdge(hist_skymap_2d_zscore.GetNbinsY()+1)
#x_axis = np.linspace(MapEdge_left,MapEdge_right,map_nbins)
#y_axis = np.linspace(MapEdge_lower,MapEdge_upper,map_nbins)
#grid_z = np.zeros((map_nbins, map_nbins))
#for ybin in range(0,len(y_axis)):
#    for xbin in range(0,len(x_axis)):
#        hist_bin_x = xbin+1
#        hist_bin_y = ybin+1
#        if hist_bin_x<1: continue
#        if hist_bin_y<1: continue
#        if hist_bin_x>hist_skymap_2d_zscore.GetNbinsX(): continue
#        if hist_bin_y>hist_skymap_2d_zscore.GetNbinsY(): continue
#        grid_z[ybin,xbin] = hist_skymap_2d_zscore.GetBinContent(hist_bin_x,hist_bin_y)
#
#axbig = fig.add_subplot()
#im = axbig.imshow(grid_z, origin='lower', extent=(x_axis.min(),x_axis.max(),y_axis.min(),y_axis.max()),zorder=0,vmin=0,vmax=5)
#for star in range(0,len(pwn_marker)):
#    print (pwn_marker[star])
#    axbig.scatter(pwn_marker[star][0], pwn_marker[star][1], c='r', marker='+', label=pwn_index[star])
#    text_offset_y = 0
#    text_offset_x = 0.05
#    plt.annotate(pwn_index[star], (pwn_marker[star][0]+text_offset_x, pwn_marker[star][1]+text_offset_y), fontsize=10, color='r')
##cbar = fig.colorbar(im,orientation="horizontal")
##cbar.set_label('Test Statistics')
#plt.savefig("output_plots/%s_zscore_hist.png"%(plotname),bbox_inches='tight')
#axbig.remove()

image_data_residual = np.full((image_data[start_pix_e, :, :].shape),0.)
image_data_residual += image_data_reduced
image_data_residual -= bkg_image_data_reduced

print ('set initial avlues and bounds')
initial_prms = []
bound_upper_prms = []
bound_lower_prms = []
for src in range(0,len(sorted_src)):
    lon = sorted_src[src][0]
    lat = sorted_src[src][1]
    if lon>180.:
        lon = lon-360.
    sigma = sorted_src[src][2]
    initial_prms += [(-lon,lat,sigma,100.)]
    bound_lower_prms += [(-lon-0.5,lat-0.5,0.,0.)]
    bound_upper_prms += [(-lon+0.5,lat+0.5,10.,1e6)]
# Flatten the initial guess parameter list.
p0 = [p for prms in initial_prms for p in prms]
p0_lower = [p for prms in bound_lower_prms for p in prms]
p0_upper = [p for prms in bound_upper_prms for p in prms]
#print ('p0')
#print (p0)
#print ('p0_lower')
#print (p0_lower)
#print ('p0_upper')
#print (p0_upper)

print ('fit 2D map')
lon_min, lat_min, energy = wcs.all_pix2world(0, 0, 0, 1)
if lon_min>180.:
    lon_min = lon_min-360.
lon_max, lat_max, energy = wcs.all_pix2world(nbins_x-1, nbins_y-1, 0, 1)
if lon_max>180.:
    lon_max = lon_max-360.
x_axis = np.linspace(lon_max,lon_min,nbins_x)
y_axis = np.linspace(lat_min,lat_max,nbins_y)
X_grid, Y_grid = np.meshgrid(x_axis, y_axis)
# We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
XY_stack = np.vstack((X_grid.ravel(), Y_grid.ravel()))
popt, pcov = curve_fit(_gaussian, XY_stack, image_data_residual.ravel(), p0, bounds=(p0_lower,p0_upper))
#popt, pcov = curve_fit(_gaussian, XY_stack, image_data_residual.ravel(), p0)
#print ('popt')
#print (popt)

n_detect_fit = 0.
for i in range(len(popt)//4):
    pwn_lon = popt[i*4+0]
    pwn_lat = popt[i*4+1]
    pwn_rad = popt[i*4+2]
    pwn_amp = popt[i*4+3]
    pwn_lon_err = pow(pcov[i*4+0][i*4+0],0.5)
    pwn_lat_err = pow(pcov[i*4+1][i*4+1],0.5)
    pwn_rad_err = pow(pcov[i*4+2][i*4+2],0.5)
    pwn_amp_err = pow(pcov[i*4+3][i*4+3],0.5)
    lines_to_write += ['++++++++++++++++++++++++++++++++++++++++++++++++++++++']
    lines_to_write += ['PWN = %s'%(i)]
    lines_to_write += ['pwn_lon = %0.5f'%(pwn_lon)]
    lines_to_write += ['pwn_lat = %0.5f'%(pwn_lat)]
    lines_to_write += ['pwn_rad = %0.5f'%(pwn_rad)]
    lines_to_write += ['pwn_amp = %0.5f'%(pwn_amp)]
    lines_to_write += ['pwn_lon_err = %0.5f'%(pwn_lon_err)]
    lines_to_write += ['pwn_lat_err = %0.5f'%(pwn_lat_err)]
    lines_to_write += ['pwn_rad_err = %0.5f'%(pwn_rad_err)]
    lines_to_write += ['pwn_amp_err = %0.5f'%(pwn_amp_err)]
    lines_to_write += ['src_gal_l = %0.5f'%(src_gal_l[i])]
    lines_to_write += ['src_gal_b = %0.5f'%(src_gal_b[i])]
    lines_to_write += ['src_sigma = %0.5f'%(src_sigma[i])]
    lines_to_write += ['src_norm = %0.5f'%(src_norm[i])]
    zscore_fit = 0.
    if pwn_amp_err>0.:
        zscore_fit = pwn_amp/pwn_amp_err
    if zscore_fit>5.:
        n_detect_fit += 1.
lines_to_write += ['++++++++++++++++++++++++++++++++++++++++++++++++++++++']
lines_to_write += ['n_detect_fit = %s'%(n_detect_fit)]

fit_2d = np.zeros(image_data_residual.shape)
for i in range(len(popt)//4):
    fit_2d += gaussian(X_grid, Y_grid, *popt[i*4:i*4+4])

print ('plot 2D map')
figsize_x = delta_sky_x*2
figsize_y = delta_sky_y*2
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)
fig.clf()
axbig = fig.add_subplot()
max_z = 5.
min_z = 0.
for binx in range(0,nbins_x):
    for biny in range(0,nbins_y):
        cnt_src = image_data_residual[biny,binx]
        cnt_bkg = bkg_image_data_reduced[biny,binx]
        zscore = 0.
        if cnt_bkg>0.:
            zscore = (cnt_src)/pow(cnt_bkg,0.5)
        image_data_zscore[biny,binx] = zscore
axbig.imshow(image_data_zscore, origin='lower', cmap='plasma', extent=(x_axis.min(), x_axis.max(), y_axis.min(), y_axis.max()))
for star in range(0,len(pwn_marker)):
    axbig.scatter(-pwn_marker[star][0], pwn_marker[star][1], c='r', marker='+', label=pwn_index[star])
    text_offset_y = 0.1
    text_offset_x = 0.1
    plt.annotate(pwn_index[star], (-pwn_marker[star][0]-text_offset_x, pwn_marker[star][1]+text_offset_y), fontsize=10, color='r')
axbig.contour(X_grid, Y_grid, fit_2d, colors='w')
plotname = 'Plot2D_%s_%s'%(folder_name,energy_threshold)
plt.savefig("output_plots/%s_zscore_fit.png"%(plotname),bbox_inches='tight')
axbig.remove()

profile_data_1d = []
profile_fit_1d = []
axis_skymap_1d = []
for binx in range(0,nbins_x):
    cnt_total = 0.
    cnt_fit = 0.
    for biny in range(0,nbins_y):
        cnt_total += image_data_residual[biny,binx]
        cnt_fit += fit_2d[biny,binx]
    lon, lat, energy = wcs.all_pix2world(binx,0 , 0, 1)
    if lon>180.:
        lon = lon-360.
    profile_data_1d += [cnt_total]
    profile_fit_1d += [cnt_fit]
    axis_skymap_1d += [-lon]
profile_pwn_1d = []
for src in range(0,len(sorted_src)):
    fit_src_2d = np.zeros(image_data_residual.shape)
    fit_src_2d += gaussian(X_grid, Y_grid, *popt[src*4:src*4+4])
    profile_onepwn_1d = []
    for binx in range(0,nbins_x):
        cnt_pwn = 0.
        for biny in range(0,nbins_y):
            cnt_pwn += fit_src_2d[biny,binx]
        profile_onepwn_1d += [cnt_pwn]
    profile_pwn_1d += [profile_onepwn_1d]

figsize_x = delta_sky_x*2
figsize_y = delta_sky_y*2
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

fig.clf()
axbig = fig.add_subplot()
axbig.plot(axis_skymap_1d,profile_data_1d,label='data',color='k')
axbig.plot(axis_skymap_1d,profile_fit_1d,label='fit',color='r')
for src in range(0,len(sorted_src)):
    axbig.plot(axis_skymap_1d,profile_pwn_1d[src],alpha=0.3)
axbig.set_xlabel('Gal. l [deg]')
axbig.set_ylabel('count')
axbig.set_ylim(bottom=1.)
axbig.set_yscale('log')
axbig.legend(loc='best')
plotname = 'ModelProfile_%s_%s'%(folder_name,energy_threshold)
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

with open('output_plots/result_%s_emin%s.txt'%(folder_name,emin), 'w') as f:
    for line in lines_to_write:
        f.write(line)
        f.write('\n')

n_bins_r = 20
radial_range = 1.0
for src in range(0,len(sorted_src)):
    print ('Making radial plot for PWN %s'%(src))
    fit_src_2d = np.zeros(image_data_residual.shape)
    fit_src_2d += gaussian(X_grid, Y_grid, *popt[src*4:src*4+4])
    radial_axis = []
    evts_in_slice = []
    pwns_in_slice = []
    mods_in_slice = []
    errs_in_slice = []
    bins_in_slice = []
    for rb in range(0,n_bins_r):
        radial_axis += [float(rb)*radial_range/float(n_bins_r)]
        evts_in_slice += [0.]
        pwns_in_slice += [0.]
        mods_in_slice += [0.]
        errs_in_slice += [0.]
        bins_in_slice += [0.]
        for binx in range(0,nbins_x):
            for biny in range(0,nbins_y):
                lon, lat, energy = wcs.all_pix2world(binx, biny, 0, 1)
                pwn_cnt = fit_src_2d[biny,binx]
                mod_cnt = fit_2d[biny,binx]
                src_cnt = image_data_residual[biny,binx]
                bkg_cnt = bkg_image_data_reduced[biny,binx]
                delta_x = abs(lon-sorted_src[src][0])
                if delta_x>180.:
                    delta_x = 360.-delta_x
                delta_y = abs(lat-sorted_src[src][1])
                dist_to_pwn = pow(pow(delta_x,2)+pow(delta_y,2),0.5)
                if dist_to_pwn>=radial_axis[rb] and dist_to_pwn<radial_axis[rb]+radial_range/float(n_bins_r):
                    pwns_in_slice[rb] += (pwn_cnt)
                    mods_in_slice[rb] += (mod_cnt)
                    evts_in_slice[rb] += (src_cnt)
                    errs_in_slice[rb] += (bkg_cnt)
                    bins_in_slice[rb] += 1.
        if bins_in_slice[rb]>0.:
            pwns_in_slice[rb] = pwns_in_slice[rb]/(bins_in_slice[rb])
            mods_in_slice[rb] = mods_in_slice[rb]/(bins_in_slice[rb])
            evts_in_slice[rb] = evts_in_slice[rb]/(bins_in_slice[rb])
            errs_in_slice[rb] = pow(errs_in_slice[rb],0.5)/(bins_in_slice[rb])

    radial_axis = np.array(radial_axis)
    pwns_in_slice = np.array(pwns_in_slice)
    mods_in_slice = np.array(mods_in_slice)
    evts_in_slice = np.array(evts_in_slice)
    errs_in_slice = np.array(errs_in_slice)

    figsize_x = 6.4
    figsize_y = 4.8
    fig.set_figheight(figsize_y)
    fig.set_figwidth(figsize_x)

    fig.clf()
    axbig = fig.add_subplot()
    axbig.errorbar(radial_axis,evts_in_slice,errs_in_slice,color='k',marker='s',ls='none',label='data')
    axbig.plot(radial_axis,mods_in_slice,label='all PWNe')
    axbig.plot(radial_axis,pwns_in_slice,label='center PWN')
    axbig.set_xlabel('radial distance to PWN [deg]')
    axbig.set_ylabel('count / pixel')
    axbig.legend(loc='best')
    plotname = 'RadialProfile_PWN%s_%s_%s'%(src,folder_name,energy_threshold)
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

