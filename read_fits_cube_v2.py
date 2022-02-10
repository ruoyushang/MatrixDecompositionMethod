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

filename = 'input_data/DHT08_Quad1_interp.fits'
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
image_data = hdu.data

# This is a three-dimensional dataset which you can check by looking at the header information by:
print (hdu.header)
print (wcs)
#Number of WCS axes: 3
#CTYPE : 'VOPT'  'GLON-CAR'  'GLAT-CAR'
#CRVAL : -90.2706  75.0  0.0
#CRPIX : 1.0  1.0  65.0
#CDELT : 0.65019  -0.125  0.125
#NAXIS : 417  497  113

#print (wcs.all_world2pix(-90.2706,75.0,0.,1))

#vel_low = 0.
#vel_up = 25.
vel_low = 25.
vel_up = 50.
#vel_low = 50.
#vel_up = 75.
center_lon = 40.5
center_lat = -0.5
#center_lon = 41.1
#center_lat = -0.3
delta_lon = 1.5
delta_lat = 1.5

pixs_start = wcs.all_world2pix(vel_low,center_lon+delta_lon,center_lat-delta_lat,1)
pixs_end = wcs.all_world2pix(vel_up,center_lon-delta_lon,center_lat+delta_lat,1)
print (pixs_start)
print (pixs_end)
vel_idx_start = int(pixs_start[0])
vel_idx_end = int(pixs_end[0])
gal_x_idx_start = int(pixs_start[1])
gal_x_idx_end = int(pixs_end[1])
gal_y_idx_start = int(pixs_start[2])
gal_y_idx_end = int(pixs_end[2])
vel_idx_center = int((pixs_start[0]+pixs_end[0])/2.)
vel_idx_size = int((pixs_end[0]-pixs_start[0])/2.)
gal_x_idx_center = int((pixs_start[1]+pixs_end[1])/2.)
gal_x_idx_size = int((pixs_end[1]-pixs_start[1])/2.)
gal_y_idx_center = int((pixs_start[2]+pixs_end[2])/2.)
gal_y_idx_size = int((pixs_end[2]-pixs_start[2])/2.)

print ((image_data[:, :, vel_idx_start].shape))
print ('vel_idx_start = %s'%(vel_idx_start))
print ('vel_idx_end = %s'%(vel_idx_end))
image_data_reduced = np.full((image_data[:, :, vel_idx_start].shape),0.)
for idx in range(vel_idx_start,vel_idx_end):
    image_data_reduced += image_data[:, :, idx]
print ('gal_x_idx_center = %s'%(gal_x_idx_center))
print ('gal_y_idx_center = %s'%(gal_y_idx_center))
print ('gal_x_idx_size = %s'%(gal_x_idx_size))
print ('gal_y_idx_size = %s'%(gal_y_idx_size))
position = (gal_x_idx_center,gal_y_idx_center)
size = (gal_x_idx_size,gal_y_idx_size)
print (wcs.dropaxis(0))
image_data_cutout = Cutout2D(image_data_reduced, position=position, size=size, wcs=wcs.dropaxis(0))
ax = plt.subplot(projection=image_data_cutout.wcs)
ax.imshow(image_data_cutout.data[:,:])
plotname = 'Plot2D_GalXY_reduced'
plt.savefig("output/%s.png"%(plotname),bbox_inches='tight')

# If we wanted to plot the spectral axes for one pixel we can do this by slicing down to one dimension.
ax = plt.subplot(projection=wcs, slices=('x', 277, 57))
print (image_data[gal_y_idx_start, gal_x_idx_start, :].shape)
image_data_reduced = np.full((image_data[gal_y_idx_start, gal_x_idx_start, :].shape),0.)
print ('gal_x_idx_start = %s'%(gal_x_idx_start))
print ('gal_x_idx_end = %s'%(gal_x_idx_end))
for x_idx in range(gal_x_idx_start,gal_x_idx_end):
    for y_idx in range(gal_y_idx_start,gal_y_idx_end):
        image_data_reduced += image_data[y_idx, x_idx, :]
ax.plot(image_data_reduced[:])
# As this is still a WCSAxes plot, we can set the display units for the x-axis
#vel, gal_lon, gal_lat = ax.coords
#vel.set_format_unit(u.m/u.s)
plotname = 'Plot1D'
plt.savefig("output/%s.png"%(plotname),bbox_inches='tight')


