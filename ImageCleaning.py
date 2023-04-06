
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

def gaussian(x, y, x0, y0, sigma, cos_theta, A):
    #return A * np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma*sigma))/(2*np.pi*sigma*sigma)
    sin_theta = pow(1. - cos_theta*cos_theta, 0.5)
    x_rot = cos_theta*(x-x0) + sin_theta*(y-y0)
    y_rot = -sin_theta*(x-x0) + cos_theta*(y-y0)
    return A * np.exp(-((x_rot)**2/(2*sigma*sigma)+(y_rot)**2/(2*0.4*sigma*0.4*sigma)))

def fft_filter(image_data,image_control):

    # Compute 2D Fourier transform
    F = fftpack.fft2(image_data)
    F_control = fftpack.fft2(image_control)
    
    # Shift the zero-frequency component to the center
    F_shifted = fftpack.fftshift(F)
    F_control_shifted = fftpack.fftshift(F_control)
    
    rows, cols = image_data.shape
    for row in range(0,rows):
        for col in range(0,cols):
            F_shifted[row,col] = F_shifted[row,col]-F_control_shifted[row,col]

    # Filter out high-frequency components
    n = 5  # Number of frequency components to keep
    crow, ccol = rows // 2, cols // 2  # Center coordinates
    F_shifted[crow-n:crow+n+1, ccol-n:ccol+n+1] = 0
    
    # Shift the zero-frequency component back to the corner
    F_filtered = fftpack.ifftshift(F_shifted)
    
    # Compute the inverse Fourier transform to get the filtered image
    image_fft_cleaned = image_data - np.real(fftpack.ifft2(F_filtered))

    return image_fft_cleaned


image_size = 25
image_signal = np.zeros((image_size, image_size))
image_data = np.zeros((image_size, image_size))
image_shape = np.shape(image_signal)

# Generate a 2D array of shape (n, m) with random numbers that follow a Poisson distribution with lambda = 3
NSB = 7
image_noise = np.random.poisson(NSB, size=(image_size, image_size))
image_control = np.random.poisson(NSB, size=(image_size, image_size))

max_z = 20.

x0 = float(image_size/2.)
y0 = float(image_size/2.)+5
sigma = 2.
A = 10.
cos_theta = np.random.uniform(0, 1)
for binx in range(0,image_shape[0]):
    for biny in range(0,image_shape[1]):
        image_signal[binx,biny] = gaussian(float(binx), float(biny), x0, y0, sigma, cos_theta, A)

for p in range (0,10):
    random_int_x = np.random.randint(0, image_size-1)
    random_int_y = np.random.randint(0, image_size-1)
    image_data[random_int_x,random_int_y] = max_z

image_data = image_data+image_signal+image_noise
avg_nsb = np.mean(image_control)
image_data_nonsb = image_data-avg_nsb


image_fft_cleaned = fft_filter(image_data_nonsb,image_control)

k = 10
U_signal, S_signal, V_signal = np.linalg.svd(image_fft_cleaned,full_matrices=False)
image_svd_cleaned = U_signal[:, :k] @ np.diag(S_signal[:k]) @ V_signal[:k, :]


fig, ax = plt.subplots()
figsize_x = 8.6
figsize_y = 6.4
fig.set_figheight(figsize_y)
fig.set_figwidth(figsize_x)

im = ax.imshow(image_signal, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_signal.png",bbox_inches='tight')

im = ax.imshow(image_data, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_data.png",bbox_inches='tight')

im = ax.imshow(image_control, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_control.png",bbox_inches='tight')

im = ax.imshow(image_data_nonsb, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_data_nonsb.png",bbox_inches='tight')

im = ax.imshow(image_fft_cleaned, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_fft_cleaned.png",bbox_inches='tight')

im = ax.imshow(image_svd_cleaned, origin='lower', vmin=0., vmax=max_z)
fig.savefig("output_plots/image_svd_cleaned.png",bbox_inches='tight')

