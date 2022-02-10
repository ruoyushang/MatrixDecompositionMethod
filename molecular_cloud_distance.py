
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def rotation_curve(r):
    # r is galactic radius in kpc, V in km/s
    # sun R0 = 7.5 kpc, sun V0 = 190 km/s
    V = 228.*math.exp(-r/50.-pow(3.6/r,2))+350.*math.exp(-r/3.25-0.1/r)
    return V
def circular_solution(r):
    # r is galactic radius in kpc, V in km/s
    V = r/R0*(V_LSR/(math.sin(gal_l*3.14/180.)*math.cos(gal_b*3.14/180.))+V0)
    return V
def difference_velocity(r):
    dV = rotation_curve(r)-circular_solution(r)
    return dV

def find_cloud_distance(r):
    Z = (r*r-R0*R0)/pow(math.cos(gal_b*3.14/180.),2)+R0*R0*pow(math.cos(gal_l*3.14/180.)/math.cos(gal_b*3.14/180.),2)
    d_plus = R0*math.cos(gal_l*3.14/180.)/math.cos(gal_b*3.14/180.) + pow(Z,0.5)
    d_minus = R0*math.cos(gal_l*3.14/180.)/math.cos(gal_b*3.14/180.) - pow(Z,0.5)
    return d_plus,d_minus

R0 = 7.5 # kpc
V0 = 190. # km/s

gal_l = 40.5
gal_b = -0.5

#V_LSR = 25.
#V_LSR = 50.
V_LSR = 70.

#gal_l = 78.2
#gal_b = 2.1
#V_LSR = 0.

gal_radius = np.linspace(0.1, 15., 100)
vectorize_rotation_curve_velocity = np.vectorize(rotation_curve)
rotation_curve_velocity = vectorize_rotation_curve_velocity(gal_radius)
vectorize_circular_velocity = np.vectorize(circular_solution)
circular_solution_velocity = vectorize_circular_velocity(gal_radius)

root_gal_radius = optimize.newton(difference_velocity, 7.0)
print ('root_gal_radius = %0.3f kpc'%(root_gal_radius))

d_plus,d_minus = find_cloud_distance(root_gal_radius)
print ('d_plus = %0.3f kpc'%(d_plus))
print ('d_minus = %0.3f kpc'%(d_minus))

fig, ax = plt.subplots()

fig.clf()
axbig = fig.add_subplot()
axbig.plot(gal_radius,rotation_curve_velocity,color='b')
axbig.plot(gal_radius,circular_solution_velocity,color='r')
axbig.set_xlabel('galactic radius (kpc)')
axbig.set_ylabel('Velocity (km/s)')
plotname = 'CO_rotation_curve'
fig.savefig("output_plots/%s.png"%(plotname),bbox_inches='tight')
axbig.remove()

