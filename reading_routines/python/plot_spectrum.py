########################################################################3
#
#   Plotting Example:  Shell_Spectra
#
#   We plot the kinetic energy spectrum from one
#   shell spectrum output , modulo a factor of 2 rho
##################################

from diagnostic_reading import ShellSpectra
import matplotlib.pyplot as plt
import numpy as np
a = ShellSpectra(filename='01330000',path='Shell_Spectra/')

# We use the lookup table to find where vr, vtheta, and vphi are stored
vr_index = a.lut[1]
vt_index = a.lut[2]
vp_index = a.lut[3]
rad_ind = 3
t_ind = 0

# Next we grab one radial index and one time instance of each variable
dims = (a.nell,a.nm)
vrc = np.reshape(a.vals[:,:,rad_ind, vr_index, t_ind],dims)
vtc = np.reshape(a.vals[:,:,rad_ind, vt_index, t_ind],dims)
vpc = np.reshape(a.vals[:,:,rad_ind, vp_index, t_ind],dims)


pspec  = np.zeros(a.nell,dtype='float64')

for l in range(a.nell):
    for m in range(a.nm):
        pspec[l] = pspec[l]+np.real(vrc[l,m])**2 +np.imag(vrc[l,m])**2
        pspec[l] = pspec[l]+np.real(vtc[l,m])**2 +np.imag(vtc[l,m])**2
        pspec[l] = pspec[l]+np.real(vpc[l,m])**2 +np.imag(vpc[l,m])**2

#time_index = 0
plt.figure(1)
plt.plot(pspec[1:])
plt.yscale('log')
plt.xscale('log')
plt.xlim(xmin=1,xmax=a.lmax-1)
plt.xlabel('Spherical Harmonic Degree '+r'$\ell$')
plt.ylabel('Kinetic Energy/2'+r"$\rho$")

plt.savefig('KE_spectrum.png')
